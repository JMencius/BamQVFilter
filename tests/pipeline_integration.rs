//! Black-box regression tests for the bounded, order-preserving pipeline.
//!
//! These tests deliberately use only the project's normal dependencies.

use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::Read;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

static NEXT_TEST_ID: AtomicU64 = AtomicU64::new(0);

#[test]
fn repeated_reader_writer_handoff_is_memory_safe() {
    let tmp = TestDir::new("record-handoff");
    let input = tmp.path("input.bam");
    let mut expected_names = Vec::new();

    {
        let mut writer = bam::Writer::from_path(&input, &header(), bam::Format::Bam).unwrap();
        for index in 0..2_000 {
            let quality = if index % 3 == 0 { 0 } else { 20 };
            write_record(&mut writer, index, [1, 31, 997, 127][index % 4], quality);
            if quality >= 10 {
                expected_names.push(format!("r{index:06}").into_bytes());
            }
        }
    }

    // rust-htslib 0.46 reader Records carry an Rc-backed header. Repeating a
    // small-batch, mixed keep/drop run makes an unsafe cross-thread hand-off
    // fail quickly under allocator checks. The production pipeline sends deep,
    // header-detached Record clones instead.
    for repetition in 0..12 {
        let output = tmp.path(&format!("output-{repetition}.bam"));
        let ratio = if repetition % 2 == 0 { "0" } else { "0.5" };
        let result = run_filter(&input, &output, 4, ratio, 10);
        assert!(
            result.status.success(),
            "repetition {repetition} failed: {}",
            String::from_utf8_lossy(&result.stderr)
        );

        let mut reader = bam::Reader::from_path(output).unwrap();
        let actual_names: Vec<_> = reader
            .records()
            .map(|record| record.unwrap().qname().to_vec())
            .collect();
        assert_eq!(actual_names, expected_names);
    }
}

struct TestDir(PathBuf);

impl TestDir {
    fn new(label: &str) -> Self {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let id = NEXT_TEST_ID.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "bamqvfilter-{label}-{}-{nonce}-{id}",
            std::process::id()
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn header() -> bam::Header {
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 10_000_000),
    );
    header
}

fn write_record(writer: &mut bam::Writer, index: usize, sequence_len: usize, quality: u8) {
    write_record_with_flags(writer, index, sequence_len, quality, 0);
}

fn write_record_with_flags(
    writer: &mut bam::Writer,
    index: usize,
    sequence_len: usize,
    quality: u8,
    flags: u16,
) {
    let name = format!("r{index:06}");
    let sequence = vec![b'A'; sequence_len];
    let qualities = vec![quality; sequence_len];
    let cigar = CigarString(vec![Cigar::Match(sequence_len as u32)]);
    let mut record = bam::Record::new();
    let cigar = if flags & 0x4 == 0 { Some(&cigar) } else { None };
    record.set(name.as_bytes(), cigar, &sequence, &qualities);
    record.set_flags(flags);
    if flags & 0x4 == 0 {
        record.set_tid(0);
        record.set_pos((index * 2) as i64);
        record.set_mapq(60);
    } else {
        record.set_tid(-1);
        record.set_pos(-1);
        record.set_mapq(0);
    }
    record.push_aux(b"XI", Aux::I32(index as i32)).unwrap();
    writer.write(&record).unwrap();
}

fn run_filter_with_limits(
    input: &Path,
    output: &Path,
    threads: usize,
    ratio: &str,
    qv: u8,
    batch_records: usize,
    batch_mib: usize,
) -> Output {
    run_filter_with_mode(
        input,
        output,
        threads,
        ratio,
        qv,
        batch_records,
        batch_mib,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_filter_with_mode(
    input: &Path,
    output: &Path,
    threads: usize,
    ratio: &str,
    qv: u8,
    batch_records: usize,
    batch_mib: usize,
    primary_option: Option<&str>,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_bamqvfilter"));
    command
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output)
        .arg("-q")
        .arg(qv.to_string())
        .arg("-t")
        .arg(threads.to_string())
        .arg("--compression-thread-ratio")
        .arg(ratio)
        .arg("--batch-records")
        .arg(batch_records.to_string())
        .arg("--batch-mib")
        .arg(batch_mib.to_string());

    if let Some(option) = primary_option {
        command.arg(option);
    }

    command.output().unwrap()
}

fn run_filter(input: &Path, output: &Path, threads: usize, ratio: &str, qv: u8) -> Output {
    // A tiny record limit forces many batch boundaries in a cheap fixture.
    run_filter_with_limits(input, output, threads, ratio, qv, 37, 1)
}

fn final_stderr_line(output: &Output) -> &str {
    std::str::from_utf8(&output.stderr)
        .unwrap()
        .lines()
        .rfind(|line| !line.trim().is_empty())
        .unwrap_or("")
}

#[test]
fn alignment_types_are_kept_by_default_and_primary_only_is_mapped_primary() {
    let tmp = TestDir::new("alignment-types");
    let input = tmp.path("input.bam");
    let all_output = tmp.path("all.bam");
    let primary_output = tmp.path("primary.bam");
    let primary_alias_output = tmp.path("primary-alias.bam");

    {
        let mut writer = bam::Writer::from_path(&input, &header(), bam::Format::Bam).unwrap();
        write_record_with_flags(&mut writer, 0, 100, 20, 0);
        write_record_with_flags(&mut writer, 1, 100, 20, 0x100);
        write_record_with_flags(&mut writer, 2, 100, 20, 0x800);
        write_record_with_flags(&mut writer, 3, 100, 20, 0x4);
    }

    let all_result = run_filter_with_mode(&input, &all_output, 4, "0.25", 10, 2, 1, None);
    assert!(
        all_result.status.success(),
        "default run failed: {}",
        String::from_utf8_lossy(&all_result.stderr)
    );
    assert_eq!(
        final_stderr_line(&all_result),
        "Kept 4 reads out of 4 reads"
    );

    let mut reader = bam::Reader::from_path(&all_output).unwrap();
    let all_records: Vec<_> = reader.records().map(|record| record.unwrap()).collect();
    assert_eq!(
        all_records
            .iter()
            .map(|record| record.flags())
            .collect::<Vec<_>>(),
        [0, 0x100, 0x800, 0x4]
    );

    let primary_result = run_filter_with_mode(
        &input,
        &primary_output,
        4,
        "0.25",
        10,
        2,
        1,
        Some("--primary-only"),
    );
    assert!(
        primary_result.status.success(),
        "primary-only run failed: {}",
        String::from_utf8_lossy(&primary_result.stderr)
    );
    assert_eq!(
        final_stderr_line(&primary_result),
        "Kept 1 reads out of 4 reads"
    );

    let mut reader = bam::Reader::from_path(primary_output).unwrap();
    let primary_records: Vec<_> = reader.records().map(|record| record.unwrap()).collect();
    assert_eq!(primary_records.len(), 1);
    assert_eq!(primary_records[0].qname(), b"r000000");
    assert_eq!(primary_records[0].flags(), 0);

    let alias_result = run_filter_with_mode(
        &input,
        &primary_alias_output,
        4,
        "0.25",
        10,
        2,
        1,
        Some("--primary"),
    );
    assert!(
        alias_result.status.success(),
        "--primary alias run failed: {}",
        String::from_utf8_lossy(&alias_result.stderr)
    );
    assert_eq!(
        final_stderr_line(&alias_result),
        "Kept 1 reads out of 4 reads"
    );

    let mut reader = bam::Reader::from_path(primary_alias_output).unwrap();
    let alias_names: Vec<_> = reader
        .records()
        .map(|record| record.unwrap().qname().to_vec())
        .collect();
    assert_eq!(alias_names, [b"r000000"]);
}

#[test]
fn output_is_the_exact_filtered_input_subsequence_across_batch_boundary() {
    let tmp = TestDir::new("ordered");
    let input = tmp.path("input.bam");
    let output = tmp.path("output.bam");
    let byte_limited_output = tmp.path("byte-limited-output.bam");
    let input_header = header();

    // More than the proposed default 4096-record batch size. Alternating
    // lengths gives worker jobs very different costs and makes a completion-
    // order writer especially likely to fail this regression test.
    let mut expected_names = Vec::new();
    {
        let mut writer = bam::Writer::from_path(&input, &input_header, bam::Format::Bam).unwrap();
        for index in 0..5_123 {
            let sequence_len = match index % 4 {
                0 => 1,
                1 => 31,
                2 => 997,
                _ => 127,
            };
            let quality = match index % 5 {
                0 => 0,
                1 => 10, // equality with threshold must pass
                2 => 20,
                3 => 9,
                _ => 40,
            };
            write_record(&mut writer, index, sequence_len, quality);
            if quality >= 10 {
                expected_names.push(format!("r{index:06}").into_bytes());
            }
        }

        // Empty sequence/QV: average is undefined and must be rejected.
        let mut empty = bam::Record::new();
        empty.set(b"empty_qual", None, b"", &[]);
        empty.set_tid(0);
        empty.set_pos(20_000);
        writer.write(&empty).unwrap();

        // 255 is BAM's missing-QV sentinel. It must not be interpreted as an
        // exceptionally good read.
        let mut missing = bam::Record::new();
        let missing_seq = vec![b'A'; 20];
        let missing_qual = vec![255; 20];
        let missing_cigar = CigarString(vec![Cigar::Match(20)]);
        missing.set(
            b"missing_qual",
            Some(&missing_cigar),
            &missing_seq,
            &missing_qual,
        );
        missing.set_tid(0);
        missing.set_pos(20_002);
        writer.write(&missing).unwrap();
    }

    let status = run_filter_with_limits(&input, &output, 2, "0", 10, 37, 64);
    assert!(
        status.status.success(),
        "filter failed: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    assert_eq!(
        final_stderr_line(&status),
        "Kept 3074 reads out of 5125 reads"
    );

    let input_reader = bam::Reader::from_path(&input).unwrap();
    let input_header_text = input_reader.header().as_bytes().to_vec();
    drop(input_reader);

    let mut output_reader = bam::Reader::from_path(&output).unwrap();
    assert_eq!(output_reader.header().as_bytes(), input_header_text);
    let records: Vec<_> = output_reader
        .records()
        .map(|result| result.expect("output BAM must be readable"))
        .collect();
    let actual_names: Vec<_> = records
        .iter()
        .map(|record| record.qname().to_vec())
        .collect();
    assert_eq!(actual_names.len(), expected_names.len());
    for (output_index, (actual, expected)) in actual_names.iter().zip(&expected_names).enumerate() {
        assert_eq!(
            actual, expected,
            "record order/content differs at retained output index {output_index}"
        );
    }

    // Verify more than names/order: retained records were not changed while
    // moving across threads.
    for record in records {
        let name = std::str::from_utf8(record.qname()).unwrap();
        let index: usize = name.strip_prefix('r').unwrap().parse().unwrap();
        assert_eq!(record.tid(), 0);
        assert_eq!(record.pos(), (index * 2) as i64);
        assert_eq!(record.mapq(), 60);
        assert_eq!(record.aux(b"XI").unwrap(), Aux::I32(index as i32));
    }

    // Run the same semantic check with a record limit too large to fire, so
    // the 1-MiB allocated-record limit is what creates the batch boundaries.
    let status = run_filter_with_limits(&input, &byte_limited_output, 8, "1", 10, 100_000, 1);
    assert!(
        status.status.success(),
        "byte-limited filter failed: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    assert_eq!(
        final_stderr_line(&status),
        "Kept 3074 reads out of 5125 reads"
    );
    let mut reader = bam::Reader::from_path(byte_limited_output).unwrap();
    let byte_limited_names: Vec<_> = reader
        .records()
        .map(|result| result.unwrap().qname().to_vec())
        .collect();
    assert_eq!(byte_limited_names, expected_names);
}

#[test]
fn probability_mean_qv_is_not_the_arithmetic_mean() {
    let tmp = TestDir::new("qv-formula");
    let input = tmp.path("input.bam");
    let output = tmp.path("output.bam");
    {
        let mut writer = bam::Writer::from_path(&input, &header(), bam::Format::Bam).unwrap();
        let mut record = bam::Record::new();
        let cigar = CigarString(vec![Cigar::Match(2)]);
        // Probability-mean QV is 3.009865683871185; arithmetic mean is 20.
        record.set(b"mixed", Some(&cigar), b"AA", &[0, 40]);
        record.set_tid(0);
        record.set_pos(1);
        writer.write(&record).unwrap();
    }

    let status = run_filter(&input, &output, 4, "0.25", 4);
    assert!(status.status.success());
    let mut reader = bam::Reader::from_path(output).unwrap();
    assert_eq!(reader.records().count(), 0);
}

#[test]
fn one_record_larger_than_byte_limit_is_emitted_once_without_hanging() {
    let tmp = TestDir::new("oversized-record");
    let input = tmp.path("input.bam");
    let output = tmp.path("output.bam");
    {
        let mut writer = bam::Writer::from_path(&input, &header(), bam::Format::Bam).unwrap();
        write_record(&mut writer, 0, 10, 20);
        // Sequence + qualities alone exceed the 1-MiB test batch target.
        write_record(&mut writer, 1, 1_000_000, 20);
        write_record(&mut writer, 2, 10, 20);
    }

    let result = run_filter_with_limits(&input, &output, 4, "0.25", 10, 100, 1);
    assert!(
        result.status.success(),
        "oversized-record run failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let mut reader = bam::Reader::from_path(output).unwrap();
    let names: Vec<_> = reader
        .records()
        .map(|record| record.unwrap().qname().to_vec())
        .collect();
    assert_eq!(names, [b"r000000", b"r000001", b"r000002"]);
}

#[test]
fn empty_input_produces_a_valid_header_only_bam() {
    let tmp = TestDir::new("empty");
    let input = tmp.path("input.bam");
    let output = tmp.path("output.bam");
    {
        let _writer = bam::Writer::from_path(&input, &header(), bam::Format::Bam).unwrap();
    }
    let status = run_filter(&input, &output, 1, "1", 10);
    assert!(
        status.status.success(),
        "filter failed: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    assert_eq!(final_stderr_line(&status), "Kept 0 reads out of 0 reads");
    let mut reader = bam::Reader::from_path(output).unwrap();
    assert!(String::from_utf8_lossy(reader.header().as_bytes()).contains("SO:coordinate"));
    assert_eq!(reader.records().count(), 0);
}

#[test]
fn invalid_thread_arguments_are_rejected_by_clap() {
    let tmp = TestDir::new("cli");
    let input = tmp.path("input.bam");
    let output = tmp.path("output.bam");
    {
        let _writer = bam::Writer::from_path(&input, &header(), bam::Format::Bam).unwrap();
    }

    for ratio in ["-0.01", "1.01", "NaN", "inf", "-inf"] {
        let result = run_filter(&input, &output, 4, ratio, 10);
        assert!(!result.status.success(), "ratio {ratio} was accepted");
    }

    let result = run_filter(&input, &output, 0, "0.25", 10);
    assert!(!result.status.success(), "-t 0 was accepted");
}

fn bgzf_blocks(bytes: &[u8]) -> Vec<(usize, usize)> {
    let mut blocks = Vec::new();
    let mut offset = 0;
    while offset + 18 <= bytes.len() {
        assert_eq!(&bytes[offset..offset + 4], &[0x1f, 0x8b, 0x08, 0x04]);
        let size = u16::from_le_bytes([bytes[offset + 16], bytes[offset + 17]]) as usize + 1;
        assert!(size >= 28);
        assert!(offset + size <= bytes.len());
        blocks.push((offset, size));
        offset += size;
    }
    assert_eq!(offset, bytes.len());
    blocks
}

fn run_filter_with_timeout(
    input: &Path,
    output: &Path,
    timeout: Duration,
) -> (std::process::ExitStatus, Vec<u8>) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_bamqvfilter"))
        .args(["-q", "10", "-t", "4", "--compression-thread-ratio", "0.25"])
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output)
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let start = Instant::now();
    loop {
        if child.try_wait().unwrap().is_some() {
            let output = child.wait_with_output().unwrap();
            return (output.status, output.stderr);
        }
        if start.elapsed() >= timeout {
            let _ = child.kill();
            let _ = child.wait();
            panic!("bamqvfilter hung after corrupt-input read failure");
        }
        thread::sleep(Duration::from_millis(10));
    }
}

#[test]
fn truncated_data_block_is_an_error_and_does_not_deadlock_writer() {
    let tmp = TestDir::new("corrupt");
    let valid = tmp.path("valid.bam");
    let corrupt = tmp.path("corrupt.bam");
    let output = tmp.path("output.bam");
    {
        let mut writer = bam::Writer::from_path(&valid, &header(), bam::Format::Bam).unwrap();
        for index in 0..20_000 {
            write_record(&mut writer, index, 100, 20);
        }
    }

    let bytes = fs::read(valid).unwrap();
    let blocks = bgzf_blocks(&bytes);
    assert!(blocks.len() > 3, "fixture needs multiple data BGZF blocks");
    // Truncate within a non-EOF data block. Merely deleting the 28-byte BGZF
    // EOF marker is allowed/warn-only in some HTSlib configurations.
    let (block_start, block_size) = blocks[blocks.len() / 2];
    fs::write(&corrupt, &bytes[..block_start + block_size / 2]).unwrap();

    let (status, stderr) = run_filter_with_timeout(&corrupt, &output, Duration::from_secs(10));
    assert!(
        !status.success(),
        "corrupt BAM was silently accepted; stderr={}",
        String::from_utf8_lossy(&stderr)
    );
    assert!(!String::from_utf8_lossy(&stderr).contains("Kept "));
}
