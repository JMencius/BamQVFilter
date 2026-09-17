use clap::Parser;
use crossbeam::channel::{bounded, Receiver, Sender};
use rayon::{prelude::*, ThreadPool};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{self, Read};
use std::mem::size_of;
use std::path::{Path, PathBuf};
use std::thread;

const DEFAULT_BATCH_RECORDS: usize = 4096;
const DEFAULT_BATCH_MIB: usize = 64;
const DEFAULT_COMPRESSION_THREAD_RATIO: f64 = 0.25;
const MIB: usize = 1024 * 1024;
const QV_COMPARISON_EPSILON: f64 = 1e-12;

type AppResult<T> = Result<T, String>;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Filtering BAM file according to read QV.",
    long_about = None
)]
struct Cli {
    /// Sets a minimum Phred average quality score
    #[arg(short = 'q', long = "quality")]
    quality: u8,

    /// Worker-pool budget: QV workers plus extra BAM compression workers
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 4,
        value_parser = parse_positive_usize
    )]
    threads: usize,

    /// Maximum fraction of --threads assigned to extra BAM compression workers
    #[arg(
        long = "compression-thread-ratio",
        value_name = "RATIO",
        default_value_t = DEFAULT_COMPRESSION_THREAD_RATIO,
        value_parser = parse_compression_thread_ratio
    )]
    compression_thread_ratio: f64,

    /// Maximum number of BAM records held in one processing batch
    #[arg(
        long = "batch-records",
        default_value_t = DEFAULT_BATCH_RECORDS,
        value_parser = parse_positive_usize
    )]
    batch_records: usize,

    /// Approximate maximum allocated BAM record memory per batch, in MiB
    #[arg(
        long = "batch-mib",
        default_value_t = DEFAULT_BATCH_MIB,
        value_parser = parse_positive_usize
    )]
    batch_mib: usize,

    /// Keep only mapped primary alignments
    #[arg(long = "primary", visible_alias = "primary-only")]
    primary_only: bool,

    /// Input BAM filename
    #[arg(short = 'i', long = "input")]
    input: String,

    /// Output BAM filename
    #[arg(short = 'o', long = "output")]
    output: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ThreadAllocation {
    qv_threads: usize,
    compression_threads: usize,
}

struct FilteredBatch {
    ordinal: u64,
    records: Vec<Record>,
}

fn main() {
    if let Err(error) = run(Cli::parse()) {
        eprintln!("Error: {error}");
        std::process::exit(1);
    }
}

fn run(args: Cli) -> AppResult<()> {
    ensure_distinct_paths(&args.input, &args.output)?;

    let allocation = allocate_threads(args.threads, args.compression_thread_ratio);
    let batch_bytes = args
        .batch_mib
        .checked_mul(MIB)
        .ok_or_else(|| "--batch-mib is too large".to_string())?;

    let qv_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(allocation.qv_threads)
        .build()
        .map_err(|error| format!("Could not create QV thread pool: {error}"))?;

    let mut reader = bam::Reader::from_path(&args.input)
        .map_err(|error| format!("Could not open input BAM '{}': {error}", args.input))?;
    let header = bam::Header::from_template(reader.header());

    let mut output = bam::Writer::from_path(&args.output, &header, bam::Format::Bam)
        .map_err(|error| format!("Could not open output BAM '{}': {error}", args.output))?;

    if allocation.compression_threads > 0 {
        output
            .set_threads(allocation.compression_threads)
            .map_err(|error| format!("Could not start BAM compression threads: {error}"))?;
    }

    eprintln!(
        "Worker allocation: {} QV thread(s), {} extra BAM compression thread(s)",
        allocation.qv_threads, allocation.compression_threads
    );

    // A zero-capacity channel creates a true two-batch pipeline: while the
    // writer owns batch N, the producer can read and filter batch N+1. If the
    // producer gets farther ahead, send() blocks and applies backpressure.
    let (sender, receiver) = bounded::<FilteredBatch>(0);
    let writer_handle = thread::Builder::new()
        .name("bam-writer".to_string())
        .spawn(move || write_batches(output, receiver))
        .map_err(|error| format!("Could not start BAM writer thread: {error}"))?;

    let producer_result = produce_batches(
        &mut reader,
        &qv_pool,
        &sender,
        args.quality,
        args.primary_only,
        args.batch_records,
        batch_bytes,
    );

    // The writer exits only after every Sender has been dropped. This must
    // happen before join(), including when reading or filtering failed.
    drop(sender);

    let writer_result = writer_handle
        .join()
        .map_err(|_| "BAM writer thread panicked".to_string())
        .and_then(|result| result);

    // If a write failed, prefer its specific error over the producer's generic
    // "writer stopped" channel error.
    match (producer_result, writer_result) {
        (_, Err(writer_error)) => Err(writer_error),
        (Err(producer_error), Ok(_)) => Err(producer_error),
        (Ok(total_reads), Ok(kept_reads)) => {
            if kept_reads > total_reads {
                return Err(format!(
                    "Internal counting error: wrote {kept_reads} records after reading {total_reads}"
                ));
            }

            eprintln!("Kept {kept_reads} reads out of {total_reads} reads");
            Ok(())
        }
    }
}

fn produce_batches(
    reader: &mut bam::Reader,
    qv_pool: &ThreadPool,
    sender: &Sender<FilteredBatch>,
    min_qv: u8,
    primary_only: bool,
    max_batch_records: usize,
    max_batch_bytes: usize,
) -> AppResult<u64> {
    let batch_capacity = initial_batch_capacity(max_batch_records, max_batch_bytes);
    let mut batch = Vec::with_capacity(batch_capacity);
    let mut batch_bytes = 0usize;
    let mut ordinal = 0u64;
    let mut total_reads = 0u64;

    for result in reader.records() {
        let record = result.map_err(|error| format!("Could not read input BAM record: {error}"))?;
        total_reads = total_reads
            .checked_add(1)
            .ok_or_else(|| "Input read count exceeded u64::MAX".to_string())?;
        let record_bytes = record_allocated_bytes(&record);

        // If this record would push a non-empty batch over the byte target,
        // dispatch the current batch first. A single oversized record is still
        // processed as a one-record batch.
        if !batch.is_empty() && batch_bytes.saturating_add(record_bytes) > max_batch_bytes {
            let completed = std::mem::replace(&mut batch, Vec::with_capacity(batch_capacity));
            dispatch_batch(completed, ordinal, qv_pool, sender, min_qv, primary_only)?;
            ordinal = next_ordinal(ordinal)?;
            batch_bytes = 0;
        }

        batch_bytes = batch_bytes.saturating_add(record_bytes);
        batch.push(record);

        if batch.len() >= max_batch_records || batch_bytes >= max_batch_bytes {
            let completed = std::mem::replace(&mut batch, Vec::with_capacity(batch_capacity));
            dispatch_batch(completed, ordinal, qv_pool, sender, min_qv, primary_only)?;
            ordinal = next_ordinal(ordinal)?;
            batch_bytes = 0;
        }
    }

    if !batch.is_empty() {
        dispatch_batch(batch, ordinal, qv_pool, sender, min_qv, primary_only)?;
    }

    Ok(total_reads)
}

fn dispatch_batch(
    batch: Vec<Record>,
    ordinal: u64,
    qv_pool: &ThreadPool,
    sender: &Sender<FilteredBatch>,
    min_qv: u8,
    primary_only: bool,
) -> AppResult<()> {
    let records = filter_batch(batch, qv_pool, min_qv, primary_only);
    sender
        .send(FilteredBatch { ordinal, records })
        .map_err(|_| "BAM writer stopped before receiving all batches".to_string())
}

fn filter_batch(
    mut batch: Vec<Record>,
    qv_pool: &ThreadPool,
    min_qv: u8,
    primary_only: bool,
) -> Vec<Record> {
    // par_iter() over a Vec is indexed, so the keep flags correspond exactly
    // to the input positions even when workers finish out of order.
    let keep: Vec<bool> = qv_pool.install(|| {
        batch
            .par_iter()
            .map(|record| {
                (!primary_only || is_mapped_primary(record)) && filter_by_quality(record, min_qv)
            })
            .collect()
    });

    // rust-htslib 0.46 reader Records contain an Rc<HeaderView>, despite that
    // version declaring Record as Send. Moving those original Records to the
    // writer while the reader keeps cloning the Rc can corrupt its non-atomic
    // reference count. Record::clone() performs a deep bam_copy1 into a fresh
    // Record with no attached header, so only these detached copies cross the
    // thread boundary. All reader-owned Records are dropped on this thread.
    clone_selected(&mut batch, &keep)
}

fn clone_selected<T: Clone>(items: &mut Vec<T>, keep: &[bool]) -> Vec<T> {
    assert_eq!(items.len(), keep.len(), "keep mask length mismatch");
    let kept_count = keep.iter().filter(|&&retain| retain).count();
    let mut selected = Vec::with_capacity(kept_count);

    for (item, retain) in items.drain(..).zip(keep.iter().copied()) {
        if retain {
            selected.push(item.clone());
        }
    }

    selected
}

fn write_batches(mut output: bam::Writer, receiver: Receiver<FilteredBatch>) -> AppResult<u64> {
    let mut expected_ordinal = 0u64;
    let mut kept_reads = 0u64;

    while let Ok(batch) = receiver.recv() {
        if batch.ordinal != expected_ordinal {
            return Err(format!(
                "Internal ordering error: expected batch {}, received batch {}",
                expected_ordinal, batch.ordinal
            ));
        }

        for record in batch.records {
            let next_kept_reads = kept_reads
                .checked_add(1)
                .ok_or_else(|| "Output read count exceeded u64::MAX".to_string())?;
            output
                .write(&record)
                .map_err(|error| format!("Could not write output BAM record: {error}"))?;
            kept_reads = next_kept_reads;
        }

        expected_ordinal = next_ordinal(expected_ordinal)?;
    }

    Ok(kept_reads)
}

fn record_allocated_bytes(record: &Record) -> usize {
    (record.inner().m_data as usize).saturating_add(size_of::<Record>())
}

fn initial_batch_capacity(max_batch_records: usize, max_batch_bytes: usize) -> usize {
    // Do not let a very large --batch-records value bypass --batch-mib merely
    // through Vec's up-front Record-slot allocation.
    max_batch_records
        .min(max_batch_bytes / size_of::<Record>())
        .max(1)
}

fn next_ordinal(ordinal: u64) -> AppResult<u64> {
    ordinal
        .checked_add(1)
        .ok_or_else(|| "Too many BAM batches to represent".to_string())
}

fn allocate_threads(total_threads: usize, compression_ratio: f64) -> ThreadAllocation {
    debug_assert!(total_threads > 0);
    debug_assert!(compression_ratio.is_finite());
    debug_assert!((0.0..=1.0).contains(&compression_ratio));

    let requested_compression_threads =
        ((total_threads as f64) * compression_ratio).floor() as usize;
    let compression_threads = requested_compression_threads.min(total_threads.saturating_sub(1));

    ThreadAllocation {
        qv_threads: total_threads - compression_threads,
        compression_threads,
    }
}

fn parse_positive_usize(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("'{value}' is not a valid positive integer"))?;

    if parsed == 0 {
        Err("value must be greater than zero".to_string())
    } else {
        Ok(parsed)
    }
}

fn parse_compression_thread_ratio(value: &str) -> Result<f64, String> {
    let ratio = value
        .parse::<f64>()
        .map_err(|_| format!("'{value}' is not a valid number"))?;

    if ratio.is_finite() && (0.0..=1.0).contains(&ratio) {
        Ok(ratio)
    } else {
        Err("compression thread ratio must be a finite number from 0.0 to 1.0".to_string())
    }
}

fn ensure_distinct_paths(input: &str, output: &str) -> AppResult<()> {
    let input_path = std::fs::canonicalize(input)
        .map_err(|error| format!("Could not resolve input BAM '{input}': {error}"))?;
    let output_path = canonical_output_path(output)?;

    if input_path == output_path {
        Err("Input and output BAM paths must be different".to_string())
    } else {
        Ok(())
    }
}

fn canonical_output_path(output: &str) -> AppResult<PathBuf> {
    let output_path = Path::new(output);

    if output_path.exists() {
        return std::fs::canonicalize(output_path)
            .map_err(|error| format!("Could not resolve output BAM '{output}': {error}"));
    }

    let file_name = output_path
        .file_name()
        .ok_or_else(|| format!("Output BAM path '{output}' has no filename"))?;
    let parent = output_path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let canonical_parent = std::fs::canonicalize(parent).map_err(|error| {
        format!(
            "Could not resolve output directory '{}': {error}",
            parent.display()
        )
    })?;

    Ok(canonical_parent.join(file_name))
}

fn filter_by_quality(record: &Record, min_qv: u8) -> bool {
    let average_qv = average_quality(record.qual());
    quality_passes(average_qv, min_qv)
}

fn is_mapped_primary(record: &Record) -> bool {
    !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary()
}

fn quality_passes(average_qv: f64, min_qv: u8) -> bool {
    average_qv.is_finite()
        && average_qv + QV_COMPARISON_EPSILON >= min_qv as f64
        && (-QV_COMPARISON_EPSILON..=90.0 + QV_COMPARISON_EPSILON).contains(&average_qv)
}

fn average_quality(quals: &[u8]) -> f64 {
    let probability_sum = quals
        .iter()
        .map(|&q| {
            let q = q as f64;
            10_f64.powf(q / -10.0)
        })
        .sum::<f64>();
    (probability_sum / quals.len() as f64).log10() * -10.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allocates_threads_without_exceeding_ratio() {
        assert_eq!(
            allocate_threads(24, 0.25),
            ThreadAllocation {
                qv_threads: 18,
                compression_threads: 6,
            }
        );
        assert_eq!(
            allocate_threads(4, 0.25),
            ThreadAllocation {
                qv_threads: 3,
                compression_threads: 1,
            }
        );
        assert_eq!(
            allocate_threads(2, 0.25),
            ThreadAllocation {
                qv_threads: 2,
                compression_threads: 0,
            }
        );
        assert_eq!(
            allocate_threads(8, 1.0),
            ThreadAllocation {
                qv_threads: 1,
                compression_threads: 7,
            }
        );
        assert_eq!(
            allocate_threads(1, 1.0),
            ThreadAllocation {
                qv_threads: 1,
                compression_threads: 0,
            }
        );
    }

    #[test]
    fn batch_preallocation_respects_the_byte_target() {
        assert_eq!(initial_batch_capacity(37, MIB), 37);

        let capacity = initial_batch_capacity(usize::MAX, MIB);
        assert!(capacity <= MIB / size_of::<Record>());
        assert!(capacity > 0);
    }

    #[test]
    fn rejects_invalid_cli_values() {
        assert!(parse_positive_usize("0").is_err());
        assert!(parse_positive_usize("4").is_ok());

        for value in ["-0.1", "1.1", "NaN", "inf", "-inf", "not-a-number"] {
            assert!(parse_compression_thread_ratio(value).is_err());
        }
        for value in ["0", "0.25", "1"] {
            assert!(parse_compression_thread_ratio(value).is_ok());
        }
    }

    #[test]
    fn cloned_selection_preserves_relative_order() {
        let mut items = vec!["A", "B", "C", "D", "E"];
        let selected = clone_selected(&mut items, &[true, false, true, false, true]);
        assert_eq!(selected, vec!["A", "C", "E"]);
        assert!(items.is_empty());
    }

    #[test]
    fn average_quality_uses_mean_error_probability() {
        for q in [0u8, 10, 20, 40, 90] {
            let quals = vec![q; 100];
            assert!((average_quality(&quals) - q as f64).abs() < 1e-9);
        }

        let mixed = average_quality(&[0, 40]);
        assert!((mixed - 3.009_865_683_871_185).abs() < 1e-12);
    }

    #[test]
    fn empty_quality_is_not_accepted() {
        assert!(average_quality(&[]).is_nan());
        assert!(!quality_passes(average_quality(&[]), 0));
    }

    #[test]
    fn exact_threshold_is_not_lost_to_rounding() {
        let average = average_quality(&[10; 31]);
        assert!(average < 10.0);
        assert!(quality_passes(average, 10));
    }
}
