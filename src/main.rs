use clap::Parser;
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::Record;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::thread;
use crossbeam::channel::{Sender, unbounded};


#[derive(Parser, Debug)]
#[clap(author, version, about="Filtering BAM file according to read QV.", long_about = None)]
#[structopt(name = "bamqvfilter", about = "Filters BAM files based on read quality values.")]
struct Cli {
    /// Sets a minimum Phred average quality score
    #[arg(short = 'q', long = "quality", value_parser)]
    quality: u8,

    /// Use N parallel threads
    #[arg(short = 't', long = "threads", value_parser, default_value_t = 4)]
    threads: usize,


    /// Input filename
    #[arg(short = 'i', long = "input", value_parser)]
    input: Option<String>,

    /// Output filename
    #[arg(short = 'o', long = "output", value_parser)]
    output: Option<String>,
}



fn main() {
    let args = Cli::parse();
    let input_file = args.input.unwrap_or_else(|| {
        eprintln!("No input file provided");
        std::process::exit(1);
    });

    let output_file = args.output.unwrap_or_else(|| {
        eprintln!("No output file provided");
        std::process::exit(1);
    });

    let min_qv: u8 = args.quality;
    let num_threads: usize = args.threads;

    // Set the number of threads for Rayon
    rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().unwrap();

    // Open the input BAM file
    let mut reader = bam::Reader::from_path(input_file).expect("Error opening BAM file");
    
    // Get the header from the input BAM file
    let header = bam::Header::from_template(reader.header());

    // Create output BAM writer
    // let output_file_clone = output_file.clone();
    // let mut output = bam::Writer::from_path(output_file, &header, bam::Format::Bam).expect("Error opening output BAM file");

    
    // Create a channel to send filtered records
    // let (sender, receiver): (Sender<Record>, Receiver<Record>) = channel();
    //

    let (sender, receiver) = unbounded::<Record>();

    // Spawn a thread to write filtered records to the output BAM file
    let output_thread = thread::spawn(move || {
        let mut output = bam::Writer::from_path(output_file, &header, bam::Format::Bam).expect("Error opening output BAM file");
        while let Ok(record) = receiver.recv() {
            output.write(&record).expect("Error writing to output BAM file");
        }
    });

    // Read and filter BAM records in parallel
    reader.records()
        .par_bridge()
        .filter_map(Result::ok)
        .filter(|record| filter_by_quality(record, min_qv))
        .for_each(|record| {
            sender.send(record).expect("Error sending record to output thread");
        
        });

    // Close the channel and wait for the output thread to finish
    drop(sender);
    output_thread.join().expect("Error joining output thread");
}

fn filter_by_quality(record: &Record, min_qv: u8) -> bool {
    // println!("{:?}", record.qual());

    let average_qv = average_quality(record.qual());
    // println!("{}", average_qv);
    average_qv >= min_qv as f64 && average_qv >= 0.0 && average_qv <= 90.0
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

