//Count fastq reads in a gzipped file
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use std::io::{self};
use flate2::read::MultiGzDecoder;

fn main() {
    let gzreader = MultiGzDecoder::new(io::stdin());
    let mut reader = fastq::Reader::new(gzreader);
    let mut record = fastq::Record::new();
    let mut record_count = 0;
    loop {
        reader.read(&mut record).expect("err");
	record_count+=1;
	if record.is_empty() {
	    break;
	    }
    };
    println!("this is the final count {}", record_count);
}