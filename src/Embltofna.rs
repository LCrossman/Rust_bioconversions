use std::io;
use std::fs;
use std::env;
use itertools::Itertools;
use std::vec::Vec;
use std::convert::AsRef;
use std::path::Path;
use std::process;
use anyhow::Context;


//const MAX_GBK_BUFFER_SIZE: usize = 512;
/// A Gbk reader.


pub trait GbkRead {
    fn read(&mut self, record: &mut Record) -> io::Result<()>;
}

#[derive(Debug)]
pub struct Reader<B> {
    reader: B,
    line_buffer: String,
}

impl Reader<io::BufReader<fs::File>> {
    /// Read Gbk from given file path in given format.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .with_context(|| format!("Failed to read Gbk from {:#?}", path))
    }
}

impl<R> Reader<io::BufReader<R>>
where
     R: io::Read,
{
    //// Create a new Gbk reader given an instance of `io::Read` in given format
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
	    line_buffer: String::new(),
        }
    }
}

impl<B> Reader<B>
where
    B: io::BufRead,
{
    pub fn from_bufread(bufreader: B) -> Self {
        Reader {
	     reader: bufreader,
	     line_buffer: String::new(),
	}
    }
    pub fn records(self) -> Records<B> {
        Records {
	   reader: self,
	   error_has_occurred: false,
	}
    }
}

impl<B> GbkRead for Reader<B>
where
    B: io::BufRead,
{
    fn read(& mut self, record: &mut Record) -> io::Result<()> {
        record.rec_clear();
	let mut sequences = String::new();
    	if self.line_buffer.is_empty() {
	    self.reader.read_line(&mut self.line_buffer)?;
	    if self.line_buffer.is_empty() {
	        return Ok(());
	        }
            }
	while !self.line_buffer.is_empty() {
	    if self.line_buffer.starts_with("ID") {
			record.rec_clear();
	            	let mut header_fields: Vec<&str> = self.line_buffer.split(";").collect();
	                let mut header_iter = header_fields.iter();
			let mut id = header_iter.next().map(|s| s.to_string()).unwrap();
			let mut rid: Vec<&str> = id.split_whitespace().collect();
			record.id = rid[1].trim().to_string();
		//	record.id = header_iter.next().map(|s| s.to_string()).unwrap();
			for _i in 0..5 {
	                    header_iter.next();
			}
	                let lens = header_iter.next().map(|s| s.to_string()).unwrap();
			let ll: Vec<&str> = lens.split_whitespace().collect();
	                record.length = ll[0].trim().parse::<u32>().unwrap();
			self.line_buffer.clear();
			}
	    if self.line_buffer.starts_with("SQ") {
	        let mut sequences = String::new();
	        let result_seq = loop {  
		     self.line_buffer.clear();
		     self.reader.read_line(&mut self.line_buffer)?;
                     if self.line_buffer.starts_with("//") {
		         break sequences;
                     } else {
	                 let mut s: Vec<&str> = self.line_buffer.split_whitespace().collect();
		         s.pop();
		         let sequence = s.iter().join("");
		         sequences.push_str(&sequence);
                         }     
	             };
		record.sequence = result_seq.to_string();
		break;
                }
	 self.line_buffer.clear();
	 self.reader.read_line(&mut self.line_buffer)?;
	 }
         Ok(())
     }
}



#[derive(Default, Clone, Debug)]
pub struct Record {
    id: String,
    length: u32,
    sequence: String
}


impl Record {
    /// Create a new instance.
    pub fn new() -> Self {
        Record {
            id: "".to_owned(),
            length: 0,
            sequence: "".to_owned()
        }
    }
    pub fn is_empty(&mut self) -> bool {
        self.id.is_empty() && self.length == 0
    }
    pub fn check(&mut self) -> Result<(), &str> {
        if self.id().is_empty() {
            return Err("Expecting id for Gbk record.");
        }
        Ok(())
    }
    pub fn id(&mut self) -> &str {
        self.id.as_ref()
    }
    pub fn length(&mut self) -> u32 {
        self.length
    }   
    pub fn sequence(&mut self) -> &str {
        self.sequence.as_ref()
    }
    fn rec_clear(&mut self) {
        self.id.clear();
        self.length = 0;
        self.sequence.clear();
    }
}

/// An iterator over the records of a Fasta file.
pub struct Records<B>
where
    B: io::BufRead,
{
    reader: Reader<B>,
    error_has_occurred: bool,
}


impl<B> Iterator for Records<B>
where
    B: io::BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
            if self.error_has_occurred {
	        None
	        } else {
                let mut record = Record::new();
                match self.reader.read(&mut record) {
	            Ok(()) if record.is_empty() => None,
                    Ok(()) => Some(Ok(record)),
		    Err(err) => {
		        self.error_has_occurred = true;
		        Some(Err(err))
		        }
                    }
                }
     }
}

struct Config {
    filename: String,
}

impl Config {
    fn new(args: &[String]) -> Result<Config, &str> {
    if args.len() < 2 {
        panic!("not enough arguments, please provide filename");
    }
    let filename = args[1].clone();

    Ok(Config { filename })
    }
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let config = Config::new(&args).unwrap_or_else(|err| {
        println!("Problem with parsing file arguments: {}", err);
	process::exit(1);
	});
    let file_gbk = fs::File::open(config.filename)?;
    let mut reader = Reader::new(file_gbk);
    for result in reader.records() {
        let record = result.expect("err");
	println!(">{}_{:?}\n{}", record.id, record.length, record.sequence);
    } 
    Ok(())
}