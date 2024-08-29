use std::io;
use std::fs;
use std::env;
use regex::Regex;
use itertools::Itertools;
use std::vec::Vec;
use std::str;
use std::convert::AsRef;
use protein_translate::translate;
use std::path::Path;
use std::process;
use bio::alphabets::dna::revcomp;
use std::collections::BTreeMap;
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
	let mut cds = BTreeMap::new();
    	if self.line_buffer.is_empty() {
	    self.reader.read_line(&mut self.line_buffer)?;
	    if self.line_buffer.is_empty() {
	        return Ok(());
	        }
            }
	'outer: while !self.line_buffer.is_empty() {
	    if self.line_buffer.starts_with("LOCUS") {
			record.rec_clear();
	            	let mut header_fields: Vec<&str> = self.line_buffer.split_whitespace().collect();
	                let mut header_iter = header_fields.iter();
	                header_iter.next();
	                record.id = header_iter.next().map(|s| s.to_string()).unwrap();
	                let lens = header_iter.next().map(|s| s.to_string()).unwrap();
	                record.length = lens.trim().parse::<u32>().unwrap();
			self.line_buffer.clear();
			}
	    if self.line_buffer.starts_with("     CDS") {
	        let re = Regex::new(r"([0-9]+)[[:punct:]]+([0-9]+)").unwrap();
		let location = re.captures(&self.line_buffer).unwrap();
		let start = &location[1];
		let end = &location[2];
		let strand: i32 = if self.line_buffer.contains("complement") {-1} else {1};
		let thestart = start.trim().parse::<u32>().unwrap();
		let thestart = thestart - 1;
		let theend = end.trim().parse::<u32>().unwrap();
                let mut locus_tag = String::new();
                let mut codon_start: u8 = 1;
		loop {
		        self.line_buffer.clear();
			self.reader.read_line(&mut self.line_buffer)?;
                        if self.line_buffer.contains("/locus_tag") {
                            let loctag: Vec<&str> = self.line_buffer.split('\"').collect();
                            locus_tag = loctag[1].to_string();
                            }
                        if self.line_buffer.contains("/codon_start") {
                            let codstart: Vec<&str> = self.line_buffer.split('=').collect();
                            let valstart = codstart[1].trim().parse::<u8>().unwrap();
                            codon_start = valstart;
                            }
			if self.line_buffer.starts_with("     CDS") {
			    cds.insert(locus_tag, (thestart, theend, strand, codon_start));
                            continue 'outer;
			    }
			if self.line_buffer.starts_with("ORIGIN") {
			    continue 'outer;
			}
		   }
                }
	    if self.line_buffer.starts_with("ORIGIN") {
	        let mut sequences = String::new();
	        let result_seq = loop {  
		     self.line_buffer.clear();
		     self.reader.read_line(&mut self.line_buffer)?;
                     if self.line_buffer.starts_with("//") {
		         break sequences;
                     } else {
	                 let s: Vec<&str> = self.line_buffer.split_whitespace().collect();
		         let s = &s[1..];
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
	for (key,val) in cds.iter() {
              let (a,b,c,d) = val;    
	      let sta = *a as usize;
	      let sto = *b as usize;
	      let cod = *d as usize - 1;
	      let mut sliced_sequence: &str = "";
	      if *c == -1 {
	          if cod > 1 {
		     sliced_sequence = &record.sequence[sta+cod..sto];
		     }
		  else {
	             sliced_sequence = &record.sequence[sta..sto];
		     }
	          let cds_char = sliced_sequence.as_bytes();
		  let prot_seq =  translate(&revcomp(cds_char));
		  let parts: Vec<&str> = prot_seq.split('*').collect();
	          println!(">{}\n{}", key,parts[0]);
	      } else {
	          if cod > 1 {
		      sliced_sequence = &record.sequence[sta+cod..sto];
		      }
		  else {
		      sliced_sequence = &record.sequence[sta..sto];
		      }
		  let cds_char = sliced_sequence.as_bytes();
		  let prot_seq = translate(cds_char);
		  let parts: Vec<&str> = prot_seq.split('*').collect();
                  println!(">{}\n{}", key, parts[0]);
                  }
	      }
        Ok(())
     }
}


#[derive(Default, Clone, Debug)]
pub struct Record {
    id: String,
    length: u32,
    sequence: String,
    start: u32,
    end: u32,
    strand: i32,
    locus_tag: String,
}


impl Record {
    /// Create a new instance.
    pub fn new() -> Self {
        Record {
            id: "".to_owned(),
            length: 0,
            sequence: "".to_owned(),
	    start: 0,
	    end: 0,
	    strand: 0,
	    locus_tag: "".to_owned(),
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
        &self.id
    }
    pub fn length(&mut self) -> u32 {
        self.length
    }   
    pub fn sequence(&mut self) -> &str {
        &self.sequence
    }
    pub fn start(&mut self) -> u32 {
        self.start
    }
    pub fn end(&mut self) -> u32 {
        self.end
    }
    pub fn strand(&mut self) -> i32 {
        self.strand
    }
    pub fn locus_tag(&mut self) -> &str {
        &self.locus_tag
    }
    fn rec_clear(&mut self) {
        self.id.clear();
	self.length = 0;
	self.sequence.clear();
	self.start = 0;
	self.end = 0;
	self.strand = 0;
	self.locus_tag.clear();
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
    } 
    Ok(())
}
