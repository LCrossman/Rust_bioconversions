# Rust_bioconversions

Build with Cargo

Gbk2fna.rs = Rust tool to convert genbank flat file to fasta DNA sequence file<br>

Embl2fna.rs = Rust tool to convert EMBL flat file to fasta DNA sequence file<br>
Both scripts need the input file to have proper genbank or embl headers.


Approximately 10 x faster than biopython tool for the same task!


Usage: prints the output to stdout

Required crates:
use std::io;
use std::fs;
use std::env;
use itertools::Itertools;
use std::vec::Vec;
use std::convert::AsRef;
use std::path::Path;
use std::process;
use anyhow::Context;

for gbk2ffn.rs also require:
use std::collections::BTreeMap;

for gbk2faa.rs you require:
use std::collections::BTreeMap;
use protein_translate::translate;
