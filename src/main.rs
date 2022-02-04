extern crate seq_io;
extern crate noodles;
use std::str;
use regex;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;
use std::io::{self, Write};

use seq_io::fasta::{Reader,Record};
use lazy_static::lazy_static;

use clap::{Parser};

// VCF
use noodles::vcf as vcf;

#[derive(Parser,Debug)]
#[clap(name = "basefinder")]
#[clap(version = "0.0.1")]
#[clap(about = "Finds corresponding base in closest species", long_about = None)]

struct Cli {
    /// Sequence alignment in fasta format
    #[clap(short,long)]
    alignment: String,
    /// #CHROM POS ALT REF HGVS.c STRAND
    #[clap(short, long)]
    tsv: String,
    /// tab-separated distance matrix (id in leftmost column only, not on top)
    /// Ex: rapidnj <fasta.fa> -o m | tail -n +2 | tr ' ' '\t'
    #[clap(short, long)]
    distance_matrix: String,
    /// hide disregarded
    #[clap(short, long)]
    no_disregarded: bool,
     /// Clean output (hide some columns)
    #[clap(short, long)]
    clean_output: bool,
    /// Hide non-cds
    #[clap(short, long)]
    remove_utr: bool,
}
struct distance_struct {
    name: String,
    distance: f32
}

fn main() {
    // Parse command line arguments
    let cli = Cli::parse(); 
    // Read files
    let mut fasta_reader = Reader::from_path(&cli.alignment).unwrap();
    let mut vcf_reader = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(File::open(&cli.tsv).unwrap());
    let mut matrix_reader = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(File::open(&cli.distance_matrix).unwrap());

    // Collect the results
    let records: Result<Vec<_>, _> = fasta_reader.records().collect();
    let csv_records: Result<Vec<csv::StringRecord>, _> = vcf_reader.records().collect();
    let matrix: Result<Vec<csv::StringRecord>, _> = matrix_reader.records().collect();
    
    // Initialise some variables
    let mut i = 0;
    if cli.clean_output {

    writeln!(io::stdout(), "{0: <10}\t{1: <5}\t{2: <1}\t{3: <5}\t{4: <20}\t{5: <5}\t{6: <5}\t{7: <5}\t{8: <20}\t{9: <5}", "VCF_CHR", "VCF_POS", "VCF_REF", "VCF_ALT", "IN_SPC", "IN_NUC", "REF_COR", "OUT_NUC", "OUT_SPC", "OUT_RNK");
    } else {
    // Print header
    writeln!(io::stdout(), "{0: <10}\t{1: <5}\t{2: <1}\t{3: <5}\t{4: <10}\t{5: <5}\t{6: <20}\t{7: <5}\t{8: <5}\t{9: <5}\t{10: <5}\t{11: <5}\t{12: <5}\t{13: <5}\t{14: <5}\t{15: <5}\t{16: <5}\t{17: <20}\t{18: <5}", "VCF_CHR", "VCF_POS", "VCF_REF", "VCF_ALT", "VCF_ANN", "GFF_STR", "IN_SPC", "IN_NUC", "REF!=IN", "REF_COR", "ALT_COR", "ALI_POS", "IN_POS", "GAP", "OUT_NUC", "DSRGD", "NOT_CDS", "OUT_SPC", "OUT_RNK");
    }
    
    // Get which species are present in distance matrix 
    let n_matrix = matrix.as_ref().unwrap().iter().len();
    let mut species_record: Vec<&str> = Vec::with_capacity(n_matrix);
                lazy_static! { 
                    static ref RE: regex::Regex = regex::Regex::new(r"\d+").unwrap();
                }
    
    for ma in matrix.as_ref().unwrap() {
        species_record.push(&ma[0].as_ref());
    }

    // For every FASTA record
    for record in records.as_ref().unwrap().iter() {
        // Save the name 
    let mut outgroup_rank = 0;
    let mut outgroup_id = 0;
        let mut current_record = str::from_utf8(record.head()).unwrap();
        
        // Go through every matrix row
        // Since we onl
        for line in matrix.as_ref().unwrap() {
            
            let mut test: Vec<&String> = Vec::new();
            let mut x = Vec::new();
            
            // Match the FASTA sequece contains the distance matrix species
            
            if current_record.contains(&line[0]) {
                let mut i = 0;
            
                // Then go through all columns
                for qa in matrix.as_ref().unwrap() {
                    
                    // Skip the first and don't go further than the last
                    if i <= n_matrix {
                       if i > 1 {
                           // Get the value at current column 
                           let number = &line[i];
                           // Put it into a vector 
                           let some = distance_struct {name: (&species_record[i-1]).to_string(), distance: number.parse::<f32>().unwrap()};
                           x.push(some);
                       }
                    }
                    i += 1;
                }
                // Sort the vector by distance
                x.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap());
                let mut i = 0;
                for element in &x {
                    if i > 0 {
                            test.push(&element.name);
                    }
                    i += 1; 
                } 
            } 

            // When distance record matches fasta record, break the loop 
            // and get the id of the fasta record
            '_outer: for ex in test {
                let mut i = 0;
                    '_inner: for record in records.as_ref().unwrap().iter() {
                        let head = str::from_utf8(record.head()).unwrap();
                        if head.contains(ex) {
                            outgroup_id = i;
                            break '_outer;
                } 
                i += 1;
            }    
            outgroup_rank += 1;
        } 
    }
    

    let outgroup_record = str::from_utf8(records.as_ref().unwrap()[outgroup_id].head()).unwrap();      // So now we go through every record in the VCF file

    for recordx in csv_records.as_ref().unwrap() {
        
        let chrom = &recordx[0];
        
        // Check first so no need to go through everything

        if current_record.contains(chrom) {
            // Get fields 
            let pos = &recordx[1];
            let reference = &recordx[2];
            let alt = &recordx[3];
            let hgvs = &recordx[4];
            let actual_strand = &recordx[5];

            let parsed_position = pos.trim().parse::<i32>().expect("Numeric!");
            // And if the current FASTA record contains the VCF Record 
            // Not sure which way to go here, might be too loose naming
            // The VCF has to be more spefic than the distance matrix
            // And it should be == to the fasta.
            // However if fasta has eval|eval_4808_1 etc... 
            // Best if they are equal and _not_ using contains
            // Then compare the sequences
            
            let current_sequence = record.seq();
            let outgroup_sequence = records.as_ref().unwrap()[outgroup_id].seq();
            

            let mut gap = 0;
            let mut nucleotide = 0;
            let mut position = 0;
            
            // Go through every base in alignment...can this go on the outside of csv-loop?
            for (base_c, base_o) in current_sequence.iter().zip(outgroup_sequence.iter()) {
                if base_c == &b'-' {
                    gap += 1;
                } else {
                    nucleotide += 1;
                }

                position += 1;
                
                
                let mut disregard = 0; 
                let mut not_cds = 0; 
                
                let bref = reference.as_bytes()[0];
                let balt = alt.as_bytes()[0];

                let mut strand_corrected = &bref;
                let mut strand_correctedalt = &balt;

                let ss = &hgvs.to_string();
                    
                if ss.contains("-") || ss.contains("*") {
                    not_cds += 1;
                }
                // only finds first digit
                //let corr_pos = ss.chars().find(|a| a.is_digit(10)).and_then(|a| a.to_digit(10)).unwrap();
                
                let cap = RE.captures(ss).unwrap();
                let mut base_vs_reference = 0;

                    let corr_pos = cap.get(0).map_or("", |m| m.as_str()).parse::<i32>().unwrap();
                //if nucleotide == parsed_position {
                // Not sure if the && ... part is ok..
                //if nucleotide == corr_pos && *base_c != b'-' {
                if nucleotide == corr_pos {

                    if *base_c != reference.as_bytes()[0] {
                        base_vs_reference = 1;
                    
                        if actual_strand.as_bytes()[0] == b'-' {
                            match bref {
                                b'A' => strand_corrected = &&b'T',
                                b'C' => strand_corrected = &&b'G',
                                b'G' => strand_corrected = &&b'C',
                                b'T' => strand_corrected = &&b'A',
                                _ => panic!("panik!"),
                            } 
                            match balt {
                                b'A' => strand_correctedalt = &&b'T',
                                b'C' => strand_correctedalt = &&b'G',
                                b'G' => strand_correctedalt = &&b'C',
                                b'T' => strand_correctedalt = &&b'A',
                                _ => panic!("panik!"),
                            }
                        }
                    
                      } 
                    
                    if strand_corrected != base_o && strand_correctedalt != base_o {
                        disregard += 1;
                    }
                                    
                    if cli.no_disregarded {
                        if disregard > 0 {
                            continue;
                        }
                    }

                    if cli.remove_utr {
                    if not_cds > 0 {
                            continue;
                        }
                    }
                    if cli.clean_output {

                        writeln!(io::stdout(), "{0: <10}\t{1: <5}\t{2: <1}\t{3: <5}\t{4: <20}\t{5: <5}\t{6: <5}\t{7: <5}\t{8: <20}\t{9: <5}", chrom, pos, reference, alt, current_record, *base_c as char, *strand_corrected, *base_o as char, outgroup_record, outgroup_rank);
                    } else {
                    
                        writeln!(io::stdout(), "{0: <10}\t{1: <5}\t{2: <5}\t{3: <5}\t{4: <10}\t{5: <5}\t{6: <20}\t{7: <5}\t{8: <5}\t{9: <5}\t{10: <5}\t{11: <5}\t{12: <5}\t{13: <5}\t{14: <5}\t{15: <5}\t{16: <5}\t{17: <20}\t{18: <5}", chrom, pos, reference, alt, hgvs.to_string(), actual_strand, current_record, *base_c as char, base_vs_reference,*strand_corrected as char, *strand_correctedalt as char, position, nucleotide, gap, *base_o as char, disregard, not_cds, outgroup_record, outgroup_rank);
                 //   if reference.as_bytes()[0] != *base_c {
                 //       panic!("Reference: {} does not match ingroup base: {}!", reference.as_bytes()[0], base_c);
                 //   }
                    }
                    // could this be here?
                    break;
                    }
            } 

        }
    }
}
}
