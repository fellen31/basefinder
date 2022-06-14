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
    /// #CHROM POS ALT REF HGVS.c AF STRAND
    #[clap(short, long)]
    tsv: String,
    /// tab-separated distance matrix (id in leftmost column only, not on top)
    /// Ex: rapidnj <fasta.fa> -o m | tail -n +2 | tr ' ' '\t'
    #[clap(short, long)]
    distance_matrix: String,
    /// hide non-cds rows
    #[clap(short, long)]
    hide_utr: bool,
     /// Hide header
    #[clap(long)]
    no_header: bool,  
    /// Verbose
    #[clap(short,long)]
    verbose: bool,

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
    
    //Create writer
    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = io::BufWriter::new(handle);
    // Write header
    if !cli.no_header {
        if cli.verbose {
    writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5: <4}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <5}\t{13: <5}\t{14: <5}\t{15: <3}\t{16: <3}\t{17: <3}\t{18: <10}\t{19: <3}\t{20: <5}\t{21: <5}", "#CHROM", "POS", "REF", "ALT", "HGVS", "AF", "STR", "IN", "IN_N", "!=", "C_REF", "C_ALT", "ALPOS", "INPOS", "GAP", "O_NUC", "!KEEP", "!CDS", "OUT_SPC", "RNK", "DAF", "V");
        }
        else {
        writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5: <4}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <10}\t{13: <3}\t{14: <5}\t{15: <5}", 
                 "#CHROM", "POS", "REF", "ALT", "HGVS", "AF", "STR", "IN", "IN_N", "C_REF", "C_ALT", "O_NUC", "OUT_SPC", "RNK", "DAF", "V");
        }
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
                for (pos,qa) in matrix.as_ref().unwrap().iter().enumerate() {
                    
                    // Skip the first and don't go further than the last
                   // if i <= n_matrix {
                       // The first record of the row is the name of the current record
                       // So we have to skip that
                       if pos > 0  {
                           // Get the value at current column 
                           let number = &line[pos];
                           // Put it into a vector 
                           // Since first record is the name, should pos-1 be correct here?
                           let some = distance_struct {name: (&species_record[pos-1]).to_string(), distance: number.parse::<f32>().unwrap()};
                           x.push(some);
                       }
                   // }
                   // i += 1;
                }
                // Sort the vector by distance
                x.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap());
                let mut i = 0;

                for element in &x {
                    //  i > 0? Is it because the first element would be 0 - the same species
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

        if current_record == chrom {
            // Get fields
            let mut pos = &recordx[1];
            let mut reference = &recordx[2];
            let mut alt = &recordx[3];
            let mut hgvs = &recordx[4]; 
            let mut af = &recordx[5];
            let mut actual_strand = &recordx[7];
            let mut v = &recordx[6];

            let mut allele_frequency_numeric = af.parse::<f32>().unwrap();


            let parsed_position = pos.trim().parse::<i32>().expect("Numeric!");
            
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
                if nucleotide == corr_pos && *base_c != b'-' {
                //if nucleotide == corr_pos {

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
                                    
                    if current_record == outgroup_record {
                        disregard += 1;
                    }

                    if cli.hide_utr {
                        if not_cds > 0 {
                                continue;
                            }
                        }
                    
                    if !cli.verbose {
                        if not_cds > 0 {
                                continue;
                            }
                        if disregard > 0 {
                                continue;
                            }
                            
                            
                    }
                    let mut derived_frequency = 0_f32;

                    if base_o == strand_corrected {
                        // 1_f32 - af är rätt
                        derived_frequency = allele_frequency_numeric;
                    } else {
                        derived_frequency = 1_f32 - allele_frequency_numeric;
                    }

                   if cli.verbose {

                        writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5:.2}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <5}\t{13: <5}\t{14: <5}\t{15: <3}\t{16: <3}\t{17: <3}\t{18: <10}\t{19: <3}\t{20: <5}\t{21: <5}", chrom, pos, reference, alt, hgvs.to_string(), allele_frequency_numeric, actual_strand, current_record, *base_c as char, base_vs_reference,*strand_corrected as char, *strand_correctedalt as char, position, nucleotide, gap, *base_o as char, disregard, not_cds, outgroup_record, outgroup_rank, derived_frequency,v);
                    } else {
                        writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5: <4}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <10}\t{13: <3}\t{14: <5}\t{15: <5}", chrom, pos, reference, alt, hgvs.to_string(), allele_frequency_numeric, actual_strand, current_record, *base_c as char, *strand_corrected as char, *strand_correctedalt as char, *base_o as char, outgroup_record, outgroup_rank, derived_frequency,v);
                    }
                    }
            } 

        }
    }
}
}
