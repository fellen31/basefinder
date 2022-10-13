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
    /// hide non-cds rows (only active in verbose mode)
    #[clap(short, long)]
    hide_utr: bool,
     /// Hide header
    #[clap(long)]
    no_header: bool,  
    /// Verbose
    #[clap(short,long)]
    verbose: bool,
    /// Groups
    #[clap(short,long)]
    groups: bool,
    /// Ingroup
    #[clap(short, long)]
    ingroup: Option<String>,
    /// Outgroup
    #[clap(short, long)]
    outgroup: Option<String>,
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
    writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5: <4}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <5}\t{13: <5}\t{14: <5}\t{15: <3}\t{16: <3}\t{17: <3}\t{18: <10}\t{19: <3}\t{20: <5}\t{21: <5}\t{22: <5}", "#CHROM", "POS", "REF", "ALT", "HGVS", "AF", "STR", "IN", "IN_N", "!=", "C_REF", "C_ALT", "ALPOS", "INPOS", "GAP", "O_NUC", "!KEEP", "!CDS", "OUT_SPC", "RNK", "DAF", "V", "DERIVED");
        }
        else {
        writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5: <4}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <10}\t{13: <3}\t{14: <5}\t{15: <5}\t{16: <5}", 
                 "#CHROM", "POS", "REF", "ALT", "HGVS", "AF", "STR", "IN", "IN_N", "C_REF", "C_ALT", "O_NUC", "OUT_SPC", "RNK", "DAF", "V", "DERIVED");
        }
    }
    
    // Get which species are present in distance matrix 
    // Number of species in matrix 
    let n_matrix = matrix.as_ref().unwrap().iter().len();
    // Make vector to contain matrix species 
    let mut species_record: Vec<&str> = Vec::with_capacity(n_matrix);
    
    // Add this regex up here (not sure if it makes a difference since its lazy, but might be
    // faster than inside loop).
    lazy_static! { 
        static ref RE: regex::Regex = regex::Regex::new(r"\d+").unwrap();
        static ref letters: regex::Regex = regex::Regex::new(r"[[:alpha:]]+").unwrap();
    }
    
    // Add every species (column 0) to species_record
    for species in matrix.as_ref().unwrap() {
        species_record.push(&species[0].as_ref());
    }
    if let Some(ref x) = cli.ingroup {
                let mut fasta_records = Vec::new();
                for record in records.as_ref().unwrap().iter() {
                    let h = str::from_utf8(record.head()).unwrap();
                    let cap = letters.captures(h).unwrap();
                    let cap = cap.get(0).map_or("", |m| m.as_str());
                    fasta_records.push(cap);
            }
                if !fasta_records.contains(&x.as_str()) {
                    panic!("Ingroup not found");
                }
            }

    // For every record in FASTA file
    for record in records.as_ref().unwrap().iter() {
 
        let mut outgroup_rank = 0;
        let mut outgroup_id = 0;
        // Name of current record
        let mut current_record = str::from_utf8(record.head()).unwrap();
            
            // If specifying ingroup, regex letters from fasta and skip records not matching CLI
            if let Some(ref x) = cli.ingroup {
                    let cap = letters.captures(current_record).unwrap();
                    let cap = cap.get(0).map_or("", |m| m.as_str());
                    
                if cap != x {
                    continue;
                }
            }
 
        // Go through every matrix row
        for line in matrix.as_ref().unwrap() {
            
            let mut sorted_distance_name: Vec<&String> = Vec::new();
            let mut distance_vector = Vec::new();
            
            // Match the FASTA sequece contains the distance matrix species
            if current_record.contains(&line[0]) {
                // Could probably be done with enumerate
                let mut i = 0;
            
                // Then go through all columns
                for (n,col) in matrix.as_ref().unwrap().iter().enumerate() {
                    
                       // The first record of the row is the name of the current record
                       // So we have to skip that
                       if n > 0  {
                           // Get the value at current column 
                           let number = &line[n];
                           // Put it into a vector 
                           // Since first record is the name, should pos-1 be correct here?
                           let species_distance = distance_struct {
                               name: (&species_record[n-1]).to_string(), 
                               distance: number.parse::<f32>().unwrap()
                           };

                           distance_vector.push(species_distance);
                       }
                   // }
                   // i += 1;
                }
                // Sort the vector by distance
                distance_vector.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap());
                
                let mut i = 0;
                
                // Might be a better way, but get the names in order sorted by distance 
                // and put in another vector
                for element in &distance_vector {
                    //  i > 0? Is it because the first element would be 0 - the same species
                    if i > 0 {
                            sorted_distance_name.push(&element.name);
                    }
                    i += 1; 
                } 
            } 

            // When distance record matches fasta record, break the loop 
            // and get the id of the fasta record
            
            // But first, if outgruop is specified as CLI arg, 
            // check if fasta contains speciefied outgrup by 
            // selecting only letters in fasta records 
            if let Some(ref x) = cli.outgroup {
                let mut fasta_records = Vec::new();
                for record in records.as_ref().unwrap().iter() {
                    let h = str::from_utf8(record.head()).unwrap();
                    let cap = letters.captures(h).unwrap();
                    let cap = cap.get(0).map_or("", |m| m.as_str());
                    fasta_records.push(cap);
            }
                if !fasta_records.contains(&x.as_str()) {
                    panic!("Outgroup not found");
                }
            }
            


            '_outer: for spec in sorted_distance_name {
                let mut i = 0;
                    '_inner: for record in records.as_ref().unwrap().iter() {
                        let head = str::from_utf8(record.head()).unwrap();
                            if let Some(ref y) = cli.outgroup {
                                if head.contains(y) {
                                    outgroup_id = i;
                                    break '_outer;
                                } 
                            } else { 
                                if head.contains(spec) {
                                    outgroup_id = i;
                                    break '_outer;
                                }
                } 
                i += 1;
            }    
            outgroup_rank += 1;
        } 
    }
    
    // Get the outgroup FASTA record from the outgroup id we found above 
    let outgroup_record = str::from_utf8(records.as_ref().unwrap()[outgroup_id].head()).unwrap();
    
    // And now we go through every record in the VCF/tsv file
    for tsv_record in csv_records.as_ref().unwrap() {
        
        let chrom = &tsv_record[0];
        // Check first so no need to go through everything
        if current_record == chrom {
            // Get fields
            let mut pos = &tsv_record[1];
            let mut reference = &tsv_record[2];
            let mut alt = &tsv_record[3];
            let mut hgvs = &tsv_record[4]; 
            let mut af = &tsv_record[5];
            let mut actual_strand = &tsv_record[7];
            let mut v = &tsv_record[6];

            let mut allele_frequency_numeric = af.parse::<f32>().unwrap();
            let parsed_position = pos.trim().parse::<i32>().expect("Numeric!");
            
            // Get FASTA sequences
            let current_sequence = record.seq();
            let outgroup_sequence = records.as_ref().unwrap()[outgroup_id].seq();
            

            let mut gap = 0;
            let mut nucleotide = 0;
            let mut position = 0;
            
            // Go through every base in alignment...can this go on the outside of csv-loop?
            // For every base in current record (base_c) and "outgroup" (nearest species) (base_o) 
            for (base_c, base_o) in current_sequence.iter().zip(outgroup_sequence.iter()) {
                
                // If current base is gap add cap count, else add nucleotide count
                if base_c == &b'-' {
                    gap += 1;
                } else {
                    nucleotide += 1;
                }
                // Add position count
                position += 1;
                
                let mut disregard = 0; 
                let mut not_cds = 0; 
                
                // Convert reference and alt to bytes
                let bref = reference.as_bytes()[0];
                let balt = alt.as_bytes()[0];

                // To corrent for strand 
                let mut strand_corrected = &bref;
                let mut strand_correctedalt = &balt;

                let hgvs_string = &hgvs.to_string();
                    
                // If the HGVS contains a - or *, then it is not CDS
                if hgvs_string.contains("-") || hgvs_string.contains("*") {
                    not_cds += 1;
                }
                
                // Do a regex for digits (the position) in HGVS
                let cap = RE.captures(hgvs_string).unwrap();
                
                let mut base_vs_reference = 0;
                // Get corrected position from the HGVS
                let corr_pos = cap.get(0).map_or("", |m| m.as_str()).parse::<i32>().unwrap();
                
                // So if the position matches the HGVS potition and it's not a gap
                if nucleotide == corr_pos && *base_c != b'-' {
                   
                    // If current base is not equal to the reference base
                    // we need to check strand 
                    if *base_c != reference.as_bytes()[0] {
                        // Keep track of if current base does not matches the reference base
                        base_vs_reference = 1;
                        
                        // If the strand is minus, then make put the corrent bases 
                        // in "strand_corrected" / "strand_correctedalt"
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
                    
                    // Then check if none (?) of the strand corrented bases 
                    // is not equal to the "outgroup" base, then we can't use it (?)

                    if strand_corrected != base_o && strand_correctedalt != base_o {
                        disregard += 1;
                    }
                    
                    // If the current FASTA record is the outgroup record, disregard.
                    if current_record == outgroup_record {
                        disregard += 1;
                    }

                    // If we want to not write UTRs / not CDS
                    if cli.hide_utr {
                        if not_cds > 0 {
                                continue;
                            }
                        }

                    // If not verbose don't print UTRs either (?)
                    if !cli.verbose {
                        if not_cds > 0 {
                                continue;
                            }
                        if disregard > 0 {
                                continue;
                            }
                    }

                    let mut derived_frequency = 0_f32;
                    
                    // If the outgroup base is the same as the strand corrected referece, then 
                    // it's ancestral (?) and the derived allele freq should be 1 - alterative
                    // allele freq (?)
                    //
                    let mut derived = b'X' as char;
                    if base_o == strand_corrected {
                        derived = *strand_correctedalt as char;
                        derived_frequency = allele_frequency_numeric;
                    } else {
                        derived_frequency = 1_f32 - allele_frequency_numeric;
                        derived = *strand_corrected as char;
                    }

                   if cli.verbose {

                        writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5:.2}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <5}\t{13: <5}\t{14: <5}\t{15: <3}\t{16: <3}\t{17: <3}\t{18: <10}\t{19: <3}\t{20: <5}\t{21: <5}\t{22: <5}", chrom, pos, reference, alt, hgvs.to_string(), allele_frequency_numeric, actual_strand, current_record, *base_c as char, base_vs_reference,*strand_corrected as char, *strand_correctedalt as char, position, nucleotide, gap, *base_o as char, disregard, not_cds, outgroup_record, outgroup_rank, derived_frequency,v,derived);
                    } else {
                        writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5: <4}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <10}\t{13: <3}\t{14: <5}\t{15: <5}\t{16: <5}", chrom, pos, reference, alt, hgvs.to_string(), allele_frequency_numeric, actual_strand, current_record, *base_c as char, *strand_corrected as char, *strand_correctedalt as char, *base_o as char, outgroup_record, outgroup_rank, derived_frequency,v,derived);
                    }
                    }
            } 

        }
    }
}
}
