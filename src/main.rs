extern crate seq_io;
use std::str;
use regex;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::io::{self, Write};

use seq_io::fasta::{Reader,Record};
use lazy_static::lazy_static;

use clap::{Parser};

#[derive(Parser,Debug)]
#[clap(name = "basefinder")]
#[clap(version = "0.0.1")]
#[clap(about = "Finds corresponding base in another species", long_about = None)]

struct Cli {
    /// Sequence alignment in FASTA-format containing one sequence per species.
    #[clap(short,long)]
    alignment: String,
    /// A _headerless_ tab-separated file containing 8 columns extracted from a VCF-file:
    ///
    /// #CHROM    POS    ALT    REF    HGVS.c    AF    VARIANT    STRAND
    ///
    /// Example:
    ///
    /// elon_1_1    166    G    A    c.-131G>A    0.0625    synonymous_variant    +
     #[clap(short,long, verbatim_doc_comment)]
    tsv: String,
    /// hide non-cds rows (only active in verbose mode)
    #[clap(short, long)]
    show_utr: bool,
     /// Hide header
    #[clap(long)]
    no_header: bool,  
    /// Ingroup
    #[clap(short, long)]
    ingroup: String,
    /// Outgroup
    #[clap(short, long)]
    outgroup: String,
}

lazy_static! { 
        static ref digits_regex: regex::Regex = regex::Regex::new(r"\d+").unwrap();
        static ref alphabetic_regex: regex::Regex = regex::Regex::new(r"[[:alpha:]]+").unwrap();
}

fn main() {
    // Parse command line arguments
    let cli = Cli::parse();

    let user_specified_ingroup = cli.ingroup.as_str();
    let user_specified_outgroup = cli.outgroup.as_str();

    // Read files
    let mut fasta_reader = Reader::from_path(&cli.alignment).unwrap();
    let mut vcf_reader = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(File::open(&cli.tsv).unwrap());
    
    // Collect results in memory
    let records: Result<Vec<_>, _> = fasta_reader.records().collect();
    
    //Create writer
    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = io::BufWriter::new(handle);
    
    // First, lets go through the FASTA file, extract only alphabetic letters from the FASTA-file
    // (the name of species), and store it.
    let mut fasta_records = Vec::new();

    let mut outgroup_id = 0;
    for (i, fasta_record) in records.as_ref().unwrap().iter().enumerate() {
       
        let record_header = get_header(fasta_record);
        let species = extract_alphabetic_letters(record_header);
        fasta_records.push(species);
        
        if species == user_specified_outgroup {
            outgroup_id = i;
        }
    }

    // Check no duplicated species, ingroup and outgroup exists
    check_fasta(&fasta_records, &user_specified_ingroup, &user_specified_outgroup);
    // Get outgroup name
    let outgroup_record = get_header(&records.as_ref().unwrap()[outgroup_id]);
    let outgroup_sequence = records.as_ref().unwrap()[outgroup_id].seq();

    // Write header
    if !cli.no_header {
        writeln!(writer, "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5: <4}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <10}\t{13: <3}\t{14: <5}\t{15: <5}\t{16: <5}", 
                 "#CHROM", "POS", "REF", "ALT", "HGVS", "AF", "STR", "IN", "IN_N", "C_REF", "C_ALT", "O_NUC", "OUT_SPC", "RNK", "DAF", "V", "DERIVED");
    }
    
    // For every record in FASTA file
    for record in records.as_ref().unwrap().iter() {
        
        let mut current_record = get_header(record);
        let current_species = extract_alphabetic_letters(get_header(record));   
       
        // Skip cases 
        if (get_header(record) == outgroup_record) || (current_species != user_specified_ingroup) {
            continue;
        }

        let mut outgroup_rank = 0;

        // And now we go through every record in the VCF/tsv file
        for result in vcf_reader.records() {
            let tsv_record = result.expect("Could not open TSV record");
            // Parse fields
            let chrom = &tsv_record[0];
            // But check first so no need to go through alignment
            if get_header(record) != chrom {
                continue;
            }
            // Parse fields
            let pos = &tsv_record[1].trim().parse::<i32>().expect("POS is not numeric.");
            let reference = &tsv_record[2].as_bytes()[0];
            let alt = &tsv_record[3].as_bytes()[0];
            let hgvs = &tsv_record[4].to_string(); 
            let allele_frequency_numeric = tsv_record[5].parse::<f32>().expect("AF is not numeric.");
            let variant = &tsv_record[6];
            let strand = &tsv_record[7].as_bytes()[0];
            
            let regex_capture = digits_regex
                .captures(hgvs)
                .expect("Failed to capture HGVS digits");
            
            let hgvs_position = regex_capture
                .get(0)
                .map_or("", |m| m.as_str())
                .parse::<i32>()
                .expect("Failed to parse HGVS position");

            let mut nucleotide_count = 0;
           
            // Go through every base 
            for (base_c, base_o) in record.seq().iter().zip(outgroup_sequence.iter()) {
                // If current base is gap add cap count, else add nucleotide count
                if !is_gap(base_c) {
                    nucleotide_count += 1;
                }
                
                // To possible corrent for strand 
                let mut strand_corrected = reference;
                let mut strand_correctedalt = alt;

                // So if the position matches the HGVS potition and it's not a gap
                if nucleotide_count == hgvs_position && !is_gap(base_c) {
                   
                    // If current base is not equal to the reference base
                    // we need to check strand 
                    if base_c != reference {
                        // If the strand is minus, then make put the corrent bases 
                        // in "strand_corrected" / "strand_correctedalt"
                        if *strand == b'-' {
                            strand_corrected = get_complement_base(&reference);
                            strand_correctedalt = get_complement_base(&alt);
                        }    
                    } 
                    
                    // Then check if none (?) of the strand corrented bases 
                    // is not equal to the "outgroup" base, then we can't use it (?)

                    if strand_corrected != base_o && strand_correctedalt != base_o {
                        continue;
                    }
                   
                    // Hide not CDS regions because they are wrong.
                    if !cli.show_utr && is_cds(&hgvs) {
                        continue;
                    }
                    
                    // If the outgroup base is the same as the strand corrected referece, then 
                    // it's ancestral (?) and the derived allele freq should be 1 - alterative
                    // allele freq (?)
                
                    let derived_allele = get_derived_base(&base_o, &strand_corrected, &strand_correctedalt);
                    let derived_frequency = get_derived_freq(&base_o, &strand_corrected, &allele_frequency_numeric);
                    
                    // Write results
                    writeln!(writer, 
                            "{0: <10}\t{1: <5}\t{2: <3}\t{3: <3}\t{4: <10}\t{5: <4}\t{6: <3}\t{7: <10}\t{8: <3}\t{9: <3}\t{10: <3}\t{11: <3}\t{12: <10}\t{13: <3}\t{14: <5}\t{15: <5}\t{16: <5}",
                            chrom, pos, *reference as char, *alt as char, hgvs, allele_frequency_numeric, *strand as char, current_record, *base_c as char, *strand_corrected as char, *strand_correctedalt as char, *base_o as char, outgroup_record, outgroup_rank, derived_frequency, variant, derived_allele as char);
                    }
                } 
            }
        }
    }

fn extract_alphabetic_letters(text: &str) -> &str {
    let captured_letters = alphabetic_regex.captures(text).unwrap();
    let captured_letters = captured_letters.get(0).map_or("", |x| x.as_str());
    return captured_letters;
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn contains_duplicates(fasta_records: &Vec<&str>) -> bool {
    let fasta_records_clone = &fasta_records;

    let mut fasta_records_copy = fasta_records.to_vec();
    fasta_records_copy.sort();
    let fasta_records_length = fasta_records_copy.len();
    fasta_records_copy.dedup();
    
    // If unequal length, contains duplicates 
    fasta_records_length != fasta_records_copy.len() 
      
}

fn get_header(fasta_record: &seq_io::fasta::OwnedRecord) -> &str {
   str::from_utf8(fasta_record.head()).expect("Couldn't get FASTA header")
}

fn get_complement_base(base: &u8) -> &u8 {
    match base {
        b'A' => &b'T',
        b'C' => &b'G',
        b'G' => &b'C',
        b'T' => &b'A',
           _ => panic!("Base not A, C, G or T."),
    }
}

fn is_cds(hgvs: &str) -> bool {
    let hgvs_string = &hgvs.to_string();
    // If the HGVS contains a - or *, then it is not CDS
    if hgvs_string.contains("-") || hgvs_string.contains("*") {
        true
    } else {
        false
    }
}

fn is_gap(base: &u8) -> bool {
    base == &b'-'
}

fn check_fasta(fasta_records: &Vec<&str>, user_specified_ingroup: &str, user_specified_outgroup: &str) {

    if contains_duplicates(&fasta_records) {
        panic!("Contains duplicate species");
    }    
    // Check that it contains the specified ingroup
    if !fasta_records.contains(&user_specified_ingroup) {
        panic!("Ingroup not found");
    }
    // Check that it contains the specified outgroup, if so save the position of the outgroup in
    // the FASTA-file
    if !fasta_records.contains(&user_specified_outgroup) {
        panic!("Outgrup not found");
    } 
}

fn get_derived_base(base_o: &u8, strand_corrected: &u8, strand_corrected_alt: &u8) -> u8 {
     if base_o == strand_corrected {
                        *strand_corrected_alt
                    } else {
                        *strand_corrected
                    }
}

fn get_derived_freq(base_o: &u8, strand_corrected: &u8, allele_frequency_numeric: &f32) -> f32 {
     if base_o == strand_corrected {
        return *allele_frequency_numeric;
    } else {
        return 1_f32 - allele_frequency_numeric;
    }

}
