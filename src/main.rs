extern crate seq_io;

use std::env;
use std::str;

use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

use seq_io::fasta::{Reader,Record};

use std::process;
use std::collections::HashMap;

#[derive(Debug)]

struct Person {
    name: String,
    age: f32
}

fn main() {
    
    let args: Vec<String> = env::args().collect();

    let mut reader = Reader::from_path(&args[1]).unwrap();

    let mut ranking: Vec<String>;
    
    let lines = lines_from_file(&args[2]);
    
//    let positions = lines_from_file(&args[3]);

    let mut rdr = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(File::open(&args[3]).unwrap());
    let mut rdr_matrix = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(File::open(&args[4]).unwrap());

    let csv_records: Result<Vec<csv::StringRecord>, _> = rdr.records().collect();
    let matrix: Result<Vec<csv::StringRecord>, _> = rdr_matrix.records().collect();
    

    // This is the slow way, might actually work the fast way now
    let records: Result<Vec<_>, _> = reader.records().collect();
    
    let length = records.as_ref().unwrap().iter().len();

            let mut outgroup_rank = 0;
            let mut outgroup_id = 0;
    let mut i = 0;
    
    //let mut outgroup_id = 0;
    // Check if highest ranking outgroup is present
    // So for each line
   // '_outer: for line in lines {
        // Go through each record 
        
     //   let mut records_id = 0;

       // '_inner: for record in records.as_ref().unwrap().iter() {
            // Save the current record 
         //   let current_record = str::from_utf8(record.head()).unwrap();
            // And if the current record matches the line 
           // if current_record.contains(&line) {
                //println!("{}", current_record);
                //println!("{}", line);
             //   i += 1;
               // outgroup_id = records_id;
               // break '_outer;

           // }
            //records_id += 1;
       // }
    //}
    //if i == 0 {
      //  eprintln!("{}", "no outgroup found!");
       // process::exit(1);
    //}
    //println!("{}", outgroup_id);

    println!("{0: <18}\t{1: <5}\t{2: <5}\t{3: <5}\t{4: <18}\t{5: <5}\t{6: <5}\t{7: <5}\t{8: <5}\t{9: <5}\t{10: <18}", "VCF_SPC", "VCF_POS", "VCF_REF", "VCF_ALT", "IN_SPC", "IN_NUC", "ALI_POS", "IN_POS", "GAP", "OUT_NUC", "OUT_SPC");
    let mut species_record: Vec<&str> = Vec::with_capacity(24);
    for ma in matrix.as_ref().unwrap() {
        species_record.push(&ma[0].as_ref());
    }
    // For every FASTA record
    for record in records.as_ref().unwrap().iter() {
        let current_record = str::from_utf8(record.head()).unwrap();
        
        // Här någonstans välja rätt 
        
    // x
    
    // Go through every matrix row
    for ma in matrix.as_ref().unwrap() {
        //Match the FASTA sequece with the Matrix 
            let mut test: Vec<&String> = Vec::new();
            let mut x = Vec::new();
            if current_record.contains(&ma[0]) {

                //println!("{} {}", current_record, &ma[0]);
            
                let mut i = 0;
            // Go through all columns
                for qa in matrix.as_ref().unwrap() {
                
                    if i <= matrix.as_ref().unwrap().len() {
                       if i > 1 {
                           let number = &ma[i];
                           //println!("{} {}", number, &species_record[i-1]);
                           let some = Person {name: (&species_record[i-1]).to_string(), age: number.parse::<f32>().unwrap()};
                           x.push(some);
                       }
                    }
            
                    i += 1;
                }

                x.sort_by(|a, b| a.age.partial_cmp(&b.age).unwrap());
                //println!("Current: {} {:?}", &ma[0], x);
                let mut i = 0;
                for element in &x {
                    if i > 0 {
                        //let string = String:
                        //println!("{} {}", current_record, element.name);
                        //if current_record.contains(&element.name) {
                            test.push(&element.name);
                        //}
                    }
                    i += 1; 
                } 
                //println!("Current: {} {:?}", &ma[0], test);
                //println!("{:?}", test);
        
                //println!("{:?}", species);
            } // If the current fasta record contains
            '_outer: for ex in test {

                //println!("{}", ex);
            let mut i = 0;
            '_inner: for record in records.as_ref().unwrap().iter() {
                let head = str::from_utf8(record.head()).unwrap();
                if head.contains(ex) {
                    outgroup_id = i;
                    //println!("{} {}", head, ex);
                    break '_outer;
                }
                i += 1;
            }    
            outgroup_rank += 1;
        }
        
    }
    

        let outgroup_record = str::from_utf8(records.as_ref().unwrap()[outgroup_id].head()).unwrap();
        
        if current_record != outgroup_record {

            let current_sequence = record.seq();
            //println!("{:?}", str::from_utf8(current_sequence));
            let outgroup_sequence = records.as_ref().unwrap()[outgroup_id].seq();
            //println!("{:?}", str::from_utf8(outgroup_sequence));
            
            let mut gap = 0;
            let mut nucleotide = 0;
            let mut position = 0;
            
            for (base_c, base_o) in current_sequence.iter().zip(outgroup_sequence.iter()) {
                if base_c == &b'-' {
                    gap += 1;
                } else {
                    nucleotide += 1;
                }

            position += 1;
            
            for recordx in csv_records.as_ref().unwrap() {
                let chrom = &recordx[0];
                let pos = &recordx[1];
                let reference = &recordx[2];
                let alt = &recordx[3];
                let strand = &recordx[4];

                let temp = pos.trim().parse::<i32>().expect("Numeric!");
                if current_record.contains(&chrom) {
                    if nucleotide == temp {
                        //println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", *base_c as char, nucleotide, gap, position, *base_o as char, temp, chrom, outgroup_record);
                        println!("{0: <18}\t{1: <5}\t{2: <5}\t{3: <5}\t{4: <18}\t{5: <5}\t{6: <5}\t{7: <5}\t{8: <5}\t{9: <5}\t{10: <18}\t{11: <5}", chrom, pos, reference, alt, current_record, *base_c as char, position, nucleotide, gap, *base_o as char, outgroup_record, outgroup_rank);
                
                    }
                } 
            }

            }
        }
    }
    
}

fn lines_from_file(filename: impl AsRef<Path>) -> Vec<String> {
    let file = File::open(filename).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}
