use rustody::singlecelldata::SingleCellData;
use rustody::singlecelldata::cell_data::GeneUmiHash;

use rustody::cellids::CellIds;
use rustody::cellids10x::CellIds10x;
use rustody::traits::CellIndex;
use rustody::mapping_info::MappingInfo;

use gtf_gene_structure::GTF::Gtf;
use gtf_gene_structure::ExonIterator::ExonIterator;
use gtf_gene_structure::GTF::QueryErrors;
use gtf_gene_structure::Gene::RegionStatus;

extern crate bam;

 use crate::bam::RecordReader;

use rustody::ofiles::{Ofiles, Fspot};

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};

use std::path::PathBuf;
use std::fs::File;
use std::path::Path;
use std::io::Write;

use std::thread;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;

use clap::Parser;


#[derive(Parser)]
#[clap(version = "1.2.5", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the bam file to quantify
    #[clap(short, long)]
    bam: String,
    /// the gtf file fitting to the Bam file
    #[clap(short, long)]
    gtf: String,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    specie: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the minimum (UMI) reads per cell (sample + genes + antibody combined)
    #[clap(short, long)]
    min_umi: usize,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    version: String,
    /// Optional: end the analysis after processing <max_reads> cell fastq entries 
    #[clap(default_value_t=usize::MAX, long)]
    max_reads: usize,
    /// this is a BD rhapsody or a 10x expression experiment? 
    #[clap(default_value="bd", long)]
    exp: String,
}

fn process_feature(
    bam_feature: &bam::Record, 
    gtf: &Gtf, 
    iterator: &mut ExonIterator, 
    gex: &mut SingleCellData, 
    report: &mut Report, 
    mapping_info: &mut MappingInfo  // Now tracks errors using a HashMap
) {

    let mut cell:Option<String> = None;
    let mut umi:Option<String>  = None;

    // Try to extract the Cell ID (CB), and report if missing
    for (tag, value) in bam_feature.tags() {
        if  tag == b"CB"{
            cb = Some(value.to_string());
        }else if tag == b"UB" {
            umi = Some(value.to_string());
        }
    }

    
    let cell_id = match cell {
        Some(id) => id, // Convert to string, or use empty string
        None => {
            mapping_info.report("missing_CellID");  // Report missing Cell ID
            return;
        }
    };

    // Try to extract the UMI (UB), and report if missing
    let umi_id = match umi {
        Some(u) => u, // Convert to string, or use empty string
        None => {
            mapping_info.report("missing_UMI");  // Report missing UMI
            return;
        }
    };

    // Extract the chromosome (reference name)
    let chr = match bam_feature.reference_name() {
        Some(name) => name.to_string(),  // Get chromosome name as a string
        None => {
            mapping_info.report("missing_Chromosome");  // Report missing chromosome
            return;
        }
    };

    // Extract the start and end positions
    let start = bam_feature.start() as u64;  // BAM is 0-based, start is inclusive
    // crap - this needs to be computed from the CIGAR!!!!
    let end = bam_feature.get_end() as u64;      // BAM is 0-based, end is exclusive

    // Find the gene ID that overlaps or matches the region
    let gene_id = match gtf.find_gene(&chr, start, end, iterator) {
        Some(gene) => gene[0].gene_id(),  // Get the gene_id if found
        None => {
            mapping_info.report("missing_Gene");  // Report missing gene match
            return;
        }
    };

    // Generate a GeneUmiHash
    let guh = GeneUmiHash(gene_id, umi.parse::<u64>().unwrap_or_default()); // Hash UMI as u64

    // Try to insert into SingleCellData (gex)
    if !gex.try_insert(
        &cell_id.parse::<u64>().unwrap_or_default(),  // Parse Cell ID as u64
        guh,
        report
    ) {
        mapping_info.report("UMI_duplicate");
        report.pcr_duplicates += 1;  // Increment duplicate count if insertion fails
    }
}


fn main() {

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    if fs::metadata(&opts.outpath).is_err() {
        if let Err(err) = fs::create_dir_all(&opts.outpath) {
            eprintln!("Error creating directory {}: {}", &opts.outpath, err);
        } else {
            println!("New output directory created successfully!");
        }
    }

    println!("Analysis will stop after having processed {} fastq entries containing a cell info\n", opts.max_reads);

    println!("reading GTF file");    

    let mut reader = bam::BamReader::from_path( &opts.bam , 1).unwrap();
    let mut gtf = Gtf::new(); 

    let _ = gtf.parse_gtf( &opts.gtf );

    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));

    let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");
    pb.set_style(spinner_style);
    pb.set_message( "" );

    let mut reads = 0;
    let mut lines:u64 = 0;
    let split = 1_000_000_u64;

    let mut mapping_info = MappingInfo::new( None, 3.0, opts.max_reads ,None );

    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }
        if lines % split == 0{
            pb.set_message( format!("{} mio reads processed", lines / split) );
            pb.inc(1);
        }
        lines +=1;

        reads += subsetter.process_record( record.clone(), tag );
    }

    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            println!("\nI have selected {reads} reads from the bam file in {milli}h {min}min {sec} sec {mil}milli sec\n");
        },
        Err(e) => {println!("Error: {e:?}");}
    }

    subsetter.print();

    
}
