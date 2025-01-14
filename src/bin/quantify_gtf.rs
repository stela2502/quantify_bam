
//use rustody::cellids::CellIds;
//use rustody::cellids10x::CellIds10x;
//use rustody::traits::CellIndex;
use rustody::mapping_info::MappingInfo;

use quantify_bam::gtf::GTF;
use quantify_bam::main_logics::{process_data};
//use quantify_bam::main_logics::bam::RecordReader;

extern crate bam;

//use rustody::ofiles::{Ofiles, Fspot};


use std::path::PathBuf;
use std::fs;
use std::fs::File;

//use std::thread;
//use rayon::prelude::*;

use std::time::SystemTime;

use clap::Parser;



#[derive(Parser)]
#[clap(version = "0.4.1", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the bam file to quantify
    #[clap(short, long)]
    bam: String,
    /// the gtf file fitting to the Bam file
    #[clap(short, long)]
    gtf: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the minimum (UMI) reads per cell (sample + genes + antibody combined)
    #[clap(short, long)]
    min_umi: usize,
    /// used processor cores (default all)
    #[clap(short, long)]
    num_proc: Option<usize>,
    /// tag name for the CELL information (default CB for velocity default - change to CR for CellRanger)
    #[clap(short, long)]
    cell_tag:Option<String>,
    /// tag name for the UMI information (default UB for velocity default - change to UR for CellRanger)
    #[clap(short, long)]
    umi_tag:Option<String>,
}


// Main function
fn main() {
    let now = SystemTime::now();
    let opts: Opts = Opts::parse();

    // Create output directory if needed
    if fs::metadata(&opts.outpath).is_err() {
        if let Err(err) = fs::create_dir_all(&opts.outpath) {
            eprintln!("Error creating directory {}: {}", &opts.outpath, err);
        } else {
            println!("New output directory created successfully!");
        }
    }

    // Set up logging
    let log_file_str = PathBuf::from(&opts.outpath).join("Mapping_log.txt");
    let log_file = File::create(log_file_str).expect("Failed to create log file");
    let umi_tag: [u8; 2] = opts.umi_tag.unwrap_or_else(|| "UB".to_string()).into_bytes().try_into().expect("umi-tag must be exactly 2 chars long");
    let cell_tag: [u8; 2] = opts.cell_tag.unwrap_or_else(|| "CB".to_string()).into_bytes().try_into().expect("umi-tag must be exactly 2 chars long");

    let num_threads = opts.num_proc.unwrap_or_else(rayon::current_num_threads);

    let mut mapping_info = MappingInfo::new(Some(log_file), 3.0, 0, None);
    mapping_info.start_counter();

    // Parse BAM and GTF
    println!("reading GTF file");
    

    let mut gtf = GTF::new(None);
    gtf.parse_gtf(&opts.gtf).unwrap();

    // Process data
    let (mut gex, genes) =  match process_data(
        &opts.bam,
        &mut mapping_info,
        &gtf,
        cell_tag,
        umi_tag,
        num_threads
    ){
        Ok(ret) => ret,
        Err(e) => {
            panic!("{e:?}");
        }
    };

    // Final reporting and cleanup

    let file_path_sp = PathBuf::from(&opts.outpath).join("BD_Rhapsody_expression");
    println!("Writing data to path {:?}", file_path_sp);

   

    let info = match gex.write_sparse_sub(file_path_sp, &genes, &genes.get_all_gene_names(), opts.min_umi){
        Ok(i) => i,
        Err(e) => panic!("While writing the data an error occures: {e:?} "),
    };
    mapping_info.write_to_log( info );
    mapping_info.log_report( );

    println!("The total issues report:\n{}", mapping_info.report_to_string());
    println!("Runtime assessment:\n{}", mapping_info.program_states_string());

    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            mapping_info.write_to_log( format!("\nfinished in {milli}h {min}min {sec} sec {mil}milli sec\n") );
            println!("\nI have analyzed the bam file in {milli}h {min}min {sec} sec {mil}milli sec\n");
        },
        Err(e) => {println!("Error: {e:?}");}
    }

    mapping_info.report_to_string();
}

