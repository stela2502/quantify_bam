use rustody::singlecelldata::SingleCellData;
use rustody::singlecelldata::cell_data::GeneUmiHash;
use rustody::int_to_str::IntToStr;
use rustody::singlecelldata::IndexedGenes;

//use rustody::cellids::CellIds;
//use rustody::cellids10x::CellIds10x;
//use rustody::traits::CellIndex;
use rustody::mapping_info::MappingInfo;

use quantify_bam::gtf::{GTF, ExonIterator, RegionStatus };
use quantify_bam::bam_helper::{ get_values, DataTuple };
use quantify_bam::main_logics::{process_data};
//use quantify_bam::main_logics::bam::RecordReader;

extern crate bam;
use crate::bam::RecordReader;

//use rustody::ofiles::{Ofiles, Fspot};

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};

use std::path::PathBuf;
use std::fs;
use std::fs::File;

//use std::thread;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;

use std::collections::HashMap;
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
    num_proc: Option<usize>,
    /// used processor cores (default all)
    #[clap(short, long)]
    min_umi: usize,
    /// tag name for the UMI information (default UB for CellRanger Bam files)
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

    let num_threads = opts.num_proc.unwrap_or_else(rayon::current_num_threads);

    let mut mapping_info = MappingInfo::new(Some(log_file), 3.0, 0, None);
    mapping_info.start_counter();

    // Parse BAM and GTF
    println!("reading GTF file");
    

    let mut gtf = GTF::new();
    gtf.parse_gtf(&opts.gtf).unwrap();

    // Process data
    let (mut gex, genes) = process_data(
        &opts.bam,
        &mut mapping_info,
        &gtf,
        umi_tag,
        num_threads
    );

    // Final reporting and cleanup

    let file_path_sp = PathBuf::from(&opts.outpath).join("BD_Rhapsody_expression");
    println!("Writing data to path {:?}", file_path_sp);

    gex.write_sparse_sub(file_path_sp, &genes, &genes.get_all_gene_names(), opts.min_umi).unwrap();
    mapping_info.log_report();

    println!("The total issues report:\n{}", mapping_info.report_to_string());
    println!("Runtime assessment:\n{}", mapping_info.program_states_string());

}


/*


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

    let log_file_str = PathBuf::from(&opts.outpath).join(
        "Mapping_log.txt"
    );

    let log_file = match File::create( log_file_str ){
            Ok(file) => file,
            Err(err) => {
            panic!("Log file creation Error: {err:#?}" );
        }
    };

    let umi_tag:[u8;2] = match opts.umi_tag{
        Some(tag) => match tag.into_bytes().as_slice() {
                [a, b] => [*a, *b],
                _ => panic!("umi-tag must be exactly 2 chars long"),
            },
            None => *b"UB"
    };

    let mut mapping_info = MappingInfo::new( Some(log_file), 3.0, 0 ,None );
    mapping_info.start_counter();


    //println!("Analysis will stop after having processed {} fastq entries containing a cell info\n", opts.max_reads);

    println!("reading GTF file");    

    let mut reader = bam::BamReader::from_path( &opts.bam , 1).unwrap();
    let header_view = reader.header().to_owned();  // Clone or use as needed

    // Create a HashMap to store reference ID to reference name
    let mut ref_id_to_name: HashMap<i32, String> = HashMap::new();

    // Populate the map with reference ID and names
    for (id, name) in header_view.reference_names().into_iter().enumerate() {
        //let name = String::from_utf8_lossy(name_bytes).to_string();
        ref_id_to_name.insert(id as i32, name.to_string());
    }


    let mut gtf = GTF::new(); 

    let _ = gtf.parse_gtf( &opts.gtf );

    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));

    let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");
    pb.set_style(spinner_style);
    pb.set_message( "" );

    let mut lines = 0_u64;
    
    let mut record = bam::Record::new();

    let num_threads = match &opts.num_proc {
        Some(n) => *n,
        None => rayon::current_num_threads()
    }; 

    let split = BUFFER_SIZE * num_threads;

    println!("Using {} processors", num_threads);

    let mut buffer = Vec::with_capacity(split);
    

    //let mut iterator = ExonIterator::new("main");
    let mut gex = SingleCellData::new(1);
    let mut genes = IndexedGenes::empty( Some(0) );

    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }

        if lines % 1_000_000 == 0 {
            let info = format!("{} mio reads processed", lines / 1_000_000 );
            pb.set_message( info.to_string() );
            pb.inc(1);
            //mapping_info.write_to_log(info);
            //mapping_info.log_report();
        }
        lines += 1;

        let data_tuple = match get_values( &record, &ref_id_to_name, &umi_tag ){
            Ok(res) => res,
            Err( "missing_Chromosome" ) => {
                // We reached the end of the mapped reads - exit the loop
                eprintln!("We have a missing chromosome for this bam entry - assuming we reached the end of the usable data!\n{:?}",record );
                break;
            }
            Err(err) => {
                //eprintln!("ERROR {:?} for this bam entry!\n{:?}",err, record );
                mapping_info.report( err );
                continue;
            }
        };

        buffer.push(data_tuple); // Store only the tuple in the buffer

        if buffer.len() == split {
            let info = format!("{} mio reads -> processing", lines / 1_000_000 );
            pb.set_message( info.to_string() );
            pb.inc(1);
            mapping_info.write_to_log(info);

            mapping_info.stop_file_io_time();

            let chunk_size = (buffer.len() / num_threads).max(1);

            // Process the buffer in parallel
            let results = buffer
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut local_iterator = ExonIterator::new("part");// Initialize a thread-local iterator here if needed
                let mut local_collector = SingleCellData::new(1);
                let mut local_report = MappingInfo::new( None, 3.0, 0 ,None );
                let mut local_genes =  IndexedGenes::empty( Some(0) );
                let mut last_chr = "unset";
                chunk.iter().for_each(|(cell_id, umi, start, cigar, chr, is_reverse_strand )| {
                    if last_chr != *chr {
                        //println!("Init chr!");
                        let _ = gtf.init_search(chr, (*start).try_into().unwrap(), &mut local_iterator);
                        last_chr = chr;
                    }

                    process_feature( 
                        cell_id,
                        umi,
                        *start,
                        cigar,
                        chr,
                        is_reverse_strand,
                        &gtf,
                        &mut local_iterator,
                        &mut local_collector,
                        &mut local_genes,
                        &mut local_report,
                        //&ref_id_to_name,
                    );

                });

                return (local_collector, local_report, local_genes )
            })
            .collect::<Vec<_>>(); // Trigger computation

            mapping_info.stop_multi_processor_time();

            for result in results{
                let translation = genes.merge( &result.2 );
                gex.merge_re_id_genes( result.0, &translation );
                mapping_info.merge( &result.1 );
            }
            mapping_info.log_report();
            buffer.clear(); // Clear the buffer after processing
        }
    }

    if buffer.len() > num_threads { 
        let info = format!("{} mio reads -> processing", lines / 1_000_000 );
        pb.set_message( info.to_string() );
        pb.inc(1);
        mapping_info.write_to_log(info);
        //mapping_info.log_report();

        mapping_info.stop_file_io_time();

        let chunk_size = (buffer.len() / num_threads).max(1);

        // Process the buffer in parallel
        let results = buffer
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local_iterator = ExonIterator::new("part");// Initialize a thread-local iterator here if needed
            let mut local_collector = SingleCellData::new(1);
            let mut local_report = MappingInfo::new( None, 3.0, 0 ,None );
            let mut local_genes =  IndexedGenes::empty( Some(0) );
            let mut last_chr = "unset";
            chunk.iter().for_each(|(cell_id, umi, start, cigar, chr, is_reverse_strand)| {
                if last_chr != *chr {
                    let _ = gtf.init_search(chr, (*start).try_into().unwrap(), &mut local_iterator);
                    last_chr = chr;
                }

                process_feature( 
                    cell_id,
                    umi,
                    *start,
                    cigar,
                    chr,
                    is_reverse_strand,
                    &gtf,
                    &mut local_iterator,
                    &mut local_collector,
                    &mut local_genes,
                    &mut local_report,
                    //&ref_id_to_name,
                );

            });

            return (local_collector, local_report, local_genes )
        })
        .collect::<Vec<_>>(); // Trigger computation

        mapping_info.stop_multi_processor_time();

        for result in results{
            let translation = genes.merge( &result.2 );
            gex.merge_re_id_genes( result.0, &translation );
            mapping_info.merge( &result.1 );
        }

    }else if buffer.len() > 0 {

        let mut local_iterator = ExonIterator::new("part");// Initialize a thread-local iterator here if needed
        let mut local_collector = SingleCellData::new(1);
        let mut local_report = MappingInfo::new( None, 3.0, 0 ,None );
        let mut local_genes =  IndexedGenes::empty( Some(0) );
        #[allow(unused_mut)]
        let mut last_chr = "unset";

        buffer.iter().for_each(|(cell_id, umi, start, cigar, chr, is_reverse_strand)| {
            if last_chr != *chr {
                //println!("Init chr!");
                let _ = gtf.init_search(chr, (*start).try_into().unwrap(), &mut local_iterator);
            }

            process_feature( 
                cell_id,
                umi,
                *start,
                cigar,
                chr,
                is_reverse_strand,
                &gtf,
                &mut local_iterator,
                &mut local_collector,
                &mut local_genes,
                &mut local_report,
                //&ref_id_to_name,
            )

        });

        let translation = genes.merge( &local_genes );
        gex.merge_re_id_genes( local_collector, &translation );
        mapping_info.merge( &local_report );

    }



    let file_path_sp = PathBuf::from(&opts.outpath).join(
        "BD_Rhapsody_expression"
        );

    println!("Writing data to path {:?}", file_path_sp);



    match gex.write_sparse_sub ( file_path_sp, &genes , &genes.get_all_gene_names(), opts.min_umi ) {
        Ok(_) => (),
        Err(err) => panic!("Error in the data write: {err}")
    };

    mapping_info.log_report();

    println!("The total issues report:\n{}", mapping_info.report_to_string());

    println!("Runtime assessment:\n{}",mapping_info.program_states_string());

    mapping_info.write_to_log( format!("Runtime assessment:\n{}",mapping_info.program_states_string()));

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

*/
