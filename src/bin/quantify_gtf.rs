use rustody::singlecelldata::SingleCellData;
use rustody::singlecelldata::cell_data::GeneUmiHash;
use rustody::int_to_str::IntToStr;
use rustody::singlecelldata::IndexedGenes;

use rustody::cellids::CellIds;
use rustody::cellids10x::CellIds10x;
use rustody::traits::CellIndex;
use rustody::mapping_info::MappingInfo;

use quantify_bam::gtf::{GTF, SplicedRead, ExonIterator, RegionStatus};

use bam::record::tags::StringType;

extern crate bam;
use bam::record::Record;
use bam::record::tags::TagViewer;
use bam::record::tags::TagValue;

use crate::bam::RecordReader;

use rustody::ofiles::{Ofiles, Fspot};

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};

use std::path::PathBuf;
use std::fs;
use std::fs::File;
use std::path::Path;
use std::io::Write;

use std::thread;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;

use std::collections::HashMap;
use std::time::SystemTime;

use clap::Parser;

const BUFFER_SIZE: usize = 1_000_000;



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


/*
pub enum TagValue<'a> {
    Char(u8),
    Int(i64, IntegerType),
    Float(f32),
    String(&'a [u8], StringType),
    IntArray(IntArrayView<'a>),
    FloatArray(FloatArrayView<'a>),
*/

/// get_tag will only return string tags - integer or floats nbeed to be implemented when needed.
fn get_tag(bam_feature: &Record, tag: &[u8;2]) -> Option<String> {
    match bam_feature.tags().get(tag) {
        Some(TagValue::String(s, StringType::String)) => {
            String::from_utf8(s.to_vec()).ok()
        }, // Handle String directly
        Some(TagValue::String(bytes, StringType::Hex )) => Some(
            bytes.iter().map(|b| format!("{:02X}", b)).collect() // Convert byte array to hex string
        ),
        _ => None, // Return None for non-string-like types
    }
}


//fn get_values( bam_feature: &bam::Record, chromosmome_mappings:&HashMap<i32, String> ) 
//        -> Result<( String, String, u32, String, String ), &str> {
/// get_values is very specififc for this usage here - read bam reads and compare them to a gtf file
/// Bam is 0-based start position and end exclusive whereas gtf is 1-based and end inclusive.
/// This means that the end is actually the same value - just the start for GTF needs to be bam start +1
/// That is what this function returns for a start!
fn get_values<'a>( bam_feature: &'a bam::Record, chromosmome_mappings:&'a HashMap<i32, String> ) 
         -> Result<( String, String, i32, String, String ), &'a str> {
    
    let cell_id = match get_tag( bam_feature, b"CB" ) {
        Some(id) => id.to_string(), // Convert to string, or use empty string
        None => {
            return Err( "missing_CellID" )
            //mapping_info.report("missing_CellID");  // Report missing Cell ID
            //return;
        }
    };

    // Try to extract the UMI (UB), and report if missing
    let umi = match get_tag( bam_feature, b"UB" ) {
        Some(u) => u.to_string(), // Convert to string, or use empty string
        None => {
            return Err( "missing_UMI");  // Report missing UMI
        }
    };

    // Extract the chromosome (reference name)
    let chr = match chromosmome_mappings.get( &bam_feature.ref_id() ){
        Some(name) => name.to_string(),
        None => {
            return Err( "missing_Chromosome");  // Report missing chromosome
        }
    };

    // Extract the start and end positions
    let start = bam_feature.start();  // BAM is 0-based, start is inclusive
    // crap - this needs to be computed from the CIGAR!!!!

    Ok( ( cell_id, umi, start+1, bam_feature.cigar().to_string(), chr) ) // convert bam to gtf notation
}


fn process_feature(
    cell_id: &str, 
    umi: &str, 
    start: i32, 
    cigar: &str, 
    chr: &str,
    gtf: &GTF, 
    iterator: &mut ExonIterator, 
    gex: &mut SingleCellData,
    genes: &mut IndexedGenes,
    mapping_info: &mut MappingInfo,  // Now tracks errors using a HashMap
    chromosmome_mappings: &HashMap<i32, String>
) {

    // Find the gene ID that overlaps or matches the region
    /*
    pub enum RegionStatus {
    InsideExon,
    SpanningBoundary,
    InsideIntron,
    AfterGene,
    BeforeGene,
    }
    */
    let gene_id = match gtf.match_cigar_to_gene(&chr, &cigar, start.try_into().unwrap() , iterator) {
        Some(read_result) => {
            match read_result.match_type{
                RegionStatus::InsideExon => {
                    read_result.gene.to_string()
                },
                RegionStatus::SpanningBoundary => {
                    read_result.gene.to_string() + "_unspliced"
                },
                RegionStatus::InsideIntron => {
                    read_result.gene.to_string() + "_unspliced"
                }
                _ => {
                    mapping_info.report("missing_Gene");
                    return;
                }
            }
        },  // Get the gene_id if found
        None => {
            mapping_info.report("missing_Gene");  // Report missing gene match
            return;
        }
    };
    
    // Generate a GeneUmiHash
    // This needs to be fixed - how do I use the genes here?! I totally forgot - but it should be rather straight forward!
    #[cfg(debug_assertions)]
    println!("Chr {chr}, cigar {cigar} and this start position {start} : CellID: {cell_id}");

    let umi_u64 = IntToStr::new( umi.into(), 32 ).into_u64();

    let gene_id_u64 = genes.get_gene_id( &gene_id );

    let guh = GeneUmiHash(gene_id_u64, umi_u64 );

    #[cfg(debug_assertions)]
    println!("\t And I got a gene: {guh}");

    // Try to insert into SingleCellData (gex)
    // the cellid is already a numeric!
    let cell_id_u64:u64 =  match cell_id.parse::<u64>() {
        Ok(number) => number,
        Err(e) => IntToStr::new( cell_id.into(), 32).into_u64(),
    };

    if !gex.try_insert(
        &cell_id_u64,
        guh,
        mapping_info
    ) {
        mapping_info.report("UMI_duplicate");
    }
    
}


fn main() {

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    let mut mapping_info = MappingInfo::new( None, 3.0, opts.max_reads ,None );
    mapping_info.start_counter();

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

    let mut lines:u64 = 0;
    let split = 1_000_000_u64;
    let mut record = bam::Record::new();

    let num_threads = rayon::current_num_threads(); 
    let mut buffer = Vec::with_capacity(BUFFER_SIZE);
    

    let mut iterator = ExonIterator::new("main");
    let mut gex = SingleCellData::new(1);
    let mut last_chr = "unset".to_string();

    let mut genes = IndexedGenes::empty( Some(0) );

    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }

        if lines % split == 0 {
            pb.set_message(format!("{} mio reads processed", lines / split));
            pb.inc(1);
        }
        lines += 1;

        let data_tuple = match get_values( &record, &ref_id_to_name ){
            Ok(res) => res,
            Err(err) => {
                mapping_info.report( err );
                return;
            }
        };

        buffer.push(data_tuple); // Store only the tuple in the buffer

        if buffer.len() == BUFFER_SIZE { 
            mapping_info.stop_file_io_time();

            let chunk_size = (buffer.len() / num_threads).max(1);

            // Process the buffer in parallel
            let results = buffer
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut local_iterator = ExonIterator::new("part");// Initialize a thread-local iterator here if needed
                let mut local_collector = SingleCellData::new(1);
                let mut local_report = MappingInfo::new( None, 3.0, opts.max_reads ,None );
                let mut local_genes =  IndexedGenes::empty( Some(0) );

                chunk.iter().for_each(|(cell_id, umi, start, cigar, chr)| {
                    if last_chr != *chr {
                        println!("Init chr!");
                        let _ = gtf.init_search(chr, (*start).try_into().unwrap(), &mut local_iterator);
                    }

                    process_feature( 
                        cell_id,
                        umi,
                        *start,
                        cigar,
                        chr,
                        &gtf,
                        &mut local_iterator,
                        &mut local_collector,
                        &mut local_genes,
                        &mut local_report,
                        &ref_id_to_name,
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

            buffer.clear(); // Clear the buffer after processing
        }
    }

    if buffer.len() > num_threads { 
        mapping_info.stop_file_io_time();

        let chunk_size = (buffer.len() / num_threads).max(1);

        // Process the buffer in parallel
        let results = buffer
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local_iterator = ExonIterator::new("part");// Initialize a thread-local iterator here if needed
            let mut local_collector = SingleCellData::new(1);
            let mut local_report = MappingInfo::new( None, 3.0, opts.max_reads ,None );
            let mut local_genes =  IndexedGenes::empty( Some(0) );
            let mut last_chr = "unset";
            chunk.iter().for_each(|(cell_id, umi, start, cigar, chr)| {
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
                    &gtf,
                    &mut local_iterator,
                    &mut local_collector,
                    &mut local_genes,
                    &mut local_report,
                    &ref_id_to_name,
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
        let mut local_report = MappingInfo::new( None, 3.0, opts.max_reads ,None );
        let mut local_genes =  IndexedGenes::empty( Some(0) );

        buffer.iter().for_each(|(cell_id, umi, start, cigar, chr)| {
            if last_chr != *chr {
                println!("Init chr!");
                let _ = gtf.init_search(chr, (*start).try_into().unwrap(), &mut local_iterator);
            }

            process_feature( 
                cell_id,
                umi,
                *start,
                cigar,
                chr,
                &gtf,
                &mut local_iterator,
                &mut local_collector,
                &mut genes,
                &mut local_report,
                &ref_id_to_name,
            )

        });

        let translation = genes.merge( &local_genes );
        gex.merge_re_id_genes( local_collector, &translation );
        mapping_info.merge( &local_report );

    }

    let file_path_sp = PathBuf::from(&opts.outpath).join(
        "BD_Rhapsody_expression"
        );

    match gex.write_sparse_sub ( file_path_sp, &genes , &genes.get_all_gene_names(), opts.min_umi ) {
        Ok(_) => (),
        Err(err) => panic!("Error in the data write: {err}")
    };

    println!("{}",mapping_info.program_states_string());

    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            println!("\nI have analyzed the bam file in {milli}h {min}min {sec} sec {mil}milli sec\n");
        },
        Err(e) => {println!("Error: {e:?}");}
    }

    mapping_info.report_to_string();

    
}
