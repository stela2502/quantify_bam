use rustody::singlecelldata::SingleCellData;
use rustody::singlecelldata::cell_data::GeneUmiHash;
use rustody::int_to_str::IntToStr;
use rustody::singlecelldata::IndexedGenes;

//use rustody::cellids::CellIds;
//use rustody::cellids10x::CellIds10x;
//use rustody::traits::CellIndex;
use rustody::mapping_info::MappingInfo;

use quantify_bam::gtf::{GTF, ExonIterator, RegionStatus };

use bam::record::tags::StringType;

extern crate bam;
use bam::record::Record;
//use bam::record::tags::TagViewer;
use bam::record::tags::TagValue;

use crate::bam::RecordReader;

//use rustody::ofiles::{Ofiles, Fspot};

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};

use std::path::PathBuf;
use std::fs;
use std::fs::File;
//use std::path::Path;
//use std::io::Write;

//use std::thread;
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
    num_proc: Option<usize>,
    /// used processor cores (default all)
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


///fn get_values( bam_feature: &bam::Record, chromosmome_mappings:&HashMap<i32, String> ) 
///        -> Result<( String, String, u32, String, String, bool ), &str> {
/// with the return values being 
/// ( CellID, UMI, start+1, Cigar, chromosome, flag().is_reverse_strand() )
/// get_values is very specififc for this usage here - read bam reads and compare them to a gtf file
/// Bam is 0-based start position and end exclusive whereas gtf is 1-based and end inclusive.
/// This means that the end is actually the same value - just the start for GTF needs to be bam start +1
/// That is what this function returns for a start!
fn get_values<'a>( bam_feature: &'a bam::Record, chromosmome_mappings:&'a HashMap<i32, String> ) 
         -> Result<( String, String, i32, String, String, bool ), &'a str> {

    // Extract the chromosome (reference name)
    let chr = match chromosmome_mappings.get( &bam_feature.ref_id() ){
        Some(name) => name.to_string(),
        None => {
            return Err( "missing_Chromosome");  // Report missing chromosome
        }
    };

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

    // Extract the start and end positions
    let start = bam_feature.start();  // BAM is 0-based, start is inclusive
    // crap - this needs to be computed from the CIGAR!!!!
    //let mut output = Vec::new();
    //bam_feature.write_sam(&mut output).expect("Failed to write SAM");
    /*println!("This bam record {:?} was parsed into that tupel {:?}", 
        bam_feature, 
        ( &cell_id, &umi, start+1, &bam_feature.cigar().to_string(), &chr, &bam_feature.flag().is_reverse_strand() ) 
    );*/

    Ok( ( cell_id, umi, start+1, bam_feature.cigar().to_string(), chr,  bam_feature.flag().is_reverse_strand() ) ) // convert bam to gtf notation
}


fn process_feature(
    cell_id: &str, 
    umi: &str, 
    start: i32, 
    cigar: &str, 
    chr: &str,
    is_reverse_strand: &bool,
    gtf: &GTF, 
    iterator: &mut ExonIterator, 
    gex: &mut SingleCellData,
    genes: &mut IndexedGenes,
    mapping_info: &mut MappingInfo,  // Now tracks errors using a HashMap
    //chromosmome_mappings: &HashMap<i32, String>
)
{

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

    /*

    let read_result = match gtf.match_cigar_to_gene(&chr, &cigar, start.try_into().unwrap(), iterator) {
        Some(result) => result,
        None => {
            mapping_info.report("missing_Gene");
            return;
        },
    };
    //
    let orientation_mismatch = read_result.sens_orientation == *is_reverse_strand;

    if read_result.sens_orientation == *is_reverse_strand
 //       && matches!(read_result.match_type, RegionStatus::ExtTag | RegionStatus::InsideExon)
    {
        mapping_info.report("Orientation mismatch");
        return;
    }

    let gene_suffix = match read_result.match_type {
        //RegionStatus::InsideExon if read_result.sens_orientation == *is_reverse_strand => "_antisense",
        RegionStatus::InsideExon => "",
        RegionStatus::SpanningBoundary => "_unspliced",
        RegionStatus::InsideIntron => "_unspliced",
        RegionStatus::ExtTag => "_ext",
        _ => {
            mapping_info.report("missing_Gene");
            return;
        },
    };

    let gene_id = format!("{}{}", read_result.gene, gene_suffix);

    #[cfg(debug_assertions)]
    println!("Chr {chr}, cigar {cigar}, start {start} : CellID: {cell_id}");

    let umi_u64 = IntToStr::new(umi.into(), 32).into_u64();
    let cell_id_u64: u64 = cell_id
        .parse::<u64>()
        .unwrap_or_else(|_| IntToStr::new(cell_id.into(), 32).into_u64());

    if !gex.try_insert(
        &cell_id_u64,
        GeneUmiHash(genes.get_gene_id(&gene_id), umi_u64),
        mapping_info,
    ) {
        mapping_info.report("UMI_duplicate");
    }
}
*/
    let gene_id = match gtf.match_cigar_to_gene(&chr, &cigar, start.try_into().unwrap() , iterator) {
        
        Some(read_result) => {
            /*
            A00681:1014:HWGVKDMXY:2:2146:22064:36855    0   chr6    70703449    255 71M *   0   0   TCCCTGCATCCAGTGAGCAGTTAACATCTGGAGGTGCCTCAGTCGTGTGCTTCTTGAACAACTTCTACCCC FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFF:FFFFFFF:FFFFFFFFF NH:i:1  HI:i:1  AS:i:65 nM:i:2  XF:Z:Igkc   TR:Z:*  TF:Z:*  CB:Z:39780979   MR:Z:AATCGAAC   CN:Z:T  ST:Z:03 UB:Z:AATCGAAC
            Checking reads orientation: false vs gtf's orientation (id true): Igkc

            Checking reads orientation: true vs gtf's orientation (id false): Ighg1
            A00681:1014:HWGVKDMXY:2:2144:12608:13495    16  chr12   113294085   255 60M11S  *   0   0   GTTAGTTTGGGCAGCAGATCCAGGGGCCAGTGGATAGACAGATGGGGGTGTCGTTTTGGCAGAGGCGACGG FF:FF:FF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF:FFF,FFFFFFFFFFFFFFFFF NH:i:1  HI:i:1  AS:i:59 nM:i:0  XF:Z:Ighg1  TR:Z:*  TF:Z:*  CB:Z:42250405   MR:Z:CTTGGATT   CN:Z:T  ST:Z:04 UB:Z:CTTGGATT

            */
            //println!("Checking reads orientation: {} vs gtf's orientation (id {}): {}", is_reverse_strand , read_result.sens_orientation, read_result.gene);
            if read_result.sens_orientation == *is_reverse_strand { // antisense read!
                //println!("Is antisense!");
                if matches!(read_result.match_type,  RegionStatus::InsideIntron) {
                    mapping_info.report(&format!("{:?} Orientation mismatch", read_result.match_type) );
                    return;
                }
            }
            let add = if read_result.sens_orientation == *is_reverse_strand{
                "_antisense"
            }else {
                ""
            };

            match read_result.match_type {
                RegionStatus::InsideExon => read_result.gene + add,
                RegionStatus::SpanningBoundary => read_result.gene+"_unspliced" +add,
                RegionStatus::InsideIntron => read_result.gene+"_unspliced" +add,
                RegionStatus::ExtTag => read_result.gene+"_ext",
                _ => {
                    mapping_info.report("missing_Gene");
                    return;
                },
            }
            //format!("{}{}", read_result.gene, gene_suffix)
            /*
            // takes too long to actually collect all this info.
            // is not necessary after the debuging.
            match read_result.match_type{
                
                RegionStatus::InsideExon => {
                    if read_result.sens_orientation != *is_reverse_strand {
                        read_result.gene.to_string() + "_antisense"
                    }
                    else {
                        mapping_info.report("InsideExon"); 
                        read_result.gene.to_string()
                    }
                },
                RegionStatus::SpanningBoundary => {
                    mapping_info.report("SpanningBoundary"); 
                    read_result.gene.to_string() + "_unspliced"
                },
                RegionStatus::InsideIntron => {
                    mapping_info.report("InsideIntron"); 
                    read_result.gene.to_string() + "_unspliced"
                },
                RegionStatus::ExtTag => {
                    mapping_info.report("ExtTag read");
                    read_result.gene.to_string() + "_ext"
                },
                _ => {
                    mapping_info.report("missing_Gene");
                    return;
                }
            }*/
            
        },  // Get the gene_id if found
        None => {
            mapping_info.report("missing_Gene");  // Report missing gene match
            return;
        },

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
        Err(_e) => IntToStr::new( cell_id.into(), 32).into_u64(),
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

    let mut mapping_info = MappingInfo::new( Some(log_file), 3.0, opts.max_reads ,None );
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

        let data_tuple = match get_values( &record, &ref_id_to_name ){
            Ok(res) => res,
            Err( "missing_Chromosome" ) => {
                // We reached the end of the mapped reads - exit the loop
                break;
            }
            Err(err) => {
                mapping_info.report( err );
                return;
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
                let mut local_report = MappingInfo::new( None, 3.0, opts.max_reads ,None );
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
            let mut local_report = MappingInfo::new( None, 3.0, opts.max_reads ,None );
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
        let mut local_report = MappingInfo::new( None, 3.0, opts.max_reads ,None );
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
