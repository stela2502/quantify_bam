
use std::collections::HashMap;

use rayon::slice::ParallelSlice;
use rayon::iter::ParallelIterator;

use indicatif::MultiProgress;
use indicatif::ProgressStyle;
use indicatif::ProgressBar;

use crate::gtf::{GTF, ExonIterator, RegionStatus };
use crate::mutation_processor::MutationProcessor;
use crate::bam_helper::{get_values, get_values_bulk, DataTuple };

use rustody::singlecelldata::{SingleCellData, IndexedGenes, cell_data::GeneUmiHash };
use rustody::mapping_info::MappingInfo;
use rustody::int_to_str::IntToStr;

use std::path::Path;

//extern crate bam;
//use bam::record::tags::TagViewer;
use bam::RecordReader;
use bam::Header;

const BUFFER_SIZE: usize = 1_000_000;

use std::env;
use lazy_static::lazy_static;

lazy_static! {
    pub static  ref PROGRAM_NAME: String = {
        if let Some(program_path) = env::args().next() {
            if let Some(program_name) = Path::new(&program_path).file_name() {
                program_name.to_string_lossy().to_string()
            } else {
                String::from("Unknown")
            }
        } else {
            String::from("Unknown")
        }
    };
}



// Function to create a mapping of reference ID to name
fn create_ref_id_to_name_map(header_view: &Header) -> HashMap<i32, String> {
    header_view
        .reference_names()
        .iter()
        .enumerate()
        .map(|(id, name)| (id as i32, name.to_string()))
        .collect()
}

// Function to process data
pub fn process_data(
    bam_file: &str,
    mapping_info: &mut MappingInfo,
    gtf: &GTF,
    cell_tag: [u8; 2],
    umi_tag: [u8; 2],
    num_threads: usize,
    mutations: &Option<MutationProcessor>,
)  -> Result<( (SingleCellData, IndexedGenes),(SingleCellData, IndexedGenes) ), String> {


    let mut reader = bam::BamReader::from_path( bam_file, 1).unwrap();
    let header_view = reader.header().to_owned();
    let ref_id_to_name = create_ref_id_to_name_map(&header_view);

    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    pb.set_style(ProgressStyle::default_bar().template("{prefix:.bold.dim} {spinner} {wide_msg}").unwrap());
    pb.set_message("");

    let mut lines = 0_u64;
    let mut record = bam::Record::new();
    let split = BUFFER_SIZE * num_threads;

    println!("Using {} processors and processinf {} reads a batch", num_threads, split);

    let mut buffer = Vec::with_capacity(split);
    let mut expr_gex = SingleCellData::new(1);
    let mut expr_idx = IndexedGenes::empty(Some(0));

    let mut mut_gex = SingleCellData::new(1);
    let mut mut_idx = IndexedGenes::empty(Some(0));

    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }

        if lines % 1_000_000 == 0 {
            pb.set_message(format!("{} mio reads processed", lines / 1_000_000));
            pb.inc(1);
        }
        lines += 1;

        let data_tuple = match get_values(&record, &ref_id_to_name, &cell_tag, &umi_tag) {
            Ok(res) => res,
            Err("missing_Chromosome") => {
                eprintln!("Missing chromosome for BAM entry - assuming end of usable data.\n{:?}", record);
                break;
            }
            Err(err) => {
                //println!("error {err:?}");
                mapping_info.report(err);
                continue;
            }
        };

        buffer.push(data_tuple);

        if buffer.len() >= split {
            pb.set_message(format!("{} mio reads - processing", lines / 1_000_000));
            process_buffer(
                &buffer,
                num_threads, 
                &mut expr_gex, 
                &mut expr_idx,
                &mut mut_gex, 
                &mut mut_idx,
                mapping_info, 
                gtf,
                mutations
            );
            pb.set_message(format!("{} mio reads - processing finished", lines / 1_000_000));
            buffer.clear();
            if expr_gex.len() == 0 {
                // after running millions of bam features we have not found a single cell?!
                return Err( format!("After analyzing {} reads I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n",buffer.len()) )
            }
            break;
        }
    }

    if !buffer.is_empty() {
        pb.set_message(format!("{} mio reads - processing", lines / 1_000_000));
        process_buffer(
            &buffer,
            num_threads, 
            &mut expr_gex, 
            &mut expr_idx,
            &mut mut_gex, 
            &mut mut_idx,
            mapping_info, 
            gtf,
            mutations
        );
        pb.set_message(format!("{} mio reads - processing finished", lines / 1_000_000));
        if expr_gex.len() == 0 {
            // after running millions of bam features we have not found a single cell?!
            return Err( format!("After analyzing {} reads I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n",buffer.len()) )
        }
    }
    if expr_gex.len() == 0 {
        // after running millions of bam features we have not found a single cell?!
        return Err( format!("After analyzing one batch I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n" ))
    }
    return Ok( (( expr_gex, expr_idx ), ( mut_gex, mut_idx )) )
}

// Function to process data
pub fn process_data_bulk(
    bam_file: &str,
    mapping_info: &mut MappingInfo,
    gtf: &GTF,
    cell_tag: [u8; 2],
    umi_tag: [u8; 2],
    num_threads: usize,
    mutations: &Option<MutationProcessor>,
)  -> Result<( (SingleCellData, IndexedGenes),(SingleCellData, IndexedGenes) ), String> {


    let mut reader = bam::BamReader::from_path( bam_file, 1).unwrap();
    let header_view = reader.header().to_owned();
    let ref_id_to_name = create_ref_id_to_name_map(&header_view);

    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    pb.set_style(ProgressStyle::default_bar().template("{prefix:.bold.dim} {spinner} {wide_msg}").unwrap());
    pb.set_message("");

    let mut lines = 0_u64;
    let mut record = bam::Record::new();
    let split = BUFFER_SIZE * num_threads;

    println!("Using {} processors and processinf {} reads a batch", num_threads, split);

    let mut buffer = Vec::with_capacity(split);
    let mut expr_gex = SingleCellData::new(1);
    let mut expr_idx = IndexedGenes::empty(Some(0));

    let mut mut_gex = SingleCellData::new(1);
    let mut mut_idx = IndexedGenes::empty(Some(0));

    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {},
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }

        if lines % 1_000_000 == 0 {
            pb.set_message(format!("{} mio reads processed", lines / 1_000_000));
            pb.inc(1);
        }
        lines += 1;
        let data_tuple = match get_values_bulk(&record, &ref_id_to_name, &cell_tag, &umi_tag, lines ) {
            Ok(res) => res,
            Err("missing_Chromosome") => {
                eprintln!("Missing chromosome for BAM entry - assuming end of usable data.\n{:?}", record);
                break;
            }
            Err(err) => {
                //println!("error {err:?}");
                mapping_info.report(err);
                continue;
            }
        };

        buffer.push(data_tuple);

        if buffer.len() >= split {
            pb.set_message(format!("{} mio reads - processing", lines / 1_000_000));
            process_buffer(
                &buffer,
                num_threads, 
                &mut expr_gex, 
                &mut expr_idx,
                &mut mut_gex, 
                &mut mut_idx,
                mapping_info, 
                gtf,
                mutations
            );
            pb.set_message(format!("{} mio reads - processing finished", lines / 1_000_000));
            if expr_gex.len() == 0 {
                // after running millions of bam features we have not found a single cell?!
                return Err( format!("After analyzing {} reads I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n",buffer.len()) )
            }
            buffer.clear();
        }
    }

    if !buffer.is_empty() {
        pb.set_message(format!("{} mio reads - processing", lines / 1_000_000));
        process_buffer(
            &buffer,
            num_threads, 
            &mut expr_gex, 
            &mut expr_idx,
            &mut mut_gex, 
            &mut mut_idx,
            mapping_info, 
            gtf,
            mutations
        );
        pb.set_message(format!("{} mio reads - processing finished", lines / 1_000_000));
        if expr_gex.len() == 0 {
            // after running millions of bam features we have not found a single cell?!
            return Err( format!("After analyzing {} reads I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n",buffer.len() ))
        }
    }
    if expr_gex.len() == 0 {
        // after running millions of bam features we have not found a single cell?!
        return Err( format!("After analyzing the first batch I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n" ))
    }
    return Ok( (( expr_gex, expr_idx ), ( mut_gex, mut_idx )) )
}


// Function to process a buffer
fn process_buffer(
    buffer: &[DataTuple],
    num_threads: usize,
    expr_gex: &mut SingleCellData,
    expt_idx: &mut IndexedGenes,
    mut_gex: &mut SingleCellData,
    mut_idx: &mut IndexedGenes,
    mapping_info: &mut MappingInfo,
    gtf: &GTF,
    mutations: &Option<MutationProcessor>,
) {
    // everything outside this function is file io
    mapping_info.stop_file_io_time();
    let chunk_size = (buffer.len() / num_threads).max(1);

    let results: Vec<_> = buffer
        .par_chunks(chunk_size)
        .map(| chunk| process_chunk(chunk, gtf, mutations ))
        .collect();

    // only the above is the multiprocessor step
    mapping_info.stop_multi_processor_time();

    for (gex_results, mut_results, local_report) in results {

        // collect expr data
        let expr_trans = expt_idx.merge( &gex_results.1 );
        expr_gex.merge_re_id_genes( gex_results.0, &expr_trans );
        
        // collect mutations
        let mut_trans = mut_idx.merge( &mut_results.1 );
        mut_gex.merge_re_id_genes( mut_results.0, &mut_trans);
        
        // fill in the report
        mapping_info.merge(&local_report);
    }
    
    // and that is the main single processor part
    mapping_info.stop_single_processor_time();

}

// Function to process a chunk
fn process_chunk(chunk: &[DataTuple], gtf: &GTF, mutations: &Option<MutationProcessor>, ) -> (( SingleCellData, IndexedGenes),  (SingleCellData,IndexedGenes), MappingInfo ) {
    let mut local_iterator = ExonIterator::new("part");
    // for the expression data
    let mut expr_gex = SingleCellData::new(1);
    let mut expr_idx = IndexedGenes::empty(Some(0)) ;
    // for the mutation counts
    let mut mut_gex = SingleCellData::new(1);
    let mut mut_idx = IndexedGenes::empty(Some(0)) ;

    let mut local_report = MappingInfo::new(None, 3.0, 0, None);
    let mut last_chr = "unset";
    //     0        1    2      3      4    5                  6    7 
    //for (data.0, umi, start, cigar, chr, is_reverse_strand, seq, qual ) in chunk {
    for data in chunk{
        if last_chr != &data.4 {
            match gtf.init_search( &data.4, (data.2).try_into().unwrap(), &mut local_iterator){
                Ok(_) => {},
                Err(e) => {
                    local_report.report( &format!("{:?}", e) );
                    continue;
                    //panic!("Does the GTF match to the bam file?! ({:?})",e)
                },
            };
            last_chr = &data.4 ;
        }

        process_feature(
            data,
            gtf,
            mutations,
            &mut local_iterator,
            &mut expr_gex,
            &mut expr_idx,
            &mut mut_gex,
            &mut mut_idx,
            &mut local_report,
        );
    }
    ( (expr_gex, expr_idx), (mut_gex, mut_idx), local_report )

}

fn process_feature(
    data: &DataTuple,
    gtf: &GTF, 
    mutations: &Option<MutationProcessor>,
    iterator: &mut ExonIterator,
    exp_gex: &mut SingleCellData,
    exp_idx: &mut IndexedGenes,
    mut_gex: &mut SingleCellData,
    mut_idx: &mut IndexedGenes,
    mapping_info: &mut MappingInfo, // Now tracks errors using a HashMap
    //chromosmome_mappings: &HashMap<i32, String>
)
{
        
    let gene_ids: Vec<String> = match gtf.match_cigar_to_gene(&data.4, &data.3, data.2.try_into().unwrap(), iterator) {
        Some(read_result) => {
            // Filter results with the opposite strand
            let good: Vec<String> = read_result
                .iter()
                .filter(|result| result.sens_orientation != data.5)
                .map(|result| match result.match_type {
                    RegionStatus::InsideExon => result.gene.clone(),
                    RegionStatus::SpanningBoundary | RegionStatus::InsideIntron => format!("{}_unspliced", result.gene),
                    RegionStatus::ExtTag => format!("{}_ext", result.gene),
                    _ => {
                        mapping_info.report("missing_Gene_Info");
                        "missing".to_string()
                    },
                })
                .collect();

            // Check the length of the filtered results
            if good.len() > 1 {
                good.into_iter().map(|name| format!("{}_ambiguous", name)).collect()
            } else if good.len() == 1 {
                good
            } else {
                // Filter results with the same strand (antisense)
                let anti: Vec<String> = read_result
                    .iter()
                    .filter(|result| result.sens_orientation == data.5)
                    .map(|result| match result.match_type {
                        RegionStatus::InsideExon => format!("{}_antisense", result.gene),
                        RegionStatus::SpanningBoundary | RegionStatus::InsideIntron => format!("{}_unspliced_antisense", result.gene),
                        RegionStatus::ExtTag => format!("{}_ext_antisense", result.gene),
                        _ => {
                            mapping_info.report("missing_Gene");
                            "missing_antisense".to_string()
                        },
                    })
                    .collect();

                if anti.len() > 1 {
                    anti.into_iter().map(|name| format!("{}_ambiguous", name)).collect()
                } else {
                    anti
                }
            }
        }
        None => {
            mapping_info.report("missing_Gene");
            return; // Exit early if there's no match
        }
    };
    

    // Generate a GeneUmiHash
    #[cfg(debug_assertions)]
    println!("Chr {}, cigar {} and this start position {} : CellID: {}", data.4, data.3, data.2, data.0);

    for gene_id in &gene_ids{
        let gene_id_u64 = exp_idx.get_gene_id( &gene_id );

        let guh = GeneUmiHash(gene_id_u64, data.1 );

        #[cfg(debug_assertions)]
        println!("\t And I got a gene: {guh}");

        // Try to insert into SingleCellData (gex)
        // the cellid is already a numeric!
        let cell_id =  match data.0.parse::<u64>() {
            Ok(number) => number,
            Err(_e) => IntToStr::new( data.0.as_bytes().to_vec(), 32).into_u64(),
        };

        if exp_gex.try_insert(
            &cell_id,
            guh,
            mapping_info
        ) {
            // Here I need to handle possible mutations, too
            
            match &mutations {
                //     0        1    2      3      4    5                  6    7 
                //for (data.0, umi, start, cigar, chr, is_reverse_strand, seq, qual ) in chunk {
                Some( processor ) => {
                    for name in processor.get_all_mutations( 
                        data.2.try_into().unwrap(),  
                        &data.3,
                        &data.6, 
                        &data.7,
                        mapping_info,
                    ){
                        let snip = gene_id.to_string() +"/" + &data.4 + "/" + &name;
                        let mut_id = mut_idx.get_gene_id( &snip );
                        let ghum = GeneUmiHash( mut_id, data.1 );

                        let _ = mut_gex.try_insert(
                            &cell_id,
                            ghum,
                            mapping_info
                        );
                    }
                },
                None => {
                    // ok no mutations - nothing to do here
                }
            }

        }else {
            mapping_info.report("UMI_duplicate");
            break // would be a UMI_duplicate for the other entry, too
        }
    };
    
}