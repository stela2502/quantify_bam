
use std::collections::HashMap;

use rayon::slice::ParallelSlice;
use rayon::iter::ParallelIterator;

use indicatif::MultiProgress;
use indicatif::ProgressStyle;
use indicatif::ProgressBar;

use crate::gtf::GTF;
use crate::gtf::ExonIterator;
use crate::gtf::RegionStatus;

use rustody::singlecelldata::SingleCellData;
use rustody::singlecelldata::IndexedGenes;
use rustody::mapping_info::MappingInfo;
use rustody::int_to_str::IntToStr;
use rustody::singlecelldata::cell_data::GeneUmiHash;

use crate::bam_helper::{get_values, get_values_bulk};
use crate::bam_helper::DataTuple;


//extern crate bam;
//use bam::record::tags::TagViewer;
use bam::RecordReader;
use bam::Header;

const BUFFER_SIZE: usize = 1_000_000;



fn process_feature(
    cell_id: &str, 
    umi: &u64, 
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

    let gene_id_u64 = genes.get_gene_id( &gene_id );

    let guh = GeneUmiHash(gene_id_u64, *umi );

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
)  -> Result<( SingleCellData, IndexedGenes ), String> {


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
    let mut gex = SingleCellData::new(1);
    let mut genes = IndexedGenes::empty(Some(0));

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
            process_buffer(&buffer, num_threads, &mut gex, &mut genes, mapping_info, gtf);
            pb.set_message(format!("{} mio reads - processing finished", lines / 1_000_000));
            buffer.clear();
            break;
        }
    }

    if !buffer.is_empty() {
        pb.set_message(format!("{} mio reads - processing", lines / 1_000_000));
        process_buffer(&buffer, num_threads, &mut gex, &mut genes, mapping_info, gtf);
        pb.set_message(format!("{} mio reads - processing finished", lines / 1_000_000));
    }
    if gex.len() == 0 {
        // after running millions of bam features we have not found a single cell?!
        return Err( format!("After analyzing one batch I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n" ))
    }
    return Ok( (gex, genes) )
}

// Function to process data
pub fn process_data_bulk(
    bam_file: &str,
    mapping_info: &mut MappingInfo,
    gtf: &GTF,
    cell_tag: [u8; 2],
    umi_tag: [u8; 2],
    num_threads: usize,
)  -> Result<( SingleCellData, IndexedGenes ), String> {


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
    let mut gex = SingleCellData::new(1);
    let mut genes = IndexedGenes::empty(Some(0));

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
            process_buffer(&buffer, num_threads, &mut gex, &mut genes, mapping_info, gtf );
            pb.set_message(format!("{} mio reads - processing finished", lines / 1_000_000));
            if gex.len() == 0 {
                // after running millions of bam features we have not found a single cell?!
                return Err( format!("After analyzing {} reads I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n",buffer.len()) )
            }
            buffer.clear();
        }
    }

    if !buffer.is_empty() {
        pb.set_message(format!("{} mio reads - processing", lines / 1_000_000));
        process_buffer(&buffer, num_threads, &mut gex, &mut genes, mapping_info, gtf);
        pb.set_message(format!("{} mio reads - processing finished", lines / 1_000_000));
        if gex.len() == 0 {
            // after running millions of bam features we have not found a single cell?!
            return Err( format!("After analyzing {} reads I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n",buffer.len() ))
        }
    }
    if gex.len() == 0 {
        // after running millions of bam features we have not found a single cell?!
        return Err( format!("After analyzing the first batch I could not detect a single cell/gene expression value?!\nSorry, but the --cell-tag and/or the --umi-tag might not be correct for this data.\n" ))
    }
    return Ok(( gex, genes ))
}


// Function to process a buffer
fn process_buffer(
    buffer: &[DataTuple],
    num_threads: usize,
    gex: &mut SingleCellData,
    genes: &mut IndexedGenes,
    mapping_info: &mut MappingInfo,
    gtf: &GTF,
) {
    let chunk_size = (buffer.len() / num_threads).max(1);

    let results: Vec<_> = buffer
        .par_chunks(chunk_size)
        .map(| chunk| process_chunk(chunk, gtf ))
        .collect();

    for (local_collector, local_report, local_genes) in results {
        let translation = genes.merge(&local_genes);
        gex.merge_re_id_genes(local_collector, &translation);
        mapping_info.merge(&local_report);
    }
}

// Function to process a chunk
fn process_chunk(chunk: &[DataTuple], gtf: &GTF ) -> (SingleCellData, MappingInfo, IndexedGenes) {
    let mut local_iterator = ExonIterator::new("part");
    let mut local_collector = SingleCellData::new(1);
    let mut local_report = MappingInfo::new(None, 3.0, 0, None);
    let mut local_genes = IndexedGenes::empty(Some(0));
    let mut last_chr = "unset";

    for (cell_id, umi, start, cigar, chr, is_reverse_strand) in chunk {
        if last_chr != *chr {
            match gtf.init_search(chr, (*start).try_into().unwrap(), &mut local_iterator){
                Ok(_) => {},
                Err(e) => {
                    local_report.report( &format!("{:?}", e) );
                    continue;
                    //panic!("Does the GTF match to the bam file?! ({:?})",e)
                },
            };
            last_chr = chr;
        }

        process_feature(
            cell_id,
            umi,
            *start,
            cigar,
            chr,
            is_reverse_strand,
            gtf,
            &mut local_iterator,
            &mut local_collector,
            &mut local_genes,
            &mut local_report,
        );
    }

    (local_collector, local_report, local_genes)
}
