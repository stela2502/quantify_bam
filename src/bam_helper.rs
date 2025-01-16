
use bam::record::tags::{StringType, TagValue};
use bam::record::Record;

use std::collections::HashMap;

use rustody::int_to_str::IntToStr;

/// The DataTuple stores the data from the Bam files:
/// (cell_id, umi, start, cigar, chr, is_reverse_strand, sequence as u8's, quality as u8's  )
pub type DataTuple = (String, u64, i32, String, String, bool, Vec<u8>, Vec<u8> );

/// get_tag will only return string tags - integer or floats nbeed to be implemented when needed.
pub fn get_tag(bam_feature: &Record, tag: &[u8;2]) -> Option<String> {
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
pub fn get_values<'a>( bam_feature: &'a bam::Record, chromosmome_mappings:&'a HashMap<i32, String>, 
    bam_cell_tag: &[u8;2], bam_umi_tag:&[u8;2], ) 
         -> Result< DataTuple , &'a str> {

    // Extract the chromosome (reference name)
    let chr = match chromosmome_mappings.get( &bam_feature.ref_id() ){
        Some(name) => name.to_string(),
        None => {
            return Err( "missing_Chromosome");  // Report missing chromosome
        }
    };

    let cell_id = match get_tag( bam_feature, bam_cell_tag ) {
        Some(id) => id.to_string(), // Convert to string, or use empty string
        None => {
            return Err( "missing_CellID" )
            //mapping_info.report("missing_CellID");  // Report missing Cell ID
            //return;
        }
    };

    // Try to extract the UMI (UB), and report if missing
    let umi = match get_tag( bam_feature, bam_umi_tag ) {
        Some(u) => u.to_string(), // Convert to string, or use empty string
        None => {
            return Err( "missing_UMI");  // Report missing UMI
        }
    };
    let umi_u64 = IntToStr::new( umi.into(), 32 ).into_u64();

    // Extract the start and end positions
    let start = bam_feature.start();  // BAM is 0-based, start is inclusive
    // crap - this needs to be computed from the CIGAR!!!!
    //let mut output = Vec::new();
    //bam_feature.write_sam(&mut output).expect("Failed to write SAM");
    /*println!("This bam record {:?} was parsed into that tupel {:?}", 
        bam_feature, 
        ( &cell_id, &umi, start+1, &bam_feature.cigar().to_string(), &chr, &bam_feature.flag().is_reverse_strand() ) 
    );*/
    
    let res: DataTuple =  (
        cell_id, 
        umi_u64,  
        start+1, 
        bam_feature.cigar().to_string(), 
        chr,  
        bam_feature.flag().is_reverse_strand(),
        bam_feature.sequence().raw().to_vec(),
        bam_feature.qualities().raw().to_vec()
        );  // convert bam to gtf notation
    Ok( res )
}

///fn get_values( bam_feature: &bam::Record, chromosmome_mappings:&HashMap<i32, String> ) 
///        -> Result<( String, String, u32, String, String, bool ), &str> {
/// with the return values being 
/// ( CellID, UMI, start+1, Cigar, chromosome, flag().is_reverse_strand() )
/// get_values is very specififc for this usage here - read bam reads and compare them to a gtf file
/// Bam is 0-based start position and end exclusive whereas gtf is 1-based and end inclusive.
/// This means that the end is actually the same value - just the start for GTF needs to be bam start +1
/// That is what this function returns for a start!
pub fn get_values_bulk<'a>( bam_feature: &'a bam::Record, chromosmome_mappings:&'a HashMap<i32, String>,
         _bam_cell_tag: &[u8;2], bam_umi_tag:&[u8;2] , pseudo_umi:u64,) 
         -> Result< DataTuple , &'a str> {

    // Extract the chromosome (reference name)
    let chr = match chromosmome_mappings.get( &bam_feature.ref_id() ){
        Some(name) => name.to_string(),
        None => {
            return Err( "missing_Chromosome");  // Report missing chromosome
        }
    };

    let cell_id = "1".to_string();

    // Try to extract the UMI (UB), and report if missing
    let umi_u64 = if bam_umi_tag == b"No" {
        pseudo_umi
    }else {
        let umi = match get_tag( bam_feature, bam_umi_tag ) {
            Some(u) => u.to_string(), // Convert to string, or use empty string
            None => {
                return Err( "missing_UMI");  // Report missing UMI
            }
        };
        IntToStr::new( umi.into(), 32 ).into_u64()
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

    let res: DataTuple =  (
        cell_id, 
        umi_u64,  
        start+1, 
        bam_feature.cigar().to_string(), 
        chr,  
        bam_feature.flag().is_reverse_strand(),
        bam_feature.sequence().raw().to_vec(),
        bam_feature.qualities().raw().to_vec()
        );  // convert bam to gtf notation
    Ok( res )
}