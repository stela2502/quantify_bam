use rustody::analysis::bam_flag::BamFlag;
use rust_htslib::bam::record::{Record};
use rust_htslib::errors::Result;
use std::collections::HashMap;
use rustody::int_to_str::IntToStr;


#[derive(Debug, Clone)]
pub struct ReadData {
    pub cell_id: String,
    pub umi: u64,
    pub start: i32,
    pub flag: BamFlag,
    pub cigar: String,
    pub chromosome: String,
    //pub is_reverse: bool,
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
}

impl ReadData {
    /// Extracts all required info from `bam_feature` and makes it thread-safe.
    pub fn new(
        cell_id: &str,
        umi: u64,
        start: i32,
        flag: u16,
        cigar: &str,
        chromosome: &str,
        //is_reverse: bool,
        sequence: Vec<u8>,
        qualities: Vec<u8>,
    ) -> Self {
        Self {
            cell_id: cell_id.to_string(),
            umi,
            start: start + 1, // Convert BAM to GTF notation
            flag: BamFlag::new( flag ),
            cigar: cigar.to_string(),
            chromosome: chromosome.to_string(),
            //is_reverse,
            sequence,
            qualities,
        }
    }


    // Replace bam::record::tags::{StringType, TagValue}; with rust_htslib's tag handling

    /// The DataTuple stores the data from the Bam files:
    /// ( 0,      1,   2,     3,     4,   5,                 6,                7,               )
    /// (cell_id, umi, start, cigar, chr, is_reverse_strand, sequence as u8's, quality as u8's  )
    /// In this case, we are using ReadData in the return values.
    fn get_tag_value(record: &Record, tag: &[u8; 2]) -> Option<String> {
        if let Some(value) = record.aux(tag).ok() {
            match value {
                rust_htslib::bam::record::Aux::String(s) => Some(s.to_string()),
                _ => None,
            }
        } else {
            None
        }
    }

    /// get_values will extract the required fields (CellID, UMI, etc.) from a BAM record.
    pub fn from_single_cell<'a>(
        bam_feature: &'a Record,
        chromosmome_mappings: &'a HashMap<i32, String>,
        bam_cell_tag: &[u8; 2],
        bam_umi_tag: &[u8; 2],
    ) -> Result<(String, Self), &'a str> {
        
        let cell_id = match Self::get_tag_value(bam_feature, bam_cell_tag) {
            Some(id) => id, // Convert to string, or use empty string
            None => return Err("missing_CellID"),
        };

        Self::from_bulk(bam_feature, chromosmome_mappings, bam_umi_tag, 1_u64, &cell_id)
    }

    /// get_values_bulk extracts the values and creates a ReadData object.
    pub fn from_bulk<'a>(
        bam_feature: &'a Record,
        chromosmome_mappings: &'a HashMap<i32, String>,
        bam_umi_tag: &[u8; 2],
        pseudo_umi: u64,
        cell_id: &str,
    ) -> Result<(String, Self), &'a str> {
        // Extract the chromosome (reference name)
        let chr = match chromosmome_mappings.get(&bam_feature.tid()) {
            Some(name) => name,
            None => return Err("missing_Chromosome"),
        };

        // Try to extract the UMI (UB), and report if missing
        let umi_u64 = if bam_umi_tag == b"No" {
            pseudo_umi
        } else {
            let umi = match Self::get_tag_value(bam_feature, bam_umi_tag) {
                Some(u) => u,
                None => return Err("missing_UMI"),
            };
            IntToStr::new(umi.into(), 32).unwrap().into_u64()
        };

        // Extract the start position (1-based for your GTF comparison)
        let start = bam_feature.pos(); // This is 0-based, needs adjustment for GTF

        // Create a new ReadData object
        let res = Self::new(
            cell_id,
            umi_u64,
            (start + 1).try_into().unwrap(), // Adjust for GTF 1-based start
            bam_feature.flags(),
            &bam_feature.cigar().to_string(),
            chr,
            bam_feature.seq().as_bytes(),
            bam_feature.qual().to_vec(),
        );

        let id = std::str::from_utf8(bam_feature.qname()).unwrap().to_string(); // Convert read name to string

        Ok((id, res))
    }


    pub fn is(&self, what:&str ) -> bool {
        self.flag.is( what )
    }
}