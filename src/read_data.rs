use rustody::analysis::bam_flag::BamFlag;

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

    pub fn is(&self, what:&str ) -> bool {
        self.flag.is( what )
    }
}