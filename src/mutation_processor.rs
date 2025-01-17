use regex::Regex;
use rustody::mapping_info::MappingInfo;

pub struct MutationProcessor {
    pub quality_cutoff: usize, // Mean quality score threshold
}

impl MutationProcessor {

    /// the main entry into the mutations.
    /// Give me the bam entries start, the Cigar,  sequence and quality scores and I give you a vector detected mutations
    pub fn get_all_mutations(&self, bam_start:usize, cigar:&str, sequence:&[u8], qual:&[u8], mapping_info: &mut MappingInfo ) -> Vec<String>{
        let re_start = Regex::new(r"(\d+)([MIDNSHPX=Z])").unwrap();
        let mut current_position = bam_start; // Start from the provided position
        let mut current_loc_pos = 0;
        let mut ret = Vec::new();


        for cap in re_start.captures_iter(cigar) {
            let count: usize = cap[1].parse().unwrap(); // Mutation length
            let operation = cap[2].chars().next().unwrap(); // Mutation type (M, I, D, etc.)

            if operation != 'N' {
                if current_loc_pos + count > qual.len() {
                    panic!("MutationProcessor::get_all_mutations - cigar {cigar} lead to a out of range issue: {current_loc_pos} + {count} ({})> {}" ,current_loc_pos + count,qual.len() );
                }
                if let Some(name) = self.mutation_name(operation, current_position, current_loc_pos, count, sequence, qual, mapping_info) {
                    ret.push(name)
                }else {
                    mapping_info.report("not mutated");
                }
                if operation != 'D' {
                    current_position += count;
                    current_loc_pos += count;
                }
            }
        }

        ret
    }

    fn mutation_name(
        &self,
        mutation_type: char,
        start: usize,
        local_start: usize,
        length: usize,
        seq: &[u8], // Accept slice of u8 directly
        qual: &[u8], // Accept slice of u8 directly
        mapping_info: &mut MappingInfo,
    ) -> Option<String> {
        // Check the quality of nucleotides involved in the mutation
        if start+length > qual.len() {
            return None
        }
        let quality_check = self.check_quality(local_start, length, qual);

        // If the mutation quality is below the cutoff, skip it
        if quality_check < self.quality_cutoff {
            mapping_info.report(&format!("mutation quality {} - failed",quality_check ) );
            return None;
        }

        // Handle each mutation type
        match mutation_type {
            'X' => self.mismatch_name(start, local_start,length, seq),
            'D' => self.deletion_name(start, length),
            'I' => self.insertion_name(start, local_start, length, seq),
            _ => None, // If the mutation type is unhandled, return None
        }
    }

    fn check_quality(&self, start: usize, length: usize, qual: &[u8]) -> usize {
        // Substring the quality scores involved in the mutation
        
        let mutation_qual = &qual[start..start + length];

        // Calculate the mean quality score for this substring
        let total_quality: usize = mutation_qual.iter().map(|&c| (c as usize) - 33).sum();
        total_quality / mutation_qual.len()
    }

    fn mismatch_name(&self, start: usize, local_start:usize, length: usize, seq: &[u8]) -> Option<String> {
        let mismatch_seq = &seq[local_start..local_start + length];
        // Convert the mutated sequence slice into a String when needed
        Some(format!(
            "snp/{}/{}",
            start + 1, // Position is 1-based in mutation naming
            String::from_utf8_lossy(mismatch_seq) // Convert slice to String only for the mutated part
        ))
    }

    fn deletion_name(&self, start: usize, length: usize) -> Option<String> {
        Some(format!(
            "ins/{}/{}",
            start + 1, // Position is 1-based in mutation naming,
            length,
        ))
    }

    fn insertion_name(&self, start: usize, local_start:usize, length: usize, seq: &[u8]) -> Option<String> {
        let inserted_seq = &seq[local_start..local_start + length];
        // Convert the inserted sequence slice into a String when needed
        Some(format!(
            "del/{}/{}",
            start + 1, // Position is 1-based in mutation naming
            String::from_utf8_lossy(inserted_seq), // Convert slice to String only for the mutated part
        ))
    }
}
