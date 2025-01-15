pub struct MutationProcessor {
    quality_cutoff: usize, // Mean quality score threshold
}

impl MutationProcessor {
    /// Function to generate a mutation description string for CIGAR operations (X, D, I).
    ///
    /// # Arguments
    /// - `mutation_type`: The type of mutation (X, D, I).
    /// - `start`: The starting position of the mutation.
    /// - `length`: The length of the mutation.
    /// - `seq`: The sequence of the read.
    /// - `qual`: The quality score string corresponding to the read sequence.
    ///
    /// # Returns
    /// A String describing the mutation with position, nucleotides, and quality check.
    fn mutation_name(
        &self,
        mutation_type: char,
        start: usize,
        length: usize,
        seq: &str,
        qual: &str,
    ) -> Option<String> {
        // Check the quality of nucleotides involved in the mutation
        let quality_check = self.check_quality(start, length, qual);

        // If the mutation quality is below the cutoff, skip it
        if quality_check < self.quality_cutoff {
            return None;
        }

        // Handle each mutation type
        match mutation_type {
            'X' => self.mismatch_name(start, length, seq),
            'D' => self.deletion_name(start, length),
            'I' => self.insertion_name(start, length, seq),
            _ => None, // If the mutation type is unhandled, return None
        }
    }

    /// Check the quality of the nucleotides involved in the mutation.
    fn check_quality(&self, start: usize, length: usize, qual: &str) -> usize {
        // Substring the quality scores involved in the mutation
        let mutation_qual = &qual[start..start + length];

        // Calculate the mean quality score for this substring
        let total_quality: usize = mutation_qual.chars().map(|c| c as usize - 33).sum();
        total_quality / mutation_qual.len()
    }

    /// Handle mismatches (X) and return mutation name with nucleotide differences.
    fn mismatch_name(&self, start: usize, length: usize, seq: &str) -> Option<String> {
        let mismatch_seq = &seq[start..start + length];
        Some(format!(
            "Mismatch at position {}: {}",
            start + 1, // Position is 1-based in mutation naming
            mismatch_seq
        ))
    }

    /// Handle deletions (D) and return mutation name with "NA" for deleted nucleotides.
    fn deletion_name(&self, start: usize, length: usize) -> Option<String> {
        Some(format!(
            "Deletion at position {}: NA",
            start + 1 // Position is 1-based in mutation naming
        ))
    }

    /// Handle insertions (I) and return mutation name with inserted nucleotides.
    fn insertion_name(&self, start: usize, length: usize, seq: &str) -> Option<String> {
        let inserted_seq = &seq[start..start + length];
        Some(format!(
            "Insertion at position {}: {}",
            start + 1, // Position is 1-based in mutation naming
            inserted_seq
        ))
    }
}