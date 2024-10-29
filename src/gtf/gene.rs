use crate::gtf::ExonIterator;


#[derive(Debug, Clone, Copy)]
pub struct Exon {
    start: usize,
    end: usize,
}



#[derive(Debug, Clone)]
pub struct Gene {
    pub gene_id: String,
    pub gene_name: String,
    pub start: usize,
    pub end: usize,
    pub sens_orientation: bool,
    exons: Vec<Exon>, // New field to store exons
}



#[derive(Debug, PartialEq)]
pub enum RegionStatus {
    InsideExon,
    InsideIntron,
    SpanningBoundary,
    BeforeGene,
    AfterGene,
}



impl Gene {
    // Constructor for Gene
    pub fn new(gene_id: &str, gene_name: &str, start: usize, end: usize, sens_orientation: bool) -> Self {
        Gene {
            start,
            end,
            sens_orientation,
            exons: Vec::new(), // Initialize the exons vector
            gene_id: gene_id.to_string(),
            gene_name: gene_name.to_string(),
        }
    }


    pub fn find_search_position(&self, pos: usize) -> Option<usize> {
        let mut low = 0;
        let mut high = self.exons.len();

        if pos < self.start {
            return Some(0)
        }
        
        // Perform a binary search for the first gene that starts after or overlaps with pos
        while low < high {
            let mid = (low + high) / 2;
            
            if self.exons[mid].start > pos {
                high = mid; // Narrow down to the left half
            } else {
                low = mid + 1; // Narrow down to the right half
            }
        }

        if low ==0 {
            return Some(low)
        }

        // At this point, `low` is the first index where the gene's start is >= pos or closest after `pos`
        if low > 0 && self.exons[low - 1].end >= pos {
            // If the gene just before the found position overlaps with pos, return it
            return Some(low - 1);
        }

        // Otherwise, check the gene at the found index (which should be the first gene starting after pos)
        if low < self.exons.len() && self.exons[low].end < pos {
            return Some(low);
        }

        // If no valid exon is found, return None
        None
    }

    // Method to add an exon to the gene
    pub fn add_exon(&mut self, start: usize, end: usize) {
        let exon = Exon { start, end };
        self.exons.push(exon); // Add the exon to the gene's exon list

        // Optionally, update the gene's start and end to reflect the new exon
        if start < self.start {
            self.start = start;
        }
        if end > self.end {
            self.end = end;
        }
    }

    pub fn check_region(&self, region_start: usize, region_end: usize, exon_id: &mut usize) -> RegionStatus {
        // Check if the region is completely before the gene
        if region_end < self.start {
            return RegionStatus::BeforeGene;
        }
        // Check if the region is completely after the gene
        if region_start > self.end {
            return RegionStatus::AfterGene;
        }

        // Check if exon_id is valid and if it points to a relevant exon
        if *exon_id < self.exons.len() && self.exons[*exon_id].end > region_start {
            // Iterate over exons starting from the current exon_id
            for i in *exon_id..self.exons.len() {
                // Check if the current exon overlaps with the region
                if self.exons[i].end > region_start {
                    *exon_id = i; // Update the exon_id to the current index
                    break; // Stop searching once we find a relevant exon
                }
            }
        }

        // Safety check: ensure that *exon_id is still within bounds
        if *exon_id >= self.exons.len() {
            return RegionStatus::AfterGene; // No relevant exons found
        }

        // Check if the region is entirely within an exon
        if region_start >= self.exons[*exon_id].start && region_end <= self.exons[*exon_id].end {
            return RegionStatus::InsideExon;
        }

        // Check if the region is entirely outside the exons (indicates it might be in an intron)
        if (region_start < self.exons[*exon_id].start && region_end < self.exons[*exon_id].start) ||
           (region_start > self.exons[*exon_id].end && region_end > self.exons[*exon_id].end) {
            return RegionStatus::InsideIntron;
        }

        // Check if the region spans the exon-intron boundary
        if (region_start < self.exons[*exon_id].start && region_end > self.exons[*exon_id].start) ||
           (region_start < self.exons[*exon_id].end && region_end > self.exons[*exon_id].end) {
            return RegionStatus::SpanningBoundary;
        }

        // If none of the above conditions matched, something went wrong.
        unreachable!("The read is neither before, after, or in the gene - that is impossible!");
    }
}