use std::collections::HashMap;
use std::cmp::Ordering;
use std::error::Error;

use std::path::Path;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

use std::fmt;

use crate::gtf::{gene::Gene, RegionStatus,ExonIterator};


#[derive(Debug, PartialEq)]
pub enum QueryErrors {
    OutOfGenes,
    ChrNotFound,
    NoMatch,
}


#[derive(Debug)]
pub struct GTF {
    pub chromosomes: HashMap<String, Vec<Gene>>, // Group genes and exons by chromosome
}

// Implement the Display trait for Gtf
impl fmt::Display for GTF {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let num_chromosomes = self.chromosomes.len();
        writeln!(f, "Number of chromosomes: {}", num_chromosomes)?;

        for (chromosome, gene_exons) in &self.chromosomes {
            let num_genes = gene_exons.len();
            writeln!(f, "Chromosome {}: {} genes", chromosome, num_genes)?;
        }

        Ok(())
    }
}


impl GTF {
    pub fn new() -> Self {
        GTF {
            chromosomes: HashMap::new(),
        }
    }

    pub fn add_exon(&mut self, gene_id: &str, gene_name: &str, start: usize, end: usize, chromosome: String, sens_orientation: bool) {
        // Get the vector of genes for the specified chromosome
        let chromosome_genes = self.chromosomes.entry(chromosome.clone()).or_insert(Vec::new());

        // Check if the gene already exists in the vector
        if let Some(gene) = chromosome_genes.iter_mut().find(|g| g.gene_id == gene_id) {
            // If it exists, add the exon to the existing gene
            gene.add_exon(start, end);
        } else {
            // If the gene does not exist, create a new gene and add the exon
            let mut new_gene = Gene::new(gene_id, gene_name, start, end, sens_orientation);
            new_gene.add_exon(start, end);
            chromosome_genes.push(new_gene);
        }
    }

    fn sort_genes(&mut self) {
        for gene_exons in self.chromosomes.values_mut() {
            // Sort by gene start position, then exon start position
            gene_exons.sort_by(|a, b | {
                let gene_cmp = a.start.cmp(&b.start);
                if gene_cmp == Ordering::Equal {
                    a.gene_id.cmp(&b.gene_id) // Sort by gene ID to maintain consistency
                } else {
                    gene_cmp
                }
            });
        }
    }


    pub fn slice_gtf(&self, chromosome: &str, start: usize, end: usize) -> Result<Vec<Gene>, &str> {
        // Check if the chromosome exists in the gtf
        if let Some(genes) = self.chromosomes.get(chromosome) {
            // Filter genes based on their start and end positions
            Ok( genes.iter()
                .filter(|gene| gene.end > start && gene.start < end) // Ensure overlapping genes
                .cloned() // Clone each gene for the resulting vector
                .collect() // Collect results into a vector
            )
        } else {
            // If the chromosome does not exist, return an empty vector
            println!("Sorry chromosome {} is not present here", chromosome);
            Err("Sorry this chromosome is not present here")
        }
    }

    // the main query functionality
    pub fn query(&self, chr:&str, start:usize, end:usize, 
            iterator: &mut ExonIterator ) -> Result<( ( String, String ) , RegionStatus ), QueryErrors > {
        match self.chromosomes.get(chr) {
            Some(genes) => {
                // we had the chromosome - assuming the ExonIterator is set - let's check that
                match genes.get(iterator.gene_id()){
                    Some(gene) => {
                        match gene.check_region( start, end, &mut iterator.exon_id() ){
                            RegionStatus::AfterGene => {
                                // that needs to be run again...
                                iterator.set_exon_id(0);
                                iterator.set_gene_id( iterator.gene_id() +1 );
                                return self.query( chr, start, end, iterator )
                            },
                            other => {
                                let this_gene = &genes[iterator.gene_id()];
                                Ok( ( ( this_gene.gene_id.to_string(), this_gene.gene_name.to_string() ), other ) )
                            }
                        }
                    },
                    None => {
                        Err( QueryErrors::OutOfGenes )
                    }
                } 
            },
            None => {
                Err( QueryErrors::ChrNotFound )
            }
        } 
    }

    fn find_search_position(&self, genes: &Vec<Gene>, pos: usize) -> Option<usize> {
        let mut low = 0;
        let mut high = genes.len();
        
        if genes[0].start < pos {
            return Some(0)
        }
        if genes[genes.len()-1].end < pos {
            return None
        }
        // Perform a binary search for the first gene that starts after or overlaps with pos
        while low < high {
            let mid = (low + high) / 2;
            
            if genes[mid].start > pos {
                high = mid; // Narrow down to the left half
            } else {
                low = mid + 1; // Narrow down to the right half
            }
        }

        // At this point, `low` is the first index where the gene's start is >= pos or closest after `pos`
        if low > 0 && genes[low - 1].end >= pos {
            // If the gene just before the found position overlaps with pos, return it
            return Some(low - 1);
        }

        // Otherwise, check the gene at the found index (which should be the first gene starting after pos)
        if low < genes.len() && genes[low].end < pos {
            return Some(low) ;
        }

        // If no valid gene is found, return None
        None
    }

    // init the search - assuming we are iterating over a certain slice of a bam file
    pub fn init_search( &self, chr: &str, start:usize, iterator: &mut ExonIterator )-> Result<(), QueryErrors>{
        match self.chromosomes.get(chr) {
            Some(genes) => {
                match self.find_search_position(genes, start){
                    Some( gene_id ) => {
                        if let Some(exon_id) = genes[gene_id].find_search_position( start ){
                            iterator.set_gene_id( gene_id );
                            iterator.set_exon_id( exon_id );
                            Ok(())
                        }else {
                            iterator.set_gene_id( gene_id-1 );
                            iterator.set_exon_id( 0 );
                            Ok(())
                        }
                    },
                    None => {
                        Err( QueryErrors::OutOfGenes )
                    }
                }
            },
            None => {
                Err( QueryErrors::ChrNotFound )
            }
        }
    }

    // Function to parse the GTF file and populate the Gtf structure
    pub fn parse_gtf(&mut self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let path = Path::new(file_path);
        let file = File::open(&path)?;
        let reader = BufReader::new(file);

        // Helper function to extract specific attributes like gene_id or gene_name
        fn extract_attribute<'a>(attributes: &'a str, key: &'a str) -> Option< &'a str > {
            attributes
                .split(';')
                .find_map(|attr| {
                    let mut parts = attr.trim().split_whitespace();
                    if let Some(attr_key) = parts.next() {
                        if attr_key == key {
                            return parts.next(); // Get the value after the key
                        }
                    }
                    None
                })
        }

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue; // Skip comment lines
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                continue; // Skip malformed lines
            }

            let chromosome = fields[0].to_string();
            let feature_type = fields[2];
            let orientation = fields[7] == "+";

            if feature_type == "exon" {
                let start: usize = fields[3].parse()?;
                let end: usize = fields[4].parse()?;

                // Extract gene_id and gene_name from the attributes
                let attributes = fields[8];
                let gene_id = extract_attribute(attributes, "gene_id").unwrap_or("unknown");
                let gene_name = extract_attribute(attributes, "gene_name").unwrap_or("unknown");

                // Call your add_exon function (you may want to add gene_name handling there)
                self.add_exon(gene_id, gene_name, start, end, chromosome, orientation);
            }
        }

        println!("I have read this:\n{}", self);
        Ok(())
    }
    
}




