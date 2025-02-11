use std::collections::HashMap;
//use std::cmp::Ordering;
use std::error::Error;

use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufReader, BufRead, BufWriter, Write, Read};
use flate2::write::GzEncoder;
use flate2::read::GzDecoder;
use flate2::Compression;

use std::fmt;

use crate::gtf::{gene::Gene, gene::RegionStatus,ExonIterator,SplicedRead};

use crate::gtf::exon_iterator::ReadResult;

#[derive(Debug, PartialEq)]
pub enum QueryErrors {
    OutOfGenes,
    ChrNotFound,
    NoMatch,
}


#[derive(Debug)]
pub struct GTF {
    pub chromosomes: HashMap<String, Vec<Gene>>, // Group genes and exons by chromosome
    exon_file: Option<BufWriter<GzEncoder<File>>>, // Buffered gzip writer


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
    pub fn new( exon_file_path:Option<PathBuf> ) -> Self {

        let exon_file = exon_file_path.map(|path| {
            let file = File::create(path).expect("Unable to create exon file");
            let encoder = GzEncoder::new(file, Compression::default());
            let mut writer = BufWriter::new(encoder);
            let line = format!("transcript_id\tgene_name\tchr\tstart\tend\n");
            writer.write_all(line.as_bytes()).expect("Failed to write to file");
            writer
        });

        GTF {
            chromosomes: HashMap::new(),
            exon_file,
        }
    }

    pub fn add_lone_exon(&mut self, gene_id: &str, gene_name: &str, start: usize, end: usize, chromosome: String, sens_orientation: bool) {
        // Get the vector of genes for the specified chromosome
        let chromosome_genes = self.chromosomes.entry(chromosome.clone()).or_insert(Vec::new());
        
        let mut new_gene = Gene::new(gene_id, gene_id, start, end, sens_orientation);
        new_gene.add_exon(start, end);
        if let Some(ref mut writer) = self.exon_file {
            let line = format!("{}\t{}\t{}\t{}\t{}\n", gene_id, gene_name, chromosome, start, end);
            writer.write_all(line.as_bytes()).expect("Failed to write to file");
        }
        chromosome_genes.push(new_gene);
        
    }

    pub fn add_exon(&mut self, gene_id: &str, gene_name: &str, start: usize, 
        end: usize, chromosome: String, sens_orientation: bool) {
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
            if let Some(ref mut writer) = self.exon_file {
                let line = format!("{}\t{}\t{}\t{}\t{}\n", gene_id, gene_name, chromosome, start, end);
                writer.write_all(line.as_bytes()).expect("Failed to write to file");
            }
            chromosome_genes.push(new_gene);
        }
    }
    /*
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
    */

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
    /*
    fn genes_overlapping(&self, chr:&str, initial_position:usize, final_position:usize, iterator: ExonIterator  )-> Vec<Gene> {
        let mut res: Vec<Gene> = Vec::new();
        
        if let Some(genes) = self.chromosomes.get( chr ) {
            let mut start_id = iterator.gene_id();
            while genes[start_id].start < final_position && genes[start_id].end >= initial_position {
                res.push( genes[start_id].clone() );
                start_id +=1;
            }
        }
        res
    }   
    */
    /// the main query functionality
    pub fn match_cigar_to_gene(&self, chr:&str, cigar: &str, initial_position: usize, 
        iterator: &mut ExonIterator ) -> Option<Vec<ReadResult>> {

        match iterator.last_result_matches( cigar, initial_position ){
            Some(result) => {
                //println!("resused old one");
                return Some(result)
            },
            None =>  {},
        }

        let (spliced_read, final_position) = SplicedRead::new(cigar, initial_position);

        let mut results = Vec::<ReadResult>::new();
        let mut best_result_id = 0;
        let mut best_result = RegionStatus::AfterGene; // worst

        if let Some(genes) = self.chromosomes.get( chr ) {
            //println!("Skipping {:?} entries", iterator.gene_id() );
            for (index, gene) in genes.iter().enumerate().skip( iterator.gene_id() ) {
                if gene.start > final_position {
                    //println!("The next gene's start position is after my end position - no match?");
                    break
                }
                if gene.end < initial_position{
                    //println!("this gene here end before the match!");
                    iterator.set_gene_id( index + 1);
                    continue;
                }
                let result = gene.match_to( &spliced_read );
                if result.match_type.is_better( &RegionStatus::AfterGene ) {
                    //println!("Result status is no more BeforeGene: {:?} - skipping {} gene info", result.match_type, index );
                    if result.match_type.is_better( &best_result ) {
                        best_result_id = results.len();
                        best_result = result.match_type.clone();
                        results.push( result.clone() );
                        //println!("{}: found a new best result: id {best_result_id}, value {best_result:?}", iterator.gene_id());
                    }
                    if result.match_type  == RegionStatus::AfterGene {
                        break;
                    }
                }
            }
            //panic!("One is enough here!");
        }else {
            //println!("ExtTag {chr} with start {initial_position} and cigar {cigar}");
            let res = ReadResult { gene: chr.to_string(), sens_orientation: true,  match_type: RegionStatus::ExtTag };
            results.push( res );
        }
        // now we have a best matching gene as 
        if ! results.is_empty() {
            //println!("We found a match!");
            let best_type = results[best_result_id].match_type.clone();

            // Filter all results that have the same match_type as the best match
            let best_matches: Vec<_> = results
                .iter()
                .filter(|result| result.match_type == best_type)
                .cloned() // Clone the matched items
                .collect();

            Some(best_matches)

        }else {
            None
        }

    }

    /// the main query functionality
    pub fn match_cigar_to_gene_overlap(&self, chr:&str, cigar: &str, initial_position: usize, 
        iterator: &mut ExonIterator ) -> Option<Vec<ReadResult>> {

        match iterator.last_result_matches( cigar, initial_position ){
            Some(result) => {
                //println!("resused old one");
                return Some(result)
            },
            None =>  {},
        }

        let (spliced_read, final_position) = SplicedRead::new(cigar, initial_position);

        let mut results = Vec::<ReadResult>::new();
        let mut best_result_id = 0;
        let mut best_result = RegionStatus::AfterGene; // worst

        if let Some(genes) = self.chromosomes.get( chr ) {
            //println!("Skipping {:?} entries", iterator.gene_id() );
            for (index, gene) in genes.iter().enumerate().skip( iterator.gene_id() ) {
                if gene.start > final_position {
                    //println!("The next gene's start position is after my end position - no match?");
                    break
                }
                if gene.end < initial_position{
                    //println!("this gene here end before the match!");
                    iterator.set_gene_id( index + 1);
                    continue;
                }
                let result = gene.match_to_overlap( &spliced_read );
                if result.match_type.is_better( &RegionStatus::AfterGene ) {
                    //println!("Result status is no more BeforeGene: {:?} - skipping {} gene info", result.match_type, index );
                    if result.match_type.is_better( &best_result ) {
                        best_result_id = results.len();
                        best_result = result.match_type.clone();
                        results.push( result.clone() );
                        //println!("{}: found a new best result: id {best_result_id}, value {best_result:?}", iterator.gene_id());
                    }
                    if result.match_type  == RegionStatus::AfterGene {
                        break;
                    }
                }
            }
            //panic!("One is enough here!");
        }else {
            //println!("ExtTag {chr} with start {initial_position} and cigar {cigar}");
            let res = ReadResult { gene: chr.to_string(), sens_orientation: true,  match_type: RegionStatus::ExtTag };
            results.push( res );
        }
        // now we have a best matching gene as 
        if ! results.is_empty() {
            //println!("We found a match!");
            let best_type = results[best_result_id].match_type.clone();

            // Filter all results that have the same match_type as the best match
            let best_matches: Vec<_> = results
                .iter()
                .filter(|result| result.match_type == best_type)
                .cloned() // Clone the matched items
                .collect();

            Some(best_matches)

        }else {
            None
        }

    }

    fn find_search_position(&self, genes: &Vec<Gene>, pos: usize) -> Option<usize> {
        let mut low = 0;
        let mut high = genes.len();
        
        if genes[0].start > pos {
            return Some(0)
        }
        if genes[genes.len()-1].end < pos {
            return None
        }
        // Perform a binary search for the first gene that starts after or overlaps with pos
        while low < high {
            let mid = (low + high) / 2;
            
            if genes[mid].start > pos {
                //println!("Right half");
                high = mid; // Narrow down to the left half
            } else {
                //println!("Left half");
                low = mid + 1; // Narrow down to the right half
            }
        }

        // At this point, `low` is the first index where the gene's start is >= pos or closest after `pos`
        if low > 0{
            if genes[low - 1].end >= pos ||  genes[low].start > pos{
                return Some(low - 1);
            }
        }

        // Otherwise, check the gene at the found index (which should be the first gene starting after pos)
        if low < genes.len() && genes[low].end < pos {
            return Some(low) ;
        }

        // panic!("I have not found a valid gene for pos {pos}! low {} -1 end {}  low end {}", low, genes[low.saturating_sub( 1 )].end, genes[low].end  );
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
                            //println!("Init set the gene_id to {gene_id} and exon_id to {exon_id}");
                            iterator.set_gene_id( gene_id );
                            iterator.set_exon_id( exon_id );
                            Ok(())
                        }else {
                            //println!("I have not found a search position!?");
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

    // This collects exons from a TE gtf
    pub fn parse_gtf_only_exons(&mut self, file_path: &str, exact_attr:&str) -> Result<(), Box<dyn Error>> {
        let path = Path::new(file_path);

        // Determine the reader: plain text or gzip
        let reader: Box<dyn Read> = if file_path.ends_with(".gz") || Self::is_gzipped(file_path)? {
            let file = File::open(&path)?;
            Box::new(GzDecoder::new(file))
        } else {
            let file = File::open(&path)?;
            Box::new(file)
        };

        let reader = BufReader::new(reader);

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
            let orientation = match fields[6]{
                "+" => true,
                "-" => false,
                _ => panic!("field[7] should not contain this value: {}",fields[6]),
            };

            if feature_type == "exon" {
                let start: usize = fields[3].parse()?;
                let end: usize = fields[4].parse()?;

                // Extract gene_id and gene_name from the attributes
                let attributes = fields[8];
                let gene_id = Self::extract_attribute(attributes, exact_attr ).unwrap_or("unknown");
                let gene_name = Self::extract_attribute(attributes, "gene_name").unwrap_or("unknown");

                let cleaned_gene_id: String = gene_id.chars()
                    .filter(|&c| c != '"' && c != '\'')
                    .collect();

                let cleaned_gene_name: String = gene_name.chars()
                    .filter(|&c| c != '"' && c != '\'')
                    .collect();

                // Call your add_exon function (you may want to add gene_name handling there)
                //println!("Adding gene {gene_id} - {gene_name}");
                self.add_lone_exon(&cleaned_gene_id, &cleaned_gene_name, start, end, chromosome, orientation);
            }
        }

        eprintln!("I have read this:\n{}", self);
        Ok(())
    }

    // Helper function to detect gzip format by inspecting magic bytes
    fn is_gzipped(file_path: &str) -> Result<bool, Box<dyn Error>> {
        let mut file = File::open(file_path)?;
        let mut magic_bytes = [0; 2];
        file.read_exact(&mut magic_bytes)?;
        Ok(magic_bytes == [0x1F, 0x8B]) // gzip magic bytes
    }
        
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

    // Function to parse the GTF file and populate the Gtf structure
    pub fn parse_gtf(&mut self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let path = Path::new(file_path);

        // Determine the reader: plain text or gzip
        let reader: Box<dyn Read> = if file_path.ends_with(".gz") || Self::is_gzipped(file_path)? {
            let file = File::open(&path)?;
            Box::new(GzDecoder::new(file))
        } else {
            let file = File::open(&path)?;
            Box::new(file)
        };

        let reader = BufReader::new(reader);

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
            let orientation = match fields[6]{
                "+" => true,
                "-" => false,
                _ => panic!("field[7] should not contain this value: {}",fields[6]),
            };

            if feature_type == "exon" {
                let start: usize = fields[3].parse()?;
                let end: usize = fields[4].parse()?;

                // Extract gene_id and gene_name from the attributes
                let attributes = fields[8];
                let gene_id = Self::extract_attribute(attributes, "gene_id").unwrap_or("unknown");
                let gene_name = Self::extract_attribute(attributes, "gene_name").unwrap_or("unknown");

                let cleaned_gene_id: String = gene_id.chars()
                    .filter(|&c| c != '"' && c != '\'')
                    .collect();

                let cleaned_gene_name: String = gene_name.chars()
                    .filter(|&c| c != '"' && c != '\'')
                    .collect();

                // Call your add_exon function (you may want to add gene_name handling there)
                //println!("Adding gene {gene_id} - {gene_name}");
                self.add_exon(&cleaned_gene_id, &cleaned_gene_name, start, end, chromosome, orientation);
            }
        }

        eprintln!("I have read this:\n{}", self);
        Ok(())
    }
    
}




