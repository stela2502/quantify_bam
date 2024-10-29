//ExonIterator.rs
use core::fmt;

use std::time::{SystemTime};

use crate::gtf::{GTF, QueryErrors};

#[derive(Debug)]
pub struct ExonIterator {
    gene_id: usize,
    exon_id: usize,
    start: SystemTime,
    name: String,
}

/// ExonInterator is the multiprocessor storage of the current search state.
/// It stores the exon_id and gene_id of the current chromosomal area we are working on.

// Implement the Display trait for ExonIterator
impl fmt::Display for ExonIterator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        let (h, min, sec, mill ) = self.since_start();
        writeln!(f, "{} processing time: {}h {} min {} sec and {} millisec", 
            self.name, h, min, sec, mill )?;

        Ok(())
    }
}


impl ExonIterator {
    // Constructor to create a new ExonIterator with std 0,0 as the ids
    pub fn new( name:&str ) -> Self {
        ExonIterator {
            gene_id: 0,
            exon_id: 0,
            start: SystemTime::now(),
            name: name.to_string(),
        }
    }

    pub fn since_start(&self ) -> (u128, u128, u128, u128){
        let mut milli = self.start.elapsed().unwrap().as_millis();

        let mil = milli % 1000;
        milli= (milli - mil) /1000;

        let sec = milli % 60;
        milli= (milli -sec) /60;

        let min = milli % 60;
        milli= (milli -min) /60;

        (milli, min, sec, mil )
    }

    pub fn exon_id(&self) -> usize{
        self.exon_id
    }

    pub fn gene_id(&self) -> usize{
        self.gene_id
    }

    pub fn set_exon_id(&mut self, id:usize ) {
        self.exon_id = id;
    }

    pub fn set_gene_id(&mut self, id:usize ) {
        self.gene_id = id;
    }

}