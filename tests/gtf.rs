
// Tests for the gtf_parser library
#[cfg(test)]
mod tests {
    use quantify_bam::gtf::GTF;
    use quantify_bam::gtf::ExonIterator;
    use quantify_bam::gtf::RegionStatus;
    use quantify_bam::gtf::QueryErrors;


    #[test]
    fn test_parse_gtf() {
        let mut gtf = GTF::new();
        let result = gtf.parse_gtf("tests/data/test.gtf"); // Use a valid test GTF file path
        
        // Check if the parsing was successful
        assert!(result.is_ok());

        // Add your assertions here to verify the contents of the GTF
        // For example:
        assert_eq!(gtf.chromosomes.len(), 2);

        // Verify the names of the chromosomes (replace "chr1" and "chr2" with your expected names)
        let expected_chromosomes = vec!["Y".to_string(), "Y_mapped_Scaffold_9_D1573".to_string() ]; // Replace with actual expected names
        let mut actual_chromosomes: Vec<String> = gtf.chromosomes.keys().cloned().collect();
        actual_chromosomes.sort();
        assert_eq!(actual_chromosomes, expected_chromosomes);
    }

    #[test]
    fn test_slice_gtf() {
        let mut gtf = GTF::new(); // Use lowercase `gtf` for the variable
        let _ = gtf.parse_gtf("tests/data/test.gtf"); // Use a valid test GTF file path

        // Assuming we have added some test data here for the gtf

        // Perform a slice operation
        let sliced = match gtf.slice_gtf("Y", 3265434, 3276434) {
            Ok(slice) => slice,
            Err(e) => panic!("{:?}", e),
        };
        
        // Check the expected result
        assert_eq!(sliced.len(), 4);
    }

    #[test]
    fn test_query() {
        let mut gtf = GTF::new();

        // Add exons on chromosome "Y"
        gtf.add_exon("gene1_id", "gene1_name", 20, 40, "Y".to_string(), true);
        gtf.add_exon("gene2_id", "gene2_name", 50, 80, "Y".to_string(), true);
        gtf.add_exon("gene3_id", "gene3_name", 100, 200, "Y".to_string(), true);


        // Perform a search on chromosome "Y" with specific start-end positions
        let mut interator = ExonIterator::new("test");
        gtf.init_search( "Y", 10, &mut interator );

        let mut result_y = match gtf.match_cigar_to_gene("Y", 10, "5M", &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::BeforeGene );

        result_y = match gtf.match_cigar_to_gene("Y", 15, "10M", &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::SpanningBoundary );

        result_y = match gtf.match_cigar_to_gene("Y", 20, "10M", &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::InsideExon );

        result_y = match gtf.match_cigar_to_gene("Y", 10, "40M", &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::SpanningBoundary );

        result_y = match gtf.match_cigar_to_gene("Y", 10, "40M", &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::SpanningBoundary );
        assert_eq!(result_y.0, ("gene1_id".to_string(), "gene1_name".to_string() ));

        result_y = match gtf.match_cigar_to_gene("Y", 60, "20M", &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::InsideExon );
        assert_eq!(result_y.0, ("gene2_id".to_string(), "gene2_name".to_string() ));

        let _ = match gtf.match_cigar_to_gene("Y", 6000, "100M", &mut interator ){
            Ok(slice) => panic!("This should fail here") ,
            Err(e) => {
                assert_eq!(e, QueryErrors::OutOfGenes )
            },
        };
    }

}
