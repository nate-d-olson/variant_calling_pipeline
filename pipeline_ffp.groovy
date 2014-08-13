////////////////////////////////////////////////////////////
// Pipeline for ffp phylogenetic analysis using fastq datasets as input
//
// By default this pipeline will attempt to use all the available cores
// on the computer it runs on. If you don't wish to do that, limit the 
// concurrency by running it with the -n flag:
//
//    bpipe run -n 4 pipeline_ffp.groovy *.fastq
// 
// Assumes: paired end reads
// Assumes: files in form  *<sample_name>*_L001_R1_001.fastq, *<sample_name>*_L001_R2_001.fastq
//
// Author: Nate Olson, NIST
//
////////////////////////////////////////////////////////////

// Create this file by copying config.groovy.template and editing
load '/Users/nolson/Desktop/variant_calling_pipeline/ffp_config.groovy'

// All the core pipeline stages in the pipeline
load '/Users/nolson/Desktop/variant_calling_pipeline/pipeline_ffp_stages.groovy'

run {
    "%_L001_R*_001" * [flash_fastq + combine_fastq + clean_fastq + convert_fastq_to_fasta] +
    [file_names, ffp_vectors] + ffp_boot + ffp_tree  
}
