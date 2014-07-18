////////////////////////////////////////////////////////////
// GATK-based variant-calling pipeline, WGS version.
//  Pipeline port from - https://github.com/ssadedin/variant_calling_pipeline which was a port from
// This pipeline is a port of the VLSCI whole genome variant
// calling pipeline from this Git repo to Bpipe:
// 
//     https://github.com/claresloggett/variant_calling_pipeline/
//
// This example is intended as an illustration of how a full pipeline
// looks in Bpipe, especially for those who might be familiar with Ruffus.
// There are a number of places where the pipeline is non-ideal, but it has
// been kept that way to reflect the original as accurately as possible.
//
// One particular point of note is that the alignment output is stored
// in SAM format by some of the intermediate stages which would be very
// unadvisable for true whole-genome data. It is easy to make it store
// output in BAM format instead, so this would be a very advisable 
// change if you intend to run this pipeline on a lot of data. 
//
// There are a number of software requirements, which you should ensure are 
// satisfied before running the pipeline. These need to be configured 
// in the file called 'config.groovy' (which is loaded below). A template
// is provided in config.groovy.template which you can use to 
// create the file and set the right paths to your tools and reference
// data.
//
// By default this pipeline will attempt to use all the available cores
// on the computer it runs on. If you don't wish to do that, limit the 
// concurrency by running it with the -n flag:
//
//    bpipe run -n 4 pipeline.groovy example_data/input_data_wgs/*.fastq.gz
// 
// Assumes: paired end reads
// Assumes: files in form  *<sample_name>*_..._R1.fastq.gz, *<sample_name>*_..._R2.fastq.gz
//
// Author: Simon Sadedin, MCRI
//
////////////////////////////////////////////////////////////

// Create this file by copying config.groovy.template and editing
load "/home/nolson/Desktop/micro_rm_dev/variant_calling_pipeline/validation.config.groovy"

// All the core pipeline stages in the pipeline
load "$PRJ_HOME/variant_calling_pipeline/pipeline_validation_stages.groovy"

run {
	"%.fasta" * [ref_index + "%_*_R*" * [ bwaMEMalign  + samToSortedBam + indexBam + dedup + indexBam  + call_variants_freebayes ]]
}
