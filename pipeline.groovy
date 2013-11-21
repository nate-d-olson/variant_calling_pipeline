////////////////////////////////////////////////////////////
// GATK-based variant-calling pipeline, WGS version.
// 
// This pipeline is a quick conversion of the pipeline at 
// 
//     https://github.com/claresloggett/variant_calling_pipeline/
// 
// It doesn't work yet, and is largely untested
//
// Assume: paired end reads
// Assume: files in form  *<sample_name>* _R1.fastq.gz, *<sample_name>*_R2.fastq.gz
//
////////////////////////////////////////////////////////////

load 'pipeline_stages_config.groovy'

// Memory in gigabytes to be used for the pipeline
mem = 6

run {
    "%_*_R*" * [
               fastqc +
               "%.gz" * [ alignBWA ]  +
               alignToSamPE +
               samToSortedBam +
               indexBam +
               dedup + 
               indexBam 
                /* TODO: add and test more as below ...
                   depthOfCoverage // + ...
                 */
                    
        ]
    }
