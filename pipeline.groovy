"""
GATK-based variant-calling pipeline, WGS version.

This pipeline is a quick conversion of the pipeline at 

    https://github.com/claresloggett/variant_calling_pipeline/

It doesn't work yet, and is largely untested.

"""
load 'pipeline_stages_config.groovy'

run {
    // Assume: paired end reads
    // Assume: files in from  *<sample_name>* _R1.fastq.gz, *<sample_name>*_R2.fastq.gz
    "%_*_R*" * [
               fastqc  +
               "%.gz" * [ alignBWA ] + 
               alignToSamPE 

                /* TODO: add and test more as below ...
                   samToSortedBam + 
                   indexBam +
                   dedup + 
                   indexBam +
                   depthOfCoverage // + ...
                 */
                    
        ]
    }
