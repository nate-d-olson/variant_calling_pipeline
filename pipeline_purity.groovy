////////////////////////////////////////////////////////////
// 	Pipeline for processing paried end illumina datasets using pathoscope
// 	mapping against whole GenBank database
//
//    bpipe run -n 4 pipeline.groovy *.fastq.gz
// 
// Assumes: paired end reads
// Assumes: files in form  *<sample_name>*_..._R1.fastq.gz, *<sample_name>*_..._R2.fastq.gz
//
// Author: Nate Olson, NIST
//
////////////////////////////////////////////////////////////

// Create this file by copying config.groovy.template and editing
//load "/home/nolson/Desktop/micro_rm_dev/variant_calling_pipeline/validation.config.groovy"

// All the core pipeline stages in the pipeline
//load "$PRJ_HOME/variant_calling_pipeline/pipeline_validation_stages.groovy"

//REF="/media/nolson/second/DATAFILES/Bacteria/all_bac_genomes.fasta"
REF="~/GMI_bioinf/all_bac/all_bac_genomes.fasta"

@Transform("bwa.sam")
bwaMEMallOut = {
	exec """
		~/bwa/bwa mem -t 8 -a $REF $input1.fastq.gz $input2.fastq.gz > $output
	"""
}

pathoscope = {
	exec """
		python /home/ubuntu/pathoscope2/pathoscope/pathoscope.py --verbose REP
		-samfile=$input
	"""
}

run {
	"%_*_R*" * [ bwaMEMallOut + pathoscope]
}
