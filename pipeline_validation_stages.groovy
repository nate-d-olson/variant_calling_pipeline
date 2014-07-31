////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////

ref_index = { 
	exec "$BIN/bwa index -a is $input.fasta"
}

@Transform("bwa.sam")
bwaMEMalign = {
	exec """
		$BIN/bwa mem 
		    -t $n
			$input.fasta
			$input1.fastq.gz $input2.fastq.gz > $input.fasta.prefix.$output
	"""
}


@Transform("bam")
samToSortedBam = {
    doc "Sort a SAM file so that it is compatible with reference order and convert to BAM file"
    exec"""
        java $JHEAP -jar $BIN/SortSam.jar 
                    VALIDATION_STRINGENCY=LENIENT 
                    INPUT=$input.fasta.prefix.$input 
                    OUTPUT=$input.fasta.prefix.$output
                    SORT_ORDER=coordinate
    """
}

readGroups = {
	// will want to work to specify values
    exec """
        java $JHEAP -jar $BIN/AddOrReplaceReadGroups.jar 
                    INPUT=$input.fasta.prefix.$input 
                    OUTPUT=$input.fasta.prefix.$output
                    RGID=1
                    RGLB=S0h_-1_S1
                    RGPL=illumina
                    RGPU=S0h-1_S1
                    RGSM=RM8375
                    RGCN=NIST
                    RGDS=MiSeq-RM8375
    """
}

indexBam = {
    // A bit of a hack to ensure the output file is expected in the
    // same directory as the input bam, no matter where it is
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "$BIN/samtools index $input.fasta.prefix.$input"
    }
    forward input
}

@Transform("dedup.bam")
dedup = {
    exec """
        java $JHEAP  -jar $BIN/MarkDuplicates.jar
             INPUT=$input.fasta.prefix.$input 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$LOG 
             OUTPUT=$input.fasta.prefix.$output
    """
}

@Transform("vcf")
call_variants_freebayes = {
	exec """
		$BIN/freebayes -p 1
		-v $input.fasta.prefix.$output -f $input.fasta -b $input.fasta.prefix.$input
	"""
}
