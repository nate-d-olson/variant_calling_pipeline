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
			$inputs.fastq.gz > $output
	"""
}


@Transform("bam")
samToSortedBam = {
    doc "Sort a SAM file so that it is compatible with reference order and convert to BAM file"
    exec"""
        java $JHEAP -jar $BIN/SortSam.jar 
                    VALIDATION_STRINGENCY=LENIENT 
                    INPUT=$input 
                    OUTPUT=$output
                    SORT_ORDER=coordinate
    """
}

readGroups = {
	// will want to work to specify values
    exec """
        java $JHEAP -jar $BIN/AddOrReplaceReadGroups.jar 
                    INPUT=$input 
                    OUTPUT=$output
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
        exec "$BIN/samtools index $input"
    }
    forward input
}

@Transform("dedup.bam")
dedup = {
    exec """
        java $JHEAP  -jar $BIN/MarkDuplicates.jar
             INPUT=$input 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$LOG 
             OUTPUT=$output
    """
}

@Transform("vcf")
call_variants_freebayes = {
	exec """
		$BIN/freebayes -p 1 -v $output -f $input.fasta -b $input
	"""
}
