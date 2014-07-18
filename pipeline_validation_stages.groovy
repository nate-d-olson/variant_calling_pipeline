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
			$input1.fastq.gz $input2.fastq.gz > $output
	"""
    }
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

@Transform("bam.bai")
indexBam = {
    exec "$BIN/samtools index $input"
    forward input $output
}

@Filter("dedup")
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
		$BIN/freebayes -p 1 -f $input.fasta -b $input -v $output
	"""
}