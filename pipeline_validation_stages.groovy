////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////

ref_index = { 
	exec "$BIN/bwa index -a is $REF"
}

@Transform("bwa.sam")
bwaMEMalign = {
	exec """
		~/bwa/bwa mem 
		    -t n
			$REF
			$input1 $input2> $output
	"""
}


@transform("bam")
samToSortedBam = {
    doc "Sort a SAM file so that it is compatible with reference order and convert to BAM file"
    output.dir=$ALIGN
    exec"""
        java $JHEAP -Djava.io.tmpdir=$TMPDIR  -jar $BIN/SortSam.jar 
                    VALIDATION_STRINGENCY=LENIENT 
                    INPUT=$input.sam 
                    OUTPUT=$output.bam 
                    SORT_ORDER=coordinate
    """
}

readGroups = {
    output.dir="align"
    exec """
        java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $BIN/AddOrReplaceReadGroups.jar 
                    INPUT=$input.bam
                    OUTPUT=$output.bam
                    RGID=1
                    RGLB=S0h_-1_S1
                    RGPL=illumina
                    RGPU=S0h-1_S1
                    RGSM=RM8375
                    RGCN=NIST
                    RGDS=MiSeq-RM8375
    """
}

@filter("")
indexBam = {
    // A bit of a hack to ensure the output file is expected in the
    // same directory as the input bam, no matter where it is
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "~/bin/samtools index $input.bam"
    }
    forward input
}


dedup = {
    output.dir="align"
    exec """
        java -Xmx2g  -jar $BIN/MarkDuplicates.jar
             INPUT=$input.bam 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$LOG 
             OUTPUT=$output.bam
    """
}

@Transform("vcf")
call_variants_gatk = {
	exec """
		$BIN/freebayes -v $output -f $REF -b $input
	"""
}
