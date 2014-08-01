////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////

//@Produce("sample001.fasta.amb")//,"${input.fasta}.ann","${input.fasta}.bwt","${input.fasta}.pac","${input.fasta}.sa")
ref_index = { 
<<<<<<< HEAD
	exec "$BIN/bwa index -a is $input.fasta"
=======
    exec "$BIN/bwa index -a is $input.fasta"
>>>>>>> 2937bf88ff6840836e222b947c9a056ad4bfa317
}


@Transform("bwa.sam")
bwaMEMalign = {
	exec """
		$BIN/bwa mem 
		    -t $n
			$input.fasta
<<<<<<< HEAD
			$input1.fastq.gz $input2.fastq.gz > $input.fasta.prefix.$output
=======
			$inputs.fastq.gz > $input.fasta.prefix.$output
>>>>>>> 2937bf88ff6840836e222b947c9a056ad4bfa317
	"""
}


@Transform("bam")
samToSortedBam = {
    doc "Sort a SAM file so that it is compatible with reference order and convert to BAM file"
    exec"""
        java $JHEAP -jar $BIN/SortSam.jar 
                    VALIDATION_STRINGENCY=LENIENT 
<<<<<<< HEAD
                    INPUT=$input.fasta.prefix.$input 
=======
                    INPUT=$input.fasta.prefix.$input
>>>>>>> 2937bf88ff6840836e222b947c9a056ad4bfa317
                    OUTPUT=$input.fasta.prefix.$output
                    SORT_ORDER=coordinate
    """
}

readGroups = {
	// will want to work to specify values
    exec """
        java $JHEAP -jar $BIN/AddOrReplaceReadGroups.jar 
<<<<<<< HEAD
                    INPUT=$input.fasta.prefix.$input 
=======
                    INPUT=$input.fasta.prefix.$input
>>>>>>> 2937bf88ff6840836e222b947c9a056ad4bfa317
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
<<<<<<< HEAD
             INPUT=$input.fasta.prefix.$input 
=======
             INPUT=$input.fasta.prefix.$input
>>>>>>> 2937bf88ff6840836e222b947c9a056ad4bfa317
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
<<<<<<< HEAD
		$BIN/freebayes -p 1
		-v $input.fasta.prefix.$output -f $input.fasta -b $input.fasta.prefix.$input
=======
		$BIN/freebayes -p 1 -@ ../../references/sim/sim_variants.vcf.gz -v $input.fasta.prefix.$output -f $input.fasta -b $input.fasta.prefix.$input
>>>>>>> 2937bf88ff6840836e222b947c9a056ad4bfa317
	"""
}
