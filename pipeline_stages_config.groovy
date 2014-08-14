////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////

fastqc = {
    doc "Run FASTQC to generate QC metrics for the reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} $inputs.gz"
    }
}


ref_index = { 
	exec "~/bwa/bwa index -a is $input.fasta"
}

@Transform("bwa.sam")
bwaMEMalign = {
	exec """
		~/bwa/bwa mem 
		    -t 8
			$REF
			$input1 $input2> $output
	"""
}


@transform("bam")
samToSortedBam = {
    doc "Sort a SAM file so that it is compatible with reference order and convert to BAM file"
    output.dir="align"
    exec"""
        java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/SortSam.jar 
                    VALIDATION_STRINGENCY=LENIENT 
                    INPUT=$input.sam 
                    OUTPUT=$output.bam 
                    SORT_ORDER=coordinate
    """
}

readGroups = {
    output.dir="align"
    exec """
        java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/AddOrReplaceReadGroups.jar 
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

indexBam = {
    // A bit of a hack to ensure the output file is expected in the
    // same directory as the input bam, no matter where it is
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "~/bin/samtools index $input.bam"
    }
    forward input
}

flagstat = {
    exec "~/bin/samtools flagstat $input.bam > $output"
}



dedup = {
    output.dir="align"
    exec """
        java -Xmx2g -Djava.io.tmpdir=$TMPDIR -jar $PICARD_HOME/MarkDuplicates.jar
             INPUT=$input.bam 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$LOG 
             OUTPUT=$output.bam
    """
}

baseQualRecalCount = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="align"
    exec "java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $input.bam -R $REF --knownSites $DBSNP -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -log $LOG -o $output.counts"
}

baseQualRecalTabulate = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="align"
    exec "java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -I $input.bam -BQSR $input.counts -R $REF -l INFO -log $LOG -o $output"
}

callSNPs = {
    doc "Call SNPs/SNVs using GATK Unified Genotyper"
    output.dir="variants"
    exec """
            java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
               -nt $threads 
               -R $REF 
               -I $input.bam 
               -stand_call_conf 50.0 -stand_emit_conf 10.0 
               -dcov 1600 
               -l INFO 
               -A AlleleBalance -A DepthOfCoverage -A FisherStrand 
               -glm SNP -log $LOG 
               -o $output.vcf
        """
}

callIndels = {
    doc "Call variants using GATK Unified Genotyper"
    output.dir="variants"
    exec """
        java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
             -nt $threads
             -R $REF 
             -I $input.bam 
             -stand_call_conf 50.0 -stand_emit_conf 10.0 
             -dcov 1600 
             -l INFO 
             -A AlleleBalance -A DepthOfCoverage -A FisherStrand 
             -glm INDEL 
             -log $LOG -o $output.vcf
    """
}

@filter("filter")
filterSNPs = {
    // Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
    output.dir="variants"
    exec """
        java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration 
             -R $REF 
             --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' 
             --filterName 'GATK_MINIMAL_FILTER'
             --variant $input.vcf 
             -log $LOG 
             -o $output.vcf
    """
}

@filter("filter")
filterIndels = {
    doc """
            Filter data using very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
            If you have 10 or more samples GATK also recommends the filter InbreedingCoeff < -0.8
        """
    output.dir="variants"
    exec """
        java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration 
                    -R $REF 
                    --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0' 
                    --filterName 'GATK_MINIMAL_FILTER' -log $LOG 
                    --variant $input.vcf 
                    -o $output.vcf
    """
}

@filter("vep")
annotateEnsembl = {
    doc "Annotate variants using VEP to add Ensemble annotations"
    output.dir="variants"
    exec """
        perl $VEP/variant_effect_predictor.pl --cache --dir ./vep_cache -i $input.vcf --vcf -o $output.vcf -species human --canonical --per_gene --protein --sift=b --polyphen=b > $LOG
    """
}

depthOfCoverage = {
    output.dir="qc"
    transform("bam") to ("sample_statistics","sample_interval_summary") {
        exec """
            java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar 
                    -T DepthOfCoverage 
                    -R $REF -I $input.bam 
                    -omitBaseOutput 
                    -ct 1 -ct 10 -ct 20 -ct 30 
                    -o $output.prefix
        """
    }
}

@Transform("vcf")
call_variants_gatk = {
	exec """
		java -Xmx2g -jar ~/bin/GenomeAnalysisTK.jar 
			-R $REF
			-I $input
			-T HaplotypeCaller 
			-o $output
	"""
}
