
BASE="/mnt/storage/shared/genomes/hg19/gatk"
REF="$BASE/gatk.ucsc.hg19.fasta"
DBSNP="$BASE/dbsnp_132.hg19.vcf"
LOG="pipeline.log"
GOLD_STANDARD_INDELS="gold_standard_indels.vcf"// ?
INDELS_100G="1000g_indels.vcf"// ?

fastqc = {
    exec "fastqc --quiet $input.gz"
}

fastqc = {
    def fastqc_outputs = inputs.gz.collect { new File(it).name }*.replaceAll('.fastq.gz','_fastqc.zip')
    produce(fastqc_outputs) {
        exec "fastqc --quiet -o . $inputs.gz"
    }
}

alignBWA = {

    var encoding_flag : "" // previous stage can set encoding_flag="-I ..."

    doc "Aligns using BWA. Note: assumes input file are gzipped"
    exec "bwa aln -t 8 $encoding_flag $REF $input.gz > $output.sai"
}

alignToSamSE = {
    // Not sure what meta is ...
    var meta : ""
    exec "bwa samse $REF $meta $input.sai $input.gz > $output.sam"
}

alignToSamPE = {
    // Not sure what meta is ...
    var meta : ""
    exec "bwa sampe $REF $meta $input1.sai $input2.sai $input2.gz $input2.gz > $output.sam"
}

samToSortedBam = {
    exec "./SortSam 6 VALIDATION_STRINGENCY=LENIENT INPUT=$input.sam OUTPUT=$output.bam SORT_ORDER=coordinate"
}

mergeBams = {
    exec "./PicardMerge 6 $inputs.bam USE_THREADING=true VALIDATION_STRINGENCY=LENIENT AS=true OUTPUT=$output.bam"
}

indexBam = {
    exec "samtools index $input.bam"
    forward input
}

flagstat = {
    exec "samtools flagstat $input.bam > $output"
}

igvcount = {
    exec "igvtools count $input.bam $output hg19"
}

indexVCF = {
    exec "./vcftools_prepare.sh $input.vcf"
}

realignIntervals = {
    // Hard-coded to take 2 known indels files right now
    exec "./GenomeAnalysisTK 1 -T RealignerTargetCreator -R $REF -I $input.bam --known $GOLD_STANDARD_INDELS --known $INDELS_100G -log $LOG -o $output.intervals"
}

realign = {
    exec "./GenomeAnalysisTK 22 -T IndelRealigner -R $REF -I $input.bam -targetIntervals $input.intervals -log $LOG -o $output"
}

dedup = {
    exec "./MarkDuplicates 6 INPUT=$input.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT AS=true METRICS_FILE=$LOG OUTPUT=$output"
}

baseQualRecalCount = {
    exec "./GenomeAnalysisTK 12 -T CountCovariates -I $input.bam -R $REF --knownSites $DBSNP -nt 8 -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -log $LOG -recalFile $output"
}

baseQualRecalTabulate = {
    exec "./GenomeAnalysisTK 4 -T TableRecalibration -I $input.bam -R $REF -recalFile $input.csv -l INFO -log $LOG -o $output"
}

callSNPs = {
    exec "./GenomeAnalysisTK 12 -T UnifiedGenotyper -nt 8 -R $REF -I $input.bam --dbsnp $DBSNP -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1600 -l INFO -A AlleleBalance -A DepthOfCoverage -A FisherStrand -glm SNP -log $LOG -o $output"
}

callIndels = {
    exec "./GenomeAnalysisTK 12 -T UnifiedGenotyper -nt 8 -R $REF -I $input.bam --dbsnp $DBSNP -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1600 -l INFO -A AlleleBalance -A DepthOfCoverage -A FisherStrand -glm INDEL -log $LOG -o $output"
}

filterSNPs = {
    // Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
    exec "./GenomeAnalysisTK 4 -T VariantFiltration -R $REF --variant $input.vcf --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'GATK_MINIMAL_FILTER' -log $LOG -o $output"
}

filterIndels = {
    // Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
    // If you have 10 or more samples GATK also recommends the filter InbreedingCoeff < -0.8
    exec "./GenomeAnalysisTK 4 -T VariantFiltration -R $REF --variant $input.vcf --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0' --filterName 'GATK_MINIMAL_FILTER' -log $LOG -o $output"
}

annotateEnsembl = {
    // This command as written assumes that VEP and its cache have been
    // downloaded in respective locations
    // ./variant_effect_predictor_2.5
    // ./variant_effect_predictor_2.5/vep_cache
    exec "perl variant_effect_predictor_2.5/variant_effect_predictor.pl --cache --dir variant_effect_predictor_2.5/vep_cache -i $input.vcf --vcf -o $output -species human --canonical --gene --protein --sift=b --polyphen=b > $LOG"
}

depthOfCoverage = {
    exec "./GenomeAnalysisTK 4 -T DepthOfCoverage -R $REF -I $input.bam -omitBaseOutput -ct 1 -ct 10 -ct 20 -ct 30 -o $output"
}

collateReadcounts = {
    output.dir = outdir
    def dir = new File(input).parentFile.name // hack - figure out directory from input file
    exec "python count_flagstat_wgs.py $dir $outdir"
}
