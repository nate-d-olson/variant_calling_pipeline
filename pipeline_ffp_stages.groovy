////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////

flash_fastq = {
    doc "Extends reads where forward and reverse reads overlap"
    produce(input.prefix+".extendedFrags.fastq",
            input.prefix+".notCombined_1.fastq",
            input.prefix+".notCombined_2.fastq",
            input.prefix+".hist",
            input.prefix+".histogram") {
        exec "/Users/nolson/Desktop/FLASH-1.2.11/flash -M $MAXOVERLAP -o $input.prefix $inputs"
    }
}

@Filter("combined")
combine_fastq = {
    doc "combine unmatched forward, reverse, and extended fastq files"
    exec """
        cat $input1.fastq 
            $input2.fastq 
            $input3.fastq
            > $output
    """
}

@Filter("clean")
clean_fastq = {
    doc "filter fastq file based on base quality"
    exec """
        $FASTX_BIN/fastq_quality_filter -v 
            -Q33
            -q $QMIN 
            -p $PQMIN 
            -i $input 
            -o $output
    """
}

@Transform("fasta")
convert_fastq_to_fasta = {
    doc "convert fastq file to fasta file"
    exec """
        $FASTX_BIN/fastq_to_fasta -v 
            -Q33
            -i $input.fastq
            -o $output
    """
}

@Transform("vec.txt")
ffp_vectors = {
    doc "create ffp vector file"
    exec """
        ffpry -l $KMER -d $inputs.fasta > $output
    """
}

@Transform("names.txt")
file_names = {
    doc "list of file names for phylogenetic analysis"
    exec """
        ls $inputs.fasta | sed 's/_L001_R1_001.extendedFrags.combinded.clean.fasta//g'  > $output
    """
    forward input
}

@Transform("mat")
ffp_boot = {
    doc "bootstrap replicates converting vectors to matrix"
    exec """
        for i in \$(seq 1 1 $REP);
            do;
                ffpboot $input.vec.txt | 
                ffprwn | 
                ffpjsd -p $input.names.txt >> $output;
            done
    """
}

@Transform("nwk")
ffp_tree = {
    doc "generate netwick tree file from matrix"
    exec "ffptree -q $input > $output"
    
}