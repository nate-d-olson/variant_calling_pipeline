////////////////////////////////////////////////////////////
// bpipe pipline for performing de-nono assembly for the 
// GMI pilot study 
// 
////////////////////////////////////////////////////////////

@transform("")
de_novo = {
	exec """
		perl /home/ubuntu/GMI_bioinf/a5_miseq_linux_20140604/bin/a5_pipeline.pl $inputs.fastq $input.prefix
	"""
}

run {
	"%_L001_R*_001" * [de_novo]
}
