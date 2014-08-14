////////////////////////////////////////////////////////////
// bpipe pipline for performing de-nono assembly for the 
// GMI pilot study 
// 
////////////////////////////////////////////////////////////

@transform("")
de_novo = {
	exec """
<<<<<<< HEAD
		perl /home/ubuntu/GMI_bioinf/a5_miseq_linux_20140604/bin/a5_pipeline.pl $inputs.fastq $input.prefix
=======
		perl /home/ubuntu/GMI_bioinf/a5_miseq_linux_20140604/bin/a5_pipeline.pl $inputs.fastq.gz
		$output
>>>>>>> c09314fcd781094c0185be785069efb15bf717a5
	"""
}

run {
	"%_L001_R*_001" * [de_novo]
}
