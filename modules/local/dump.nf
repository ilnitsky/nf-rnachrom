// Check for fastq compressed inputs

    shell:
    ext_file = f_ext.getExtension()
    base_name_file = f_ext.getBaseName()
    if (ext_file == "fasta" | ext_file == "fa"){
	    '''
	    cp -n !{f_ext} !{base_name_file}.fasta
	    '''
    }else if(ext_file == "zip"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
	    '''
      gunzip -f -S .zip !{f_ext}
	    cp -n !{old_base_name_file} !{base_name_file}.fasta
	    '''
    }else if(ext_file == "gz"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
	    '''
	    gunzip -f !{f_ext}
	    cp -n !{old_base_name_file} !{base_name_file}.fasta
	    '''
    }else {
      '''
      echo "Your pathogen genome files appear to have the wrong extension. \n Currently, the pipeline only supports .fasta or .fa, or compressed files with .zip or .gz extensions."
      '''
    }
}




