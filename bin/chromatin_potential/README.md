# nf-rnachrom-chromatin_potential

>cd /home/snap/projects/lncRNA_app/chromatin_potential

>conda activate pipline

Usage: sh count_contacts_all.sh -d <distance> -i <input_path> -o <output_path> -u <uu_file> -m <um_file> -n <N_contacts_min> -f <fdr_threshold> -r <input_path_RNAseq>

Example:
/usr/bin/time -v sh count_contacts_all.sh \
            -d 500000 \
            -i /home/snap/projects/lncRNA_app/voting/output_SRR17331253_UU_UM \
            -o /home/snap/projects/lncRNA_app/voting/output_SRR17331253_UU_UM/chP \
            -u contacts.voting.UU.bed \
            -m contacts.voting.UM.bed \
            -n 100 \
            -f 0.05 \
            -r /home/snap/Downloads