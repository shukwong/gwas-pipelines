version development


workflow run_metal_workflow {
    input {
        File association_summary_files_list
    	String prefix 
		String FREQLABEL
    	String? MARKERLABEL
    }    

    Array[File] association_summary_files  = read_lines(association_summary_files_list)

    call run_metal {
        input: 
            association_summary_files = association_summary_files,
     		prefix = prefix,
			FREQLABEL = FREQLABEL,
    		MARKERLABEL = MARKERLABEL
    }
}	


#TODO: add allele frequency for saige
task run_metal {
    input {
    Array[File] association_summary_files
    String prefix 

    String FREQLABEL
    String? MARKERLABEL = "SNP"
   

    Int? memory = 32
    Int? disk = 20
    Int? threads = 2
    Int? preemptible_tries = 3
    }

     command <<<
        set -euo pipefail

        wget https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/scripts/process_metal_outputs.R

        echo -e "MARKERLABEL ~{MARKERLABEL} \n \
ALLELELABELS Tested_Allele Other_Allele \n \
EFFECTLABEL BETA \n \
STDERRLABEL SE \n \
CUSTOMVARIABLE TotalSampleSize \n \
LABEL TotalSampleSize as N \n \
SCHEME STDERR \n \
GENOMICCONTROL ON \n \
PROCESSFILE ~{sep=',PROCESSFILE ' association_summary_files}  \n\
OUTFILE ~{prefix}.metal .tsv \n\
ANALYZE HETEROGENEITY \n\
QUIT" > metal_command_temp

        cat  metal_command_temp | tr ',' '\n' >metal_command

        metal metal_command

        echo -e "MarkerName\tCHR\tPOS\tAllele1\tAllele2\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tTotalSampleSize">~{prefix}.metal.tsv
        grep -v MarkerName ~{prefix}.metal1.tsv | tr ':' '\t' | awk '{print $1":"$2":"$3":"$4"\t"$0}' | cut -f1-3,6-16 >>~{prefix}.metal.tsv

        gzip ~{prefix}.metal.tsv
        mv ~{prefix}.metal1.tsv.info ~{prefix}.metal.tsv.info
        
    >>>

    runtime {
		docker: "quay.io/shukwong/metal:2011-03-25"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: "${preemptible_tries}"
	}

    output {
        File metal_output_file =  prefix + ".metal.tsv.gz"
        File metal_info_file = prefix + ".metal.tsv.info"
    }
}
