nextflow.enable.dsl=2
ANSI_RESET = "\u001B[0m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";

def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }


//TODO Prepare boost, julia, etc
//TODO Singularity support
process PrepareSoftware {
    log.info print_cyan("Installing https://github.com/AndrewSigorskih/RNAChrom_ALLvsALL_data_processing.git ")
    conda "${projectDir}/envs/rnachromprocessing.yaml"
    output:
    path "binaries_ready.txt" 
    
    script:
    def rnachrom = new File("${projectDir}/bin/RnaChromATA/setup.py")
    if (!rnachrom.exists()) {
        """
        git clone https://github.com/AndrewSigorskih/RNAChrom_ALLvsALL_data_processing.git ${projectDir}/bin/RnaChromATA
        pip install ${projectDir}/bin/RnaChromATA/
        touch binaries_ready.txt
        """
    } else {
        """
        python3 ${projectDir}/bin/RnaChromATA/src/RnaChromProcessing/Processing/AllStages.py
        touch binaries_ready.txt
        """
    }

}

// pip install ${projectDir}/bin/RNAChrom_ALLvsALL_data_processing 