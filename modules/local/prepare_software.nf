nextflow.enable.dsl=2
ANSI_RESET = "\u001B[0m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";

def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }


//TODO Prepare boost, julia, etc
//TODO Singularity support
process PrepareSoftware {
    
    conda "${projectDir}/envs/rnachromprocessing.yaml"
    output:
    path "binaries_ready.txt" 
    
    script:
    def rnachrom = new File("${projectDir}/bin/RnaChromATA/setup.py")
    def bitap = new File("${projectDir}/bin/bitap")

    """
    set -e

    install_rnachrom() {
        git clone https://github.com/AndrewSigorskih/RNAChrom_ALLvsALL_data_processing.git ${projectDir}/bin/RnaChromATA
        pip install ${projectDir}/bin/RnaChromATA/
    }

    install_bitap() {
        git clone https://github.com/ikm4rkov/RawReadsProcessor.git ${projectDir}/bin/bitap_temp
        unzip ${projectDir}/bin/bitap_temp/main.zip -d ${projectDir}/bin/bitap_temp/
        g++ -o ${projectDir}/bin/bitap ${projectDir}/bin/bitap_temp/bitap.cpp
        mv ${projectDir}/bin/bitap_temp/bitap ${projectDir}/bin/
    }

    install_stereogene() {
        git clone https://github.com/favorov/stereogene.git 
        cd stereogene/src; make
        find . -type f -executable -print0 | xargs -0 -I {} mv {} ${projectDir}/bin/
    }



    if [ ! -f "${projectDir}/bin/RnaChromATA/setup.py" ]; then
        install_rnachrom
        touch binaries_ready.txt
    else
        echo "RnaChromATA is already installed."
    fi

    if [ ! -f "${projectDir}/bin/bitap" ]; then
        install_bitap
        touch binaries_ready.txt
    else
        echo "Bitap is already installed."
    fi

    if [ ! -f "${projectDir}/bin/stereogene" ]; then
        touch binaries_ready.txt
        install_stereogene
    else
        echo "Stereogene is already installed."
    fi
    """
    // stereogene installation



}

// pip install ${projectDir}/bin/RNAChrom_ALLvsALL_data_processing 

    // if (!rnachrom.exists()) {
    //     log.info print_cyan("Installing https://github.com/AndrewSigorskih/RNAChrom_ALLvsALL_data_processing.git ")
    //     """
    //     git clone https://github.com/AndrewSigorskih/RNAChrom_ALLvsALL_data_processing.git ${projectDir}/bin/RnaChromATA
    //     pip install ${projectDir}/bin/RnaChromATA/
    //     touch binaries_ready.txt
    //     """
    // } 
 
    // if (!bitap.exists()) {
    //     log.info print_cyan("Installing https://github.com/ikm4rkov/RawReadsProcessor.git ")
    //     """
    //     git clone https://github.com/ikm4rkov/RawReadsProcessor.git ${projectDir}/bin/bitap_temp
    //     unzip ${projectDir}/bin/bitap_temp/main.zip
    //     g++ -o bitap  ${projectDir}/bin/bitap_temp/bitap.cpp
    //     mv ${projectDir}/bin/bitap_temp/bitap  ${projectDir}/bin/
    //     """
    // } 

    // unzip ${projectDir}/bin/bitap_temp/main.zip -d ${projectDir}/bin/bitap_temp/