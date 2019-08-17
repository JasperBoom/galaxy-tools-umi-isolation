#!/usr/bin/env bash
# Copyright ©2007 Free Software Foundation, Inc.

# Author: Jasper Boom

# Prequisites:
# - sudo apt-get install python3
# - sudo apt-get install python3-pip
# - sudo pip3 install pandas
# - sudo apt-get install libargtable2-dev
# - Download VSEARCH from GitHub
# - Unpack downloaded file
# - sudo ./autogen.sh
# - sudo ./configure && sudo make && sudo make install

# The getFormatFlow function.
# This function creates two temporary storage directories in the output directory.
# It then calls the getUmiIsolation.py script with the correct input values.
# The output tabular file is copied to the expected location and removed.
# The output BLAST file is copies to the expected lacation and removed.
# A zip file is created of the generated preVSEARCH fastA file and copied to
# the expected location. After that the temporary storage directories are removed.
getFormatFlow() {
    strScriptDir=$(dirname "$(readlink -f "$0")")
    strDirectory=$(mktemp -d /media/GalaxyData/database/files/XXXXXX)
    mkdir -p "${strDirectory}_temp"
    mkdir -p "${strDirectory}_clusterCheck"
    python3 $strScriptDir"/getUmiIsolation.py" -i ${fisInput} \
                                               -o ${strDirectory}_temp/flTempCsv.csv \
                                               -z ${strDirectory}_temp/ \
                                               -q ${strDirectory}_temp/flTempBlast.fasta \
                                               -f ${disFormat} -p ${disProcess} \
                                               -l ${disUmiLength} -s ${disSearch} \
                                               -a ${disForward} -b ${disReverse} \
                                               -c ${strDirectory}_clusterCheck/ \
                                               -d ${disIdentity} \
                                               -u ${disAbundance}
    cat ${strDirectory}_temp/flTempCsv.csv > ${fosOutputTabular}
    rm ${strDirectory}_temp/flTempCsv.csv
    cat ${strDirectory}_temp/flTempBlast.fasta > ${fosBlastFile}
    rm ${strDirectory}_temp/flTempBlast.fasta
    find ${strDirectory}_temp/ -name "derep*" -delete
    find ${strDirectory}_temp/ -name "sorted*" -delete
    find ${strDirectory}_temp/ -name "UMI#*" -print | zip -jqr ${strDirectory}_temp/flTempZip.zip -@
    cat ${strDirectory}_temp/flTempZip.zip > ${fosOutputZip}
    rm -rf ${strDirectory}_temp
    rm -rf ${strDirectory}_clusterCheck
}

# The main function.
main() {
    getFormatFlow
}

# The getopts function.
while getopts ":i:o:z:q:p:f:l:s:a:b:d:u:vh" opt; do
    case ${opt} in
        i)
            fisInput=${OPTARG}
            ;;
        o)
            fosOutputTabular=${OPTARG}
            ;;
        z)
            fosOutputZip=${OPTARG}
            ;;
        q)  
            fosBlastFile=${OPTARG}
            ;;
        p)
            disProcess=${OPTARG}
            ;;
        f)
            disFormat=${OPTARG}
            ;;
        l)
            disUmiLength=${OPTARG}
            ;;
        s)
            disSearch=${OPTARG}
            ;;
        a)
            disForward=${OPTARG}
            ;;
        b)
            disReverse=${OPTARG}
            ;;
        d)
            disIdentity=${OPTARG}
            ;;
        u)
            disAbundance=${OPTARG}
            ;;
        v)
            echo ""
            echo "getUmiIsolation.sh [0.1.0]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: getUmiIsolation.sh [-h] [-v] [-i INPUT] [-o TABULAR]"
            echo "                          [-z ZIP] [-q BLAST] [-p PROCESS]"
            echo "                          [-f FORMAT] [-l LENGTH] [-s SEARCH]"
            echo "                          [-a FORWARD] [-b REVERSE]"
            echo ""
            echo "Optional arguments:"
            echo " -h                    Show this help page and exit"
            echo " -v                    Show the software's version number"
            echo "                       and exit"
            echo " -i                    The location of the input file(s)"
            echo " -o                    The location of the tabular output file(s)"
            echo " -z                    The location of the pre-vsearch zip output file(s)"
            echo " -q                    The location of the BLAST output file(s)"
            echo " -p                    The UMI search approach"
            echo "                       [primer/scaffold(adapter)/zero]"
            echo " -f                    The format of the input file(s)"
            echo "                       [fasta/fastq]"
            echo " -l                    The length of the UMI sequences"
            echo " -s                    Search UMIs at 5'-end [umi5],"
            echo "                       3'-end [umi3] or at 5'-end and" 
            echo "                       3'-end [umidouble]"
            echo " -a                    The 5'-end search nucleotides"
            echo " -b                    The 3'-end search nucleotides"
            echo " -d                    The identity percentage with which to"
            echo "                       perform the final VSEARCH check"
            echo " -u                    The minimum abundance a read has to be"
            echo "                       order to be part of the final VSEARCH"
            echo "                       check present in"
            echo ""
            echo "Use a python script to accumulate all UMIs and output a"
            echo "tabular file, a BLAST file and a zip file. The tabular file"
            echo "will contain all unique UMI nucleotides, a count of the"
            echo "number of reads that are associated with that UMI and a"
            echo "unique identifier for every UMI."
            echo "The BLAST file can be used to identify all UMI clusters."
            echo "The zip file will contain fastA files for every unique UMI"
            echo "and contain reads associated with that UMI, this zip file"
            echo "is created before VSEARCH is used for a final check."
            echo ""

            exit
            ;;
        \?)
            echo ""
            echo "You've entered an invalid option: -${OPTARG}."
            echo "Please use the -h option for correct formatting information."
            echo ""

            exit
            ;;
        :)
            echo ""
            echo "You've entered an invalid option: -${OPTARG}."
            echo "Please use the -h option for correct formatting information."
            echo ""

            exit
            ;;
    esac
done

main

# Additional information:
# =======================
#
# Files in fastA format should always have a .fasta extension.
# Files in fastQ format should always have a .fastq extension.
# Every read in a fastA/fastQ file should be on one line. For instance they can
# not be "human readable" and have a \n after every 80 characters.
