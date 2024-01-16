#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# Naturalis internship repository for UMI isolation tool.
# Copyright (C) 2019 Jasper Boom

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Contact information: info@jboom.org.
# -----------------------------------------------------------------------------

format_flow() {
    # The format_flow function.
    #     This function creates two temporary storage directories in the output
    #     directory. It then calls the umi-isolation.py script with the correct
    #     input values. The output tabular file is copied to the expected
    #     location and removed. The output blast file is copies to the expected
    #     location and removed. A zip file is created of the generated
    #     prevsearch fasta file and copied to the expected location. After that
    #     the temporary storage directories are removed.
    script_directory=$(dirname "$(readlink -f "$0")")
    directory_name=$(mktemp -d /media/GalaxyData/database/files/XXXXXX)
    mkdir \
        -p "${directory_name}_temp"
    mkdir \
        -p "${directory_name}_cluster_check"
    python3 \
        $script_directory"/umi-isolation.py" \
            -i ${input_file} \
            -o ${directory_name}_temp/csv_temp_file.csv \
            -z ${directory_name}_temp/ \
            -q ${directory_name}_temp/blast_temp_file.fasta \
            -f ${format} -p ${process} \
            -l ${umi_length} -s ${search_method} \
            -a ${forward} -b ${reverse} \
            -c ${directory_name}_cluster_check/ \
            -d ${identity_score} \
            -u ${abundance}
    cat ${directory_name}_temp/csv_temp_file.csv \
        > ${output_tabular_file}
    rm ${directory_name}_temp/csv_temp_file.csv
    cat ${directory_name}_temp/blast_temp_file.fasta \
        > ${output_blast_file}
    rm ${directory_name}_temp/blast_temp_file.fasta
    find ${directory_name}_temp/ \
        -name "derep*" \
        -delete
    find ${directory_name}_temp/ \
        -name "sorted*" \
        -delete
    find ${directory_name}_temp/ \
        -name "UMI#*" \
        -print \
        | zip \
              -jqr ${directory_name}_temp/zip_temp_file.zip -@
    cat ${directory_name}_temp/zip_temp_file.zip \
        > ${output_zip_file}
    rm \
        -rf ${directory_name}_temp
    rm \
        -rf ${directory_name}_cluster_check
}

main() {
    # The main function.
    #     This function calls all processing functions in correct order.
    format_flow
}

# The getopts function.
# https://kodekloud.com/blog/bash-getopts/
OPT_STRING=":i:o:z:q:p:f:l:s:a:b:d:u:vh"
while getopts ${OPT_STRING} option;
do
    case ${option} in
        i)
            input_file=${OPTARG}
            ;;
        o)
            output_tabular_file=${OPTARG}
            ;;
        z)
            output_zip_file=${OPTARG}
            ;;
        q)  
            output_blast_file=${OPTARG}
            ;;
        p)
            process=${OPTARG}
            ;;
        f)
            format=${OPTARG}
            ;;
        l)
            umi_length=${OPTARG}
            ;;
        s)
            search_method=${OPTARG}
            ;;
        a)
            forward=${OPTARG}
            ;;
        b)
            reverse=${OPTARG}
            ;;
        d)
            identity_score=${OPTARG}
            ;;
        u)
            abundance=${OPTARG}
            ;;
        v)
            echo ""
            echo "umi-isolation.sh [0.1.0]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: umi-isolation.sh [-h] [-v] [-i INPUT] [-o TABULAR]"
            echo "                        [-z ZIP] [-q BLAST] [-p PROCESS]"
            echo "                        [-f FORMAT] [-l LENGTH] [-s SEARCH]"
            echo "                        [-a FORWARD] [-b REVERSE]"
            echo ""
            echo "Optional arguments:"
            echo " -h          Show this help page and exit"
            echo " -v          Show the software's version number and exit"
            echo " -i          The location of the input file(s)"
            echo " -o          The location of the tabular output file(s)"
            echo " -z          The location of the pre-vsearch zip output"
            echo "             file(s)"
            echo " -q          The location of the blast output file(s)"
            echo " -p          The umi search approach"
            echo "             [primer/scaffold(adapter)/zero]"
            echo " -f          The format of the input file(s) [fasta/fastq]"
            echo " -l          The length of the UMI sequences"
            echo " -s          Search umis at 5'-end [umi5], 3'-end [umi3] or"
            echo "             at 5'-end and 3'-end [umidouble]"
            echo " -a          The 5'-end search nucleotides"
            echo " -b          The 3'-end search nucleotides"
            echo " -d          The identity percentage with which to perform"
            echo "             the final vsearch check"
            echo " -u          The minimum abundance a read has to be order to"
            echo "             be part of the final vsearch check present in"
            echo ""
            echo "Use a python script to accumulate all umis and output a"
            echo "tabular file, a blast file and a zip file. The tabular file"
            echo "will contain all unique UMI nucleotides, a count of the"
            echo "number of reads that are associated with that umi and a"
            echo "unique identifier for every umi."
            echo "The blast file can be used to identify all umi clusters."
            echo "The zip file will contain fasta files for every unique umi"
            echo "and contain reads associated with that umi, this zip file"
            echo "is created before vsearch is used for a final check."
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
# Files in fastq format should always have a .fasta extension.
# Files in fastq format should always have a .fastq extension.
# Every read in a fasta/fastq file should be on one line. For instance they can
# not be "human readable" and have a \n after every 80 characters.
# Prequisites:
# - sudo apt-get install python3
# - sudo apt-get install python3-pip
# - sudo pip3 install pandas
# - sudo apt-get install libargtable2-dev
# - Download vsearch from GitHub
# - Unpack downloaded file
# - sudo ./autogen.sh
# - sudo ./configure && sudo make && sudo make install