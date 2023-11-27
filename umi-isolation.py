#!/usr/bin/env python3

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

# Imports:
import os
import sys
import argparse
import re
import pandas as pd
import subprocess as sp


def create_output_files(cluster_directory, blast_file, tabular_file):
    """
    The create_output_files function:
        This function creates the tabular and BLAST output files.
    """
    output = pd.DataFrame(
        columns=["UMI ID", "UMI SEQ", "Read Count", "Centroid Read"]
    )
    count = 0
    for file_name in os.listdir(cluster_directory):
        umi_number = file_name.split("_")[0]
        umi_string = file_name.split("_")[1][:-6]
        line_count = 0
        with open(cluster_directory + file_name) as cluster_file:
            for line in cluster_file:
                line_count += 1
        with open(cluster_directory + file_name) as umi_file:
            if line_count == 2:
                for line in umi_file:
                    if line.startswith(">"):
                        header = line.split("=")[1].strip("\n")
                        read = next(umi_file)
                        output.loc[count] = [
                            umi_number,
                            umi_string,
                            header.strip("\n"),
                            read.strip("\n").upper(),
                        ]
                        with open(blast_file, "a") as output_file:
                            output_file.write(">" + umi_number + "\n")
                            output_file.write(read.strip("\n").upper() + "\n")
                    else:
                        pass
            elif line_count > 2:
                version_count = 1
                for line in umi_file:
                    if line.startswith(">"):
                        header = line.split("=")[1].strip("\n")
                        read = next(umi_file)
                        umi_version = umi_number + "." + str(version_count)
                        output.loc[count] = [
                            umi_version,
                            umi_string,
                            header.strip("\n"),
                            read.strip("\n").upper(),
                        ]
                        with open(blast_file, "a") as output_file:
                            output_file.write(">" + umi_version + "\n")
                            output_file.write(read.strip("\n").upper() + "\n")
                        version_count += 1
                        count += 1
                    else:
                        pass
            else:
                pass
        count += 1
    output = output.set_index("UMI ID")
    output.to_csv(tabular_file, sep="\t", encoding="utf-8")


def get_vsearch_cluster_size(zip, cluster_directory, identity_score):
    """
    The get_vsearch_cluster_size function:
        This function controls the VSEARCH clustering. Every fasta file created
        by get_vsearch_sort_by_size is clustered using VSEARCH. The expected
        result is a single centroid sequence. This is checked in the
        create_output_files function.
    """
    for file_name in os.listdir(zip):
        if file_name.startswith("sorted"):
            input_command = zip + file_name
            output_command = cluster_directory + file_name[11:]
            vsearch_cluster_command = sp.Popen(
                [
                    "vsearch",
                    "--cluster_size",
                    input_command,
                    "--fasta_width",
                    "0",
                    "--id",
                    identity_score,
                    "--sizein",
                    "--minseqlength",
                    "1",
                    "--centroids",
                    output_command,
                    "--sizeout",
                ],
                stdout=sp.PIPE,
                stderr=sp.PIPE,
            )
            out, error = vsearch_cluster_command.communicate()
        else:
            pass


def get_vsearch_sort_by_size(zip, minimal_size_abundance):
    """
    The get_vsearch_sort_by_size function:
        This function controls the VSEARCH sorting. Every fasta file created by
        get_vsearch_derep is sorted based on abundance. Any reads with a
        abundance lower than minimal_size_abundance will be discarded.
    """
    for file_name in os.listdir(zip):
        if file_name.startswith("derep"):
            input_command = zip + file_name
            output_command = zip + "sorted" + file_name
            vsearch_sort_command = sp.Popen(
                [
                    "vsearch",
                    "--sortbysize",
                    input_command,
                    "--output",
                    output_command,
                    "--minseqlength",
                    "1",
                    "--minsize",
                    minimal_size_abundance,
                ],
                stdout=sp.PIPE,
                stderr=sp.PIPE,
            )
            out, error = vsearch_sort_command.communicate()
        else:
            pass


def get_vsearch_derep(zip):
    """
    The get_vsearch_derep function:
        This function controls the VSEARCH dereplication. Every fasta file
        created by get_fasta_files is dereplicated. This step is necessary for
        the sorting step to work.
    """
    for file_name in os.listdir(zip):
        if file_name.endswith(".fasta"):
            input_command = zip + file_name
            output_command = zip + "derep" + file_name
            vsearch_derep_command = sp.Popen(
                [
                    "vsearch",
                    "--derep_fulllength",
                    input_command,
                    "--output",
                    output_command,
                    "--minseqlength",
                    "1",
                    "--sizeout",
                ],
                stdout=sp.PIPE,
                stderr=sp.PIPE,
            )

            out, error = vsearch_derep_command.communicate()
        else:
            pass


def get_fasta_files(header, read, umi_code, unique_umi_dictionary, zip):
    """
    The get_fasta_files function:
        This function creates separate fasta files for every unique UMI. The
        function creates a unique name for every UMI file and combines that
        with the desired output path. A file is opened or created based on this
        combination. The read header and the read itself are appended to it.
    """
    file_identifier = (
        "UMI#"
        + str(unique_umi_dictionary[umi_code])
        + "_"
        + umi_code
        + ".fasta"
    )
    file_name = zip + file_identifier
    with open(file_name, "a") as output_file:
        output_file.write(header)
        output_file.write(read)


def get_target_zero(read, umi_length, search_method, forward, reverse):
    """
    The get_target_zero function:
        This function checks if both the forward and reverse primer can be
        found, if that succeeds, the forward or reverse (when working with
        single UMIs) or forward and reverse (when working with double UMIs)
        nucleotides are isolated based on the length of the UMI. This isolation
        is done from the first nucleotide at the 5'-end of a read and the last
        nucleotide at the 3'-end of a read.
    """
    check_forward_primer = re.search(forward, read)
    if check_forward_primer != None:
        check_reverse_primer = re.search(reverse, read)
        if check_reverse_primer != None:
            if search_method == "umi5":
                return read[0 : int(umi_length)]
            elif search_method == "umidouble":
                return (read[0 : int(umi_length)], read[-int(umi_length) :])
            elif search_method == "umi3":
                return read[-int(umi_length) :]
            else:
                pass
        else:
            pass
    else:
        pass


# The getTargetFront function.
# This function searches for a regex string in the provided read. It will
# isolate either a forward or reverse UMI or double UMIs. The isolation is based on
# this read structure SCAFFOLDF-UMI-PRIMERF-PRODUCT-PRIMERR-UMI-SCAFFOLDR.
# When looking for the forward UMI, the last position of SCAFFOLDF is used, when
# looking for the reverse UMI, the first position of SCAFFOLDR is used, when
# looking for double UMIs both positions are used.
# The mentioned positions + or - the UMI length result in a UMI code.
# When not searching for both UMIs a check needs to be passed, this check
# makes sure the reverse scaffold (in the case of umi5) or the forward scaffold
# (in the case of umi3) are present.
def getTargetFront(read, umi_length, search_method, forward, reverse):
    if search_method == "umi5" or search_method == "umidouble":
        intPositionForward = re.search(forward, read).end()
        intPositionUmiForward = intPositionForward + int(umi_length)
        strUmiCodeForward = read[intPositionForward:intPositionUmiForward]
        if search_method == "umi5":
            check_reverse_primer = re.search(reverse, read)
            if check_reverse_primer != None:
                return strUmiCodeForward
            else:
                pass
        elif search_method == "umidouble":
            intPositionReverse = re.search(reverse, read).start()
            intPositionUmiReverse = intPositionReverse - int(umi_length)
            strUmiCodeReverse = read[intPositionUmiReverse:intPositionReverse]
            return strUmiCodeForward, strUmiCodeReverse
        else:
            pass
    elif search_method == "umi3":
        check_forward_primer = re.search(forward, read)
        if check_forward_primer != None:
            intPositionReverse = re.search(reverse, read).start()
            intPositionUmiReverse = intPositionReverse - int(umi_length)
            strUmiCodeReverse = read[intPositionUmiReverse:intPositionReverse]
            return strUmiCodeReverse
        else:
            pass
    else:
        pass


# The getTargetBehind function.
# This function searches for a regex string in the provided read. It will
# isolate either a forward or reverse UMI or double UMIs. The isolation is based on
# this read structure UMI-PRIMERF-PRODUCT-PRIMERR-UMI.
# When looking for the forward UMI, the first position of PRIMERF is used, when
# looking for the reverse UMI, the last position of PRIMERR is used, when
# looking for double UMIs both positions are used.
# The mentioned positions + or - the UMI length result in a UMI code.
# When not searching for both UMIs a check needs to be passed, this check
# makes sure the reverse primer (in the case of umi5) or the forward primer
# (in the case of umi3) are present.
def getTargetBehind(read, umi_length, search_method, forward, reverse):
    if search_method == "umi5" or search_method == "umidouble":
        intPositionForward = re.search(forward, read).start()
        intPositionUmiForward = intPositionForward - int(umi_length)
        strUmiCodeForward = read[intPositionUmiForward:intPositionForward]
        if search_method == "umi5":
            check_reverse_primer = re.search(reverse, read)
            if check_reverse_primer != None:
                return strUmiCodeForward
            else:
                pass
        elif search_method == "umidouble":
            intPositionReverse = re.search(reverse, read).end()
            intPositionUmiReverse = intPositionReverse + int(umi_length)
            strUmiCodeReverse = read[intPositionReverse:intPositionUmiReverse]
            return strUmiCodeForward, strUmiCodeReverse
        else:
            pass
    elif search_method == "umi3":
        check_forward_primer = re.search(forward, read)
        if check_forward_primer != None:
            intPositionReverse = re.search(reverse, read).end()
            intPositionUmiReverse = intPositionReverse + int(umi_length)
            strUmiCodeReverse = read[intPositionReverse:intPositionUmiReverse]
            return strUmiCodeReverse
        else:
            pass
    else:
        pass


# The getReverseComplement function.
# This function creates a complementary string using a sequence as input. The
# function loops through a list version of the sequence and checks and changes
# every character. It then returns a joined string.
def getReverseComplement(line):
    dicComplementCodes = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "M": "K",
        "R": "Y",
        "W": "W",
        "S": "S",
        "Y": "R",
        "K": "M",
        "V": "B",
        "H": "D",
        "D": "H",
        "B": "V",
        "N": "N",
    }
    lstLine = list(line)
    for intPosition in range(len(lstLine)):
        lstLine[intPosition] = dicComplementCodes[lstLine[intPosition]]
    return "".join(lstLine)


# The getRegex function.
# This function creates a regex string using a sequence as input. This regex
# string is based on the IUPAC ambiguity codes. The function loops through
# a list version of the sequence and checks per character if it is a
# ambiguous character, and changes it into regex code if so. It then returns
# a joined string.
def getRegex(line):
    dicAmbiguityCodes = {
        "M": "[AC]",
        "R": "[AG]",
        "W": "[AT]",
        "S": "[CG]",
        "Y": "[CT]",
        "K": "[GT]",
        "V": "[ACG]",
        "H": "[ACT]",
        "D": "[AGT]",
        "B": "[CGT]",
        "N": "[GATC]",
    }
    lstLine = list(line)
    for intPosition in range(len(lstLine)):
        if (
            lstLine[intPosition] != "A"
            and lstLine[intPosition] != "T"
            and lstLine[intPosition] != "G"
            and lstLine[intPosition] != "C"
        ):
            lstLine[intPosition] = dicAmbiguityCodes[lstLine[intPosition]]
        else:
            pass
    return "".join(lstLine)


# The getUmiCode function.
# This function controls the umi searching process. It first uses the functions
# getRegex and getReverseComplement to create regex strings of both the
# forward primer/scaffold and the reverse complement primer/scaffold. These regex
# strings are then directed to the desired functions, this depends on the
# search method choice [primer/scaffold/zero].
def getUmiCode(read, strProcess, umi_length, search_method, forward, reverse):
    read = read.strip("\n")
    strRegexForward = getRegex(forward)
    strRegexComplementReverse = getRegex(getReverseComplement(reverse[::-1]))
    if strProcess == "primer":
        try:
            return getTargetBehind(
                read,
                umi_length,
                search_method,
                strRegexForward,
                strRegexComplementReverse,
            )
        except AttributeError:
            pass
    elif strProcess == "scaffold":
        try:
            return getTargetFront(
                read,
                umi_length,
                search_method,
                strRegexForward,
                strRegexComplementReverse,
            )
        except AttributeError:
            pass
    elif strProcess == "zero":
        try:
            return get_target_zero(
                read,
                umi_length,
                search_method,
                strRegexForward,
                strRegexComplementReverse,
            )
        except AttributeError:
            pass
    else:
        pass


# The getUmiCollection function.
# This function opens the input file and loops through it. It isolates the
# read headers and reads. For every read the getUmiCode function is called
# which outputs one or two UMIs. In the case of a double UMI search [umidouble]
# the two UMIs are put together. For every read that contains a UMI the
# get_fasta_files function is called. After all reads have been processed, the
# getVSEARCH and create_output_files functions are called.
def getUmiCollection(
    flInput,
    cluster_directory,
    tabular_file,
    zip,
    blast_file,
    strProcess,
    umi_length,
    search_method,
    forward,
    reverse,
    strFormat,
    strOperand,
    identity_score,
    minimal_size_abundance,
):
    unique_umi_dictionary = {}
    intNoUmiInReadCount = 0
    intUniqueUmi = 1
    with open(flInput) as oisInput:
        for line in oisInput:
            if (
                line[0] == strOperand
                and bool(re.match("[A-Za-z0-9]", line[1])) == True
            ):
                header = line
                read = next(oisInput)
                try:
                    strUmiCode = getUmiCode(
                        read.upper(),
                        strProcess,
                        umi_length,
                        search_method,
                        forward.upper(),
                        reverse.upper(),
                    )
                except UnboundLocalError:
                    intNoUmiInReadCount += 1
                try:
                    if strUmiCode != None:
                        if search_method == "umi5" or search_method == "umi3":
                            umi_code = strUmiCode
                        elif search_method == "umidouble":
                            umi_code = strUmiCode[0] + strUmiCode[1]
                        else:
                            pass
                        if umi_code not in unique_umi_dictionary:
                            unique_umi_dictionary[umi_code] = intUniqueUmi
                            intUniqueUmi += 1
                        else:
                            pass
                    else:
                        pass
                except UnboundLocalError:
                    pass
                try:
                    if umi_code != None:
                        get_fasta_files(
                            header, read, umi_code, unique_umi_dictionary, zip
                        )
                    else:
                        pass
                except UnboundLocalError:
                    pass
            strUmiCode = None
            umi_code = None
    get_vsearch_derep(zip)
    get_vsearch_sort_by_size(zip, minimal_size_abundance)
    get_vsearch_cluster_size(zip, cluster_directory, identity_score)
    create_output_files(cluster_directory, blast_file, tabular_file)


# The setFormat function.
# This function specifies the first character of the read headers based on the
# input file format. It then calls the getUmiCollection function.
def setFormat(
    flInput,
    cluster_directory,
    tabular_file,
    zip,
    blast_file,
    strProcess,
    strFormat,
    umi_length,
    search_method,
    forward,
    reverse,
    identity_score,
    minimal_size_abundance,
):
    if strFormat == "fasta":
        strOperand = ">"
    elif strFormat == "fastq":
        strOperand = "@"
    else:
        pass
    getUmiCollection(
        flInput,
        cluster_directory,
        tabular_file,
        zip,
        blast_file,
        strProcess,
        umi_length,
        search_method,
        forward,
        reverse,
        strFormat,
        strOperand,
        str(identity_score),
        str(minimal_size_abundance),
    )


# The argvs function.
def parseArgvs():
    parser = argparse.ArgumentParser(
        description="Use a python script to\
                                                  accumulate all UMIs and\
                                                  output a tabular file, a\
                                                  BLAST file and a zip file."
    )
    parser.add_argument("-v", action="version", version="%(prog)s [0.1.0]")
    parser.add_argument(
        "-i",
        action="store",
        dest="fisInput",
        help="The location of the input file(s)",
    )
    parser.add_argument(
        "-c",
        action="store",
        dest="fosClusterDirectory",
        help="The location of the clustering output file(s)",
    )
    parser.add_argument(
        "-o",
        action="store",
        dest="fosOutput",
        help="The location of the tabular output file(s)",
    )
    parser.add_argument(
        "-z",
        action="store",
        dest="fosOutputZip",
        help="The location of the pre-vsearch zip output file(s)",
    )
    parser.add_argument(
        "-q",
        action="store",
        dest="fosBlastFile",
        help="The location of the BLAST output file(s)",
    )
    parser.add_argument(
        "-p",
        action="store",
        dest="disProcess",
        help="The UMI search approach [primer/scaffold(adapter)/zero]",
    )
    parser.add_argument(
        "-f",
        action="store",
        dest="disFormat",
        help="The format of the input file(s) [fasta/fastq]",
    )
    parser.add_argument(
        "-l",
        action="store",
        dest="disUmiLength",
        help="The length of the UMI sequences",
    )
    parser.add_argument(
        "-s",
        action="store",
        dest="disSearch",
        help="Search UMIs at 5'-end [umi5], 3'-end [umi3] or at\
                              5'-end and 3'-end [umidouble]",
    )
    parser.add_argument(
        "-a",
        action="store",
        dest="disForward",
        help="The 5'-end search nucleotides",
    )
    parser.add_argument(
        "-b",
        action="store",
        dest="disReverse",
        help="The 3'-end search nucleotides",
    )
    parser.add_argument(
        "-d",
        action="store",
        dest="disIdentity",
        help="The identity percentage with which to perform the\
                              final VSEARCH check",
    )
    parser.add_argument(
        "-u",
        action="store",
        dest="disAbundance",
        help="The minimum abundance a read has to be present in\
                              order to be part of the final VSEARCH check",
    )
    argvs = parser.parse_args()
    return argvs


def main():
    """
    The main function:
        This function handles the arguments parsed to the script and calls
        the first function setFormat.
    """
    argvs = parseArgvs()
    setFormat(
        argvs.fisInput,
        argvs.fosClusterDirectory,
        argvs.fosOutput,
        argvs.fosOutputZip,
        argvs.fosBlastFile,
        argvs.disProcess,
        argvs.disFormat,
        argvs.disUmiLength,
        argvs.disSearch,
        argvs.disForward,
        argvs.disReverse,
        argvs.disIdentity,
        argvs.disAbundance,
    )


if __name__ == "__main__":
    main()

# Additional information:
# =======================
# Prequisites:
#     sudo apt-get install python3
#     sudo apt-get install python3-pip
#     sudo pip3 install pandas
#     sudo apt-get install libargtable2-dev
#     Download VSEARCH from GitHub
#     Unpack downloaded file
#     sudo ./autogen.sh
#     sudo ./configure && sudo make && sudo make install
