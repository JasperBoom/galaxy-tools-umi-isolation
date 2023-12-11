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
import argparse
import re
import pandas as pd
import subprocess as sp


def create_output_files(cluster_directory, output_blast_file, tabular_file):
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
                        with open(output_blast_file, "a") as output_file:
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
                        with open(output_blast_file, "a") as output_file:
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


def get_vsearch_cluster_size(zip_file, cluster_directory, identity_score):
    """
    The get_vsearch_cluster_size function:
        This function controls the VSEARCH clustering. Every fasta file created
        by get_vsearch_sort_by_size is clustered using VSEARCH. The expected
        result is a single centroid sequence. This is checked in the
        create_output_files function.
    """
    for file_name in os.listdir(zip_file):
        if file_name.startswith("sorted"):
            input_command = zip_file + file_name
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


def get_vsearch_sort_by_size(zip_file, minimal_size_abundance):
    """
    The get_vsearch_sort_by_size function:
        This function controls the VSEARCH sorting. Every fasta file created by
        get_vsearch_derep is sorted based on abundance. Any reads with a
        abundance lower than minimal_size_abundance will be discarded.
    """
    for file_name in os.listdir(zip_file):
        if file_name.startswith("derep"):
            input_command = zip_file + file_name
            output_command = zip_file + "sorted" + file_name
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


def get_vsearch_derep(zip_file):
    """
    The get_vsearch_derep function:
        This function controls the VSEARCH dereplication. Every fasta file
        created by get_fasta_files is dereplicated. This step is necessary for
        the sorting step to work.
    """
    for file_name in os.listdir(zip_file):
        if file_name.endswith(".fasta"):
            input_command = zip_file + file_name
            output_command = zip_file + "derep" + file_name
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


def get_fasta_files(header, read, umi_code, unique_umi_dictionary, zip_file):
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
    file_name = zip_file + file_identifier
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


def get_target_front(read, umi_length, search_method, forward, reverse):
    """
    The get_target_front function:
        This function searches for a regex string in the provided read. It will
        isolate either a forward or reverse UMI or double UMIs. The isolation
        is based on this read structure
        SCAFFOLDF-UMI-PRIMERF-PRODUCT-PRIMERR-UMI-SCAFFOLDR. When looking for
        the forward UMI, the last position of SCAFFOLDF is used, when looking
        for the reverse UMI, the first position of SCAFFOLDR is used, when
        looking for double UMIs both positions are used. The mentioned
        positions + or - the UMI length result in a UMI code. When not
        searching for both UMIs a check needs to be passed, this check makes
        sure the reverse scaffold (in the case of umi5) or the forward
        scaffold (in the case of umi3) are present.
    """
    if search_method == "umi5" or search_method == "umidouble":
        forward_position = re.search(forward, read).end()
        umi_forward_position = forward_position + int(umi_length)
        forward_umi_code = read[forward_position:umi_forward_position]
        if search_method == "umi5":
            check_reverse_primer = re.search(reverse, read)
            if check_reverse_primer != None:
                return forward_umi_code
            else:
                pass
        elif search_method == "umidouble":
            reverse_position = re.search(reverse, read).start()
            umi_reverse_position = reverse_position - int(umi_length)
            reverse_umi_code = read[umi_reverse_position:reverse_position]
            return forward_umi_code, reverse_umi_code
        else:
            pass
    elif search_method == "umi3":
        check_forward_primer = re.search(forward, read)
        if check_forward_primer != None:
            reverse_position = re.search(reverse, read).start()
            umi_reverse_position = reverse_position - int(umi_length)
            reverse_umi_code = read[umi_reverse_position:reverse_position]
            return reverse_umi_code
        else:
            pass
    else:
        pass


def get_target_behind(read, umi_length, search_method, forward, reverse):
    """
    The get_target_behind function:
        This function searches for a regex string in the provided read. It will
        isolate either a forward or reverse UMI or double UMIs. The isolation
        is based on this read structure UMI-PRIMERF-PRODUCT-PRIMERR-UMI. When
        looking for the forward UMI, the first position of PRIMERF is used,
        when looking for the reverse UMI, the last position of PRIMERR is used,
        when looking for double UMIs both positions are used. The mentioned
        positions + or - the UMI length result in a UMI code. When not
        searching for both UMIs a check needs to be passed, this check makes
        sure the reverse primer (in the case of umi5) or the forward primer
        (in the case of umi3) are present.
    """
    if search_method == "umi5" or search_method == "umidouble":
        forward_position = re.search(forward, read).start()
        umi_forward_position = forward_position - int(umi_length)
        forward_umi_code = read[umi_forward_position:forward_position]
        if search_method == "umi5":
            check_reverse_primer = re.search(reverse, read)
            if check_reverse_primer != None:
                return forward_umi_code
            else:
                pass
        elif search_method == "umidouble":
            reverse_position = re.search(reverse, read).end()
            umi_reverse_position = reverse_position + int(umi_length)
            reverse_umi_code = read[reverse_position:umi_reverse_position]
            return forward_umi_code, reverse_umi_code
        else:
            pass
    elif search_method == "umi3":
        check_forward_primer = re.search(forward, read)
        if check_forward_primer != None:
            reverse_position = re.search(reverse, read).end()
            umi_reverse_position = reverse_position + int(umi_length)
            reverse_umi_code = read[reverse_position:umi_reverse_position]
            return reverse_umi_code
        else:
            pass
    else:
        pass


def create_reverse_complement(line):
    """
    The create_reverse_complement function:
        This function creates a complementary string using a sequence as input.
        The function loops through a list version of the sequence and checks
        and changes every character. It then returns a joined string.
    """
    complement_codes_dictionary = {
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
    line_list = list(line)
    for position in range(len(line_list)):
        line_list[position] = complement_codes_dictionary[line_list[position]]
    return "".join(line_list)


def generate_regex(line):
    """
    The generate_regex function:
        This function creates a regex string using a sequence as input. This
        regex string is based on the IUPAC ambiguity codes. The function loops
        through a list version of the sequence and checks per character if it
        is a ambiguous character, and changes it into regex code if so. It then
        returns a joined string.
    """
    ambiguity_codes = {
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
    line_list = list(line)
    for position in range(len(line_list)):
        if (
            line_list[position] != "A"
            and line_list[position] != "T"
            and line_list[position] != "G"
            and line_list[position] != "C"
        ):
            line_list[position] = ambiguity_codes[line_list[position]]
        else:
            pass
    return "".join(line_list)


def get_umi_code(read, process, umi_length, search_method, forward, reverse):
    """
    The get_umi_code function:
        This function controls the umi searching process. It first uses the
        functions generate_regex and create_reverse_complement to create regex
        strings of both the forward primer/scaffold and the reverse complement
        primer/scaffold. These regex strings are then directed to the desired
        functions, this depends on the search method
        choice [primer/scaffold/zero].
    """
    read = read.strip("\n")
    forward_regex = generate_regex(forward)
    reverse_complement_regex = generate_regex(
        create_reverse_complement(reverse[::-1])
    )
    if process == "primer":
        try:
            return get_target_behind(
                read,
                umi_length,
                search_method,
                forward_regex,
                reverse_complement_regex,
            )
        except AttributeError:
            pass
    elif process == "scaffold":
        try:
            return get_target_front(
                read,
                umi_length,
                search_method,
                forward_regex,
                reverse_complement_regex,
            )
        except AttributeError:
            pass
    elif process == "zero":
        try:
            return get_target_zero(
                read,
                umi_length,
                search_method,
                forward_regex,
                reverse_complement_regex,
            )
        except AttributeError:
            pass
    else:
        pass


def get_umi_collection(
    input_file,
    cluster_directory,
    tabular_file,
    zip_file,
    output_blast_file,
    process,
    umi_length,
    search_method,
    forward,
    reverse,
    format_string,
    operand,
    identity_score,
    minimal_size_abundance,
):
    """
    The get_umi_collection function:
        This function opens the input file and loops through it. It isolates the
        read headers and reads. For every read the get_umi_code function is
        called which outputs one or two UMIs. In the case of a double UMI
        search [umidouble] the two UMIs are put together. For every read that
        contains a UMI the get_fasta_files function is called. After all reads
        have been processed, the VSEARCH and create_output_files functions
        are called.
    """
    unique_umi_dictionary = {}
    count_reads_without_umi = 0
    count_unique_umis = 1
    with open(input_file) as input:
        for line in input:
            if (
                line[0] == operand
                and bool(re.match("[A-Za-z0-9]", line[1])) == True
            ):
                header = line
                read = next(input)
                try:
                    umi_code = get_umi_code(
                        read.upper(),
                        process,
                        umi_length,
                        search_method,
                        forward.upper(),
                        reverse.upper(),
                    )
                except UnboundLocalError:
                    count_reads_without_umi += 1
                try:
                    if umi_code != None:
                        if search_method == "umi5" or search_method == "umi3":
                            umi_code = umi_code
                        elif search_method == "umidouble":
                            umi_code = umi_code[0] + umi_code[1]
                        else:
                            pass
                        if umi_code not in unique_umi_dictionary:
                            unique_umi_dictionary[umi_code] = count_unique_umis
                            count_unique_umis += 1
                        else:
                            pass
                    else:
                        pass
                except UnboundLocalError:
                    pass
                try:
                    if umi_code != None:
                        get_fasta_files(
                            header,
                            read,
                            umi_code,
                            unique_umi_dictionary,
                            zip_file,
                        )
                    else:
                        pass
                except UnboundLocalError:
                    pass
            umi_code = None
            umi_code = None
    get_vsearch_derep(zip_file)
    get_vsearch_sort_by_size(zip_file, minimal_size_abundance)
    get_vsearch_cluster_size(zip_file, cluster_directory, identity_score)
    create_output_files(cluster_directory, output_blast_file, tabular_file)


def set_format(
    input_file,
    cluster_directory,
    tabular_file,
    zip_file,
    output_blast_file,
    process,
    format_string,
    umi_length,
    search_method,
    forward,
    reverse,
    identity_score,
    minimal_size_abundance,
):
    """
    The set_format function:
        This function specifies the first character of the read headers based
        on the input file format. It then calls the get_umi_collection function.
    """
    if format_string == "fasta":
        operand = ">"
    elif format_string == "fastq":
        operand = "@"
    else:
        pass
    get_umi_collection(
        input_file,
        cluster_directory,
        tabular_file,
        zip_file,
        output_blast_file,
        process,
        umi_length,
        search_method,
        forward,
        reverse,
        format_string,
        operand,
        str(identity_score),
        str(minimal_size_abundance),
    )


def parse_argvs():
    """
    The parse_argvs function:
        This function handles all positional arguments that the script accepts,
        including version and help pages.
    """
    parser = argparse.ArgumentParser(
        description="Use a python script to accumulate all UMIs and output a\
                     tabular file, a BLAST file and a zip file."
    )
    parser.add_argument("-v", action="version", version="%(prog)s [0.1.0]")
    parser.add_argument(
        "-i",
        action="store",
        dest="input_file",
        help="The location of the input file(s)",
    )
    parser.add_argument(
        "-c",
        action="store",
        dest="cluster_directory",
        help="The location of the clustering output file(s)",
    )
    parser.add_argument(
        "-o",
        action="store",
        dest="output_tabular_file",
        help="The location of the tabular output file(s)",
    )
    parser.add_argument(
        "-z",
        action="store",
        dest="output_zip_file",
        help="The location of the pre-vsearch zip output file(s)",
    )
    parser.add_argument(
        "-q",
        action="store",
        dest="output_blast_file",
        help="The location of the BLAST output file(s)",
    )
    parser.add_argument(
        "-p",
        action="store",
        dest="process",
        help="The UMI search approach [primer/scaffold(adapter)/zero]",
    )
    parser.add_argument(
        "-f",
        action="store",
        dest="format",
        help="The format of the input file(s) [fasta/fastq]",
    )
    parser.add_argument(
        "-l",
        action="store",
        dest="umi_length",
        help="The length of the UMI sequences",
    )
    parser.add_argument(
        "-s",
        action="store",
        dest="search_method",
        help="Search UMIs at 5'-end [umi5], 3'-end [umi3] or at 5'-end and\
              3'-end [umidouble]",
    )
    parser.add_argument(
        "-a",
        action="store",
        dest="forward",
        help="The 5'-end search nucleotides",
    )
    parser.add_argument(
        "-b",
        action="store",
        dest="reverse",
        help="The 3'-end search nucleotides",
    )
    parser.add_argument(
        "-d",
        action="store",
        dest="identity_score",
        help="The identity percentage with which to perform the final VSEARCH\
              check",
    )
    parser.add_argument(
        "-u",
        action="store",
        dest="abundance",
        help="The minimum abundance a read has to be present in order to be\
              part of the final VSEARCH check",
    )
    argvs = parser.parse_args()
    return argvs


def main():
    """
    The main function:
        This function handles the arguments parsed to the script and calls
        the first function set_format.
    """
    argvs = parseArgvs()
    set_format(
        argvs.input_file,
        argvs.cluster_directory,
        argvs.output_tabular_file,
        argvs.output_zip_file,
        argvs.output_blast_file,
        argvs.process,
        argvs.format,
        argvs.umi_length,
        argvs.search_method,
        argvs.forward,
        argvs.reverse,
        argvs.identity_score,
        argvs.abundance,
    )


if __name__ == "__main__":
    main()

# Additional information:
# =======================
# Prequisites:
# * sudo apt-get install python3
# * sudo apt-get install python3-pip
# * sudo pip3 install pandas
# * sudo apt-get install libargtable2-dev
# * Download VSEARCH from GitHub
# * Unpack downloaded file
# * sudo ./autogen.sh
# * sudo ./configure && sudo make && sudo make install
