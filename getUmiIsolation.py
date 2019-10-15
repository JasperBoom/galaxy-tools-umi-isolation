#!/usr/bin/env python3
# Copyright ©2001-2019 Python Software Foundation

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

# Imports
import os
import sys
import argparse
import re
import pandas as pd
import subprocess as sp

# The setOutputFiles function.
# This function creates the tabular and BLAST output files.
def setOutputFiles(strClusterDir, flBlast, flTabular):
    dfOutput = pd.DataFrame(columns=["UMI ID", "UMI SEQ", "Read Count",
                                     "Centroid Read"])
    intCount = 0
    for strFileName in os.listdir(strClusterDir):
        strUmiNumber = strFileName.split("_")[0]
        strUmiString = strFileName.split("_")[1][:-6]
        intLineCount = 0
        with open(strClusterDir + strFileName) as oisClusterFile:
            for strLine in oisClusterFile:
                intLineCount += 1
        with open(strClusterDir + strFileName) as oisUmiFile:
            if intLineCount == 2:
                for strLine in oisUmiFile:
                    if strLine.startswith(">"):
                        strHeader = strLine.split("=")[1].strip("\n")
                        strRead = next(oisUmiFile)
                        dfOutput.loc[intCount] = [strUmiNumber, strUmiString,
                                                  strHeader.strip("\n"),
                                                  strRead.strip("\n").upper()]
                        with open(flBlast, "a") as flOutput:
                            flOutput.write(">" + strUmiNumber + "\n")
                            flOutput.write(strRead.strip("\n").upper() + "\n")
                    else:
                        pass
            elif intLineCount > 2:
                intVersionCount = 1
                for strLine in oisUmiFile:
                    if strLine.startswith(">"):
                        strHeader = strLine.split("=")[1].strip("\n")
                        strRead = next(oisUmiFile)
                        strUmiVersion = strUmiNumber + "." + str(intVersionCount)
                        dfOutput.loc[intCount] = [strUmiVersion, strUmiString,
                                                  strHeader.strip("\n"),
                                                  strRead.strip("\n").upper()]
                        with open(flBlast, "a") as flOutput:
                            flOutput.write(">" + strUmiVersion + "\n")
                            flOutput.write(strRead.strip("\n").upper() + "\n")
                        intVersionCount += 1
                        intCount += 1
                    else:
                        pass
            else:
                pass
        intCount += 1
    dfOutput = dfOutput.set_index("UMI ID")
    dfOutput.to_csv(flTabular, sep="\t", encoding="utf-8")

# The getVSEARCHclusterSize function.
# This function controls the VSEARCH clustering. Every fasta file created by
# getVSEARCHsortBySize is clustered using VSEARCH. The expected result is a single
# centroid sequence. This is checked in the setOutputFiles function.
def getVSEARCHclusterSize(flZip, strClusterDir, strIdentity):
    for strFileName in os.listdir(flZip):
        if strFileName.startswith("sorted"):
            strInputCommand = flZip + strFileName
            strOutputCommand = strClusterDir + strFileName[11:]
            rafVSEARCHcluster = sp.Popen(["vsearch", "--cluster_size", strInputCommand,
                                          "--fasta_width", "0", "--id", strIdentity,
                                          "--sizein", "--minseqlength", "1",
                                          "--centroids", strOutputCommand,
                                          "--sizeout"], stdout=sp.PIPE, stderr=sp.PIPE)
            strOut, strError = rafVSEARCHcluster.communicate()
        else:
            pass
                                          
# The getVSEARCHsortBySize function.
# This function controls the VSEARCH sorting. Every fasta file created by
# getVSEARCHderep is sorted based on abundance. Any reads with a abundance lower
# than strMinSizeAbundance will be discarded.
def getVSEARCHsortBySize(flZip, strMinSizeAbundance):
    for strFileName in os.listdir(flZip):
        if strFileName.startswith("derep"):
            strInputCommand = flZip + strFileName
            strOutputCommand = flZip + "sorted" + strFileName
            rafVSEARCHsort = sp.Popen(["vsearch", "--sortbysize", strInputCommand,
                                       "--output", strOutputCommand, "--minseqlength",
                                       "1", "--minsize", strMinSizeAbundance],
                                       stdout=sp.PIPE, stderr=sp.PIPE)
            strOut, strError = rafVSEARCHsort.communicate()
        else:
            pass

# The getVSEARCHderep function.
# This function controls the VSEARCH dereplication. Every fasta file created by
# getFastaFiles is dereplicated. This step is necessary for the sorting step to
# work.
def getVSEARCHderep(flZip):
    for strFileName in os.listdir(flZip):
        if strFileName.endswith(".fasta"):
            strInputCommand = flZip + strFileName
            strOutputCommand = flZip + "derep" + strFileName
            rafVSEARCHderep = sp.Popen(["vsearch", "--derep_fulllength",
                                        strInputCommand, "--output", strOutputCommand,
                                        "--minseqlength", "1", "--sizeout"],
                                        stdout=sp.PIPE, stderr=sp.PIPE)

            strOut, strError = rafVSEARCHderep.communicate()
        else:
            pass

# The getFastaFiles function.
# This function creates separate fasta files for every unique UMI. The function
# creates a unique name for every UMI file and combines that with the desired
# output path. A file is opened or created based on this combination. The
# read header and the read itself are appended to it. 
def getFastaFiles(strHeader, strRead, strCode, dicUniqueUmi, flZip):
    strFileIdentifier = "UMI#" + str(dicUniqueUmi[strCode]) + "_" + strCode + ".fasta"
    strFileName = flZip + strFileIdentifier
    with open(strFileName, "a") as flOutput:
        flOutput.write(strHeader)
        flOutput.write(strRead)

# The getTargetZero function.
# This function checks if both the forward and reverse primer can be found, if
# that succeeds, the forward or reverse (when working with single UMIs) or
# forward and reverse (when working with double UMIs) nucleotides are isolated
# based on the length of the UMI. This isolation is done from the first
# nucleotide at the 5'-end of a read and the last nucleotide at the 3'-end of a
# read.
def getTargetZero(strRead, intUmiLength, strSearch, strForward, strReverse):
    tplCheckForward = re.search(strForward, strRead)
    if tplCheckForward != None:
        tplCheckReverse = re.search(strReverse, strRead)
        if tplCheckReverse != None:
            if strSearch == "umi5":
                return strRead[0:int(intUmiLength)]
            elif strSearch == "umidouble":
                return (strRead[0:int(intUmiLength)],
                        strRead[-int(intUmiLength):])
            elif strSearch == "umi3":
                return strRead[-int(intUmiLength):]
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
def getTargetFront(strRead, intUmiLength, strSearch, strForward, strReverse):
    if strSearch == "umi5" or strSearch == "umidouble":
        intPositionForward = re.search(strForward, strRead).end()
        intPositionUmiForward = intPositionForward + int(intUmiLength)
        strUmiCodeForward = strRead[intPositionForward:intPositionUmiForward]
        if strSearch == "umi5":
            tplCheckReverse = re.search(strReverse, strRead)
            if tplCheckReverse != None:
                return strUmiCodeForward
            else:
                pass
        elif strSearch == "umidouble":
            intPositionReverse = re.search(strReverse, strRead).start()
            intPositionUmiReverse = intPositionReverse - int(intUmiLength)
            strUmiCodeReverse = strRead[intPositionUmiReverse:intPositionReverse]
            return strUmiCodeForward, strUmiCodeReverse
        else:
            pass
    elif strSearch == "umi3":
        tplCheckForward = re.search(strForward, strRead)
        if tplCheckForward != None:
            intPositionReverse = re.search(strReverse, strRead).start()
            intPositionUmiReverse = intPositionReverse - int(intUmiLength)
            strUmiCodeReverse = strRead[intPositionUmiReverse:intPositionReverse]
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
def getTargetBehind(strRead, intUmiLength, strSearch, strForward, strReverse):
    if strSearch == "umi5" or strSearch == "umidouble":
        intPositionForward = re.search(strForward, strRead).start()
        intPositionUmiForward = intPositionForward - int(intUmiLength)
        strUmiCodeForward = strRead[intPositionUmiForward:intPositionForward]
        if strSearch == "umi5":
            tplCheckReverse = re.search(strReverse, strRead)
            if tplCheckReverse != None:
                return strUmiCodeForward
            else:
                pass
        elif strSearch == "umidouble":
            intPositionReverse = re.search(strReverse, strRead).end()
            intPositionUmiReverse = intPositionReverse + int(intUmiLength)
            strUmiCodeReverse = strRead[intPositionReverse:intPositionUmiReverse]
            return strUmiCodeForward, strUmiCodeReverse
        else:
            pass
    elif strSearch == "umi3":
        tplCheckForward = re.search(strForward, strRead)
        if tplCheckForward != None:
            intPositionReverse = re.search(strReverse, strRead).end()
            intPositionUmiReverse = intPositionReverse + int(intUmiLength)
            strUmiCodeReverse = strRead[intPositionReverse:intPositionUmiReverse]
            return strUmiCodeReverse
        else:
            pass
    else:
        pass

# The getReverseComplement function.
# This function creates a complementary string using a sequence as input. The
# function loops through a list version of the sequence and checks and changes
# every character. It then returns a joined string.
def getReverseComplement(strLine):
    dicComplementCodes = {"A": "T", "T": "A", "G": "C", "C": "G", "M": "K",
                          "R": "Y", "W": "W", "S": "S", "Y": "R", "K": "M",
                          "V": "B", "H": "D", "D": "H", "B": "V", "N": "N"}
    lstLine = list(strLine)
    for intPosition in range(len(lstLine)):
        lstLine[intPosition] = dicComplementCodes[lstLine[intPosition]]
    return "".join(lstLine)

# The getRegex function.
# This function creates a regex string using a sequence as input. This regex
# string is based on the IUPAC ambiguity codes. The function loops through
# a list version of the sequence and checks per character if it is a
# ambiguous character, and changes it into regex code if so. It then returns
# a joined string.
def getRegex(strLine):
    dicAmbiguityCodes = {"M": "[AC]", "R": "[AG]", "W": "[AT]", "S": "[CG]",
                         "Y": "[CT]", "K": "[GT]", "V": "[ACG]", "H": "[ACT]",
                         "D": "[AGT]", "B": "[CGT]", "N": "[GATC]"}
    lstLine = list(strLine)
    for intPosition in range(len(lstLine)):
        if lstLine[intPosition] != "A" and lstLine[intPosition] != "T" and\
           lstLine[intPosition] != "G" and lstLine[intPosition] != "C":
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
def getUmiCode(strRead, strProcess, intUmiLength, strSearch, strForward,
               strReverse):
    strRead = strRead.strip("\n")
    strRegexForward = getRegex(strForward)
    strRegexComplementReverse = getRegex(getReverseComplement(\
                                         strReverse[::-1]))
    if strProcess == "primer":
        try:
            return getTargetBehind(strRead, intUmiLength, strSearch,
                                   strRegexForward, strRegexComplementReverse)
        except AttributeError:
            pass
    elif strProcess == "scaffold":
        try:
            return getTargetFront(strRead, intUmiLength, strSearch,
                                  strRegexForward, strRegexComplementReverse)
        except AttributeError:
            pass
    elif strProcess == "zero":
        try:
            return getTargetZero(strRead, intUmiLength, strSearch,
                                 strRegexForward, strRegexComplementReverse)
        except AttributeError:
            pass
    else:
        pass

# The getUmiCollection function.
# This function opens the input file and loops through it. It isolates the
# read headers and reads. For every read the getUmiCode function is called
# which outputs one or two UMIs. In the case of a double UMI search [umidouble]
# the two UMIs are put together. For every read that contains a UMI the
# getFastaFiles function is called. After all reads have been processed, the
# getVSEARCH and setOutputFiles functions are called.
def getUmiCollection(flInput, strClusterDir, flTabular, flZip, flBlast,
                     strProcess, intUmiLength, strSearch, strForward,
                     strReverse, strFormat, strOperand, strIdentity,
                     strMinSizeAbundance):
    dicUniqueUmi = {}
    intNoUmiInReadCount = 0
    intUniqueUmi = 1
    with open(flInput) as oisInput:
        for strLine in oisInput:
            if strLine[0] == strOperand and bool(re.match("[A-Za-z0-9]",
               strLine[1])) == True:
                strHeader = strLine
                strRead = next(oisInput)
                try:
                    strUmiCode = getUmiCode(strRead.upper(), strProcess,
                                            intUmiLength, strSearch,
                                            strForward.upper(),
                                            strReverse.upper())
                except UnboundLocalError:
                    intNoUmiInReadCount += 1
                try:
                    if strUmiCode != None:
                        if strSearch == "umi5" or strSearch == "umi3":
                            strCode = strUmiCode
                        elif strSearch == "umidouble":
                            strCode = strUmiCode[0] + strUmiCode[1]
                        else:
                            pass
                        if strCode not in dicUniqueUmi:
                            dicUniqueUmi[strCode] = intUniqueUmi
                            intUniqueUmi += 1
                        else:
                            pass
                    else:
                        pass
                except UnboundLocalError:
                    pass
                try:
                    if strCode != None:
                        getFastaFiles(strHeader, strRead, strCode, dicUniqueUmi,
                                  flZip)
                    else:
                        pass
                except UnboundLocalError:
                    pass
            strUmiCode = None
            strCode = None
    getVSEARCHderep(flZip)
    getVSEARCHsortBySize(flZip, strMinSizeAbundance)
    getVSEARCHclusterSize(flZip, strClusterDir, strIdentity)
    setOutputFiles(strClusterDir, flBlast, flTabular)

# The setFormat function.
# This function specifies the first character of the read headers based on the
# input file format. It then calls the getUmiCollection function.
def setFormat(flInput, strClusterDir, flTabular, flZip, flBlast, strProcess,
              strFormat, intUmiLength, strSearch, strForward, strReverse,
              strIdentity, strMinSizeAbundance):
    if strFormat == "fasta":
        strOperand = ">"
    elif strFormat == "fastq":
        strOperand = "@"
    else:
        pass
    getUmiCollection(flInput, strClusterDir, flTabular, flZip, flBlast,
                     strProcess, intUmiLength, strSearch, strForward,
                     strReverse, strFormat, strOperand, str(strIdentity),
                     str(strMinSizeAbundance))

# The argvs function.
def parseArgvs():
    parser = argparse.ArgumentParser(description="Use a python script to\
                                                  accumulate all UMIs and\
                                                  output a tabular file, a\
                                                  BLAST file and a zip file.")
    parser.add_argument("-v", action="version", version="%(prog)s [0.1.0]")
    parser.add_argument("-i", action="store", dest="fisInput",
                        help="The location of the input file(s)")
    parser.add_argument("-c", action="store", dest="fosClusterDirectory",
                        help="The location of the clustering output file(s)")
    parser.add_argument("-o", action="store", dest="fosOutput",
                        help="The location of the tabular output file(s)")
    parser.add_argument("-z", action="store", dest="fosOutputZip",
                        help="The location of the pre-vsearch zip output file(s)")
    parser.add_argument("-q", action="store", dest="fosBlastFile",
                        help="The location of the BLAST output file(s)")
    parser.add_argument("-p", action="store", dest="disProcess",
                        help="The UMI search approach [primer/scaffold(adapter)/zero]")
    parser.add_argument("-f", action="store", dest="disFormat",
                        help="The format of the input file(s) [fasta/fastq]")
    parser.add_argument("-l", action="store", dest="disUmiLength",
                        help="The length of the UMI sequences")
    parser.add_argument("-s", action="store", dest="disSearch",
                        help="Search UMIs at 5'-end [umi5], 3'-end [umi3] or at\
                              5'-end and 3'-end [umidouble]")
    parser.add_argument("-a", action="store", dest="disForward",
                        help="The 5'-end search nucleotides")
    parser.add_argument("-b", action="store", dest="disReverse",
                        help="The 3'-end search nucleotides")
    parser.add_argument("-d", action="store", dest="disIdentity",
                        help="The identity percentage with which to perform the\
                              final VSEARCH check")
    parser.add_argument("-u", action="store", dest="disAbundance",
                        help="The minimum abundance a read has to be present in\
                              order to be part of the final VSEARCH check")
    argvs = parser.parse_args()
    return argvs

# The main function.
def main():
    argvs = parseArgvs()
    setFormat(argvs.fisInput, argvs.fosClusterDirectory, argvs.fosOutput,
              argvs.fosOutputZip, argvs.fosBlastFile, argvs.disProcess,
              argvs.disFormat, argvs.disUmiLength, argvs.disSearch,
              argvs.disForward, argvs.disReverse, argvs.disIdentity,
              argvs.disAbundance)

if __name__ == "__main__":
    main()
