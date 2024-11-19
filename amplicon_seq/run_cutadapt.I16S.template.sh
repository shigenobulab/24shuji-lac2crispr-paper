#!/bin/sh
#$ -S /bin/bash
#$ -cwd

#===
# run_cutadapt.PE.template.sh
# 
# A template script for pre-processing Illumina TruSeq paried-end reads, 
# using cutadapt software.
#
# Validated cutadapt versions: 1.9.1
#
# author: Shuji Shigenobu <shige@nibb.ac.jp>
# date:   Jan 22, 2016
#

OVERHANG_COMMON_REV=CTGTCTCTTATACACATCT


QV=20
O=8
MINCUT=25
INPUT1=%READ1%
INPUT2=%READ2%

OUTPUT1=`basename $INPUT1 .fastq.gz`.cln.Q${QV}L${MINCUT}.fastq.gz
OUTPUT2=`basename $INPUT2 .fastq.gz`.cln.Q${QV}L${MINCUT}.fastq.gz

CUTADAPT=cutadapt

$CUTADAPT \
  -q $QV \
  -O $O \
  -a $OVERHANG_COMMON_REV \
  -A $OVERHANG_COMMON_REV \
  --minimum-length $MINCUT \
  -o $OUTPUT1 \
  -p $OUTPUT2 \
  $INPUT1 \
  $INPUT2

# Other popular options:
# --quality-base=33 or 64
#   usually 33, but 64 for old Illumina

# Reference
# - cutadapt manuals / section on illumina-truseq: 
#    http://cutadapt.readthedocs.org/en/stable/guide.html#illumina-truseq
#
# Acknowledgements
# - K Yamaguchi.
