SEQ1=%READ1%
SEQ2=%READ2%
NAME=%NAME%

REF=../../../Data/Acyrthosiphon_pisum/NIBB_ApisBuc1/bowtie2/ApisBuc1.genome
OUTBAM=${NAME}.on.`basename $REF`.bowtie2.bam
NCPU=12

bowtie2 -p $NCPU \
  -x $REF \
  -1 $SEQ1 \
  -2 $SEQ2 \
  --end-to-end \
  | samtools view -bS - \
  > $OUTBAM

## bam sort

#=== config
MAX_MEMORY=8G
NCPU=8
#===
OUTBAM_SORTED=`basename $OUTBAM .bam`.sorted.bam

samtools  sort -m $MAX_MEMORY -o $OUTBAM_SORTED -@ $NCPU $OUTBAM
samtools  index $OUTBAM_SORTED
