#!/bin/bash

# Invoked by:
# ORDER=Mammalia GENOME=Pan.troglodytes NUM=500 sbatch randomSampling.sh

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH --time=1-00:00 
#SBATCH --mem=32GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@adelaide.edu.au      

# Load the necessary modules
module load BEDTools/2.25.0-foss-2015b 
module load SAMtools/1.2-foss-2015b
module load Biopython/1.66-foss-2016uofa-Python-2.7.11

# make a new dir for each genome
mkdir -p "$GENOME"
cd $GENOME/

# Copy all L1 regions to this directory
cp /data/rc003/atma/LASTZ_extraction/testResults/$ORDER/$GENOME/VERIFIED/"$GENOME"_L1_verified.fasta .

# Extract all L1s between 4.5kb to 5.5kb in length
# (Check that this includes the bat-chimp putative transfer)
usearch -sortbylength "$GENOME"_L1_verified.fasta -minseqlength 4500 -maxseqlength 5500 -fastaout "$GENOME"_L1_4.5to5.5kb.fasta
#6380 L1s for human

# Make a 3-column bed file from the FASTA headers of verified L1 sequences
cat "$GENOME"_L1_verified.fasta \
| grep '>' \
| sed 's/^>//g' \
| awk -F ":" '{print $1 " " $2}' \
| sed 's/(-)/.rev/g' \
| sed 's/(+)/.fwd/g' \
| sed 's/-/ /g' \
| sed 's/.fwd/ +/g' \
| sed 's/.rev/ -/g' \
| awk '{print $1 "\t" $2 "\t" $3 }' \
> "$GENOME"_L1intervals.bed

# Index the genome
samtools faidx /data/rc003/atma/LASTZ_extraction/Genomes/$ORDER/$GENOME/*.fa 

# Only need the first two columns of the faidx file (chr name and length)
cut -f 1,2 /data/rc003/atma/LASTZ_extraction/Genomes/$ORDER/$GENOME/*.fa.fai > "$GENOME"_chrom.sizes

# Now, randomly grab 100 L1s from the 4.5-5.5kb fasta file
cat "$GENOME"_L1_4.5to5.5kb.fasta \
| awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' \
| shuf \
| head -n "$NUM" \
| awk '{printf("%s\n%s\n",$1,$2)}' \
> "$GENOME"_"$NUM"randomL1s.fasta

# make a 6-column bed file (strand included) of these 100 random L1s
cat "$GENOME"_"$NUM"randomL1s.fasta | \
grep '>' |\
sed 's/^>//g' |\
awk -F ":" '{print $1 " " $2}' |\
sed 's/(-)/.rev/g' |\
sed 's/(+)/.fwd/g' |\
sed 's/-/ /g' |\
sed 's/.fwd/ +/g' |\
sed 's/.rev/ -/g' |\
awk '{print $1 "\t" $2 "\t" $3 "\t" "randomL1" "\t" "1" "\t" $4}' \
> "$GENOME"_"$NUM"randomL1s.bed

# sort all bed files
bedtools sort -i "$GENOME"_L1intervals.bed > "$GENOME"_L1intervals.bed.sorted
bedtools sort -i "$GENOME"_"$NUM"randomL1s.bed > "$GENOME"_"$NUM"randomL1s.bed.sorted
sort "$GENOME"_chrom.sizes > "$GENOME"_chrom.sizes.sorted

# shuffle the 100 sample L1s to non-L1 locations in the genome
bedtools shuffle -i "$GENOME"_"$NUM"randomL1s.bed.sorted -excl "$GENOME"_L1intervals.bed.sorted -g "$GENOME"_chrom.sizes.sorted > "$GENOME"_"$NUM"randomNonL1s.bed

# sort the nonL1s
bedtools sort -i "$GENOME"_"$NUM"randomNonL1s.bed > "$GENOME"_"$NUM"randomNonL1s.bed.sorted

# extract FASTA sequences of 100 randoms nonL1s 
bedtools getfasta -s -fi /data/rc003/atma/LASTZ_extraction/Genomes/$ORDER/$GENOME/*.fa -bed "$GENOME"_"$NUM"randomNonL1s.bed.sorted -fo "$GENOME"_"$NUM"randomNonL1s.fasta

# convert both nonL1s and L1s into uppercase
cat "$GENOME"_"$NUM"randomNonL1s.fasta \
| awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' \
> "$GENOME"_"$NUM"randomNonL1s_uppercase.fasta

cat "$GENOME"_"$NUM"randomL1s.fasta \
| awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' \
> "$GENOME"_"$NUM"randomL1s_uppercase.fasta

# remove sequences that have Ns
removeNfromfas.py "$GENOME"_"$NUM"randomNonL1s_uppercase.fasta #output is N_removed.fasta
mv N_removed.fasta "$GENOME"_"$NUM"randomNonL1s_noNs.fasta

removeNfromfas.py "$GENOME"_"$NUM"randomL1s_uppercase.fasta #output is N_removed.fasta
mv N_removed.fasta "$GENOME"_"$NUM"randomL1s_noNs.fasta

# resample a subset of 100 from the seqs without Ns
cat "$GENOME"_"$NUM"randomNonL1s_noNs.fasta \
| awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' \
| shuf \
| head -n 100 \
| awk '{printf("%s\n%s\n",$1,$2)}' \
> "$GENOME"_100_nonL1s.fasta

cat "$GENOME"_"$NUM"randomL1s_noNs.fasta \
| awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' \
| shuf \
| head -n 100 \
| awk '{printf("%s\n%s\n",$1,$2)}' \
> "$GENOME"_100_L1s.fasta
