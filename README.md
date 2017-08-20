# FSFix
FSFix uses Illumina short reads to fix frameshift (indel) errors in a Pacbio assembly. It fixes many of the remaining insertions or deletions (indels) in homopolymer regions, and indels caused by reads that were derived from heterozygous haplotigs, even after assembly polishing by Arrow, Quiver, and Pilon. This program was used in a xxx genome assembly project (link to publication xxx), where it fixed 1255 indels in a 720Mb assembly.
This program was written by Chern Han Yong. Please contact me at chernycherny@hotmail.com if you have any questions.


CITATION
If you use this program, please cite: xxx. 


REQUIREMENTS
1. Samtools
2. Perl interpreter


INPUT
1. unmasked.bed - bed file of unmasked region (ie. non-repeat region) of the genome to fix. Fixing only the unmasked region will require much less time and space. Generate this file eg. from RepeatMasker's output. If you want to fix the entire genome, then create a bed file covering the entire genome.
2. alignment.bam - alignment of Illumina short reads onto Pacbio assembly.
3. genome.fasta - Pacbio assembly


USAGE
1. Extract alignments of short reads on the unmasked region, with map quality > 20
samtools view alignment.bam -f 0x2 -q 20 -L unmasked.bed > dna_reads_unmasked_tmp.sam

2. Remove short read alignments that are not uniquely mapped
grep -v "XA:Z" dna_reads_unmasked_tmp.sam > dna_reads_unmasked.sam

3. Analyze DNA read alignments for frameshifts
perl find_frameshifts.pl -s dna_reads_unmasked.sam > out_find_frameshifts

4. Fix frameshifts
perl fix_frameshifts.pl -i out_find_frameshifts -s genome.fasta -o genome_FSfix.fasta > out_fix_frameshifts


DETAILS
FSFix corrects many of the remaining insertions or deletions (indels) in homopolymer regions, and indels caused by reads that were derived from heterozygous haplotigs, from a Pacbio assembly, even after assembly polishing by Arrow, Quiver, and Pilon. First, FSFix aligns the Illumina short reads onto the assembly, keeping only unique alignments with mapping quality >20. Next, FSFix looks indels in the alignment supported by at least 90% of the aligned reads, and with at least 5 supporting reads. Reads for which the indel under consideration falls within the first 5 or last 5 bases are ignored in calculating this support. Finally, FSFix corrects the assembly using these supported indels: supported deletions in the aligned reads are used to delete base(s) from the assembly, while supported insertions in the aligned reads are used to insert base(s) into the assembly. FSFix also looks for a special case of indels likely derived from reads assembled from heterozygous haplotigs, where the assembly contains both alleles from the heterozygous haplotigs. This is manifested as two nearby deletions in the aligned reads. This is not detected by FSFix in the general case, because each deletion would only be supported by about 50% of reads. Thus, if FSFix finds two deletions, X1 and X2, that are are within 3bp of each other, and X1’s deletion-supporting reads are distinct from X2’s deletion-supporting reads, FSFix considers X1 and X2 as the same deletion event, taking its supporting reads as the union of X1 and X2’s supporting reads, and calculates if the deletion is supported as before. If so, it deletes only one of either X1 or X2. Please see (link to publication xxx) for full details.
