# DMDuper
Workflow to detect duplication in a gene and assign tandem/non-tandem status

This is a beta version of the dmd duplication caller. It checks for duplication signatures based on softclip reads, performs haplotype-aware assembly of the soft-clipped reads and then maps it to the genome to identify the supplementary alignments to call the tandem/non-tandem status. 

Dependencies:
* Minimap2
* Seqkit
* Flye
* Output of sniffles2 vcf

![alt text](https://github.com/jagadhesh89/dmduper/blob/main/DMD.png)


Usage:
```
python dmduper.py -b sample_hg38.bam -r hg38.fasta -sv sniffles.vcf.gz -o <output_folder>
```
