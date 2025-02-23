
# Large files 
Notice that in the code snippets below, the parameters, file names may be different or directories may be missing.  
This is (currently) just a data dump (with original file names as in local directories)...
Whenever not mentioned, files are smaller than 10MB . 


File [06.campy_bins_metawrap.txz](06.campy_bins_metawrap.txz) contains the bins generated by metawrap for the Campylobacter genomes.
File [01.gff3_from_19_bakta.txz](01.gff3_from_19_bakta.txz) contains the gff3 files generated from the bakta directories for the MAG bins generated by metawrap.

Files
[SPAdes_on_collection_786__contigs__fasta.txz](SPAdes_on_collection_786__contigs__fasta.txz)
[SPAdes_on_collection_786__scaffolds__fasta.txz](SPAdes_on_collection_786__scaffolds__fasta.txz)
(55MB each) contain the SPADES results from the (non-MAG) isolates. All of them. 

File [bakta_gff3_isolates.txz](bakta_gff3_isolates.txz) (78MB) contains the gff3 files generated from the bakta directories for the (non-MAG) isolates. All of them.

The files below are for selected isolates (based on phylonium distances), not all of them.


[core_gene_alignment.aln.xz](core_gene_alignment.aln.xz) was generated with 
```
nohup panaroo -i 00.60_genomes_bakta/PID_* ../12.all_19_metagenomes/*gff3 -o panaroo_MAGS_and_60genomes-core --remove-invalid-genes -t 38 --core_threshold 0.2 -a core  --clean-mode strict &
```
IQTREE file [trimal99.aln.treefile](trimal99.aln.treefile) and alignment [trimal99.aln.xz](trimal99.aln.xz) were generated with 

```
trimal -in core_gene_alignment.aln -out trimal_alignment.aln  -gt 0.99
nohup sh -c 'iqtree -s triaml_alignment.aln .aln  -m HKY+G -nt 12' > nohup.out  &
```
