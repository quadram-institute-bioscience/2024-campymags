
# Large files 
Notice that in the code snippets below, the file names may be different or directories may be missing.  

[core_gene_alignment.aln.xz](core_gene_alignment.aln.xz) was generated with 
```
nohup panaroo -i 00.60_genomes_bakta/PID_* ../12.all_19_metagenomes/*gff3 -o panaroo_MAGS_and_60genomes-core --remove-invalid-genes -t 38 --core_threshold 0.2 -a core  --clean-mode strict &
```
IQTREE file [trimal99.aln.treefile](trimal99.aln.treefile) and alignment [trimal99.aln.xz](trimal99.aln.xz) were generated with 

```
trimal -in core_gene_alignment.aln -out trimal_alignment.aln  -gt 0.9
nohup sh -c 'iqtree -s triaml_alignment.aln .aln  -m HKY+G -nt 12' > nohup.out  &
```
