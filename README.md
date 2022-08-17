# vulgar2gff3

vulgar2gff converts exonerate vulgar output to hierarchical gff3 format that is compatible with many R packages.

### Usage
Get some data:
```
curl 'http://rest.ensembl.org/sequence/region/mouse/5:150445000-150495000' -H 'Content-type:text/x-fasta' > mouse_brca2_gene.nuc
curl 'http://rest.ensembl.org/sequence/id/ENSP00000369497' -H 'Content-type:text/x-fasta' > human_brca2_protein.pep
```
Run exonerate:
```
exonerate -m p2g -t mouse_brca2_gene.nuc -q human_brca2_protein.pep --showalignment F --showtargetgff F --useaatla F > exonerate_brca2.vulgar
```

Convert full output (six alignments) to gff3:
```
perl vulgar2gff3.pl exonerate_brca2.vulgar > mouse_brca2.gff3
```

Manually select an alignment before converting to gff3:
```
head -3 exonerate_brca2.vulgar > exonerate_brca2_first.vulgar
perl vulgar2gff3.pl exonerate_brca2_first.vulgar > mouse_brca2_first.gff3
```
