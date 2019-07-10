---
title: "Sequence based analysis "
output: 
  html_document:
    keep_md: true
---



## Codon bias

This function allows you to map the codon bias in given multi-fasta file.
The codon bias calculation can be done on:
1.  RSCU: Relative Synonymous Codon Usage, default
1.   eff: codon counts
1.  freq: codon relative frequencies

You can summarise the occurrence of codon based on mean or median values.

### Run this function as:


```r
source("get_codon_usage_index.R")
multifasta_file  <-  "SequenceSample.fa"
get_codon_usage_index(multifasta_file, summariseBy = 1,condon_usage_by = "freq", mylist=NULL,outfile = "sample")
```

### Output of the function is pie-chart displaying distribution of codons in given multifasta file.

![](CodonBiasPlot.png)

## Amino acid or nucleotide frequency distirbution
Frequency of nucleotides or amino acids can be determined from multi-fasta DNA or protein file.
Frequency can also be calculated for specific genelist by specifying mylist or the enitre multi-fasta file. Subset of the fasta file will be generated for given list from reference multi-fasta file. Further the function checks first and alerts if the uploaded file is protein or dna. 
### Run this function as:


```r
source("get_freq_distribution.R")
multifasta_file  <-  "SequenceSample.fa"
get_freq_distribution(multifasta_file = multifasta_file, mylist = NULL, outfile = "SequenceSample")
```

### 

The function gives three out put

1. subset of the multi-fasta file based on the uploaded list
1. Average frequency file
1. Plot of gene wise frequency
![Genewise frequency plot](SequenceSample_genewise_freq.png)
1. Plot of average frequency
![Average frequency plot](SequenceSample_average_freq.png)







