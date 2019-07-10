#---- date: 18.04.2019
#---- define word size, 1: AA, 3: nucl codons
#---- number_of_codons, how many first codons you need
#---- run as
get_first_n_codons(fasta_file = "SequenceSample.fa",word_size = 1,number_of_codons = 5,output_name = "Sample")



#--- function
get_first_n_codons <- function(fasta_file, word_size, number_of_codons, output_name){
          library(tidyverse)
          library(seqinr)
          library(purrr)
          
          
          #--- read fasta file
          fa_seq <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE)
          
          #--- make tidy form
          fa_tidy <- fa_seq %>% 
                    as_tibble() %>% 
                    gather(key="gene", value="sequence")
          
          #---- get first n codons
          #---- define word size, 1: AA, 3: nucl codons
          
          fa_codon <- fa_tidy %>% mutate(fa_split= map(sequence, function(i){
                    seqinr::splitseq(s2c(i), word=word_size)
          })) %>%
                    mutate(selected_codons=map(fa_split, function(ii){
                              ii[1:number_of_codons]
                    })) %>%
                    dplyr::select(-c(fa_split, sequence)) %>% 
                    unnest(selected_codons) %>%
                    group_by(gene) %>%
                    mutate(col=seq_along(gene)) %>%
                    spread(col, selected_codons)
          
          #--- write output file
          write_delim(fa_codon, paste(output_name, number_of_codons,"tabulated.txt", sep="_"), delim="\t",col_names = TRUE)
          
          
}


