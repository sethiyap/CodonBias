#' remove_first_n_nucleotides
#'
#' @param multifasta_file 
#' @param remove_first_n 
#' @param outfile 
#'
#' @return
#' @export
#'
#' @examples
#' run as: remove_first_n_nucleotides(multifasta_file,remove_first_n = 1,outfile = "an_prot_")


remove_first_n_nucleotides <- function(multifasta_file, remove_first_n, outfile){
          library(Biostrings)
          library(GenomicRanges)
          library(tidyverse)
          library(data.table)
          
          
          check_seq = seqinr::read.fasta(multifasta_file,as.string = TRUE,forceDNAtolower = FALSE)
          
          seq_characters <- seqinr::getSequence(check_seq) %>% tibble(data=.)  %>% unnest() %>% unique() %>% as.matrix()
          match_dna <- c("A|C|G|T")
          check <- seq_characters %like% match_dna
          
          if(FALSE %in% check){
                    message("Input multi-fasta is protein file...")
                    seq = Biostrings::readAAStringSet(multifasta_file)
                   
          }
          
          else{
                    message("Input multi-fasta is DNA file...")
                    seq = Biostrings::readDNAStringSet(multifasta_file)
                    
          }

          print(seq)
          
          #---- subseting the input file
          
          message("subseting input file")
          get_fasta <- Biostrings::subseq(seq, start=(remove_first_n+1))
          print(get_fasta)
          
          
          message("writing multi-fasta file")
          
          # --- write the sequences
          Biostrings::writeXStringSet(get_fasta, paste(outfile,"_remove_first", remove_first_n,".fa", sep=""))
          
          
}
