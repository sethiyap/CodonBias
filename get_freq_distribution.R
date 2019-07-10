multifasta_file = "/Users/Pooja/Documents/Data-Analysis/Others/Reference_Annotation/An/A_nidulans_FGSC_A4_version_s10-m04-r12_orf_trans_all.fasta"
mylist=clipr::read_clip()

#' get_freq_distribution
#' @author Pooja Sethiya
#' @date 09/07/2019
#'
#' @param multifasta_file 
#' @param mylist 
#' @param outfile 
#'
#' @return
#' @export
#'
#' @examples

#run as: get_freq_distribution(multifasta_file, mylist = mylist)
# --- function
get_freq_distribution <- function(multifasta_file, mylist=NULL, outfile="sample"){
          
          # ---- Load packages
          library(Biostrings)
          library(seqinr)
          library(tidyverse)
          library(dplyr)
          library(data.table)
          
          # ---- get subset of the sequences if mylist is uploaded ----    
          check_seq <-  read.fasta(multifasta_file,as.string = TRUE,forceDNAtolower = FALSE)
          
          if(is.null(mylist) == FALSE){
                    names(check_seq) <- gsub(' .*', '', names(check_seq))
                    ID = as.matrix(mylist)
                    subset_fasta=check_seq[ID[,1]]
                    
                    message(paste(outfile,length(names(subset_fasta)),".fasta", sep=""), " contains ", length(names(subset_fasta)), " sequences")
                    
                    write.fasta(subset_fasta, names = names(subset_fasta),as.string = TRUE,
                                file.out = paste(outfile,length(names(subset_fasta)),".fasta", sep=""),nbchar = nchar(subset_fasta[1]))
                    
                    subset_seq <- paste(outfile,length(names(subset_fasta)),".fasta", sep="")
                    
          }
          else{
                    subset_seq <- multifasta_file
                    message(basename(multifasta_file), " contains ", length(names(check_seq)), " sequences")
                    
          }
          
          # ---- Check if the given file is nucleotide or protein fasta ----
          seq_characters <- getSequence(check_seq) %>% tibble(data=.)  %>% unnest() %>% unique() %>% as.matrix()
          match_dna <- c("A|C|G|T")
          check <- seq_characters %like% match_dna
          
          if(FALSE %in% check){
                    message("Input multi-fasta is protein file...")
                    seq = Biostrings::readAAStringSet(subset_seq)
                    x1="Amino acids"
                    y1="Amino acid frequency"
                    y2="average (amino acid frequency)"
                    select_column = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
          }
          
          else{
                    message("Input multi-fasta is DNA file...")
                    seq = Biostrings::readDNAStringSet(subset_seq)
                    x1="Nucleotides"
                    y1="Nucleotide frequency"
                    y2="average (nucleotide frequency)"
                    select_column = c("A", "T", "G", "C")
          }
          
          # ---- calculate the frequency ----
          message(" Calculating percent frequency ...")
          aa_count <- alphabetFrequency(seq) %>% as_tibble() %>%
                    dplyr::select(select_column)  %>%
                    mutate(names=names(seq))
          
          
          aa_freq <- aa_count %>% gather(AA, count, -names) %>% 
                    group_by(names) %>%
                    mutate(total_count=sum(count), freq_count=100*(count/total_count)) %>% 
                    dplyr::select(-c("total_count", "count")) 
          
          # --- plot heatmap of the frequency gene wise ----
          message(" Plotting raw count...")
          gg_1 <-  aa_freq %>% spread(AA, freq_count) %>% ungroup() %>% 
                    #dplyr::top_n(20) %>%
                    gather(AA, freq_count,-names) %>% 
                    ggplot(aes(AA, names, fill=freq_count))+
                    geom_tile()+theme_bw()+scale_fill_gradient(high = "#7a0177",low="#feebe2")+
                    xlab(x1)+ylab("")+theme_bw()+
                    theme(axis.text.x= element_text(face="bold", colour="black", size=10,angle=90,vjust=0.8),
                          axis.text.y = element_text(color="black",size=4.5),
                          axis.title.y=element_text(face="bold", color="black",size=14))
          
          print(gg_1)
          ggsave(gg_1, filename = paste(outfile,"_genewise_freq.pdf", sep=""))
          
          # --- plot boxplot of the frequency gene wise ----
          gg_2 <-  aa_freq %>% 
                    ggplot(aes(AA, freq_count, fill=AA))+
                    geom_violin()+theme_bw()+
                    xlab("")+ylab(y1)+theme_bw()+
                    theme(axis.text.x= element_text(face="bold", colour="black", size=10,angle=90,vjust=0.8),
                          axis.text.y = element_text(color="black",size=10),
                          axis.title.y=element_text(face="bold", color="black",size=14),
                          legend.position="none")
          
          print(gg_2)
          ggsave(gg_2, filename = paste(outfile,"_boxplot.pdf", sep=""))
          
          
          # --- plot bar chart of the average frequency ----
          
          message(" Computing and writing average ...")
          avg_dat <- aa_freq %>% spread(AA, freq_count) %>% ungroup() %>% 
                    dplyr::select(-c(names)) %>%
                    colMeans(.) %>% data.frame(mean=.) %>% rownames_to_column("AA") 
          
          write_delim(avg_dat, path =paste(outfile,"_average_freq.tab", sep=""), delim="\t" )
              
          message(" plotting average count ...")
          gg_3<-avg_dat  %>%
                              ggplot(aes(AA, mean, fill=AA))+geom_col()+xlab("")+ylab(y2)+theme_bw()+
                              #coord_polar()+
                              theme(axis.text.x= element_text(face="bold", colour="black", size=10,angle=90,vjust=0.8),
                                    axis.text.y = element_text(face="bold", color="black",size=12),
                                    axis.title.y=element_text(face="bold", color="black",size=14),
                                    legend.position="none")
          
          print(gg_3)
          ggsave(gg_3, filename = paste(outfile,"_average_freq.pdf", sep=""))
}
