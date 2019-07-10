
   
dir_fasta_files <- (".")

optimal_codon_distribution <- function(dir_fasta_files,output_name ){
          library(tidyverse)
          library(rlist)
          library(Biostrings)
          library(tidyr)
          library(purrr) 
          library(seqinr)
          
          
          cg_optimal_codon_table <- read_delim("/Users/Pooja/Documents/Data-Analysis/PhD/Pol2-RNASeq/CG-RNASeq-PolI/Cg-Pol2-All-Time-RNASeq2hr-rawdata/Genome-Reference-Files/cg_optimal_codons.txt", 
                                               col_names = TRUE, delim="\t")
          
          head(cg_optimal_codon_table) # O= Optimal , N=Non-optimal
          # codon Cgla 
          # 1 TTT   N    
          # 2 TTC   O    
          # 3 TTA   N 
 
          fa_files <- list.files(dir_fasta_files, pattern = "*.fa$", recursive = T, full.names = T)
          names(fa_files) <- gsub(pattern = "*.fa$",replacement = "", fa_files)
          names(fa_files) <- gsub(pattern = "./",replacement = "", names(fa_files))
          
          
          fa_files <- broom::tidy(fa_files)
          print(fa_files)
          
          xx= fa_files %>% mutate(cond_freq=map(x, function(ii){
                    read.fasta(ii,as.string = TRUE,forceDNAtolower = FALSE) %>% 
                              map(~uco(getSequence(.x, as.string = FALSE),frame = 0,index="eff"))  %>%
                              do.call("cbind",.) %>%  
                              as_tibble(rownames = "codon") %>%
                              mutate(codon=toupper(codon)) %>%
                              left_join(cg_optimal_codon_table, by = c("codon")) %>%
                              gather("gene","freq", -codon, -Cgla) %>%
                              dplyr::group_by(gene) %>%  
                              dplyr::summarize(Total = sum(freq),OnlyOptimal = sum(freq[Cgla=="O"])) %>% 
                              mutate(percent=100*(OnlyOptimal/Total))
          })) %>% dplyr::select(-c(x)) %>% tidyr::unnest()
          
          
          gg <- ggplot(xx,aes(names,percent,fill=names))+
                    geom_boxplot(notch = TRUE, notchwidth = 0.4, alpha=0.9) +
                    scale_fill_manual(values=c("#FC4E07", "#00AFBB" ))+
                    scale_y_log10()+
                    #geom_text(data = median, aes(label = med, y = med + 170))+
                    theme_bw()+
                    labs(x="", y="optimal codons (in %)")+
                    scale_x_discrete(position="top")+
                    theme(axis.text.x = element_text(color="black",size=12),
                          axis.text.y = element_text(color="black",size=12),
                          axis.title.y = element_text(color="black",size=12, face="bold"),
                          strip.text = element_text(color="black",size=12,face="bold"),
                          strip.background  = element_rect(fill = "white"),
                          legend.position  = "none"
                    )
          
          print(gg)
          ggsave(paste(output_name,"_optimalCodon_percentage.pdf",sep=""), path="./" ,dpi=300, device = "pdf")
          write_delim(xx,paste(output_name,"_optimalCodon_percentage.tab",sep=""),delim="\t",col_names = TRUE)
          
          
}

#---- rndom genes
library(rlist)
codon_table <- read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)

all_orf <- read.fasta("../../../Genome-Reference-Files/C_glabrata_CBS138_version_s02-m07-r27_orf_genomic.fasta", forceDNAtolower = FALSE,as.string = TRUE)
head(all_orf)



rand_a <-  list.sample(all_orf, size = 186,replace = FALSE,weight = 1, prob = NULL) %>% 
          map(~uco(getSequence(.x, as.string = FALSE),frame = 0,index="eff"))  %>%
          do.call("cbind",.) %>% data.frame() %>%
          mutate(optimal = codon_table$Cgla[match(toupper(row.names(.)), codon_table$codon)]) %>%
          melt(., id.vars = c("optimal")) %>% group_by(variable) %>%  
          summarize(Total = sum(value),OnlyOptimal = sum(value[optimal=="O"])) %>% 
          mutate(percent=100*(OnlyOptimal/Total), label=rep( "random_1", length(Total)))

rand_b <- list.sample(all_orf, size = 86,replace = FALSE,weight = 1, prob = NULL) %>% 
          map(~uco(getSequence(.x, as.string = FALSE),frame = 0,index="eff"))  %>%
          do.call("cbind",.) %>% data.frame() %>%
          mutate(optimal = codon_table$Cgla[match(toupper(row.names(.)), codon_table$codon)]) %>%
          melt(., id.vars = c("optimal")) %>% group_by(variable) %>%  
          summarize(Total = sum(value),OnlyOptimal = sum(value[optimal=="O"])) %>% 
          mutate(percent=100*(OnlyOptimal/Total), label=rep( "random_2", length(Total)))

master_rand <- rbind(rand_a,rand_b)    
ggplot(master_rand,aes(label,percent,fill=label))+geom_boxplot(alpha=0.6)+theme_bw()+
          theme(axis.text.x= element_text(face="bold", colour="black", size=12,angle=0,vjust=0.8),
                axis.text.y = element_text(face="bold", color="black",size=12),
                axis.title.y=element_text(face="bold", color="black",size=14),
                legend.position="none")+ylab("percent optimal codons")+xlab("")
t.test(rand_a$percent, rand_b$percent)


