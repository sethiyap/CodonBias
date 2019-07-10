###############
## Author : Pooja Sethiya
## Institute : Chris Lab Faculty of Health Science / University of Macau. 
## Email : yb57662@umac.mo
## Date : 13/06/2018
###############
## Get Codon Usage Index using UCO for multi-fasta file
## By Default Median is calculated of the Relative Synonymous Codon Usage (RSCU)
## You can specify if you want 1. eff: codon counts 2. freq: codon relative frequencies 
## Provide file path of the multifasta file
## codon_usage_by: "rscu", "eff", "freq"
## summariseBy:1=mean, 2:=median
## RunScript as : get_codon_usage_index(multifasta_file, condon_usage_by="rscu", summariseBy=1)
##
mylist=clipr::read_clip()

#' get_codon_usage_index
#'
#' @param multifasta_file 
#' @param condon_usage_by 
#' @param summariseBy 
#'
#' @return 
#' @export
#'
#' @examples
get_codon_usage_index = function(multifasta_file, condon_usage_by="rscu", summariseBy, mylist=NULL, outfile){
#-- Prerequisite for summary          
          summariseByValue=c("mean","median")
#--- Load Packages          
          require(tidyverse)
          require(dplyr)
          require(seqinr)
          require(data.table)

#--- Load Sequence file
          Codon_Usage=data_frame()  
          seq_check = read.fasta(multifasta_file,as.string = TRUE,forceDNAtolower = FALSE)
          
          seq_characters <- getSequence(seq_check) %>% tibble(data=.)  %>% unnest() %>% unique() %>% as.matrix()
          match_dna <- c("A|C|G|T")
          check <- seq_characters %like% match_dna
          
          if(FALSE %in% check){
                    stop(" Not a standard DNA input ")
          }
          
          else{
                    message("DNA input file uploaded")
                    
                    if(is.null(mylist) == FALSE){
                              names(seq_check ) <- gsub(' .*', '', names(seq_check ))
                              ID = as.matrix(mylist)
                              subset_fasta=seq_check[ID[,1]]
                              
                              message(paste(outfile,length(names(subset_fasta)),"_nucl.fasta", sep=""), " contains ", length(names(subset_fasta)), " sequences")
                              
                              write.fasta(subset_fasta, names = names(subset_fasta),as.string = TRUE,
                                          file.out = paste(outfile,length(names(subset_fasta)),"_nucl.fasta", sep=""),nbchar = nchar(subset_fasta[1]))
                              
                              seq <- subset_fasta
                              
                    }
                    else{
                              seq <- seq_check
                              message(basename(multifasta_file), " contains ", length(names(seq)), " sequences")
                              
                    }
                    
          } 
          
#--- for rscu         
          if(condon_usage_by=="rscu"){
                    
                    Codon_Usage =  lapply(seq_along(seq),function(i){
                              
                              name = getName(seq[[i]])
                              CBI = uco(getSequence(seq[[i]], as.string = FALSE),frame = 0,index="rscu",NA.rscu = 0)
                              
                              return(data.frame(names(CBI),name=name,CBI)) 
                              
                    }) 
                    
#-- Combine Data                     
                    Codon_Usage_df <- do.call("rbind",Codon_Usage)

                    rownames(Codon_Usage_df) <- NULL
                    
                    df_CBI = spread(Codon_Usage_df,key = names.CBI.,value = CBI) 
                    print(df_CBI)
#-- Plot                     
                    Codon_Usage_data = df_CBI %>% dplyr::select(matches("^(a|t|g|c)")) %>% 
                              summarise_all(funs(!!(summariseByValue[summariseBy]))) %>% 
                              gather(key=names) 
                    Codon_Usage_Plot <- Codon_Usage_data %>% 
                                        ggplot(aes(names, value,fill=names))+geom_col()+xlab("")+ylab("Codon Usage Index")+theme_bw()+
                                        theme(axis.text.x= element_text(face="bold", colour="black", size=10,angle=45,vjust=0.8),
                                              axis.text.y = element_text(face="bold", color="black",size=12),
                                              axis.title.y=element_text(face="bold", color="black",size=14),
                                              legend.position="none")+coord_polar()
                    
                    print(Codon_Usage_Plot)
                    
                    
                    
          }
          
#--- for eff and freq
          
          if(condon_usage_by=="eff"|condon_usage_by=="freq"){
                    
                    Codon_Usage =  lapply(seq_along(seq),function(i){
                              
                              name = getName(seq[[i]])
                              CBI = uco(getSequence(seq[[i]], as.string = FALSE),frame = 0,index=condon_usage_by,NA.rscu = 0)
                              
                              return(data.frame(name=name,CBI)) 
                              
                    }) 
                    
                    Codon_Usage_df <- do.call("rbind",Codon_Usage)
                    
    
                    df_CBI = spread(Codon_Usage_df,key = Var1,value = Freq) 
                    print(df_CBI)
                    
                    Codon_Usage_data = df_CBI %>% dplyr::select(matches("^(a|t|g|c)")) %>% 
                              summarise_all(funs(!!(summariseByValue[summariseBy]))) %>% 
                              gather(key=names) 
                    
                    Codon_Usage_Plot <-  Codon_Usage_data %>% 
                                                  ggplot(aes(names, value,fill=names))+geom_col()+xlab("")+ylab("Codon Counts Frequency")+theme_bw()+
                                                  theme(axis.text.x= element_text(face="bold", colour="black", size=10,angle=90,vjust=0.8),
                                                        axis.text.y = element_text(face="bold", color="black",size=12),
                                                        axis.title.y=element_text(face="bold", color="black",size=14),
                                                        legend.position="none")+coord_polar()
                     
                    print(Codon_Usage_Plot)
                    
                    
                    
          }
          
          write_delim(Codon_Usage_data, path=paste(basename(outfile),condon_usage_by,"codondistri.tab", sep="_"), delim="\t")
          ggsave(paste(basename(outfile),condon_usage_by,"codondistri.pdf", sep="_"),device = "pdf", width=8, height = 7)
}





