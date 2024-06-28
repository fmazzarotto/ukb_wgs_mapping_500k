#this function takes a single input, slightly different from the input used to create the map for the previous 200k WGS release. 
# In this case the input is a table with more colums, some of which needed for initial QC but not here. 

# Column 1 contains the chromosome (with chr prefix).
# Column 2 contains the total n. of files in the relative chr folder in DNAnexus.
# Column 3 contains the number of pVCF chunks for that chr (=tot n. files / 2).
# Column 4 contains chr size (as from the UCSC genome browser for hg38).
# Column 5 contains the average chunk size (chr size / number of chunks).
# Column 6 contains the exact n of 20kb chunks needed to cover the chr in question.
# Column 7 contains the number of whole 20kb chunks that fit into chr.
# Column 8 contains the n of bases covered by the last chunk.

# Cross-checks between columns 3, 6, 7 and 8 are sufficient to infer that each pVCF chunk in this release covers exactly 20kb.
# The only exception to this is the last chunk for each chromosome. Differently from the previous WGS release (200k), in this one
# there are no exceptions to this within specific chromosomal regions.

# NOTE/1: modify input/output directories accordingly
# NOTE/2: chunk_index_filepath = file "n_chunks_chr_size.tsv"

create_wgs_block_map_500k <- function(chunk_index_filepath){
  
  options(scipen=999)
  
  index = read.table(chunk_index_filepath,header=T,sep="\t")
  index[,"LAST_CHUNK"] = as.numeric(as.character(index[,"N_CHUNKS"]))-1
  
  start = c()
  end = c()
  chr = c()
  wgs_block = c()
  file_name = c()
  
  for (c in c(1:22,"X")){
    print(paste("Processing chromosome ",c,"...",sep=""))
    
    nchunks = index[which(index["CHR"]==paste0("chr",c)),"LAST_CHUNK"] #0-based
    clength = index[which(index["CHR"]==paste0("chr",c)),"CHR_SIZE"]
    
    chr = c(chr,rep(c,nchunks+1)) #the +1 is needed because blocks are numbered in a 0-based manner
    for(b in 0:nchunks){
      
      #define start and end coordinates for each block
      s = (b*20000)+1
      e = (b+1)*20000

      #for the last block, the last base will be the last base of the entire chr
      if(e >= clength){
        e = clength
      }
      else{
        #do nothing
      }
      
      #define block name for each block
      wb = paste0("b",b)
      
      #define file name for each block
      fn = paste0("ukb23374_c",c,"_",wb,"_v1.vcf.gz")
      
      start = c(start,s)
      end = c(end,e)
      wgs_block = c(wgs_block,wb)
      file_name = c(file_name,fn)
    }
  }
  block_map = cbind(chr,start,end,wgs_block,file_name)
  block_map = as.data.frame(block_map)
  block_map$chr = paste0("chr",block_map$chr)
  write.table(block_map,"WGS_500k_block_map.tsv",col.names=F,row.names=F,sep="\t",quote=F)
  }