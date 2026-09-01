options(bedtools.path="/bedtools2/bin")
library(bedtoolsr)
library(dplyr)
library(parallel)

# function: delete the overlap of top_nucleosome from nucleosome_list
delete_overlap<-function(nucleosome_list,top_nucleosome){
  delete_list<-bt.intersect(a=nucleosome_list,b=top_nucleosome,wo=T)[,c(1:4)]
  colnames(delete_list)<-c("chrom","start","end","score")
  nucleosomes<-anti_join(nucleosome_list,delete_list,by=c("chrom","start","end","score"))
  return(nucleosomes)
}

# function: extract nucleosome center to facilitate subsequent processing
extract_center_nucleosome<-function(nucleosome){
  sort_nucleosome<-bt.sort(nucleosome)
  center_nucleosome<-data.frame(chrom=sort_nucleosome$V1,center=sort_nucleosome$V2+73,score=sort_nucleosome$V4)
  down_distance<-c()
  up_distance<-c()

  for(chrome in unique(center_nucleosome$chrom)){
    chr_center<-center_nucleosome[center_nucleosome$chrom==chrome,]$center
    if(length(chr_center)>=2){
      subtract_center<-c(chr_center[2:length(chr_center)],chr_center[length(chr_center)]+73)
      down<-subtract_center-chr_center
      down_distance<-c(down_distance,down)
      up<-c(73,down[1:(length(down)-1)])
      up_distance<-c(up_distance,up)
    }else{
      down_distance<-c(down_distance,73)
      up_distance<-c(up_distance,73)
    }
  }
  center_nucleosome$d_distance<-down_distance
  center_nucleosome$u_distance<-up_distance
  return(center_nucleosome)
}

# function: From the input nucleosome(1-chr,2-start,3-end,4-p), get the nucleosomes that do not overlap with any other nucleosomes
non_overlap_nucleosome_center<-function(nucleosome){
  sort_nucleosome<-bt.sort(nucleosome)
  center_nucleosome<-data.frame(chrom=sort_nucleosome[,1],center=sort_nucleosome[,2]+73,score=sort_nucleosome[,4])
  down_distance<-c()
  up_distance<-c()

  for(chrome in unique(center_nucleosome$chrom)){
    chr_center<-center_nucleosome[center_nucleosome$chrom==chrome,]$center
    if(length(chr_center)>=2){
      subtract_center<-c(chr_center[2:length(chr_center)],chr_center[length(chr_center)]+73)
      down<-subtract_center-chr_center
      down_distance<-c(down_distance,down)
      up<-c(73,down[1:(length(down)-1)])
      up_distance<-c(up_distance,up)
    }else{
      down_distance<-c(down_distance,73)
      up_distance<-c(up_distance,73)
    }
  }

  center_nucleosome$d_distance<-down_distance
  center_nucleosome$u_distance<-up_distance
  nonoverlap_center<-center_nucleosome[center_nucleosome$d_distance>=146&center_nucleosome$u_distance>=146,]
  return(nonoverlap_center)
}

# function: parallel by chromatin
parallel_by_chr<-function(chromatin_chr){
  bed<-subset(all_bed,chrom==chromatin_chr)
  center_nucleosome<-data.frame(chrom=bed$chrom,center=bed$center,score=bed$score)
  cancer_chr_nucleosome<-data.frame(chrom=bed$chrom,start=bed$center-73,end=bed$center+73,score=bed$score)

  nonoverlap_center<-non_overlap_nucleosome_center(cancer_chr_nucleosome)
  nonoverlap_nucleosome<-data.frame(chrom=nonoverlap_center$chrom,start=nonoverlap_center$center-73,end=nonoverlap_center$center+73,score=nonoverlap_center$score)

  center_nucleosome<-extract_center_nucleosome(cancer_chr_nucleosome)
  if(nrow(nonoverlap_nucleosome)!=0){
    deoverlap_center<-anti_join(center_nucleosome,nonoverlap_center,by=c("chrom","center","score","d_distance","u_distance"))
  }else{
    deoverlap_center<-center_nucleosome
  }

  deoverlap_nucleosome<-data.frame(chrom=deoverlap_center$chrom,start=deoverlap_center$center-73,end=deoverlap_center$center+73,score=deoverlap_center$score)
  cancer_sort_nucleosome<-deoverlap_nucleosome[order(deoverlap_nucleosome$score,decreasing=TRUE),]

  if(length(cancer_sort_nucleosome[,1])>10000){
    span<-10000
  }else if(length(cancer_sort_nucleosome[,1])>1000){
    span<-1000
  }else if(length(cancer_sort_nucleosome[,1])>100){
    span<-100
  }else{
    span<-10
  }

  top_nucleosome_list=data.frame()
  top_nucleosome_list<-rbind(top_nucleosome_list,nonoverlap_nucleosome)

  while(length(cancer_sort_nucleosome[,1])>(span+1)){
    step=span
    cancer_sort_nucleosome<-cancer_sort_nucleosome[order(cancer_sort_nucleosome$score,decreasing=TRUE),]
    top_span_nucleosome<-cancer_sort_nucleosome[c(1:step),]
    top_nucleosome_list_step=data.frame()
    nonoverlap_center_top_span<-non_overlap_nucleosome_center(top_span_nucleosome)
    nonoverlap_nucleosome_top_span<-data.frame(chrom=nonoverlap_center_top_span$chrom,start=nonoverlap_center_top_span$center-73,end=nonoverlap_center_top_span$center+73,score=nonoverlap_center_top_span$score)

    if(nrow(nonoverlap_nucleosome_top_span)!=0){
      deoverlap_span_nucleosome<-anti_join(top_span_nucleosome,nonoverlap_nucleosome_top_span,by=c("chrom","start","end","score"))
    }else{
      deoverlap_span_nucleosome<-top_span_nucleosome
    }

    while(length(deoverlap_span_nucleosome[,1])!=0){
      top_nucleosome<-deoverlap_span_nucleosome[which.max(deoverlap_span_nucleosome$score),]
      top_nucleosome_list_step<-rbind(top_nucleosome_list_step,top_nucleosome)
      deoverlap_span_nucleosome<-delete_overlap(deoverlap_span_nucleosome,top_nucleosome)
      print(paste("Left Nucleosomes:",length(deoverlap_span_nucleosome[,1])))
    }

    top_nucleosome_list_step<-rbind(top_nucleosome_list_step,nonoverlap_nucleosome_top_span)
    top_nucleosome_list<-rbind(top_nucleosome_list,top_nucleosome_list_step)
    cancer_sort_nucleosome<-cancer_sort_nucleosome[c((step+1):(length(cancer_sort_nucleosome[,1]))),]

    delete_nucleosome<-bt.intersect(a=cancer_sort_nucleosome,b=top_nucleosome_list_step,wo=T)
    if(nrow(delete_nucleosome)>0){
      delete_nucleosome<-unique(delete_nucleosome[,c(1:4)])
      colnames(delete_nucleosome)<-c("chrom","start","end","score")
      cancer_sort_nucleosome<-anti_join(cancer_sort_nucleosome,delete_nucleosome,by=c("chrom","start","end","score"))
      print(length(cancer_sort_nucleosome[,1]))
    }

    if(nrow(cancer_sort_nucleosome)!=0){
      sub_nonoverlap_center<-non_overlap_nucleosome_center(cancer_sort_nucleosome)
      if(nrow(sub_nonoverlap_center)!=0){
        sub_nonoverlap_nucleosome<-data.frame(chrom=sub_nonoverlap_center$chrom,start=sub_nonoverlap_center$center-73,end=sub_nonoverlap_center$center+73,score=sub_nonoverlap_center$score)
        sub_center_nucleosome<-extract_center_nucleosome(cancer_sort_nucleosome)
        deoverlap_center<-anti_join(sub_center_nucleosome,sub_nonoverlap_center,by=c("chrom","center","score","d_distance","u_distance"))
        sub_deoverlap_nucleosome<-data.frame(chrom=deoverlap_center$chrom,start=deoverlap_center$center-73,end=deoverlap_center$center+73,score=deoverlap_center$score)
        top_nucleosome_list<-rbind(top_nucleosome_list,sub_nonoverlap_nucleosome)
        cancer_sort_nucleosome<-sub_deoverlap_nucleosome
      }
    }
    print(paste("Nucleosome left",length(cancer_sort_nucleosome[,1])))
  }

  while(length(cancer_sort_nucleosome[,1])!=0){
    top_nucleosome<-cancer_sort_nucleosome[which.max(cancer_sort_nucleosome$score),]
    top_nucleosome_list<-rbind(top_nucleosome_list,top_nucleosome)
    cancer_sort_nucleosome<-delete_overlap(cancer_sort_nucleosome,top_nucleosome)
    print(paste("Left Nucleosomes:",length(cancer_sort_nucleosome[,1])))
  }

  cancer_specific_nucleosome<-top_nucleosome_list
  write.table(cancer_specific_nucleosome,file=paste(out_de_overlap_path,paste(chromatin_chr,"_deNOPA_cancer_nucleosome.txt",sep=""),sep="/"),sep="\t",quote=FALSE,row.names=FALSE)
}

combine_bed<-function(chromatin_list){
  cancer_nucleosome_all_merge=data.frame()
  for(chromatin_sub in chromatin_list){
    chr_nucleosome_sub<-read.table(paste(out_de_overlap_path,paste(chromatin_sub,"_deNOPA_cancer_nucleosome.txt",sep=""),sep="/"),sep="\t",header=T)
    cancer_nucleosome_all_merge<-rbind(cancer_nucleosome_all_merge,chr_nucleosome_sub)
    file.remove(paste(out_de_overlap_path,paste(chromatin_sub,"_deNOPA_cancer_nucleosome.txt",sep=""),sep="/"))
  }
  write.table(cancer_nucleosome_all_merge,file=paste(out_de_overlap_path,"deNOPA_cancer_nucleosome.txt",sep="/"),sep="\t",quote=FALSE,row.names=FALSE)
}

# read in the cancer-type path
out_overlap_path="/denopa/cancer"
cancer_path="/denopa/sample/"

cancer_info<-read.table("cancer_info.txt",header=T,sep="\t")
cancer_info_df<-data.frame(sample=gsub("_atacseq_gdc_realn.bam","",cancer_info$File.Name),cancer=cancer_info$Project)

for(cancer_type_info in sort(unique(cancer_info_df$cancer))[seq_start:seq_end]){
  print(cancer_type_info)
  out_de_overlap_path<-paste(c(out_overlap_path,"/",cancer_type_info),collapse="")
  sample_list_files<-subset(cancer_info_df,cancer==cancer_type_info)$sample
  cancer_nucleosome<-data.frame()

  if(!dir.exists(out_de_overlap_path)){
    dir.create(out_de_overlap_path)
  }

  for(sample in sample_list_files){
    print(sample)
    sample_nucleosome<-read.csv(paste(cancer_path,paste(sample,"/deNOPA_sample_nucleosome.txt",sep=""),sep=""),sep="\t",header=T)
    cancer_nucleosome<-rbind(cancer_nucleosome,sample_nucleosome)
  }

  cancer_nucleosome<-unique(cancer_nucleosome)
  cancer_nucleosome_sort_bychr<-bt.sort(cancer_nucleosome)
  colnames(cancer_nucleosome_sort_bychr)<-c("chrom","start","end","score")
  chromatin<-unique(cancer_nucleosome_sort_bychr$chrom)

  all_bed<-data.frame(chrom=cancer_nucleosome_sort_bychr$chrom,center=cancer_nucleosome_sort_bychr$start+73,score=cancer_nucleosome_sort_bychr$score)

  print(date())
  cl<-makeCluster(30)
  clusterExport(cl,varlist=c("all_bed","chromatin","non_overlap_nucleosome_center","delete_overlap","extract_center_nucleosome","out_de_overlap_path","combine_bed"))
  clusterEvalQ(cl,{
    options(bedtools.path="/bedtools2/bin")
    library(bedtoolsr)
    library(dplyr)
  })
  clusterApply(cl=cl,x=chromatin,fun=parallel_by_chr)
  combine_bed(chromatin)
  stopCluster(cl)
  print(date())

  # only include nucleosomes that overlap with 2 or more samples
  merge_path<-paste(c(out_overlap_path,"/",cancer_type_info,"/deNOPA_cancer_nucleosome.txt"),collapse="")
  merge_nucleosome<-read.csv(merge_path,sep="\t",header=T)
  merge_nucleosome_sort<-bt.sort(merge_nucleosome)
  count<-rep(0,length(merge_nucleosome_sort[,1]))
  sample_list_files<-subset(cancer_info_df,cancer==cancer_type_info)$sample

  for(sample in sample_list_files){
    sample_nucleosome<-read.csv(paste(cancer_path,paste(sample,"/deNOPA_sample_nucleosome.txt",sep=""),sep=""),sep="\t",header=T)
    overlap_nucleosome<-unique(bt.intersect(a=merge_nucleosome_sort,b=sample_nucleosome,wo=T)[,c(1:4)])
    nonoverlap_nucleosome<-anti_join(merge_nucleosome_sort,overlap_nucleosome,by=c("V1","V2","V3","V4"))
    overlap_nucleosome$V5<-rep(1,length(overlap_nucleosome[,1]))
    nonoverlap_nucleosome$V5<-rep(0,length(nonoverlap_nucleosome[,1]))
    count_nucleosome<-rbind(overlap_nucleosome,nonoverlap_nucleosome)
    count_nucleosome_sort<-bt.sort(count_nucleosome)
    colnames(count_nucleosome_sort)<-c("chrom","start","end","score","count")
    count<-count+count_nucleosome_sort$count
  }

  colnames(merge_nucleosome_sort)<-c("chrom","start","end","score")
  merge_nucleosome_sort$V5<-count
  colnames(merge_nucleosome_sort)<-c("chrom","start","end","score","count")

  merge_nucleosome_sort_extract<-subset(merge_nucleosome_sort,count>=2)
  out_dir_overlap<-paste(c(out_overlap_path,"/",cancer_type_info,"/cancer_nucleosome_danpos.txt"),collapse="")
  write.table(merge_nucleosome_sort_extract,file=out_dir_overlap,sep="\t",quote=FALSE,row.names=FALSE)
  print("done1")
}
