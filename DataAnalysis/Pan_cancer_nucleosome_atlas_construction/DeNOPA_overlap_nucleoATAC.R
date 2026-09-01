options(bedtools.path="/bedtools2/bin")
options(scipen=999)

library(bedtoolsr)
library(dplyr)

# Integrate nucleosome calls from deNOPA and NucleoATAC
# for each cancer type. Nucleosomes supported by both methods
# (overlap >=73 bp) are retained and converted to 146-bp intervals.

input_path<-"/nucleosome_path"
output_path<-"/nucleosome_path/overlap"

tools_list<-c("denopa","nucleoatac_map")
cancer_list<-list.files(file.path(input_path,"nucleoatac_map","cancer"))


for(ct in cancer_list){
  print(ct)

  # Read nucleosome calls from denopa and NucleoATAC
  denopa<-read.table(
    file.path(input_path,tools_list[1],"cancer",ct,
              "cancer_nucleosome_danpos.txt"),
    sep="\t",header=T
  )[,1:3]

  nucleoatac<-read.table(
    file.path(input_path,tools_list[2],"cancer",ct,
              "cancer_nucleosome_danpos.txt"),
    sep="\t",header=T
  )[,1:3]

  # Retain nucleosomes supported by both methods
  overlap<-bt.intersect(a=denopa,b=nucleoatac)
  overlap<-subset(overlap,V3-V2>=73)

  # Convert the overlap to a 146-bp nucleosome interval
  ct_nucleosome<-data.frame(
    chr=overlap[,1],
    start=round((overlap[,2]+overlap[,3])/2-73),
    end=round((overlap[,2]+overlap[,3])/2+73)
  )

  # Save the integrated nucleosome set for each cancer type
  write.table(
    ct_nucleosome,
    file=file.path(
      output_path,
      paste0(gsub("TCGA-","",ct),"_nucleosome.txt")
    ),
    col.names=F,row.names=F,quote=F,sep="\t"
  )
}
