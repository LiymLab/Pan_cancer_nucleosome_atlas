options(bedtools.path="/bedtools2/bin")
options(scipen=999)

library(bedtoolsr)
library(dplyr)
library(parallel)

# Identify pan-cancer conserved nucleosomes supported by both deNOPA and nucleoATAC
# and NucleoATAC across all 23 cancer types.

input_path<-"/cancer_type_nucleosome_path"
tools_list<-c("denopa","nucleoatac_map")
cancer_list<-list.files(file.path(input_path,"nucleoatac_map","cancer"))

#-----------------------------#
# 1. Initialize with the first cancer type
#-----------------------------#

denopa<-read.table(
  file.path(input_path,tools_list[1],"cancer",cancer_list[1],
            "cancer_nucleosome_danpos.txt"),
  sep="\t",header=T
)[,1:3]

nucleoatac<-read.table(
  file.path(input_path,tools_list[2],"cancer",cancer_list[1],
            "cancer_nucleosome_danpos.txt"),
  sep="\t",header=T
)[,1:3]

# Retain nucleosomes supported by both methods.
# An overlap is retained when the intersection length is >=73 bp.
overlap<-bt.intersect(a=denopa,b=nucleoatac)
overlap<-subset(overlap,V3-V2>=73)

pan_cancer<-data.frame(
  chr=overlap[,1],
  start=round((overlap[,2]+overlap[,3])/2-73),
  end=round((overlap[,2]+overlap[,3])/2+73)
)

#-----------------------------#
# 2. Intersect across all 23 cancer types
#-----------------------------#

for(ct in cancer_list[2:23]){
  print(ct)

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

  # Identify nucleosomes supported by both methods
  overlap<-bt.intersect(a=denopa,b=nucleoatac)
  overlap<-subset(overlap,V3-V2>=73)

  ct_nucleosome<-data.frame(
    chr=overlap[,1],
    start=round((overlap[,2]+overlap[,3])/2-73),
    end=round((overlap[,2]+overlap[,3])/2+73)
  )

  # Retain nucleosomes shared with all previously processed cancer types
  pan_cancer<-bt.merge(
    bt.sort(
      bt.intersect(
        a=pan_cancer,
        b=ct_nucleosome
      )
    )
  )

  print(nrow(pan_cancer))
}

#-----------------------------#
# 3. Save the pan-cancer conserved nucleosome set
#-----------------------------#

output_dir<-"/pan_cancer_nucleosome/file"

write.table(
  pan_cancer,
  file=file.path(output_dir,"pan_cancer_conserved_nucleosome_overlap.txt"),
  col.names=F,row.names=F,quote=F,sep="\t"
)

# Convert merged regions to fixed 146-bp nucleosome intervals
pan_cancer_nucleosome<-data.frame(
  chr=pan_cancer[,1],
  start=round((pan_cancer[,2]+pan_cancer[,3])/2-73),
  end=round((pan_cancer[,2]+pan_cancer[,3])/2+73)
)

write.table(
  pan_cancer_nucleosome,
  file=file.path(output_dir,"pan_cancer_conserved_nucleosome.txt"),
  col.names=F,row.names=F,quote=F,sep="\t"
)
