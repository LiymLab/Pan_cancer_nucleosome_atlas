#ys-changed
import argparse
from multiprocessing import Pool
from pybedtools import BedTool
import numpy as np
import pandas as pd
import os

def chain_reverse(chain_bed,count_list):
    chain = chain_bed.tolist()
    index=0
    for flag in chain:
        if flag=="-":
            count_list[index]=np.flipud(count_list[index])
        index=index+1
    return count_list

def middle_bed(bed):
    bed_np = np.array(bed)
    bed_np = bed_np[:,0:4]  #change:3>4
    bed_np[:,1] = np.around((bed_np[:,1].astype(int)+bed_np[:,2].astype(int))/2).astype(int)
    bed_np[:,2] = bed_np[:,1].astype(int)+1
    bed_str = ''
    for i in bed_np:
        bed_str = bed_str + str(i).replace(" ","").replace("''","\t").replace("['","").replace("']","\n")
    return BedTool(bed_str,from_string  =True)

def pattern_count(array):
    d1 = array[2]-array[0]
    d2 = array[3]-array[1]
    out_array = np.zeros(2*slop+1).astype(int)
    if d1 < 0:
        d1 = 0
    if d2 >= 0:
        d2 = 2*slop
    else:
        d2 = 2*slop + d2
    out_array[d1:d2+1] = 1
    return out_array

def select_col(bed):
    bed_np=np.array(bed)
    bed_np = bed_np[:,[0,1,2]]
    bed_str = ''
    for i in bed_np:
        bed_str = bed_str + str(i).replace(" ","").replace("''","\t").replace("['","").replace("']","\n")
    return BedTool(bed_str,from_string  =True)

slop = 1000
denopa_array = np.zeros(2*slop+1).astype(int)
CT_list = ['ACC',"BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","PCPG","PRAD","SKCM","STAD","TGCT","THCA","UCEC"]
cancer_nuc_path = '/cancer_nucleosome_qc/script/file'
genome_path = "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
threads = 1

for CT in CT_list:
    print("####################"+CT+"######################")
    target_path='hg38_tss.txt'
    input_path=cancer_nuc_path+'/TCGA-'+CT+'.txt'
    input_bed = select_col(BedTool(input_path,from_string = False)) #Bedtools>>select_col
    target_bed = middle_bed(BedTool(target_path,from_string = False))
    
    intersect_bed = np.array(target_bed.slop(b = slop,g = genome_path).intersect(b = input_bed ,loj = True))
    intersect_bed = intersect_bed[intersect_bed[:,4] != '.',:]
    chain_bed=intersect_bed[:,3]
    intersect_bed = intersect_bed[:,[1,2,5,6]].astype(int)       #change:4>5;5>6
    
    if __name__ == '__main__':
        pool = Pool(threads)
        count_list = pool.map(pattern_count,intersect_bed)
        count_list = chain_reverse(chain_bed,count_list)
        count_np = np.array(count_list,dtype=np.int8)
        pool.close()

    count_np = np.sum(count_np,axis = 0)
    normalized_np = count_np/len(target_bed)
    np.savetxt('/cancer_nucleosome_qc/output/TSS/'+CT+'.txt', normalized_np,fmt='%.8f',delimiter='\t')        
    print(CT)
