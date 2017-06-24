import sys
import numpy as np
#import csv as csv
#import h5py
sys.path.append('../')
from Public_modules.Utilities import *
from Public_modules.Utilities import *
import pandas
import matplotlib.pyplot as plt
 
def interesect_protein_mrna (folder_data_protein='', file_protein_quantification='', file_protein_meta='',\
                       folder_data_rna='',file_rna_quantification='',file_rna_meta='',intersect_lines=None,intersect_genes=None):

#    folder_data_protein='/Users/mirauta/Data/MS/hipsci/phenotypes/'
#    file_protein_quantification='hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517_protein_Reporter intensity corrected_ranknorm_regbatch_genes_filtered'
#    file_protein_meta='hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517_protein_metadata_genes_filtered'
##    hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517_lines_metadata'
#    folder_data_rna='/Users/mirauta/Data/RNA/hipsci/'
#    
#    file_rna_quantification='IPSc-Full.featureCounts.genes.counts.unique.tsv_counts.CorrExpressedGenes.tsv.gz'
#    file_rna_meta='annotationFiles/Ensembl_75_Limix_Annotation_FC_Gene.txt'
#    file_rna_quantification='HipSci.ISR.Salmon_transcripts_TPM_Corrected.tsv.gz'
#    file_rna_meta='annotationFiles/Ensembl_75_Limix_Annotation_Salmon_Tx.txt'
#    
    data={}
    data['protein']=pandas.read_table(folder_data_protein+file_protein_quantification+".txt",sep='\t',index_col=0).transpose()
    data['peptide']=pandas.read_table(folder_data_protein+file_protein_quantification.replace('_protein','_peptide')+".txt",sep='\t',index_col=0).transpose()
    data['peptide_meta']=pandas.read_table(folder_data_protein+file_protein_meta.replace('_protein','_peptide')+".txt",sep='\t').set_index('ensembl_gene_id',drop=False).transpose() 
    data['protein_meta']=pandas.read_table(folder_data_protein+file_protein_meta+".txt",sep='\t').set_index('ensembl_gene_id',drop=False).transpose()
#   
    data['rna']=pandas.read_table(folder_data_rna+file_rna_quantification,sep='\t',index_col=0).transpose()
    try:
        data['rna_meta']=pandas.read_table(folder_data_rna+file_rna_meta,sep='\t').set_index('feature_id',drop=False).transpose()
        data['rna_meta']=data['rna_meta'].rename(index={'gene_id':'ensembl_gene_id'}, columns=str)
    except:
        data['rna_meta']=pandas.read_table(folder_data_rna+file_rna_meta,sep='\t').set_index('ProbeName',drop=False).transpose()
        data['rna_meta']=data['rna_meta'].rename(index={'ProbeName':'feature_id','gene_id':'ensembl_gene_id'}, columns=str).transpose()
        data['rna_meta']['ensembl_gene_id']=data['rna_meta']['feature_id']
        data['rna_meta']=data['rna_meta'].transpose()
    
    try: data['rna'].index=np.array([l.split('/')[3].split('.')[0]for l in data['rna'].index])
    except:
        data['rna'].index=np.array(['-'.join(l.split('.')[:-1]) for l in data['rna'].index])
        
    data['rna']=data['rna'].transpose()
    data['rna'].index=np.array([ g.split('_')[0] for g in data['rna'].index])
    data['rna']=data['rna'].transpose()
    
    
    common_lines=np.sort(np.intersect1d(data['protein'].index,data['rna'].index)); print (common_lines.shape)
    diff_lines=np.setdiff1d(data['protein'].index,data['rna'].index); print (diff_lines.shape)
    common_genes=np.sort(np.intersect1d(data['peptide_meta'].transpose()['ensembl_gene_id'],\
                                        np.intersect1d(data['protein_meta'].transpose()['ensembl_gene_id'],\
                                                       data['rna_meta'].transpose()['ensembl_gene_id']))); print (common_genes.shape)
    ## select common lines                       
    if intersect_lines:
        for x in ['protein','peptide',  'rna']:
            data[x]=data[x].transpose()[common_lines].transpose()
            print (data[x].shape)

    
    for x in ['protein','peptide',  'rna']:
        print (data[x+'_meta'].shape)
        
    ## select common genes
    if intersect_genes:
        for x in ['protein','peptide',  'rna']:
            data[x+'_meta']=data[x+'_meta'].transpose().set_index('ensembl_gene_id',drop=False).transpose()[common_genes].transpose().set_index('feature_id',drop=False).transpose()
            print (data[x+'_meta'].shape)
            data[x]= data[x][np.intersect1d(data[x+'_meta'].columns.values,data[x].columns.values)]
            print (data[x].shape)
            
    data['protein'].transpose().to_csv(folder_data_protein+file_protein_quantification+"_intersect_prot_rna.txt",mode='w', sep='\t', columns=None, header=True, index=True)
    data['peptide'].transpose().to_csv(folder_data_protein+file_protein_quantification.replace('_protein','_peptide')+"_intersect_prot_rna.txt",mode='w', sep='\t', columns=None, header=True, index=True)
    data['peptide_meta'].transpose().to_csv(folder_data_protein+file_protein_meta.replace('_protein','_peptide')+"_intersect_prot_rna.txt",mode='w', sep='\t', columns=None, header=True, index=True)
    data['protein_meta'].transpose().to_csv(folder_data_protein+file_protein_meta+"_intersect_prot_rna.txt",mode='w', sep='\t', columns=None, header=True, index=True)
#   
    data['rna'].transpose().to_csv(folder_data_rna+file_rna_quantification+'_intersect_prot_rna.txt',mode='w', sep='\t', columns=None, header=True, index=True)
    data['rna_meta'].transpose().to_csv(folder_data_rna+file_rna_meta+'_intersect_prot_rna.txt',mode='w', sep='\t', columns=None, header=True, index=True)
    

    return data


folder_data_protein='/Users/mirauta/Data/MS/hipsci/phenotypes/'
file_protein_quantification='hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517_protein_Reporter intensity corrected_ranknorm_regbatch_genes_filtered'
file_protein_meta='hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517_protein_metadata_genes_filtered'

folder_data_rna='/Users/mirauta/Data/RNA/hipsci/'
file_rna_quantification='IPSc-Full.featureCounts.genes.counts.unique.tsv_counts.CorrExpressedGenes.tsv.gz'
file_rna_meta='annotationFiles/Ensembl_75_Limix_Annotation_FC_Gene.txt'
#file_rna_quantification='HipSci.ISR.Salmon_transcripts_TPM_Corrected.tsv.gz'
#file_rna_meta='annotationFiles/Ensembl_75_Limix_Annotation_Salmon_Tx.txt'
# 
interesect_protein_mrna (folder_data_protein, file_protein_quantification, file_protein_meta, folder_data_rna,file_rna_quantification,file_rna_meta,intersect_lines=True,intersect_genes=True)

       