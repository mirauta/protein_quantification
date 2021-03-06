import sys
import numpy as np
import csv as csv
import h5py
sys.path.append('../')
from Public_modules.Utilities import *
import pandas

######################################################
## load folder_name, batch names
######################################################
print ('Transforms the maxquant peptide files into the input for QTL processing')

#parameters_file=sys.argv[1]
#norm_method=sys.argv[2]
parameters_file='../parameters_quantification'
norm_method='qnorm'
folder_dest,folder_data,batch_name_tmt,file_protein_quantification,line_keyword_exclude_lines,line_keyword_exclude_batch=read_parameters(parameters_file)

field_data='Reporter intensity corrected'
field_data2='Reporter intensity corrected'+'_regbatch'
field_meta_pep=np.array(['Leading razor protein','Sequence','Proteins','Gene names', 'Unique (Groups)','PEP','Score','Reverse', 'Potential contaminant','Start position'])
field_meta_pro=np.array(['Majority protein IDs', 'Protein IDs','Peptide IDs','Peptide counts (razor+unique)','Gene names', 'Peptide is razor','Mol. weight [kDa]','Only identified by site'])

#==============================================================================
#==============================================================================
# # ### read data in dictionary
#==============================================================================
#==============================================================================
prodata=load_data_from_maxquant_output(folder_data=folder_data,name=batch_name_tmt[0],data_file='.proteinGroups.txt', field_data=field_data, field_metadata=field_meta_pro)
pepdata=load_data_from_maxquant_output(folder_data=folder_data,name=batch_name_tmt[0],data_file='.peptides.txt',      field_data=field_data, field_metadata=field_meta_pep)

'''
filter proteins and peptides by 'Reverse', 'Potential contaminant'
'''
pepdata['meta']['valid']=(pepdata['meta']['Potential contaminant']=='')&(pepdata['meta']['Reverse']=='')& np.array(['REV_' not in a for a in pepdata['meta']['Leading razor protein']])& np.array(['CON_' not in a for a in pepdata['meta']['Leading razor protein']])
prodata['meta']['valid']=np.array(['REV_' not in a.split(";")[0] for a in prodata['meta']['Majority protein IDs']])&np.array(['CON_' not in a.split(";")[0] for a in prodata['meta']['Majority protein IDs']])
prodata['meta']['Lead protein']=np.array([p.split(";")[0] for p in prodata['meta']['Majority protein IDs']])
prodata['meta']['valid']=prodata['meta']['valid']&(np.in1d(prodata['meta']['Lead protein'],pepdata['meta']['Leading razor protein'][pepdata['meta']['valid']]))


'''
remove non valid peptides/ proteins
'''
for x in [prodata,pepdata]:
    index=x['meta']['valid']
    for key in x['meta'].keys():
        x['meta'][key]=x['meta'][key][index]
    x['data'][field_data]=x['data'][field_data][:,index]
     
prodata['data']['lines']=np.array([l.replace(field_data,'').split(' ')[1]for l in prodata['data']['fields']])
prodata['data']['batch']=np.copy(prodata['data']['lines']);#np.zeros(len(prodata['data']['lines']),dtype='U32')
for il,l in enumerate(prodata['data']['fields']):
    temp=l.replace(field_data,'').replace(prodata['data']['lines'][il],'').replace(' ','')
    if temp!='':prodata['data']['batch'][il]=temp
  
'''
filter lines by median number of peptides detected
'''
temp=(np.nansum(pepdata['data'][field_data]>0,1))
#temp2=(np.nansum(prodata['data'][field_data],1))
prodata['data']['valid_lines']=(temp<=np.nanmedian(temp)*1.25)&(temp>=np.nanmedian(temp)*0.75) &\
                               np.all(np.array([np.array([field_filter not in l for field_filter in line_keyword_exclude_lines]) for l in prodata['data']['lines']]),1)&\
                               np.all(np.array([np.array([field_filter not in l for field_filter in line_keyword_exclude_batch]) for l in prodata['data']['batch']]),1)&\
                               np.array(['HPSI' in l for l in prodata['data']['lines']])
    
             
#==============================================================================
#==============================================================================
# # 
#==============================================================================
#==============================================================================
'''
bring annotation merge on protein id
From Uniprot ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/
or 
From Ensembl Homo_sapiens.GRCh37.75_gene.gtf
'''
def mergeid(key1=None,key2=None):
    index=[np.where(key2==pr)[0] for ipr,pr in enumerate(key1)]
    index2=np.zeros(len(key2))+np.nan
    for i,ind in enumerate(index):
              index2[ind]=i
    return index2

#ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/
#temp=load_data(file_name="/Users/mirauta/Data/Annotation/Uniprot/UP000005640_9606_proteome.bed",delimiter='\t', quotechar='"')[1]
#annotation={}
#annotation['chromosome']=temp[0];annotation['annotation_protein']=temp[3];annotation['start']=temp[1];annotation['end']=temp[2];annotation['strand']=temp[5];
#temp=np.copy(annotation['annotation_protein'])
#for key in annotation.keys():    annotation[key]=annotation[key][np.unique(temp,return_index=1)[1]]
#temp=np.copy(annotation['annotation_protein'])
#for key in annotation.keys():   annotation[key]=annotation[key][np.in1d(temp,prodata['meta']['Lead protein'])]
#prodata['meta']['index_annotation']=mergeid(key1=annotation['annotation_protein'],key2=prodata['meta']['Lead protein'])
#pepdata['meta']['index_annotation']=mergeid(key1=annotation['annotation_protein'],key2=pepdata['meta']['Leading razor protein'])

tempi=load_data(file_name="/Users/mirauta/Data/Annotation/Ensembl_37.75/Homo_sapiens.GRCh37.75.gtf",delimiter='\t', quotechar='"',skip_before=5)[1]
temp=np.array([t[tempi[2]=='gene']for t in tempi])
annotation={}
annotation['chromosome']=temp[0];
annotation['annotation_gene']=np.array([g.split(";")[1].replace('gene_name','').replace(' ','').replace('"','') for g in temp[8]])
annotation['annotation_ensembl_id']=np.array([g.split(";")[0].replace('gene_id','').replace(' ','').replace('"','') for g in temp[8]])
annotation['start']=temp[3].astype(int);
annotation['end']=temp[4].astype(int);
annotation['strand']=temp[6];


temp=np.copy(annotation['annotation_gene'])
for key in annotation.keys():  
    annotation[key]=annotation[key][np.unique(temp,return_index=1)[1]]
temp=np.copy(annotation['annotation_gene'])
for key in annotation.keys():  
    annotation[key]=annotation[key][np.in1d(temp,prodata['meta']['Gene names'])]

prodata['meta']['index_annotation']=mergeid(key1=annotation['annotation_gene'],key2=prodata['meta']['Gene names'])
pepdata['meta']['index_annotation']=mergeid(key1=annotation['annotation_gene'],key2=pepdata['meta']['Gene names'])

for x in [prodata,pepdata]:
    for key in annotation.keys():
        x['meta'][key]=np.zeros(len(x['meta']['index_annotation']),dtype=annotation[key].dtype)
        x['meta'][key][x['meta']['index_annotation']==x['meta']['index_annotation']]= annotation[key][x['meta']['index_annotation'][x['meta']['index_annotation']==x['meta']['index_annotation']].astype(int)]

for x in [prodata,pepdata]:
    x['meta']['chromosome']=np.array([xx.replace('chr','') for xx in x['meta']['chromosome']])
    
for x in [prodata,pepdata]:
    x['meta']['valid']=x['meta']['valid']& (x['meta']['chromosome']!='')&(x['meta']['start']!='')&(x['meta']['end']!='')

'''
keep data with available chromosome and start/ end positions
'''
for x in [prodata,pepdata]:
    index=x['meta']['valid']
    for key in x['meta'].keys():
        x['meta'][key]=x['meta'][key][index]
    x['data'][field_data]=x['data'][field_data][:,index]
#==============================================================================
# remove the second gene name
#==============================================================================
for x in [prodata,pepdata]:
    x['meta']['Gene names']=np.array([xx.split(';')[0] for xx in  x['meta']['Gene names']])
#==============================================================================
#==============================================================================
# # 
#==============================================================================
'''
rank normalise lines
'''
#field_data3=field_data2+'_qnorm'
#for x in [prodata,pepdata]:
#    x['data'][field_data3]=np.array([transform_vector_qnorm(y, method='qnorm') for y in     x['data'][field_data2].T]).T
#    
import matplotlib.pyplot as plt
#
#L=35
for temp in [pepdata['data'],prodata['data']]:
    xref=temp['Reporter intensity corrected'][prodata['data']['valid_lines']]\
             [np.argmax(np.nansum(temp['Reporter intensity corrected'][prodata['data']['valid_lines']]>0,1))]
    temp['Reporter intensity corrected_ranknorm']=\
        np.array([rank_transform_new(x=temp['Reporter intensity corrected'][i],\
        ref=xref)    for i in np.arange(temp['Reporter intensity corrected'].shape[0])])
    
#plt.figure(figsize=(7,14))
#for i in np.arange(L):
#    plt.subplot(L,1,i+1)
#    plt.hist([np.log(temp[i][temp[i]>0]),np.log(temp2[i][temp2[i]>0])],bins=20)
#    plt.xlim(5,20)
#    
#plt.show()

#plt.xlim(4,25)
#sys.exit()

'''
remove batch effect
'''
for field_data in ['Reporter intensity corrected','Reporter intensity corrected_ranknorm']:
    for x in [prodata,pepdata]:
        x['data'][field_data][x['data'][field_data]==0]=np.nan
       
    prodata['data'][field_data+'_regbatch']=np.copy(prodata['data'][field_data])
    pepdata['data'][field_data+'_regbatch']=np.copy(pepdata['data'][field_data])
    for x in [prodata,pepdata]:
        medianpeptides=np.nanmedian(1+x['data'][field_data],0)
        for b in np.unique(prodata['data']['batch']):
            index=np.where(prodata['data']['batch']==b)[0]
            if index.shape[0]>3:
                print (b)
                x['data'][field_data+'_regbatch'][index]=x['data'][field_data][index]/\
                 np.nanmedian(1+x['data'][field_data][index],0)[None]*\
                 medianpeptides[None]
    

for field_data in ['Reporter intensity corrected'+'_regbatch','Reporter intensity corrected_ranknorm'+'_regbatch']:
    temp=pepdata['data'][field_data]
    plt.hist([np.log(temp[i][temp[i]>0])for i in [1,5,10,15,20,25]],bins=10,normed=1)
    plt.show()
'''
generate qnorm phenotype
'''
 
for field_data in ['Reporter intensity corrected','Reporter intensity corrected_ranknorm']:
    field_data3=field_data+'_regbatch'+'_qnorm'
    for x in [prodata,pepdata]:
        x['data'][field_data3]=np.array([transform_vector_qnorm(y, method='qnorm') for y in     x['data'][field_data+'_regbatch'].T]).T
    
#sys.exit()

'''
select independent samples
'''


prodata['data']['donor_valid_lines']=np.array([l.split('-')[-1].split('_')[0] for l in prodata['data']['lines'][prodata['data']['valid_lines']]])
prodata['data']['countpep_valid_lines']=np.nansum(pepdata['data'][field_data+'_regbatch'][prodata['data']['valid_lines']]>0,1)
prodata['data']['select_valid_lines']=np.zeros(len(prodata['data']['donor_valid_lines']),dtype='bool')
for d in np.unique(prodata['data']['donor_valid_lines']):
    prodata['data']['select_valid_lines'][prodata['data']['donor_valid_lines']==d]=\
           np.arange((prodata['data']['donor_valid_lines']==d).sum())==\
                  np.argmax(prodata['data']['countpep_valid_lines'][prodata['data']['donor_valid_lines']==d])

#==============================================================================
kin=np.array(load_data(file_name="/Users/mirauta/Data/MS/hipsci/TMT/genotype_maf10/plink.rel",delimiter='\t', quotechar='"',skip_before=0)[1]).astype(float)
names=np.array(load_data(file_name="/Users/mirauta/Data/MS/hipsci/TMT/genotype_maf10/plink.rel.id",delimiter='\t', quotechar='"',skip_before=0)[1])[1]
inphenotype=np.in1d(names, prodata['data']['lines'][prodata['data']['valid_lines']][prodata['data']['select_valid_lines']])
kinship=kin[inphenotype][:,inphenotype]
names=names[inphenotype]
#prodata['data']['select_valid_lines']=np.zeros(len(prodata['data']['donor_valid_lines']),dtype='bool')
relatedindividuals=names[np.triu(kinship,k=1).max(1)>0.05]

prodata['data']['valid_lines_independent']=np.copy(prodata['data']['valid_lines'])
prodata['data']['valid_lines_independent'][prodata['data']['valid_lines']]= prodata['data']['select_valid_lines']
prodata['data']['valid_lines_independent'][np.in1d(prodata['data']['lines'],relatedindividuals)]=False


       

'''
select genes by # detected lines and not overlaping snps
'''
pepdata['data']['peptide_overlap_snp']=pandas.read_table(folder_data+\
"phenotypes/hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517_peptide_lines_filtered_unique_overlap_peptide_snp.txt",sep='\t').set_index('feature_id',drop=0)
    
prodata['data']['at_least_x_lines']=np.nansum(prodata['data'][field_data+'_regbatch'][prodata['data']['valid_lines_independent']]>0,0)>(prodata['data']['valid_lines_independent'].sum()*0.66)

pepdata['data']['at_least_x_lines']=np.in1d(pepdata['meta']['Leading razor protein'],prodata['meta']['Lead protein'][prodata['data']['at_least_x_lines']])&\
                                   (np.nansum(pepdata['data'][field_data+'_regbatch'][prodata['data']['valid_lines_independent']]>0,0)>prodata['data']['valid_lines_independent'].sum()*0.33)&\
                                   np.in1d(pepdata['meta']['Sequence'],pepdata['data']['peptide_overlap_snp']['feature_id'][pepdata['data']['peptide_overlap_snp']['overlap_snp']==0])
                                            
'''
write data
'''
#==============================================================================
# # write all
#==============================================================================
metadata_fields=                  ['feature_id','chromosome','start','end','ensembl_gene_id','gene_name','feature_strand','superior_feature_id']
pro_metadata_corresponding_fields=['Lead protein','chromosome',  'start',         'end', 'annotation_ensembl_id','Gene names',       'strand','Gene names']
pep_metadata_corresponding_fields=['Sequence','chromosome',  'start',         'end', 'annotation_ensembl_id','Gene names',       'strand','Leading razor protein']
line_corresponding_fields=['lines','batch','valid_lines_independent']

write_data(folder_dest+file_protein_quantification+"_protein_metadata.txt", \
           mat=np.array([prodata['meta'][key]for key in pro_metadata_corresponding_fields]).T, header= metadata_fields ,delim='\t')
write_data(folder_dest+file_protein_quantification+"_peptide_metadata.txt", \
           mat=np.array([pepdata['meta'][key]for key in pep_metadata_corresponding_fields]).T, header= metadata_fields ,delim='\t')
write_data(folder_dest+file_protein_quantification+"_lines_metadata.txt", \
           mat=np.array([prodata['data'][key] for key in line_corresponding_fields]).T, header= line_corresponding_fields ,delim='\t')


for field_data in ['Reporter intensity corrected','Reporter intensity corrected_ranknorm']:
    write_data(folder_dest+file_protein_quantification+"_protein_"+field_data+'_regbatch'+".txt",\
          mat0=prodata['meta']['Lead protein'][:,None],  mat=prodata['data'][field_data+'_regbatch'].T,\
                      header=np.hstack(["feature_id",prodata['data']['lines']]),delim='\t')

    write_data(folder_dest+file_protein_quantification+"_peptide_"+field_data+'_regbatch'+".txt",\
          mat0=pepdata['meta']['Sequence'][:,None], mat=pepdata['data'][field_data+'_regbatch'].T,header=np.hstack(["feature_id",prodata['data']['lines']]),delim='\t')

#==============================================================================
#==============================================================================
# # write filtered by 
#==============================================================================
#==============================================================================

#==============================================================================

write_data(folder_dest+file_protein_quantification+"_protein_metadata_genes_filtered.txt", \
           mat=np.array([prodata['meta'][key][prodata['data']['at_least_x_lines']] for key in pro_metadata_corresponding_fields]).T, header= metadata_fields ,delim='\t')
write_data(folder_dest+file_protein_quantification+"_peptide_metadata_genes_filtered.txt", \
           mat=np.array([pepdata['meta'][key][pepdata['data']['at_least_x_lines']]for key in pep_metadata_corresponding_fields]).T, header= metadata_fields ,delim='\t')
write_data(folder_dest+file_protein_quantification+"_lines_metadata_lines_filtered_unique.txt", \
           mat=np.array([prodata['data'][key][prodata['data']['valid_lines_independent']] for key in line_corresponding_fields]).T, header= line_corresponding_fields ,delim='\t')


#==============================================================================
# 
#==============================================================================
''' indepednent lines and all genes'''
for field_data in ['Reporter intensity corrected','Reporter intensity corrected_ranknorm']:

    #no processing
    write_data(folder_dest+file_protein_quantification+"_protein_"+field_data+".txt",\
                       mat0=prodata['meta']['Lead protein'][:,None],\
                       mat=prodata['data'][field_data][prodata['data']['valid_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines']]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_peptide_"+field_data+".txt",\
                       mat0=pepdata['meta']['Sequence'][:,None],\
                       mat=pepdata['data'][field_data][prodata['data']['valid_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines']]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_protein_"+field_data+"_lines_filtered_unique.txt",\
                       mat0=prodata['meta']['Lead protein'][:,None],\
                       mat=prodata['data'][field_data][prodata['data']['valid_lines_independent']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines_independent']]]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_peptide_"+field_data+"_lines_filtered_unique.txt",\
                       mat0=pepdata['meta']['Sequence'][:,None],\
                       mat=pepdata['data'][field_data][prodata['data']['valid_lines_independent']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines_independent']]]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_protein_"+field_data+"_genes_filtered.txt",\
                       mat0=prodata['meta']['Lead protein'][prodata['data']['at_least_x_lines']][:,None],\
                       mat=prodata['data'][field_data+'_regbatch'][prodata['data']['valid_lines']][:,prodata['data']['at_least_x_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines']]]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_peptide_"+field_data+"_genes_filtered.txt",\
                       mat0=pepdata['meta']['Sequence'][pepdata['data']['at_least_x_lines']][:,None],\
                       mat=pepdata['data'][field_data+'_regbatch'][prodata['data']['valid_lines']][:,pepdata['data']['at_least_x_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines']]]),delim='\t')
      
    
for field_data in ['Reporter intensity corrected','Reporter intensity corrected_ranknorm']:
    
    #regbatch
    ''' selected genes'''
    write_data(folder_dest+file_protein_quantification+"_protein_"+field_data+'_regbatch'+"_genes_filtered.txt",\
                       mat0=prodata['meta']['Lead protein'][prodata['data']['at_least_x_lines']][:,None],\
                       mat=prodata['data'][field_data+'_regbatch'][prodata['data']['valid_lines']] [:,prodata['data']['at_least_x_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines']]]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_peptide_"+field_data+'_regbatch'+"_genes_filtered.txt",\
                       mat0=pepdata['meta']['Sequence'][pepdata['data']['at_least_x_lines']][:,None],\
                       mat=pepdata['data'][field_data+'_regbatch'][prodata['data']['valid_lines']] [:,pepdata['data']['at_least_x_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines']] ]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_protein_"+field_data+'_regbatch'+"_lines_filtered_unique.txt",\
                       mat0=prodata['meta']['Lead protein'][:,None],\
                       mat=prodata['data'][field_data+'_regbatch'][prodata['data']['valid_lines_independent']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines_independent']]]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_peptide_"+field_data+'_regbatch'+"_lines_filtered_unique.txt",\
                       mat0=pepdata['meta']['Sequence'][:,None],\
                       mat=pepdata['data'][field_data+'_regbatch'][prodata['data']['valid_lines_independent']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines_independent']]]),delim='\t')
    

    ''' indepednent lines and selected genes'''
    write_data(folder_dest+file_protein_quantification+"_protein_"+field_data+'_regbatch'+"_lines_filtered_unique_genes_filtered.txt",\
                       mat0=prodata['meta']['Lead protein'][prodata['data']['at_least_x_lines']][:,None],\
                       mat=prodata['data'][field_data+'_regbatch'][prodata['data']['valid_lines_independent']][:,prodata['data']['at_least_x_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines_independent']]]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_peptide_"+field_data+'_regbatch'+"_lines_filtered_unique_genes_filtered.txt",\
                       mat0=pepdata['meta']['Sequence'][pepdata['data']['at_least_x_lines']][:,None],\
                       mat=pepdata['data'][field_data+'_regbatch'][prodata['data']['valid_lines_independent']][:,pepdata['data']['at_least_x_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines_independent']]]),delim='\t')
    
    
    write_data(folder_dest+file_protein_quantification+"_protein_"+field_data3+"_lines_filtered_unique_genes_filtered.txt",\
                       mat0=prodata['meta']['Lead protein'][prodata['data']['at_least_x_lines']][:,None],\
                       mat=prodata['data'][field_data3][prodata['data']['valid_lines_independent']] [:,prodata['data']['at_least_x_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines_independent']] ]),delim='\t')
    
    write_data(folder_dest+file_protein_quantification+"_peptide_"+field_data3+"_lines_filtered_unique_genes_filtered.txt",\
                       mat0=pepdata['meta']['Sequence'][pepdata['data']['at_least_x_lines']][:,None],\
                       mat=pepdata['data'][field_data3][prodata['data']['valid_lines_independent']] [:,pepdata['data']['at_least_x_lines']].T,\
                       header=np.hstack(["feature_id",prodata['data']['lines'][prodata['data']['valid_lines_independent']] ]),delim='\t')

