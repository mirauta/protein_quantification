import sys
import numpy as np
#import csv as csv
#import h5py
sys.path.append('../')
from Public_modules.Utilities import *
from Public_modules.Utilities import *
import pandas
import matplotlib.pyplot as plt


folder_data='/Users/mirauta/data/MS/hipsci/phenotypes/'
file_protein_quantification='hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517'
field_data='Reporter intensity corrected_regbatch'
filter_lines='_lines_filtered_unique'
filter_genes='_genes_filtered'


data={}

data['protein_intensity']=pandas.read_table(folder_data+file_protein_quantification+"_protein_"+field_data+filter_genes+".txt",sep='\t',index_col=0).transpose()
data['peptide_intensity']=pandas.read_table(folder_data+file_protein_quantification+"_peptide_"+field_data+filter_genes+".txt",sep='\t',index_col=0).transpose()
data['peptide_meta']=pandas.read_table(folder_data+file_protein_quantification+"_peptide_metadata"+filter_genes+".txt",sep='\t',index_col=0).transpose() 
data['peptide_protein']=pandas.read_table(folder_data+file_protein_quantification+"_peptide_metadata"+filter_genes+".txt",sep='\t').set_index('superior_feature_id').transpose()
data['protein_meta']=pandas.read_table(folder_data+file_protein_quantification+"_protein_metadata"+filter_genes+".txt",sep='\t',index_col=0).transpose()
data['line_meta']=pandas.read_table(folder_data+file_protein_quantification+"_lines_metadata.txt",sep='\t')
data['line_meta']['linesi']=np.copy(data['line_meta']['lines'])#np.array([l.split('.')[0]for l in data['line_meta']['lines']])
data['line_meta']['lines']=data['peptide_intensity'].index

data['line_meta']['donor']=np.array([l.split('-')[-1].split('_')[0] for l in data['protein_intensity'].index])
#data['line_meta']['valid_lines']=np.array(['HPSI' in l for l in data['line_meta']['lines']])& np.array(['Calcium' not in l for l in data['line_meta']['batch']])
#
#temp=(np.nansum(data['peptide_intensity']>0,1))
#data['line_meta']['valid_lines']=(temp<=np.nanmedian(temp[np.where(data['line_meta']['valid_lines'])[0]])*1.25)&\
#         (temp>=np.nanmedian(temp[np.where(data['line_meta']['valid_lines'])[0]])*0.75) 
#plt.plot(np.nansum(data['peptide_intensity']>0,1))
#plt.plot(data['line_meta']['valid_lines'])
#
#data['peptide_intensity']=pandas.read_table(folder_data+\
#    file_protein_quantification+"_peptide_"+'Reporter intensity corrected_ranknorm_regbatch'+filter_genes+".txt",sep='\t',index_col=0).transpose()
#data['peptide_intensity2']=pandas.read_table(folder_data+\
#    file_protein_quantification+"_peptide_"+'Reporter intensity corrected_regbatch'+filter_genes+".txt",sep='\t',index_col=0).transpose()
#temp=data['peptide_intensity'].values[data['line_meta']['valid_lines'].values]
#temp=data['peptide_intensity2'].values[data['line_meta']['valid_lines'].values]
#plt.hist([np.log(temp[i][temp[i]>0])for i in [20,25,35,45,47,48,58]],bins=10,normed=1)
#plt.xlim(4,25)
#
#sys.exit()
#==============================================================================
# select peptides for a protein
#==============================================================================

#data['peptide_protein_intensity']=[data['peptide_intensity'][data['peptide_protein'][g].transpose()['feature_id']].T for g in np.unique(data['peptide_protein'].columns.values)]
#data['peptide_protein']=[data['peptide_intensity'][data['peptide_protein'][g].transpose()['feature_id']].T for g in np.unique(data['peptide_protein'].columns.values)]
#data['protein_intensity_scaled']=np.zeros([data['line_meta'].shape[0],len(np.unique(data['peptide_protein'].columns.values))])    +np.nan
#data['peptide_protein_intensity_scaled']=np.zeros(len(data['peptide_protein_intensity']),dtype='object')
data['peptide_scalar']=pandas.DataFrame(data=np.zeros(data['peptide_intensity'].shape[1]),index=data['peptide_intensity'].columns.values,columns=['scalar']).transpose()
data['peptide_intensity_scalarcorrected']= pandas.DataFrame(data= data['peptide_intensity'].values,index=data['peptide_intensity'].index,columns=data['peptide_intensity'].columns.values)
data['protein_intensity_scalarcorrected']= pandas.DataFrame(data= data['protein_intensity'].values,index=data['protein_intensity'].index,columns=data['protein_intensity'].columns.values)

### 

valid_replicated_lines=data['line_meta'][(data['line_meta']['donor']=='bubh')]['lines']
for ig,g in enumerate(np.unique(data['peptide_protein'].columns.values)):
#    peptidesinprotein=data['peptide_intensity'][data['peptide_protein'][g].transpose()['feature_id']].T
    peptides_in_protein= data['peptide_protein'][g].transpose()['feature_id']
    peptides_intensity=data['peptide_intensity'][peptides_in_protein]

    temp=peptides_intensity.transpose()[valid_replicated_lines]
    try:
        scalar=np.nanmedian(temp/np.nanmedian(temp),1)
        data['peptide_intensity_scalarcorrected'] [peptides_in_protein]= peptides_intensity/scalar
        data['peptide_scalar'][peptides_in_protein]= scalar
        data['protein_intensity_scalarcorrected'][g]=    np.nanmedian(peptides_intensity/scalar,1)
    except:
        data['peptide_scalar'][peptides_in_protein]= np.nan
        
        
        
#for ig,g in enumerate(np.unique(data['peptide_protein'].columns.values)[:1]): 
##plt.plot(np.log(data['peptide_protein_intensity'][ig].values.T))
##plt.plot(np.log(data['peptide_protein_intensity_scaled'][ig].T))
#
#    plt.plot(np.log(data['protein_intensity'][g].values.T))
#    plt.plot(np.log(data['protein_intensity_scaled'][:,ig].T),'r')



data['protein_intensity_scalarcorrected_qnorm']=np.array([transform_vector_qnorm(y, method='qnorm') for y in     data['protein_intensity_scalarcorrected'].values.T]).T
    

write_data(folder_data+file_protein_quantification+"_protein_intensity_scalarcorrected_"+field_data+filter_lines+filter_genes+".txt",\
                   mat0=data['protein_intensity_scalarcorrected'].columns.values[:,None],\
                   mat=data['protein_intensity_scalarcorrected'][ data['line_meta']['valid_lines_independent'].values].T.values,\
                   header=np.hstack(["feature_id", data['protein_intensity_scalarcorrected'].index[ data['line_meta']['valid_lines_independent'].values] ]),delim='\t')

write_data(folder_data+file_protein_quantification+"_protein_intensity_scalarcorrected_"+field_data+filter_lines+filter_genes+"_qnorm.txt",\
                   mat0=data['protein_intensity_scalarcorrected'].columns.values[:,None],\
                   mat=data['protein_intensity_scalarcorrected_qnorm'][ data['line_meta']['valid_lines_independent'].values].T,\
                   header=np.hstack(["feature_id", data['protein_intensity_scalarcorrected'].index[ data['line_meta']['valid_lines_independent'].values] ]),delim='\t')

write_data(folder_data+file_protein_quantification+"_peptide_intensity_scalarcorrected_"+field_data+filter_lines+filter_genes+".txt",\
                   mat0=data['peptide_intensity_scalarcorrected'].columns.values[:,None],\
                   mat=data['peptide_intensity_scalarcorrected'][ data['line_meta']['valid_lines_independent'].values].T.values,\
                   header=np.hstack(["feature_id", data['peptide_intensity_scalarcorrected'].index[ data['line_meta']['valid_lines_independent'].values] ]),delim='\t')


write_data(folder_data+file_protein_quantification+"_protein_intensity_scalarcorrected_"+field_data+filter_genes+".txt",\
                   mat0=data['protein_intensity_scalarcorrected'].columns.values[:,None],\
                   mat=data['protein_intensity_scalarcorrected'].T.values,\
                   header=np.hstack(["feature_id", data['protein_intensity_scalarcorrected'].index  ]),delim='\t')

write_data(folder_data+file_protein_quantification+"_protein_intensity_scalarcorrected_"+field_data+filter_genes+"_qnorm.txt",\
                   mat0=data['protein_intensity_scalarcorrected'].columns.values[:,None],\
                   mat=data['protein_intensity_scalarcorrected_qnorm'].T,\
                   header=np.hstack(["feature_id", data['protein_intensity_scalarcorrected'].index ]),delim='\t')

write_data(folder_data+file_protein_quantification+"_peptide_intensity_scalarcorrected_"+field_data+filter_genes+".txt",\
                   mat0=data['peptide_intensity_scalarcorrected'].columns.values[:,None],\
                   mat=data['peptide_intensity_scalarcorrected'].T.values,\
                   header=np.hstack(["feature_id", data['peptide_intensity_scalarcorrected'].index ]),delim='\t')

