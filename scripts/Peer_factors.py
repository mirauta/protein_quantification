
import sys
sys.path.append('../Public_modules/')
from Utilities import *
from Load_protein_mrna import *
import scipy.cluster.hierarchy as hier
import peer
from limix.utils import *
######################################################
## load folder_name, batch names
######################################################

folder_data='/nfs/research2/hipsci/processed_data/proteomics/'
file_protein_quantification='hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517'
filter_genes=''
filter_lines="_lines_filtered_unique"

data=load_protein_mrna( filter_genes=filter_genes,filter_lines=filter_lines,file_protein_quantification=file_protein_quantification,\
        folder_data=folder_data,field_data="Reporter intensity corrected_regbatch",only_protein=True)

## Run PEER ### Input matrix: N rows and G columns, where N is the number of samples, and G is the number of genes.

#Cov=np.vstack([np.array([f5p['metadata/lines/'][covariate][:] ==c for c in np.setdiff1d(f5p['metadata/lines/'][covariate][:],'')]) for covariate in ['gender','batch_2_ID']]).astype(float).T
model = peer.PEER()
#model.setCovariates(Cov) # N x C matrix

model.setPhenoMean(data['protein_intensity'].values[:,np.isfinite(np.sum(data['protein_intensity'].values,0))][:,:1000])
model.setAdd_mean(True)
# To infer K hidden confounders, define K and number of iterations
K=10;model.setNk(int(K)) # or PEER_setNk(model,number_of_covs)
model.setNmax_iterations(1000)
model.update()
factors = model.getX(); weights = model.getW(); precision = model.getAlpha(); expr_peer = model.getResiduals()

write_data( folder_data+file_protein_quantification+"_lines_metadata2"+filter_lines+".txt",        mat=np.hstack([data['line_meta'].values,factors]),\
        header= np.hstack([data['line_meta'].columns.values,['peer_'+str(i) for i in np.arange(factors.shape[1])]]),delim='\t')

sys.exit()
import matplotlib.pyplot as plt
plt.figure(figsize=(12,8))
plt.subplot(2,1,1)
plt.plot(factors1);
plt.subplot(2,1,2)
plt.plot(factors2);
plt.show()
