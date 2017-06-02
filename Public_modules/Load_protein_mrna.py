import sys
import numpy as np
import csv as csv
import h5py
import pandas
from Utilities import *

def load_protein_mrna(file_protein_quantification='',folder_data='', filter_genes="_genes_filtered",filter_lines="_lines_filtered_unique",field_data="Reporter intensity corrected_regbatch",only_protein=False):

    data={}
    data['protein_intensity']=pandas.read_table(folder_data+file_protein_quantification+"_protein_"+field_data+filter_lines+filter_genes+".txt",sep='\t',index_col=0).transpose()
    data['peptide_intensity']=pandas.read_table(folder_data+file_protein_quantification+"_peptide_"+field_data+filter_lines+filter_genes+".txt",sep='\t',index_col=0).transpose()
    data['peptide_meta']=pandas.read_table(folder_data+file_protein_quantification+"_peptide_metadata"+filter_genes+".txt",sep='\t').set_index('gene_name',drop=False).transpose() 
    data['peptide_protein']=pandas.read_table(folder_data+file_protein_quantification+"_peptide_metadata"+filter_genes+".txt",sep='\t').set_index('gene_name',drop=False).transpose()
    data['protein_meta']=pandas.read_table(folder_data+file_protein_quantification+"_protein_metadata"+filter_genes+".txt",sep='\t').set_index('gene_name',drop=False).transpose()
    data['line_meta']=pandas.read_table(folder_data+file_protein_quantification+"_lines_metadata"+filter_lines+".txt",sep='\t').set_index('lines',drop=False)
    data['batch_mat']=pandas.DataFrame(data=np.vstack([data['line_meta']['batch']==tmt for tmt in np.unique(data['line_meta']['batch'])]).astype(float).T,\
        columns=np.unique(data['line_meta']['batch']),  index=data['line_meta']['lines']) 
    if only_protein:
        return data

    rna_marc = h5py.File('/Users/mirauta/Data/RNA/hipsci/rna_seq_counts_salmon_genes_2304.h5','r')
    lines=np.unique(rna_marc['meta/lines'][:].astype('U'))
    
    rna_lines=rna_marc['meta/lines'][:].astype('U')
    index=np.unique(rna_lines,return_index=1)[1]
    temp_rna_counts=rna_marc['data/counts'][:][index]
    rna_lines=rna_lines[index]
    rna_all_genes=np.array([temp[1] for temp   in rna_marc['meta/genes'][:]]).astype('U')
    temp=np.unique(rna_all_genes,return_counts=1);temp
     
    common_lines=np.intersect1d(data['protein_intensity'].index,rna_lines)
    common_genes=np.intersect1d(data['peptide_meta'].transpose()['gene_name'],np.intersect1d(data['protein_meta'].transpose()['gene_name'],rna_all_genes)); print (common_genes.shape)
    
    temp_rna_counts=np.array([temp_rna_counts[rna_lines ==ll][0] for ll in common_lines]).T; print (temp_rna_counts.shape) 
    temp_rna_counts=np.array([d/d.sum()*temp_rna_counts[:,0].sum() for d in temp_rna_counts.T]).T; print (temp_rna_counts.shape) 
    #data['rna_counts']=np.array([temp_rna_counts[ rna_all_genes ==g].sum(0) for g in  common_genes])
    data['rna_counts']=pandas.DataFrame(data=np.array([temp_rna_counts[ rna_all_genes ==g].sum(0) for g in  common_genes]).T,\
        index=common_lines,  columns=common_genes)
     
    ## select common lines                       
    for x in ['protein_intensity','peptide_intensity', 'line_meta','batch_mat' ]:
        temp=data[x].transpose()
        data[x]=pandas.DataFrame(data=np.array([temp[ll] for ll in common_lines]), index=common_lines,  columns=data[x].columns.values)                   
    
    for x in ['protein_meta','peptide_meta', 'peptide_protein']:
        data[x]=data[x].transpose()[np.in1d(data[x].transpose().index,common_genes)].transpose()
    
    return data

#sys.exit()
#corr={}
#for key in ['regbatch','reg0']:    corr[key]=np.zeros(len(common_genes))+np.nan
#for ig,g in enumerate(common_genes[:5000]):
#    try :
#        yp=data['protein_intensity'][data['protein_meta'][g]['feature_id']]
#        yr=data['rna_counts'][g]
##        yp1=limix.util.preprocess.regressOut(Y= yp,X=data['batch_mat'])
##        corr['regbatch'][ig]=scst.spearmanr(yp1[yp==yp],yr[yp==yp])[0 ]
#        corr['reg0'][ig]=scst.spearmanr(yp[yp==yp],yr[yp==yp])[0 ]
#    except: 1 
#    
#plt.hist([corr[key][np.isfinite(corr[key])]for key in corr.keys()],bins=50,label=corr.keys())    
#plt.legend(loc=2)
#
##get only genes with same number of transcirpts
#temp=np.unique(rna_all_genes,return_counts=1);
##listgenes=temp[0][np.where((temp[1]<8)&(temp[1]>5))[0]]    
#data['rna_counts']=np.array([temp_rna_counts[ rna_all_genes ==g].sum(0) for g in  common_genes])
#
#data['rna_counts']=np.array([temp_rna_counts[ rna_all_genes ==g].sum(0) for g in  common_genes])
#
#
#
#
#indexg=np.where(((data['maxquant_scaled']>0).sum(1)==data['maxquant_scaled'].shape[1])&\
#                (np.array([(pep.shape[0]>5)&(pep.shape[0]<100) for pep in data['peptide']]))&\
#                (np.array([(pep>0).mean()>0.8 for pep in data['peptide_intensity']])))[0]
#
#print (indexg.shape)
#data1={}
#for key in data.keys():
#    data1[key]=data[key][indexg]
    
#
#selgenes=np.unique(data['peptide_protein'].columns.values,return_counts=1)[0][np.unique(data['peptide_protein'].columns.values,return_counts=1)[1]>5]
#print (selgenes.shape)
#
#
#data1={}
#data1['peptide_intensity']=[data['peptide_intensity'][data['peptide_meta'][g].transpose()['feature_id']].T for g in selgenes]
#data1['protein_intensity']= [data['protein_intensity'][data['protein_meta'][g].transpose()['feature_id']].T for g in selgenes] 
#data1['rna_counts']= [data['rna_counts'][g].T for g in selgenes] 
#data1['ensembl_id']=selgenes
#
#ig=0
#y=(data1['peptide_intensity'][ig]);print (y.shape)
#y=y[np.argsort(np.nansum((y>0),1))]
#plt.plot(y.values);plt.show()
#
#data['protein_meta'][g]['feature_id']
def get_peptide_scalar_average (y,Iter=None,beta_model=True,return_all=False,given_a=None,given_u=None,given_tauil=None,given_betai=None):
 
    
    I=y.shape[0];    L=y.shape[1]

    a=np.nanmedian(y ,1)/(1+np.max(np.nanmedian(y ,1)))
    u= np.array([np.nanmean(y[:,l]/a) for  l in np.arange(y.shape[1])])
    likelihood= np.nansum(scst.norm.logpdf(y,np.dot(a[:, None],u[None,:]), scale=1))

    return [a,u,np.ones([I,L]),np.ones([I]),likelihood, abs(y-np.dot(a[:,None],u[None,:]))/y]

#
#
#
#for key in ['all','pred','pred_err','pred_err2']:
#    data1[key]=np.zeros(len(data1['ensembl_id']),dtype='object')
# 
#for ig,g in enumerate(data1['ensembl_ID']):
#    print (ig)
#    y=(data1['peptide_intensity'][ig].values)
##    y[data1['missing'][ig]]=np.nan
#    y=y[np.argsort(np.nansum((y>0),1))[::-1][:20]]
#    y=y/np.nanmean(y,1)[:,None]
#    
#    data1['all'][ig]=get_peptide_scalar_average(y,Iter=2, beta_model=True)
#    a,u,tau,beta,like1, err=data1['all'][ig]
#
#    data1['pred'][ig]=np.log(u)
#    data1['pred_err'][ig]=1/np.sqrt(np.nansum(tau,0))
#    data1['pred_err2'][ig]=np.nansum(err,0)
#    
#
#
#data1['deviation_protein']={}
#data1['deviation_protein_batch']={}
##data1['deviation_protein_batch2']={}
#data1['pred_err2_batch']={}
#data1['sum_pep']={}
#for ig,g in enumerate(data1['ensembl_id']):
#    data1['sum_pep'][ig]=np.nansum(data1['peptide_intensity'][ig],0)
#    data1['deviation_protein'][ig]=(data1['pred'][ig]-data1['pred'][ig].mean())/data1['pred'][ig]
#
#    data1['deviation_protein_batch'][ig]=limix.util.preprocess.regressOut(Y= data1['deviation_protein'][ig],X=data['batch_mat'])
#    data1['pred_err2_batch'][ig]=limix.util.preprocess.regressOut(Y= data1['pred_err'][ig],X=data['batch_mat'])
#
#
#keys=['RTP_batch','RTP','RTP_sum','cor','cor_sum']
#for key in keys: data1[key]={}
#    
#for ig,g in enumerate(data1['ensembl_ID']):
#    r=np.log(1+data1['rna_counts'][ig]);r=(r-np.nanmean(r))/np.nanstd(r)
#    y=data1['pred'][ig];y=(y-np.nanmean(y))/np.nanstd(y)
#    data1['RTP_batch'][ig]=y-r; data1['cor'][ig]=scst.spearmanr(y,r)[0]
#    y=data1['sum_pep'][ig];y=(y-np.nanmean(y))/np.nanstd(y)
#    data1['RTP_sum'][ig]=y-r; data1['cor_sum'][ig]=scst.spearmanr(y,r)[0]
#
#
#rez={}
#rez['sd_RTP']=np.array([scst.spearmanr(abs(data1['RTP_batch'][ig]-data1['RTP_batch'][ig].mean()), data1['pred_err2'][ig])[0 ] for ig,g in enumerate(data1['ensembl_ID'])])
#rez['diver_RTP']=np.array([scst.spearmanr(abs(data1['RTP_batch'][ig]-data1['RTP_batch'][ig].mean()), abs(data1['deviation_protein_batch'][ig]))[0 ] for ig,g in enumerate(data1['ensembl_ID'])])
#rez['diver_sd']=np.array([scst.spearmanr( data1['pred_err2_batch'][ig], abs(data1['deviation_protein_batch'][ig]))[0 ] for ig,g in enumerate(data1['ensembl_ID'])])
#
#
##rez['pro.cv_prot_Spearman']=np.array([scst.spearmanr(abs(data1['pred2'][ig]),                                             data1['pred2_err'][ig]*data1['pred2'][ig]**2)[0 ] for ig,g in enumerate(data1['ensembl_ID'])])
#print (np.nanpercentile(rez['sd_RTP'],(25,50,75)))
#print (np.nanpercentile(rez['diver_RTP'],(25,50,75)))
#print (np.nanpercentile(rez['diver_sd'],(25,50,75)))
#plt.figure(figsize=(12,8))
#plt.hist([rez[key][rez[key]==rez[key]] for key in rez.keys()],label=rez.keys())
#plt.legend(loc=1)
#plt.show()
