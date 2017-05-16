
import sys
sys.path.append('/nfs/research2/stegle/users/mirauta/scripts/')
sys.path.append('/nfs/research2/stegle/users/mirauta/scripts/Public_modules/')
from Libraries import *
from Utilities import *
from Utilities_proteomics import *
import scipy.cluster.hierarchy as hier

######################################################
## load folder_name, batch names
######################################################

parameters_file=sys.argv[1]
norm_method=sys.argv[2]
with open(parameters_file, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='=', quotechar='"')
        parameters = np.array([row for row in reader])

try:    folder_dest=parameters[np.array(['folder_dest' in p for p in parameters[:,0]]),1][0]
except: print 'no folder dest'
try:    folder_data_lf=parameters[np.array(['folder_data_lf' in p for p in parameters[:,0]]),1][0]
except: print 'no folder data lf'
try:    folder_data_tmt=parameters[np.array(['folder_data_tmt' in p for p in parameters[:,0]]),1][0]
except: print 'no folder data tmt'
try:	batch_names_lf=np.array([parameters[i,1] for i,p in enumerate(parameters[:,0]) if 'batch_name_lf' in p])
except: print 'no lf data'
try:	batch_names_tmt=np.array([parameters[i,1] for i,p in enumerate(parameters[:,0]) if 'batch_name_tmt' in p])
except: print 'no tmt data'
try:    file_protein_quantification=parameters[np.array(['file_protein_quantification' in p for p in parameters[:,0]]),1][0]
except: print 'no tmt data'
try:    folder_sample_metadata=parameters[np.array(['folder_sample_metadata' in p for p in parameters[:,0]]),1][0]
except: print 'geno_name'
try:    file_sample_metadata=parameters[np.array(['file_sample_metadata' in p for p in parameters[:,0]]),1][0]
except: print 'geno_name'
try:   ensembl_file=parameters[np.array(['ensembl_annotation_file' in p for p in parameters[:,0]]),1][0]
except: print 'geno_name'
try:   gencode_file=parameters[np.array(['gencode_annotation_file' in p for p in parameters[:,0]]),1][0]
except: print 'geno_name'

print batch_names_lf; print batch_names_tmt
#short_names=[b.split('/')[0].replace('_uniprot','') for b in np.hstack([batch_names_lf,batch_names_tmt])]
batches=np.hstack([batch_names_lf,batch_names_tmt])
batch_lf_tmt={}
batch_lf_tmt['batch']=np.hstack([batch_names_lf,batch_names_tmt])
batch_lf_tmt['batch_short']=np.array([b.split('/')[0].replace('_uniprot','') for b in  batch_lf_tmt['batch']])
batch_lf_tmt['type']=np.hstack([np.repeat('lf',len(batch_names_lf)),np.repeat('tmt',len(batch_names_tmt))])
batch_lf_tmt['intensity']=np.hstack([np.repeat('Intensity ',len(batch_names_lf)),np.repeat('Reporter intensity corrected ',len(batch_names_tmt))])
batch_lf_tmt['identification']=np.hstack([np.repeat('Identification type ',len(batch_names_lf)),np.repeat('Identification type ',len(batch_names_tmt))])
################################


pep_name=np.array(['Leading razor protein','Sequence','Proteins', 'Unique (Groups)','PEP','Score','Reverse', 'Potential contaminant'])
pro_name=np.array(['Majority protein IDs','Peptide IDs','Peptide counts (razor+unique)','Gene names', 'Peptide is razor','Mol. weight [kDa]','Only identified by site'])
print "reading data"

field_intensity_pro='Intensity '

#sys.exit()
allpro=[load_data_from_maxquant_output(folder_data=folder_data_lf,name=batch_lf_tmt['batch'][ib],data_file='.proteinGroups.txt', field_data=batch_lf_tmt['intensity'][ib], field_metadata=pro_name) for ib,b in enumerate(batch_lf_tmt['batch'])]
#sys.exit()

print "read protein data. start reading peptide data"
allpep=[load_data_from_maxquant_output(folder_data=folder_data_lf,name=batch_lf_tmt['batch'][ib],data_file='.peptides.txt',          field_data=batch_lf_tmt['intensity'][ib], field_metadata=pep_name) for ib,b in enumerate(batch_lf_tmt['batch'])]

### allpep[ib][1]['data'] is a list [fields][peptides]
for ib,b in enumerate(batch_lf_tmt['batch']):
    temp=allpep[ib][2]['match']=load_data_from_maxquant_output(folder_data=folder_data_lf,name=batch_lf_tmt['batch'][ib],float_data=0,data_file='.peptides.txt',  field_data=batch_lf_tmt['identification'][ib], field_metadata=pep_name)[2]
    allpep[ib][2]['match']=np.zeros(allpep[ib][2]['data'].shape)
    if batch_lf_tmt['type'][ib]=='lf':
        allpep[ib][2]['match']=np.vstack(temp['data'])=='By MS/MS'
    else:
        tempfieldlines=np.array([' '.join(f.replace(batch_lf_tmt['intensity'][ib],'').split(' ')[1:]) for f in allpep[ib][2]['fields']])
        for itt, tt in enumerate(temp['fields']):
            allpep[ib][2]['match'][np.where(tempfieldlines==tt.replace(batch_lf_tmt['identification'][ib],''))[0]]=temp['data'][itt]=='By MS/MS'

print "data read finished"
# 1) get all data + modify the Mjority protein field
# 2) get common proteins
# 3) get protein abundance file 4) get peptide per protein
# 5) get protein x peptide file  and peptide x lines file
# 6) get metadatafile



################
# check if different fields

fields=[list(allpro[i][1]['fields'])  for i, b in enumerate(batches)]
if np.array([np.in1d(fields[0],fields[i]).sum()==len(fields[0])for i, b in enumerate(batches)]).sum()!=len(batches):
	print 'different fields in batches'
	sys.exit()






### 0) create joindata container
### 0.1) bring lines
#>>>>>>>>>>>>>
joindata={'fields':{'fields': allpro[0][1]['fields']}}
joindata['lines']={}
joindata['lines']['batch_ID']= np.hstack([np.repeat(b,len(allpro[i][2]['fields'])) for i, b in       enumerate(batch_lf_tmt['batch_short'])])
joindata['lines']['original_line_ID']=np.hstack([[a.replace(batch_lf_tmt['intensity'][ib],'') for a in       allpro[ib][2]['fields']] for ib, b in enumerate(batches)])
joindata['lines']['line_ID']=np.copy( joindata['lines']['original_line_ID'])
for il, l in enumerate( joindata['lines']['line_ID']):
    temp=np.array(l.split(' '))[np.array(['HPSI' in aa for aa in l.split(' ')])]
    if len(temp)>0:
        joindata['lines']['line_ID'][il]=temp[0]

joindata['lines']['batch_2_ID']=np.hstack([l.replace( joindata['lines']['line_ID'][il],'') for il,l    in enumerate( joindata['lines']['original_line_ID'])])
joindata['lines']['batch_3_ID']=np.copy( joindata['lines']['line_ID'])
for il, l in enumerate( joindata['lines']['original_line_ID']):
    temp=l.replace( joindata['lines']['line_ID'][il],'')
    if len(temp)>0:
        joindata['lines']['batch_3_ID'][il]=' '.join(temp.split(' ')[1:])

joindata['lines']['protocol_ms']=np.repeat('tmt',len( joindata['lines']['batch_ID']))
if len(batch_names_lf)>0:  joindata['lines']['protocol_ms'][:sum([len(allpro[i][2]['fields']) for i in range(len(batch_names_lf))])]='lf'

####################### 1)
field_pro=np.squeeze(np.where(allpro[0][1]['fields']=='Majority protein IDs')[0])
field_pep=np.squeeze(np.where(allpep[0][1]['fields']=='Leading razor protein')[0])
for i, b in enumerate(batches):
	allpro[i][1]['protein']=np.array([a.split(';')[0] for a in allpro[i][1]['data'][field_pro]])
        allpep[i][1]['protein']=np.array([a.split(';')[0] for a in allpep[i][1]['data'][field_pep]])



####################### 2)
# 2) get common proteins


compr=np.array(list(set(allpro[0][1]['protein']).intersection(*[allpro[i][1]['protein']for i, b in enumerate(batches)])))
comindex=np.array([np.where(np.in1d( allpro[i][1]['protein'],compr))[0][np.argsort( allpro[i][1]['protein'][np.in1d( allpro[i][1]['protein'],compr)])] for i, b in enumerate(batches)])

 #check if all proteins are in the same order
proteinstest=np.vstack([allpro[ib][1]['protein'][comindex[ib]] for i, b in enumerate(batches)])
if np.array([np.in1d(proteinstest[0],proteinstest[i]).sum()==len(proteinstest[0])for i, b in enumerate(batches)]).sum()!=len(batches):
    print 'different proteins in batches'
    print np.array([np.in1d(proteinstest[0],proteins[i]).sum()==len(proteinstest[0])for i, b in enumerate(batches)])
    sys.exit()



#>>>>>>>>>>>>>
joindata['data']={'protein_intensity_maxquant': np.hstack([allpro[i][2]['data'][:,comindex[i]].T for i, b in enumerate(batches)]).T}
joindata['data']['protein'] = allpro[0][1]['protein'][comindex[0]]
joindata  ['data']['protein_info'] =[[f[comindex[i]]for f in allpro[i][1]['data']] for i, b in enumerate(batches)]
joindata['fields']['protein_info'] = [allpro[i][1]['fields'] for i, b in enumerate(batches)]
joindata['lines']['sum_pro']=np.hstack([(allpro[i][2]['data']).sum(1) for i, b in enumerate(batches)]).T
joindata['lines']['count_pro']=np.hstack([(allpro[i][2]['data']>0).sum(1) for i, b in enumerate(batches)]).T
joindata['lines']['count_pep']=np.hstack([(allpep[i][2]['data']>0).sum(1) for i, b in enumerate(batches)]).T
joindata['lines']['count_pep_direct']=np.hstack([np.array(allpep[i][2]['match']).sum(1) for i, b in enumerate(batches)]).T
 #


#sys.exit()
# 3) for each protein recompute the abundance
# make the choice to normalise peptides
# then obtain peptides
# make the filtering choices:
# include only razor
# direct
# filter on,
# filter out PEP,  'Reverse', 'Potential contaminant'


# make the normalisation choice
#sys.exit()
referencebatch=np.argmax(np.array([np.max((allpep[i][2]['data']>0).sum(1)) for i, b in enumerate(batches)]))
referenceline=np.argmax((allpep[referencebatch][2]['data']>0).sum(1))
for ib, b in enumerate(batches):
	allpep[ib][2]['data_norm']=normalize_reference(data=allpep[ib][2]['data'].T, ntype=norm_method,referencedata=allpep[referencebatch][2]['data'][referenceline]).T


''' here correct intensy: divide it by the first  '''
# for each protein obtain peptides
#
print "computind protein_data quantifications"
#sys.exit()
protein_data={}
testing=False
field_peppro=np.squeeze(np.where(allpro[0][1]['fields']=='Peptide IDs')[0])
for ipr,pr in enumerate(joindata['data']['protein']):
	#print joindata['data']['protein_info'][ib][0][ipr]
        indices=[np.array(joindata['data']['protein_info'][ib][field_peppro][ipr].split(';')).astype(int) for ib, b in enumerate(batches)]
        #
        peppro={'indices':indices}
        peppro['peptide_filter']=[np.array(allpep[ib][1]['data'][np.squeeze(np.where(allpep[ib][1]['fields']=='Reverse')[0])][indices[ib]]=='') for ib, b in enumerate(batches)]
        peppro['peptide_filter']=[peppro['peptide_filter'][ib]*np.array(allpep[ib][1]['data'][np.squeeze(np.where(allpep[ib][1]['fields']=='Potential contaminant')[0])][indices[ib]]=='') for ib, b in enumerate(batches)]
        peppro['peptide_razor']=[np.array(joindata['data']['protein_info'][ib][np.squeeze(np.where(joindata['fields']['protein_info'][ib]=='Peptide is razor')[0])][ipr].split(';'))=='True' for ib, b in enumerate(batches)]
        #
        peppro['peptide']=[[p[indices[ib]]for p in allpep[ib][1]['data']] for ib, b in enumerate(batches)]
        peppro['peptide_intensity']=[allpep[ib][2]['data'][:,indices[ib]] for ib, b in enumerate(batches)]
        peppro['peptide_intensity_norm']=[allpep[ib][2]['data_norm'][:,indices[ib]] for ib, b in enumerate(batches)]
        peppro['peptide_direct']=[np.array([p[indices[ib]] for p in allpep[ib][2]['match']])  for ib, b in enumerate(batches)]
        #
        peppro['peptide_direct_intensity']=[peppro['peptide_intensity'][ib]*peppro['peptide_direct'][ib]/peppro['peptide_direct'][ib] for ib, b in enumerate(batches)]
        peppro['peptide_direct_intensity_norm']=[peppro['peptide_intensity_norm'][ib]*peppro['peptide_direct'][ib]/peppro['peptide_direct'][ib] for ib, b in enumerate(batches)]





        protein_data[pr]={'intensity':{'maxquant':joindata['data']['protein_intensity_maxquant'][:,ipr]}}
        protein_data[pr]['intensity']['maxquant_scaled']=protein_data[pr]['intensity']['maxquant']*joindata['lines']['sum_pro'].min()/joindata['lines']['sum_pro']

        protein_data[pr]['intensity']['direct']=np.hstack([peppro['peptide_intensity'][ib].sum(1) for ib, b in enumerate(batches)] )
        #inlude only razor proteins
        protein_data[pr]['intensity']['razor']=np.hstack([peppro['peptide_intensity'][ib][:,peppro['peptide_razor'][ib]].sum(1) for ib, b in enumerate(batches)] )
        protein_data[pr]['intensity']['razor_pepnorm']=np.hstack([peppro['peptide_intensity_norm'][ib][:,peppro['peptide_razor'][ib]].sum(1) for ib, b in enumerate(batches)] )

        protein_data[pr]['intensity']['razor_direct']=np.hstack([np.nansum(peppro['peptide_direct_intensity'][ib][:,peppro['peptide_razor'][ib]],1) for ib, b in enumerate(batches)] )
        protein_data[pr]['intensity']['razor_filter_direct']=np.hstack([np.nansum(peppro['peptide_direct_intensity'][ib][:,peppro['peptide_razor'][ib]*peppro['peptide_filter'][ib]],1) for ib, b in enumerate(batches)] )
        protein_data[pr]['intensity']['razor_filter']=np.hstack([np.nansum(peppro['peptide_intensity'][ib][:,peppro['peptide_razor'][ib]*peppro['peptide_filter'][ib]],1) for ib, b in enumerate(batches)] )
        protein_data[pr]['intensity']['peptide_count_razor_filter']=np.hstack([np.nansum(peppro['peptide_intensity'][ib][:,peppro['peptide_razor'][ib]*peppro['peptide_filter'][ib]]>0,1) for ib, b in   enumerate(batches)] )
        protein_data[pr]['intensity']['razor_pepnorm_filter_direct']=np.hstack([np.nansum(peppro['peptide_direct_intensity'][ib][:,peppro['peptide_razor'][ib]*peppro['peptide_filter'][ib]],1) for ib, b in enumerate(batches)]  )
        protein_data[pr]['intensity']['razor_pepnorm_filter']=np.hstack([np.nansum(peppro['peptide_intensity_norm'][ib][:,peppro['peptide_razor'][ib]*peppro['peptide_filter'][ib]],1) for ib, b in enumerate(batches)] )


        uniquepep=np.unique(np.hstack([b[0][peppro['peptide_razor'][ib]*peppro['peptide_filter'][ib]*( peppro['peptide_intensity'][ib].sum(0)>0)] for ib, b in  enumerate(peppro['peptide'])]))

        temp=np.zeros([joindata ['lines']['batch_ID'].shape[0],uniquepep.shape[0]])+np.nan
        for ib, b in  enumerate(batch_lf_tmt['batch_short']):
            for ipep, pep in enumerate(uniquepep):
                temp2=np.squeeze(peppro['peptide_intensity_norm'][ib][:,peppro['peptide'][ib][0]==pep])
                if np.isscalar(temp2[0]): temp[joindata ['lines']['batch_ID']==b,ipep]=temp2
        protein_data[pr]['peptide_intensity_norm']=temp
        protein_data[pr]['peptide']=uniquepep
        temp0=temp[:,(temp>=0).sum(0)==temp.shape[0]]
        protein_data[pr]['intensity']['razor_pepnorm_filter_allbatch']=np.nansum(temp0,1)
        protein_data[pr]['intensity']['razor_pepnorm_filter_allbatch_corpep']=np.nansum(temp0,1)
        if temp0.shape[1]>4:
            Z=hier.linkage(scst.spearmanr(temp0)[0], 'ward');
            cl=np.squeeze(hier.cut_tree(Z,height=3))
            cl0=np.unique(cl,return_counts=1)
            if len(cl0[0])>1:
#                sys.exit()
                protein_data[pr]['intensity']['razor_pepnorm_filter_allbatch_corpep']=np.nansum(temp0[:,cl==cl0[0][np.argmax(cl0[1])]],1)


        if testing:
            dif=(abs((protein_data[pr]['maxquant']-protein_data[pr]['razor'])/protein_data[pr]['maxquant']))
            if np.nanmax(dif) > 0.001:
                print dif
                temp=np.where(dif>0.001)[0]
                print pr
                print protein_data[pr]['maxquant'][temp]
                print protein_data[pr]['direct'][temp]
                print  '>>>>>>>'
                sys.exit()
print 'finished computing protein data; launching h5'

#peptides=[np.array([a.split(';')[0] for a in allpep[i][1]['data'][:,['fields']=='Majority protein IDs']]) for i, b in enumerate(batches)]
#compep=np.array(list(set(peptides[0]).intersection(*peptides)))
#comindexpep=np.array([np.where(np.in1d(peptides[i],compep))[0][np.argsort(peptides[i][np.in1d(peptides[i],compep)])] for i, b in enumerate(batches)])





#### reading isample metadata
with open(folder_sample_metadata+file_sample_metadata, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
    metadata = np.array([row for row in reader])



  #### bring additional metadata
    #### reorganise metadata frmo the rna file
#with open(rna_sample_metadata_file, 'rb') as csvfile:
 #   reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
  #  rnametadata = np.array([row for row in reader])

   # rnameta=h5py.File('metadata_sample',mode='w')
    #for i in np.arange(rnametadata.shape[1]): rnameta[rnametadata[0,i]]=rnametadata[1:,i]

#    f = h5py.File(datafile,mode='r')
 #    add_sample_cov=np.array([rnameta['growing_conditions_rnaseq'][:][np.squeeze(np.where(rnameta['friendly'][:]==fd)[0])]=='Feeder-dependent'  for fd in f['sample_meta/friendlyID'][:]])+0
  #   donor_add=np.array([rnameta['donor'][:][np.squeeze(np.where(rnameta['friendly'][:]==fd)[0])] for fd in f['sample_meta/friendlyID'][:]])
# additional_covariates=np.array([ add_sample_cov[donor_add==fd].mean() for fd in f['meanpheno/sample_meta/donorID'][:]])[:,None]
 #       rnameta.close()



#meta_fields=np.array(['pluri_novelty', 'reprogramming','gender','friendlyID','disease', 'derived_from_cell_type', 'covariates','growing_conditions', 'in_flagship','donorID','pluri_raw', 'sampleID','mq_batch'])
#metadata=metadata[:,np.in1d(metadata[0],meta_fields)]

for im,m in enumerate(metadata[0]):
    joindata ['lines'][m]=np.zeros( joindata ['lines']['line_ID'].shape[0],dtype='|S60')
    for il,l in enumerate(joindata ['lines']['line_ID']):
        temp=np.where(metadata[:,0]==l)[0]
        if len(temp)>0:
            joindata ['lines'][m][il]=metadata[temp[0],im]


###  here choose lines for each donor
joindata ['lines']['select_line_donor']=np.copy(joindata ['lines']['donor']); joindata ['lines']['select_line_donor'][:]=''

for donor in np.unique(joindata ['lines']['donor']):
    donor_lines=np.where(joindata ['lines']['donor']==donor)[0]
    line_selected=donor_lines[np.argmax(joindata ['lines']['count_pep'][donor_lines])]
    joindata ['lines']['select_line_donor'][line_selected]=donor
    print line_selected


#sys.exit()


### here put data in a matrix datatypes x genes x lines
matrix_protein_data={'intensity':{}}
for key in protein_data.values()[0]['intensity'].keys():
    matrix_protein_data['intensity'][key] = np.array([protein_data[pr]['intensity'][key] for ipr,pr in enumerate(joindata['data']['protein'])])
#matrix_protein_data['peptide_count_razor_filter']=np.array([ protein_data[pr]['peptide_count_razor_filter']for ipr,pr in enumerate(joindata['data']['protein'])])


f5 = h5py.File(folder_dest+file_protein_quantification+'.h5',mode='w')

# sample metadata
meta=f5.create_group('metadata')
dumpDictHdf5( joindata['lines'],meta.create_group('lines'))

#protein metadata
metap=f5.create_group('metadata/protein_info')
for indf,f in enumerate(joindata ['fields']['protein_info'][0]):
    metap[f]=joindata ['data']['protein_info'][0][indf]
metap['Protein']=joindata['data']['protein']
### here bring additionla informatyion for protein from ensembl and gencode
metap['Protein_main_isoform']=np.array([pr.split('-')[0] for pr in metap['Protein']])

metapdict2=get_protein_info(metap['Protein_main_isoform'],ensembl_file,gencode_file)
dumpDictHdf5(metapdict2,metap)


# protein and peptide abundances
f5.create_group("Protein_Intensity")
dumpDictHdf5( matrix_protein_data['intensity'],f5["Protein_Intensity"])
f5.close()



############ peptide files

f5p = h5py.File(folder_dest+file_protein_quantification+'_peptides.h5',mode='w')

 # sample metadata
meta=f5p.create_group('metadata')
dumpDictHdf5( joindata['lines'],meta.create_group('lines'))

     #protein metadata
metap=f5p.create_group('metadata/protein_info')
for indf,f in enumerate(joindata ['fields']['protein_info'][0]):
    metap[f]=joindata ['data']['protein_info'][0][indf]
metap['Protein']=joindata['data']['protein']
    ### here bring additionla informatyion for protein from ensembl and gencode
metap['Protein_main_isoform']=np.array([pr.split('-')[0] for pr in metap['Protein']])
metapdict2=get_protein_info(metap['Protein_main_isoform'],ensembl_file,gencode_file)
dumpDictHdf5(metapdict2,metap)

f5p.create_group("proteins")
for ipr,pr in enumerate(joindata['data']['protein']):
    gpr=f5p.create_group("proteins/"+pr)
    gpr.create_dataset('peptide_intensity',data=protein_data[pr]['peptide_intensity_norm'])
    gpr.create_dataset('peptide',data=protein_data[pr]['peptide'])
f5p.close()

#f5 = h5py.File(folder_dest+file_protein_quantification,mode='r')
#f5p = h5py.File(folder_dest+file_protein_quantification+'_peptides.h5',mode='r')















