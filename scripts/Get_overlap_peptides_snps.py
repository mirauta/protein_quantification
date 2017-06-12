import sys
import numpy as np
import csv as csv
import pandas
import vcf as vcf

release=sys.argv[1]
folder_annotation=sys.argv[2]
folder_data=sys.argv[3]
name=sys.argv[4]

#release='Ensembl_37.75'
#folder_annotation='/Users/mirauta/Data/Annotation/'
#folder_data='/Users/mirauta/Data/MS/hipsci/TMT/phenotypes/'
#name="hipsci.proteomics.maxquant.uniprot.TMT_batch_14.20170517_peptide_lines_filtered_unique"

snp_VE=vcf.Reader(filename=folder_annotation+release+'/Homo_sapiens_incl_consequences.vcf.gz')
snp_1000=vcf.Reader(filename=folder_annotation+release+'/1000GENOMES-phase_1_EUR.vcf.gz')

toinclude=np.array(['coding_sequence_variant', 'frameshift_variant', 'inframe_deletion', 'inframe_insertion',       'missense_variant',                     'protein_altering_variant',       'splice_acceptor_variant', 'splice_donor_variant',       'splice_region_variant', 'start_lost', 'stop_gained', 'stop_lost'])


peptide_genome_mapi=pandas.read_table(folder_data+name+"_list_"+release+".bed",sep='\t',header=None)
peptide_genome_map=pandas.DataFrame(data=np.array(peptide_genome_mapi)[:,[0,3,1,2,5,1]], index=peptide_genome_mapi.index, columns=["chromosome", "feature_id","start", "end","strand","overlap"])
peptide_genome_map['chromosome']=np.array([c.replace('chr','')for c in peptide_genome_map['chromosome']])
peptide_genome_map=peptide_genome_map.set_index('feature_id',drop=False)

peptides1= peptide_genome_map['feature_id']
print (peptides1.shape)

overlap=pandas.DataFrame(data=np.zeros(len(peptides1)),index=peptides1,columns=['overlap_snp'])
 
for pepname in peptides1:
    overlap['overlap_snp'][pepname]=0
 
    try:pep=peptide_genome_map.transpose()[pepname].icol(0)
    except:pep=peptide_genome_map.transpose()[pepname]
    print (pepname)
    try:
        for reader in snp_1000.fetch(pep['chromosome'], pep['start']-1,pep['start']+3*len(pepname)):
           if reader.INFO['AF'][0]>0.1:
               try:
                   for reader2 in snp_VE.fetch(reader.CHROM, reader.POS-1,reader.POS):
                        temp=np.array([fi  in '_'.join(reader2.INFO['VE'])for fi in toinclude])
                        print (toinclude[temp])
                        overlap['overlap_snp'][pepname]+=np.nansum(temp)
                   print(reader.INFO['AF'][0])
               except:
                   print ('missing VE')
        if pep['end']>(pep['start']+3*len(pepname)):
            for reader in snp_1000.fetch(pep['chromosome'],pep['end']-3*len(pepname),pep['end']):
                if reader.INFO['AF'][0]>0.1:
                    try:
                        for reader2 in snp_VE.fetch(reader.CHROM, reader.POS-1,reader.POS):
                            temp=np.array([fi  in '_'.join(reader2.INFO['VE'])for fi in toinclude])
                            print (toinclude[temp])
                            overlap['overlap_snp'][pepname]+=np.nansum(temp)
                        print(reader.INFO['AF'][0])
                    except:
                        print ('missing VE')
    except:
        1



overlap.to_csv(path_or_buf=folder_data+name+'_overlap_peptide_snp.txt',mode='w', sep='\t', columns=None, header=True, index=True)
