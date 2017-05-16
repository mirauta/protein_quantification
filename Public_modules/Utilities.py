
import csv as csv
import numpy as np



def write_data_csv(file_name,data,delimiter=','):
    dim=len(data)
    with open(file_name, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=delimiter, quotechar='"')
        for n in range(data[0].shape[0]):
            writer.writerow([data[c][n] for c in range (dim)])


def load_csv_old(file_name,delimiter=',', quotechar='"', floatstart=0, int_flag=False,fields_to_load=None, skip_before=0,header=0,only_header=0,skip_after=0):
	with open(file_name, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
        	for ss in range(skip_before):
           		next(reader, None)
                if header:
			headers = next(reader, None)
		if only_header:
			if fields_to_load!=None:
				return [headers[fields_to_load]]
			return [headers]
		for ss in range(skip_after):
                        next(reader, None)
          	if floatstart>0:
            		if int_flag:
                		names, dat = np.array([[row[:floatstart],np.array([int(row[j]) for j in range(floatstart, len(row))])] for row in reader]).T
            		else:
                		names, dat = np.array([[row[:floatstart],np.array([float(row[j]) for j in range(floatstart, len(row))])] for row in reader]).T
            		dat=np.vstack(dat)
            		return [headers,names,dat]
        	else:
            		dat = np.array([row for row in reader])
            		return [headers,dat]



def load_csv(file_name,delimiter=',', quotechar='"',fields_to_load=None, skip_before=0,header=0,only_header=0,skip_after=0):
        with open(file_name, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
                for ss in range(skip_before):
                        next(reader, None)
                if header:
                        headers = next(reader, None)
                else:
                    headers=[]
                if only_header:
                        if fields_to_load!=None:
                                return [np.array(headers)[fields_to_load]]
                        return [np.array(headers)]
                for ss in range(skip_after):
                        next(reader, None)
                if fields_to_load==None:
			dat=[row for row in reader]
		else:
			dat=[[row[fnr]  for fnr in  fields_to_load] for row in reader]

                dat=[np.array([dat[i][j] for i in np.arange(len(dat))]) for j in np.arange(len(dat[0]))]
                #dat = np.array([np.array(row)[fields_to_load] for row in reader])
                        #dat=np.vstack(dat)
	return [np.array(headers)[fields_to_load],dat]
                       #dat=np.vstack(dat)
                       #rdat=np.vstack(dat)
                       #dat=np.vstack(dat)
                       #rray(dat=np.vstack(dat)

def load_data_from_maxquant_output(folder_data='',  name='',data_file='.proteinGroups.txt', field_data='',float_data=1, field_metadata=[]):
    data_path=folder_data+name+data_file
    allfields=np.array(load_csv(file_name=data_path,delimiter='\t', header=1,only_header=1)[0])

    metadata={'fields': allfields[np.in1d(allfields,field_metadata)]}
    data={'fields': allfields[np.where([field_data in field for field in allfields])[0]]}

    alldata=load_csv(file_name=data_path,fields_to_load=np.where(np.in1d(allfields,data['fields']))[0],delimiter='\t', header=1,       only_header=0)
    allmetadata=load_csv(file_name=data_path,fields_to_load=np.where(np.in1d(allfields,metadata['fields']))[0],delimiter='\t',         header=1,only_header=0)

                                                                     #if float_data:
                                                                               #                  data['data']=np.squeeze(np.array(alldata[1])).astype(float)
                                                                                       #                             else:
                                                                                                #                      data['data']=np.squeeze(np.array(alldata[1])).astype('S3')
                                                                                                        # metadata['data']=np.squeeze(np.array(allmetadata[1])).astype('S60')
    if float_data: data['data']=np.squeeze(np.array(alldata[1])).astype(float)
    else:    data['data']=alldata[1]
    metadata['data']=allmetadata[1]
    return [allfields,metadata,data]




def dumpDictHdf5(RV,f5):
        """ Dump a dictionary where each page is a list or an array """
        for key in RV.keys():
                f5.create_dataset(name=key,data=np.array(RV[key]),chunks=True,compression='gzip')




def rank_transform(x=None, ref=None):
	_doc='rank transform x into ref;keep the range'
        refn0=ref[ref>0]
        xn0=x[x>0]
        sref=np.sort(refn0)
        index= np.linspace(0,len(refn0)-0.1, num=len(xn0)).astype(int);

        orderx=np.argsort(xn0)

        xtr=x*0
        xtr[x>0]=sref[index][np.argsort(orderx)]
        return xtr




def normalize_reference(data, ntype='',referencedata=None):
	if ntype=='no':
		return data

	elif ntype=='qnorm':
        	return np.array([rank_transform(x=data[:,l],ref=referencedata) for l in range(data.shape[1])]).T

	elif ntype=='snorm':
		return np.array([data[:,l]/np.nansum(data[:,l])*np.nansum(referencedata) for l in range(data.shape[1])]).T


def matchID(id1,id2):
    """ match id1 and id2 """
    idx1 = []
    idx2 = []
    for i in range(id1.shape[0]):
        if id1[i] in id2.tolist():
            idx1.append(i)
            idx2.append(np.where(id2==id1[i])[0][0])
    idx1 = np.array(idx1)
    idx2 = np.array(idx2)
    print (id1[idx1]==id2[idx2]).all()  # assert that's right
    return idx1, idx2

def matchID_order(id1,id2):
    idx1=np.where(np.in1d(id1,id2))[0]
    idx1=idx1[np.argsort(id1[idx1])]
    return 1

