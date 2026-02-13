import numpy as np
import json
from .MEDOC_tools import *
import time
import requests
class Seq_Prop_Manager:
    def __init__(self,filename):
        self.filename=filename
        self.SeqProp=self._load_data()
    def _load_data(self):
        try:
            with open(self.filename,'r') as f:
                return json.load(f)
        except (FileNotFoundError,json.JSONDecodeError):
            # Return empty dict if file is missing or broken
            return {}
    def get_bounds(self,seq_id,label):
        bounds=[]
        if not label in self.SeqProp[seq_id]:
            return None,None
        for i in list(self.SeqProp[seq_id][label].keys()):
            bounds+=[[int(i.split('_')[0]),int(i.split('_')[1])]]
        bounds=np.array(bounds)
        return bounds,list(self.SeqProp[seq_id][label].keys())

    def save(self):
        with open(self.filename,'w') as f:
            json.dump(self.SeqProp,f)
    def add_pH_information(self,seq_id,seq,pH_ref,pH,label,bounds_label):
        if seq_id not in self.SeqProp:
            self.SeqProp[seq_id]={label:{}}
        if label not in self.SeqProp[seq_id]:
            self.SeqProp[seq_id][label]={}
        if bounds_label not in self.SeqProp[seq_id][label]:
            self.SeqProp[seq_id][label][bounds_label]={}
        layers_q,Meso_G_all=MEDOC_wrapper(seq)
        q=get_q_profile_from_F(pH,Meso_G_all,layers_q,298)
        iso_q_ind=np.argwhere(abs(q)==min(abs(q)))[0][0]
        self.SeqProp[seq_id][label][bounds_label]['pI']=str(pH[iso_q_ind])

        for i in range(len(pH_ref)):
            # Find index closest to reference pH
            pH_ind=np.argwhere(abs(pH-pH_ref[i])==min(abs(pH-pH_ref[i])))[0][0]
            # Key is the reference pH, Value is the charge
            self.SeqProp[seq_id][label][bounds_label]['q_at_pH_'+str(pH_ref[i])]=str(q[pH_ind])

    def add_sequence_properties(self,seq_id,seq,label,bounds_label):
        if seq_id not in self.SeqProp:
            self.SeqProp[seq_id]={label:{}}
        if label not in self.SeqProp[seq_id]:
            self.SeqProp[seq_id][label]={}
        if bounds_label not in self.SeqProp[seq_id][label]:
            self.SeqProp[seq_id][label][bounds_label]={}
        fn,fp,NCPR,FCR,kappa=get_seq_properties(seq,window_size=5)
        self.SeqProp[seq_id][label][bounds_label]['NCPR']=str(NCPR)
        self.SeqProp[seq_id][label][bounds_label]['FCR']=str(FCR)
        self.SeqProp[seq_id][label][bounds_label]['f+']=str(fp)
        self.SeqProp[seq_id][label][bounds_label]['f-']=str(fn)
        self.SeqProp[seq_id][label][bounds_label]['Kappa']=str(kappa)

    def add_ensemble_properties(self,seq_id,seq,label,bounds_label):
        from sparrow.predictors import batch_predict
        if seq_id not in self.SeqProp:
            self.SeqProp[seq_id]={label:{}}
        if label not in self.SeqProp[seq_id]:
            self.SeqProp[seq_id][label]={}
        if bounds_label not in self.SeqProp[seq_id][label]:
            self.SeqProp[seq_id][label][bounds_label]={}
        temp={}
        temp[bounds_label]=seq[int(bounds_label.split('_')[0]):int(bounds_label.split('_')[1])]
        self.SeqProp[seq_id][label][bounds_label]['rg']=str(batch_predict.batch_predict(temp,network='scaled_rg')[bounds_label][1])
        self.SeqProp[seq_id][label][bounds_label]['re']=str(batch_predict.batch_predict(temp,network='scaled_re')[bounds_label][1])
        self.SeqProp[seq_id][label][bounds_label]['asph']=str(batch_predict.batch_predict(temp,network='asphericity')[bounds_label][1])
        self.SeqProp[seq_id][label][bounds_label]['nu']=str(batch_predict.batch_predict(temp,network='scaling_exponent')[bounds_label][1])
        self.SeqProp[seq_id][label][bounds_label]['pref']=str(batch_predict.batch_predict(temp,network='prefactor')[bounds_label][1])
    def add_percent_folded(self,seq_id,seq,bounds_fd):
        if len(bounds_fd)!=0:
            SUM=np.sum(bounds_fd[:,1]-bounds_fd[:,0])
        else :
            SUM=0
        self.SeqProp[seq_id]["all"]["0_"+str(len(seq))]["pct_folded"]=str(100*SUM/len(seq))
    def get_values(self,seq_id,label,bounds_label,val_type):
        if not val_type in self.SeqProp[seq_id][label][bounds_label]:
            return float('NaN')
        else :
            return float(self.SeqProp[seq_id][label][bounds_label][val_type])
def query_profiler(organism_org,organism_target,gene_name,max_wait=30):
    wait=0
    while wait<max_wait:
        try:
            r2=requests.post(url='https://biit.cs.ut.ee/gprofiler/api/orth/orth/',
                             json={'organism':organism_org,'target':organism_target,'query':gene_name})
            out=r2.json()
            break
        except:
            time.sleep(1)
            wait+=1
            print('Waited '+str(wait)+' seconds')

    if wait>=max_wait:
        print("Maxed out on server waiting time. Either your connection if unstable, or you have been temporarly blocked to many recent connections.")
        quit()

    return out

def get_all_ids(Orthology_all):
    ids=[]
    for orths in Orthology_all:
        for orga in Orthology_all[orths]:
            for gene_id in Orthology_all[orths][orga]:
                ids+=[gene_id]
    return ids

def clean_seq(seq) :
    seq = seq.replace('*', '')
    seq = seq.replace('U', '')
    seq = seq.replace('X', '')
    return seq


def get_seq_properties(seq,window_size=5):
    seq_q=seq_1_let_to_Q(seq)
    fn=get_f(seq_q,-1)
    fp=get_f(seq_q,+1)
    NCPR=fp-fn
    FCR=fp+fn
    d=get_delta(seq_q,window_size)
    dmax=get_delta(make_dmax_sequence(seq_q),window_size)
    if dmax!=0:
        kappa=float(d)/dmax
    else:
        kappa=0

    return fn,fp,NCPR,FCR,kappa

def get_f(window,val):
    num=0
    for i in range(len(window)):
        if window[i]==val:
            num+=1
    return float(num)/len(window)


def get_sigma(window):
    fp=get_f(window,1)
    fn=get_f(window,-1)

    if (fp+fn)!=0:
        out=(fp-fn)**2/(fp+fn)
    else:
        out=0

    return out

def make_dmax_sequence(seq_q):
    seq_q=np.ndarray.tolist(np.sort(seq_q))
    seq_out=seq_q

    return seq_out

def get_kappa(seq,window_size=5):
    #seq_q=[0 for i in range(int(np.ceil(window_size/2.)))]+seq_1_let_to_Q(seq)+[0 for i in range(int(np.ceil(window_size/2.)))]
    #temp=[0 for i in range(int(np.ceil(window_size/2.)))]+np.ndarray.tolist(np.sort(seq_1_let_to_Q(seq)))+[0 for i in range(int(np.ceil(window_size/2.)))]
    seq_q=seq_1_let_to_Q(seq)
    temp=make_dmax_sequence(seq_q)
    d=get_delta(seq_q,window_size)
    # The definition of dmax is not jsut sorting.
    dmax=get_delta(temp,window_size)

    if dmax!=0:
        out=float(d)/dmax
    else:
        out=0
    return out

def get_delta(seq_q,window_sz):
    # this version has no "window border effect"
    tot_sig=get_sigma(seq_q)
    delta_temp=0
    n=0
    for i in range(len(seq_q)):
        if i-np.floor(window_sz/2.)>=0 and i+np.floor(window_sz/2.)<len(seq_q):

            offset_p=int(np.floor(window_sz/2.))
            offset_n=int(np.ceil(window_sz/2.))
        elif i-np.floor(window_sz/2.)<0:

            offset_n=int(np.ceil(window_sz/2.))
            for j in range(int(np.floor(window_sz/2.)+1)):
                if i-j>=0:
                    offset_p=j
        else:

            for j in range(int(np.ceil(window_sz/2.))):
                if i+j<=len(seq_q):
                    offset_n=j
            offset_p=int(np.floor(window_sz/2.))

        window_curr=[seq_q[j] for j in range(i-offset_p,i+offset_n)]
        delta_temp+=(get_sigma(window_curr)-tot_sig)**2
        n+=1

    return delta_temp/n

def seq_1_let_to_Q(seq,mode=1,silent=True):
    if mode==1 : #Corresponds to not counting the salt
        dico = {"A" : 0,
                "C" : 0,
                "E" : -1,
                "D" : -1 ,
                "K" : 1 ,
                "C" : 0 ,
                "Q" : 0 ,
                "H" : 0 ,
                "G" : 0 ,
                "P" : 0 ,
                "F" : 0 ,
                "S" : 0 ,
                "M" : 0 ,
                "I" : 0 ,
                "L" : 0 ,
                "R" : 1 ,
                "T" : 0 ,
                "W" : 0 ,
                "V" : 0 ,
                "Y" : 0 ,
                "N" : 0 ,
                "d" : 0 ,
                "e" : 0 ,
                "h" : 0 ,
                "k" : 0 }
    seq_new=[]
    for i in range(len(seq)) :
        try :
            seq_new+=[dico[seq[i]]]
        except :
            if silent==False :
                print("Did not find a key for "+seq[i]+", ignoring ")
    return seq_new
