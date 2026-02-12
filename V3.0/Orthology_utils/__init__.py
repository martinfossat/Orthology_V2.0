import json
import Fossat_utils as FU
from .MEDOC_tools import *
from .File_tools import *
from .Organism_dict import *
from .Homology_tools import *
from .Plot_tools import *

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
        fn,fp,NCPR,FCR,kappa=FU.get_seq_properties(seq,window_size=5)
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
    import requests
    import time
    wait=0
    while wait<max_wait:
        r2=requests.post(url='https://biit.cs.ut.ee/gprofiler/api/orth/orth/',
                         json={'organism':organism_org,'target':organism_target,'query':gene_name})
        out=r2.json()

        try:
            r2=requests.post(url='https://biit.cs.ut.ee/gprofiler/api/orth/orth/',
                             json={'organism':organism_org,'target':organism_target,'query':gene_name})
            out=r2.json()
            break
        except:
            time.sleep(1)
            wait+=1
            print('Waited '+str(wait)+' seconds')
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
