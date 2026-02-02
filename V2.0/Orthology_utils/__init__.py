import json
import qcan_utils as QC
import Fossat_utils as FU
import numpy as np
from sparrow.predictors import batch_predict
import skbio
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
        layers_q,Meso_G_all=QC.MEDOC_wrapper(seq)
        q=QC.get_q_profile_from_F(pH,Meso_G_all,layers_q,298)
        iso_q_ind=np.argwhere(abs(q)==min(abs(q)))[0][0]
        self.SeqProp[seq_id][label][bounds_label]['pI']=str(pH[iso_q_ind])
        self.SeqProp[seq_id][label][bounds_label]['q_at_pH']={}

        for i in range(len(pH_ref)):
            # Find index closest to reference pH
            pH_ind=np.argwhere(abs(pH-pH_ref[i])==min(abs(pH-pH_ref[i])))[0][0]
            # Key is the reference pH, Value is the charge
            self.SeqProp[seq_id][label][bounds_label]['q_at_pH'][str(pH_ref[i])]=str(q[pH_ind])

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
def find_matching_folded_domains(bounds,bounds_ref,seq_oth,seq_ref,length_max_ratio=0.8,homology_min=0.4):
    bound_match=[]
    match_score=[]
    bounds_all_ref=[]
    bonds_all_oth=[]
    for i in range(len(bounds)):
        for j in range(len(bounds_ref)):
            len_ref=bounds_ref[j,1]-bounds_ref[j,0]
            len_oth=bounds[i,1]-bounds[i,0]
            ratio=min(len_ref,len_oth)/max(len_ref,len_oth)
            if ratio>length_max_ratio:
                seq_oth_sub=seq_oth[bounds[i,0]:bounds[i,1]]
                seq_ref_sub=seq_ref[bounds_ref[j,0]:bounds_ref[j,1]]
                homo,normed_homo=get_homology_score(seq_oth_sub,seq_ref_sub,local=True)
                if normed_homo>homology_min:
                    bound_match+=[str(bounds_ref[j,0])+'_'+str(bounds_ref[j,1])+'_&_'+str(bounds[i,0])+'_'+str(bounds[i,1])]
                    match_score+=[[homo,normed_homo]]

                    bounds_all_ref+=[[bounds_ref[j,0],bounds_ref[j,1]]]
                    bonds_all_oth+=[[bounds[i,0],bounds[i,1]]]

    bounds_not_folded_ref=get_bounds_inverted(bounds_all_ref,seq_ref)
    bounds_not_folded=get_bounds_inverted(bonds_all_oth,seq_oth)
    # This is a correction : if the first residue is 0 in one and not the other, there is a disodered domain mismatch
    # which must be removed
    if (bounds_not_folded_ref[0,0]==0 and not bounds_not_folded[0,0]==0):
        bounds_not_folded_ref=np.array([bounds_not_folded_ref[i] for i in range(1,len(bounds_not_folded_ref[i]))])
    elif (bounds_not_folded[0,0]==0 and not bounds_not_folded_ref[0,0]==0):
        bounds_not_folded=np.array([bounds_not_folded[i] for i in range(1,len(bounds_not_folded[i]))])

    return bound_match,match_score,bounds_not_folded,bounds_not_folded_ref

def get_homology_score(seq_ref,seq,local):
    seq_ref=skbio.sequence.Protein(seq_ref)
    seq=skbio.sequence.Protein(seq)
    if local:
        d=skbio.alignment.local_pairwise_align_protein(seq,seq_ref)
    else:
        d=skbio.alignment.global_pairwise_align_protein(seq,seq_ref)

    norm=get_self_alignement(seq_ref,AA_type,AA_scores)
    return d[1],d[1]/norm

def get_self_alignement(seq_test,AA_type,AA_scores) :
    count=0
    for aa in range(len(AA_type)):
        for s in range(len(seq_test)):
            count+=(str(seq_test[s])==AA_type[aa])*AA_scores[aa]
    return count
AA_type,AA_scores=FU.get_self_homology_score()

def get_bounds_inverted(bounds,seq):
    if len(bounds)==0:
        return np.array([[0,len(seq)]])
    if bounds[0][0]==0:
        bounds_inv=[]
    else :
        bounds_inv=[[0]]
    for b in range(len(bounds)):
        if bounds[b][0]!=0:
            bounds_inv[-1]+=[bounds[b][0]]
        if bounds[b][1]!=len(seq):
            bounds_inv+=[[bounds[b][1]]]
    if bounds[b][1]!=len(seq):
        bounds_inv[-1]+=[len(seq)]

    return np.array(bounds_inv)
