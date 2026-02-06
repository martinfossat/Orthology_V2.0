import json
import qcan_utils as QC
import Fossat_utils as FU
import numpy as np

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
        #self.SeqProp[seq_id][label][bounds_label]['q_at_pH']={}

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

def find_matching_folded_domains(bounds,bounds_ref,seq_oth,seq_ref,length_max_ratio=0.8,homology_min=0.4):
    bound_match=[]
    match_score=[]
    bounds_all_ref=[]
    bounds_all_oth=[]
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
                    bounds_all_oth+=[[bounds[i,0],bounds[i,1]]]

    bounds_not_folded_ref=get_bounds_inverted(bounds_all_ref,seq_ref)
    bounds_not_folded=get_bounds_inverted(bounds_all_oth,seq_oth)
    # if len(bounds_not_folded_ref)!=len(bounds_not_folded):
    #     print("Problem")
    temp_ref=[]
    temp=[]
    for i in range(len(bounds_not_folded_ref)):
        if not (bounds_not_folded_ref[i,0]==bounds_not_folded_ref[i,1] or bounds_not_folded[i,0]==bounds_not_folded[i,1]):
            temp+=[bounds_not_folded[i]]
            temp_ref+=[bounds_not_folded_ref[i]]
    bounds_not_folded_ref=np.array(temp_ref)
    bounds_not_folded=np.array(temp)
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
    # if bounds[0][0]==0:
    #     # bounds_inv=[]
    #     bounds_inv=[[0,0]]
    # else :
    bounds_inv=[[0]]
    for b in range(len(bounds)):

        bounds_inv[-1]+=[bounds[b][0]]
        # if bounds[b][1]!=len(seq):
        bounds_inv+=[[bounds[b][1]]]
    # if bounds[b][1]!=len(seq):
    bounds_inv[-1]+=[len(seq)]

    return np.array(bounds_inv)
def get_top_x_pct(score_list,top_fraction,names=[]):
    score_list=np.array(score_list)
    nans_bool=np.invert(np.isnan(score_list))
    temp_orth_top=np.array(score_list)[nans_bool]
    top_ind_top=int(np.ceil(len(temp_orth_top)*(1-top_fraction)))
    args_top=np.argsort(score_list)[min(top_ind_top,len(temp_orth_top)-1):]
    if names!=[]:
        out_name=names[args_top[-1]]
    else:
        out_name=None
    return score_list[args_top],out_name

def clean_seq(seq) :
    seq = seq.replace('*', '')
    seq = seq.replace('U', '')
    seq = seq.replace('X', '')
    return seq
def get_all_homologies(ref_orga,orth_ref,label,specie,save_homo,min_len_ratio,top_iso_fraction):
    save_temp_ids=[]
    save_temp=[]
    for orth in save_homo[ref_orga][specie][orth_ref]:
        temp=[]
        for prot_ref in save_homo[ref_orga][specie][orth_ref][orth]:
            for prot in save_homo[ref_orga][specie][orth_ref][orth][prot_ref]:
                temp_region=[]
                # names_iso_temp=[]

                for region in save_homo[ref_orga][specie][orth_ref][orth][prot_ref][prot][label]:
                    reg1=np.array(region.split('_&_')[0].split('_'),dtype=int)
                    len1=reg1[1]-reg1[0]
                    reg2=np.array(region.split('_&_')[1].split('_'),dtype=int)
                    len2=reg2[1]-reg2[0]

                    if min(len1,len2)/max(len1,len2)<min_len_ratio:
                        print(len1,len2)
                        continue
                    temp_region+=[float(
                        save_homo[ref_orga][specie][orth_ref][orth][prot_ref][prot][label][region]['Homology_ratio'])]
                if len(temp_region)!=0:
                    temp+=[np.mean(temp_region)]
        save_temp+=[np.mean(get_top_x_pct(temp,top_iso_fraction)[0])]
        save_temp_ids+=[orth]
    return save_temp,save_temp_ids

def plot_2d_hist(x,y,xlabel,ylabel,x_ticks,y_ticks,name,binwidth):
    import matplotlib
    matplotlib.use("pgf")
    from matplotlib.backends.backend_pgf import FigureCanvasPgf
    matplotlib.backend_bases.register_backend('pdf',FigureCanvasPgf)
    pgf_with_latex={
        "text.usetex":True,
        "pgf.preamble":
            r'\usepackage{color}',
        "font.family":"Times New Roman"}
    import matplotlib
    matplotlib.rcParams.update(pgf_with_latex)
    import matplotlib.pyplot as plt

    plt.figure(figsize=(4,4))
    plt.subplot2grid((8,8),(1,0),colspan=7,rowspan=7)

    bins_x=np.arange(0,np.nanmax(x)+binwidth,binwidth)
    bins_y=np.arange(0,np.nanmax(y)+binwidth,binwidth)

    hist=np.histogram2d(x,y,bins=[bins_x,bins_y])[0]

    plt.xticks(x_ticks,x_ticks)
    plt.yticks(y_ticks,y_ticks)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.plot([0,1],[0,1],linestyle='--',color='k')
    plt.imshow(hist,aspect='auto',origin='lower',interpolation='None',cmap='Spectral',extent=(0.,np.nanmax(y),0.,np.nanmax(x)))

    plt.subplot2grid((8,8),(1,7),colspan=1,rowspan=7)
    bins=np.arange(0,np.nanmax(x)+binwidth,binwidth)
    plt.ylim(0,np.nanmax(x))
    plt.yticks([],[])
    plt.xticks([],[])
    plt.hist(x,bins=bins,orientation='horizontal',histtype='step',color='g')

    plt.subplot2grid((8,8),(0,0),colspan=7,rowspan=1)
    bins=np.arange(0,np.nanmax(y)+binwidth,binwidth)
    plt.xlim(0,np.nanmax(y))
    plt.yticks([],[])
    plt.xticks([],[])
    plt.hist(y,bins=bins,histtype='step',color='g')
    plt.savefig(name)
    plt.close()


def plot_hist(save_all_prop,region_types,species,label_dic,bin_properties):
    import matplotlib
    matplotlib.use("pgf")
    from matplotlib.backends.backend_pgf import FigureCanvasPgf
    matplotlib.backend_bases.register_backend('pdf',FigureCanvasPgf)
    pgf_with_latex={
        "text.usetex":True,
        "pgf.preamble":
            r'\usepackage{color}',
        "font.family":"Times New Roman"}
    import matplotlib
    matplotlib.rcParams.update(pgf_with_latex)
    import matplotlib.pyplot as plt

    colors=[plt.cm.gist_rainbow(i/(len(species))) for i in range(len(species))]
    type_ensemble=[]
    for region_type in region_types:
        FU.check_and_create_rep('Plots/'+region_type)
        for prop in save_all_prop:
            if not prop in bin_properties:
                continue
            if not prop in type_ensemble :
                type_ensemble+=[prop]
            j=0
            if not region_type in save_all_prop[prop]:
                continue
            for orga in save_all_prop[prop][region_type]:
                # if prop in bin_properties:
                temp_tuple=()
                #Need to remove and do explicitly, this nis not necessary
                for i in range(2):
                    if bin_properties[prop][i] is None :

                        if i==0:
                            modulo=np.amin(save_all_prop[prop][region_type][orga])%bin_properties[prop][2]
                            temp_tuple+=(np.amin(save_all_prop[prop][region_type][orga])-modulo-bin_properties[prop][2],)

                        else :
                            modulo=np.amax(save_all_prop[prop][region_type][orga])%bin_properties[prop][2]
                            temp_tuple+=(np.amax(save_all_prop[prop][region_type][orga])+(1-modulo)+bin_properties[prop][2],)
                    else :
                        if i==0:
                            temp_tuple+=(bin_properties[prop][i]-bin_properties[prop][2],)
                        else :
                            temp_tuple+=(bin_properties[prop][i]+bin_properties[prop][2],)
                temp_tuple+=(bin_properties[prop][2],)
                bins=np.arange(temp_tuple[0],temp_tuple[1],temp_tuple[2])

                bins_center=(bins[1:]+bins[:-1])/2.
                hist,bins_temp=np.histogram(save_all_prop[prop][region_type][orga],bins=bins)
                if len(hist)>0:
                    hist=hist/np.amax(hist)
                    plt.step(bins_center,hist,linewidth=1,color=colors[j],label=orga)
                    j+=1
            if prop in label_dic:
                xlabel=label_dic[prop]
            else :
                xlabel=prop
            plt.xlabel(xlabel)
            plt.xlim(bin_properties[prop][0],bin_properties[prop][1])
            plt.legend()
            plt.savefig('Plots/'+region_type+"/Hist_ensemble_"+prop+'_'+region_type+'.pdf')
            plt.close()

def plot_bar_charts(save_all_prop,region_types,species,ensembles,name,label_dic_short):
    import matplotlib
    matplotlib.use("pgf")
    from matplotlib.backends.backend_pgf import FigureCanvasPgf
    matplotlib.backend_bases.register_backend('pdf',FigureCanvasPgf)
    pgf_with_latex={
        "text.usetex":True,
        "pgf.preamble":
            r'\usepackage{color}',
        "font.family":"Times New Roman"}
    import matplotlib
    matplotlib.rcParams.update(pgf_with_latex)
    import matplotlib.pyplot as plt

    colors=[plt.cm.gist_rainbow(i/(len(species))) for i in range(len(species))]

    type_ensemble=[]
    for region_type in region_types:
        i=0
        FU.check_and_create_rep('Plots/'+region_type)
        N_offset=1
        for prop in save_all_prop:
            if not prop in ensembles:
                continue
            if not prop in type_ensemble:
                type_ensemble+=[prop]
            j=0
            if not region_type in save_all_prop[prop]:
                continue
            for orga in save_all_prop[prop][region_type]:

                if i!=0:
                    pre='_'
                else:
                    pre=''
                pos=i-0.5+(j+1)/(len(species)+N_offset)
                if prop=='pct_folded':
                    norm=100.
                else :
                    norm=1.
                value=np.mean(np.array(save_all_prop[prop][region_type][orga])/norm)
                val_std=np.std(np.array(save_all_prop[prop][region_type][orga])/norm)
                plt.bar(pos,value,color=colors[j],width=1/(len(species)+N_offset),label=pre+orga,yerr=val_std,capsize=3)
                j+=1
            i+=1
        plt.title('Predicted ensemble properties')
        plt.legend()
        xticks=[i for i in range(len(type_ensemble))]
        x_labels=[]
        for t in range(len(type_ensemble)):
            if type_ensemble[t] in label_dic_short:
                x_labels+=[label_dic_short[type_ensemble[t]]]
            else :
                x_labels+=[type_ensemble[t].replace("_"," ")]
        plt.xticks(xticks,x_labels,rotation=45,ha='right',va='top')
        plt.tight_layout()
        plt.savefig('Plots/'+region_type+"/Ensemble_properties_"+name+".pdf")
        plt.close()