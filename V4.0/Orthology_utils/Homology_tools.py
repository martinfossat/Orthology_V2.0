import skbio
import numpy as np
global AAs_type,AAs_scores,annotation,annotation,Gene_annotation

AAs=np.array(["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])
AA_scores=np.array([5,13,8,6,8,8,10,5,6,5,7,7,10,7,7,5,5,5,15,8])

annotation=[]
Gene_annotation=[]

# The initialize function are only passing the variable in the scope of this file
def initialize_gene_annotation(pass_Gene_annotation):
    global Gene_annotation
    Gene_annotation=pass_Gene_annotation

def initialize_gene_annotation_names(pass_annotation):
    global annotation
    annotation=pass_annotation

class Homology_Compare_Manager():
    # This class should have a sequence labels, and a gene label name and a gene labels type
    def __init__(self):
        global annotation
        self.save_all_compare=np.array([])
        self.save_all_ids_ref=np.array([])
        self.save_all_ids_top=np.array([])
        self.save_all_ids_norm=np.array([])
        self.save_all_compare_norm=np.array([])
        self.save_all_compare_top=np.array([])
        self.save_all_len_ratio=np.array([])
        # I think the best is to make a mask of the different gene labels
        self.masks_all={}
        for at in annotation:
            self.masks_all[at]={}
            for an in annotation[at]:
                self.masks_all[at][an]=np.array([],dtype=bool)

    def add_entry(self,best_top_id,best_norm_id,orth_ref,save_temp_top,save_temp_norm,temp_len_ratio_top,temp_len_ratio_norm):
        global Gene_annotation,annotation
        self.save_all_ids_top=np.append(self.save_all_ids_top,[best_top_id])
        self.save_all_ids_norm=np.append(self.save_all_ids_norm,[best_norm_id])
        self.save_all_ids_ref=np.append(self.save_all_ids_ref,[orth_ref])
        self.save_all_compare=np.append(self.save_all_compare,[np.mean(save_temp_top)/np.mean(save_temp_norm)])
        self.save_all_compare_top=np.append(self.save_all_compare_top,[np.mean(save_temp_top)])
        self.save_all_compare_norm=np.append(self.save_all_compare_norm,[np.mean(save_temp_norm)])
        self.save_all_len_ratio=np.append(self.save_all_len_ratio,[np.mean(temp_len_ratio_top)/np.mean(temp_len_ratio_norm)])
        for at in annotation:
            for an in range(len(annotation[at])):
                if orth_ref in Gene_annotation:
                    if Gene_annotation[orth_ref][at] is None:
                        self.masks_all[at][annotation[at][an]]=np.append(self.masks_all[at][annotation[at][an]],[False])
                        continue
                else:
                    self.masks_all[at][annotation[at][an]]=np.append(self.masks_all[at][annotation[at][an]],[False])
                    continue

                if annotation[at][an] in Gene_annotation[orth_ref][at]:
                    self.masks_all[at][annotation[at][an]]=np.append(self.masks_all[at][annotation[at][an]],[True])
                else:
                    self.masks_all[at][annotation[at][an]]=np.append(self.masks_all[at][annotation[at][an]],[False])

    def sort_all(self):
        indices=np.argsort(self.save_all_compare)[::-1]
        self.save_all_compare=self.save_all_compare[indices]
        self.save_all_ids_top=self.save_all_ids_top[indices]
        self.save_all_ids_norm=self.save_all_ids_norm[indices]
        self.save_all_ids_ref=self.save_all_ids_ref[indices]
        self.save_all_compare_top=self.save_all_compare_top[indices]
        self.save_all_compare_norm=self.save_all_compare_norm[indices]
        self.save_all_len_ratio=self.save_all_len_ratio[indices]

    def get_compare_masked(self,at_type,an_type):
        compare=self.save_all_compare[self.masks_all[at_type][an_type]]
        compare_top=self.save_all_compare_top[self.masks_all[at_type][an_type]]
        compare_norm=self.save_all_compare_norm[self.masks_all[at_type][an_type]]
        return compare,compare_top,compare_norm
    def get_all(self):
        return self.save_all_compare,self.save_all_ids_top,self.save_all_ids_norm,self.save_all_ids_ref,self.save_all_compare_top,self.save_all_compare_norm,self.save_all_len_ratio


    def get_compare(self):
        return self.save_all_compare,self.save_all_compare_top,self.save_all_compare_norm
def find_matching_folded_domains(bounds,bounds_ref,seq_oth,seq_ref,length_max_ratio=0.8,homology_min=0.4):
    bound_match=[]
    match_score=[]
    bounds_all_ref=[]
    bounds_all_oth=[]
    for i in range(len(bounds)):
        for j in range(len(bounds_ref)):
            len_ref=bounds_ref[j,1]-bounds_ref[j,0]
            len_oth=bounds[i,1]-bounds[i,0]
            if min(len_ref,len_oth)==0:
                continue
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

    norm=get_self_alignement(seq_ref)
    return d[1],d[1]/norm

def get_self_alignement(seq_test) :
    count=0
    for aa in range(len(AAs)):
        for s in range(len(seq_test)):
            count+=(str(seq_test[s])==AAs[aa])*AA_scores[aa]
    return count

def get_bounds_inverted(bounds,seq):
    if len(bounds)==0:
        return np.array([[0,len(seq)]])
    bounds_inv=[[0]]
    for b in range(len(bounds)):
        bounds_inv[-1]+=[bounds[b][0]]
        bounds_inv+=[[bounds[b][1]]]
    bounds_inv[-1]+=[len(seq)]
    return np.array(bounds_inv)

def get_top_x_pct(score_list,len_ratio,top_fraction,names=[]):
    score_list=np.array(score_list)
    len_ratio=np.array(len_ratio)
    nans_bool=np.invert(np.isnan(score_list))
    score_list=score_list[nans_bool]
    len_ratio=len_ratio[nans_bool]
    if len(score_list)==0:
        return [],[],None
    top_ind_top=int(np.ceil(len(score_list)*(1-top_fraction)))
    args_top=np.argsort(score_list)[min(top_ind_top,len(score_list)-1):]
    out_name=names[args_top[-1]]
    return score_list[args_top],len_ratio[args_top],out_name


def get_all_homologies(ref_orga,orth_ref,label,specie,save_homo,min_len_ratio,top_iso_fraction,factor_length_ratio=False,MLO_only=False):
    # Use MLF is applied on the
    save_temp_ids=[]
    save_temp_score=[]
    save_len_ratio=[]

    for orth in save_homo[ref_orga][specie][orth_ref]:
        temp_score=[]
        temp_len_ratio=[]
        temp_iso_id_pair=[]
        # For thr MLO : we are going to have to loop over all of the combinations, then find the longest to longest for
        # the all label, then save the id and skip until id is reached
        # This measn that the
        if MLO_only:
            max_len_tmp1=0
            id_tmp_1=''
            max_len_tmp2=0
            id_tmp_2=''
            for prot_ref in save_homo[ref_orga][specie][orth_ref][orth]:
                for prot in save_homo[ref_orga][specie][orth_ref][orth][prot_ref]:
                    for region in save_homo[ref_orga][specie][orth_ref][orth][prot_ref][prot]['all']:
                        reg1=np.array(region.split('_&_')[0].split('_'),dtype=int)
                        len1=reg1[1]-reg1[0]
                        reg2=np.array(region.split('_&_')[1].split('_'),dtype=int)
                        len2=reg2[1]-reg2[0]
                        if len1>max_len_tmp1:
                            max_len_tmp1=len1
                            id_tmp_1=prot_ref

                        if len2>max_len_tmp2:
                            max_len_tmp2=len2
                            id_tmp_2=prot


        for prot_ref in save_homo[ref_orga][specie][orth_ref][orth]:
            if MLO_only and prot_ref!=id_tmp_1 :
                continue
            for prot in save_homo[ref_orga][specie][orth_ref][orth][prot_ref]:
                if MLO_only and prot!=id_tmp_2:
                    continue

                temp_score_region=[]
                temp_len_ratio_region=[]

                if not label in save_homo[ref_orga][specie][orth_ref][orth][prot_ref][prot]:
                    continue
                for region in save_homo[ref_orga][specie][orth_ref][orth][prot_ref][prot][label]:
                    reg1=np.array(region.split('_&_')[0].split('_'),dtype=int)
                    len1=reg1[1]-reg1[0]
                    reg2=np.array(region.split('_&_')[1].split('_'),dtype=int)
                    len2=reg2[1]-reg2[0]
                    len_ratio=min(len1,len2)/max(len1,len2)
                    if len_ratio<min_len_ratio:
                        continue
                    temp_len_ratio_region+=[len_ratio]
                    if factor_length_ratio:
                        factor=len_ratio
                    else :
                        factor=1.
                    temp_score_region+=[factor*float(save_homo[ref_orga][specie][orth_ref][orth][prot_ref][prot][label][region]['Homology_ratio'])]

                if len(temp_score_region)!=0:
                    temp_score+=[np.mean(temp_score_region)]
                    temp_len_ratio+=[np.mean(temp_len_ratio_region)]
                    temp_iso_id_pair+=[prot_ref+'_'+prot]

        if len(temp_score)==0:
            save_temp_score+=[float('NaN')]
            save_len_ratio+=[float('NaN')]
            save_temp_ids+=[orth+'__N/A_N/A']

        else :
            temp_score,temp_len_ratio,top_iso_pair=get_top_x_pct(temp_score,temp_len_ratio,top_iso_fraction,names=temp_iso_id_pair)
            save_temp_score+=[np.mean(temp_score)]
            save_len_ratio+=[np.mean(temp_len_ratio)]
            save_temp_ids+=[orth+'__'+top_iso_pair]

    return save_temp_score,save_len_ratio,save_temp_ids
