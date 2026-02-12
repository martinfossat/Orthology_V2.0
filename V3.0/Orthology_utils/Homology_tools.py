import skbio
import numpy as np

global AAs_type,AAs_scores

AA_types=np.array(["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])
AA_scores=np.array([5, 13, 8, 6, 8, 8, 10, 5, 6, 5, 7, 7, 10, 7, 7, 5, 5, 5, 15, 8])

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
    for aa in range(len(AA_types)):
        for s in range(len(seq_test)):
            count+=(str(seq_test[s])==AA_types[aa])*AA_scores[aa]
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
    temp_orth_top=score_list[nans_bool]
    top_ind_top=int(np.ceil(len(temp_orth_top)*(1-top_fraction)))
    args_top=np.argsort(score_list)[min(top_ind_top,len(temp_orth_top)-1):]
    if names!=[]:
        out_name=names[args_top[-1]]
    else:
        out_name=None
    return score_list[args_top],len_ratio[args_top],out_name

def get_all_homologies(ref_orga,orth_ref,label,specie,save_homo,min_len_ratio,top_iso_fraction,factor_length_ratio=False):
    save_temp_ids=[]
    save_temp_score=[]
    save_len_ratio=[]

    for orth in save_homo[ref_orga][specie][orth_ref]:
        temp_score=[]
        temp_len_ratio=[]
        temp_iso_id_pair=[]

        for prot_ref in save_homo[ref_orga][specie][orth_ref][orth]:
            for prot in save_homo[ref_orga][specie][orth_ref][orth][prot_ref]:
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
