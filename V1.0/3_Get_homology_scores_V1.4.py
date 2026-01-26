import skbio
import json
import matplotlib
import Fossat_utils as FU
import os
matplotlib.use("pgf")
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
pgf_with_latex = {
    "text.usetex": True,
    "pgf.preamble":
        r'\usepackage{color}',
    "font.family": "Times New Roman" }
import matplotlib
matplotlib.rcParams.update(pgf_with_latex)
import argparse
import numpy as np

def clean_seq(seq) :
    seq = seq.replace('*', '')
    seq = seq.replace('U', '')
    seq = seq.replace('X', '')
    return seq

def get_self_alignement(seq_test,AA_type,AA_scores) :
    count=0
    for aa in range(len(AA_type)):
        for s in range(len(seq_test)):
            count+=(str(seq_test[s])==AA_type[aa])*AA_scores[aa]
    return count

def test_overlap(bounds,bounds_ref,fct_thr):
    bound1=bounds[0]
    bound2=bounds[1]
    bound1_ref=bounds_ref[0]
    bound2_ref=bounds_ref[1]

    bool_test=(((bound1>bound1_ref and bound1<bound2_ref) or
             (bound2>bound1_ref and bound2<bound2_ref)) or
            ((bound1_ref>bound1 and bound1_ref<bound2) or
             (bound2_ref>bound1 and bound2_ref<bound2)))

    range=np.arange(bound1,bound2)
    range_ref=np.arange(bound1_ref,bound2_ref)
    full_range=np.arange(min(bound1,bound1_ref),max(bound2,bound2_ref))
    if len(full_range)!=0:
        fct_overlap=min(len(np.intersect1d(range,full_range)),len(np.intersect1d(range_ref,full_range)))/len(full_range)
        bool_test=fct_thr<fct_overlap
    else :
        bool_test=False

    return bool_test


AA_type,AA_scores=FU.get_self_homology_score()

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Checks the orthologs of a gene name list.\n
    The user needs to specify an original species name, corresponding to the gene name list, and a reference specie, which may be used for 
    homology comparison, in other steps. Additional species may be given, but those two first are required.\n
    Species name must follow the id in gProfiler (https://biit.cs.ut.ee/gprofiler/page/organism-list).\n
    The output is a file containing the gene name and gene IDs for all species, and that is required to use subsequent programs.""")
    parser.add_argument("--file_orthologs","-fo",help='File name where the file is a json file containing created by the Ortholog program')
    parser.add_argument("--file_IDRs","-fi",help='File name where the file is a json file containing created by the Get_IDR program, contaning IDR domain boundaries')
    parser.add_argument("--output_file_name","-of",help='Name of the outputfile. Default is Gene_names.txt')
    parser.add_argument("--min_overlap_fraction","-mof",help='Minimum overlap in the offset corrected domain boundaries for two regions to be compared')
    # parser.add_argument("--original_specie","-os",help="Original specie")
    parser.add_argument("--reference_specie","-rs",help="Reference specie")
    parser.add_argument("--additional_species","-as",help="Additional species",default=[],nargs='+')
    parser.add_argument("--max_size_factor","-msf",help="Factor for the length of the otrholog array compared to the input gene list size. Must be integer, bigger number is slower, but if many orhtolog exists, may be necessary")
    parser.add_argument("--delete_cross_refs","-dcr",help="Whether gene name that share an ortholog should be delete, so they don't appear twice.  1 is delete, 0 is keep.")
    args = parser.parse_args()

    if args.min_overlap_fraction :
        try :
            min_overlap_fraction=float(args.min_overlap_fraction)
        except :
            print("Wrong format for mof")
        if min_overlap_fraction<0 or min_overlap_fraction>1.:
            print("Wrong value for mof. Please use a number between 0 and 1")

    else :
        min_overlap_fraction=0.8
        print("No argument for mof, defaulting to "+str(min_overlap_fraction))

    if args.file_orthologs :
        filename=args.file_orthologs
    else :
        filename='Gene_orthology.json'

    if args.file_IDRs :
        filename_IDRs=args.file_IDRs
    else :
        filename_IDRs='Proteins_IDR_properties.json'

    if args.output_file_name :
        output_file=args.output_file_name

    else :
        output_file='Gene_homologies.json'

    if args.reference_specie :
        ref_org=args.reference_specie
    else :
        print("You must specify a reference specie")
        quit(1)

    if args.additional_species :
        other_organism=args.additional_species
    else :
        other_organism=[]

    file_path=os.path.realpath(__file__)

    with open(os.path.dirname(file_path)+"/g_Profiler_Organisms_names_dic.json", "r") as fp:
        dictionary_organisms_gprofiler = json.load(fp)

    clean_sequence=True

    dictionary_organisms_gprofiler_inv=dict(zip(dictionary_organisms_gprofiler.values(), dictionary_organisms_gprofiler.keys()))

    all_orga=[ref_org]+other_organism

    f=open(filename)
    Orthology_all=json.load(f)

    local=True
    N1=0
    cross_check=False

    f=open(filename_IDRs)
    IDRs=json.load(f)

    Sequence_homology_all={}

    Sequence_homology_all[ref_org]={oth:{} for oth in all_orga}
    Sequence_homology_all[ref_org]={oth:{} for oth in all_orga}

    for orths in Orthology_all:
        print(100*N1/len(Orthology_all))
        N1+=1
        # I should make an orthology score matrix here
        # I need a reference organism, which will carry the scores for the others
        # They will be in order
        # I need an all organisms AND orthologs AND all versions AND all IDRs comparison (Yeah that's a lot)
        if not ref_org in (Orthology_all[orths].keys()):
            continue
        for orga in Orthology_all[orths]:
            #  dimension for each orga
            for N_orth in Orthology_all[orths][orga]:
                # dimension for each ortolog
                for N_orth_ref in Orthology_all[orths][ref_org]:
                    orth_ref_id=Orthology_all[orths][ref_org][N_orth_ref]['gene_id']
                    orth_id=Orthology_all[orths][orga][N_orth]['gene_id']
                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id]={}
                     # These are version of the saem gene (i.e. alleles)
                    for vers in Orthology_all[orths][orga][N_orth]['version']:
                        for vers_ref in Orthology_all[orths][ref_org][N_orth_ref]['version']:
                            ref_id=Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq_id']
                            seq_id=Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']
                            # I will need the ortholog ID too
                            # Actually I would like the self, commenting out
                            #Sequence_homology_all[ref_org][orga][ref_id]=seq_id
                            seq_ref_raw = Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq']
                            seq_raw = Orthology_all[orths][orga][N_orth]['version'][vers]['seq']
                            if clean_sequence :
                                seq_raw=clean_seq(seq_raw)
                                seq_ref_raw= clean_seq(seq_ref_raw)

                            seq_ref=skbio.sequence.Protein(seq_ref_raw)
                            seq=skbio.sequence.Protein(seq_raw)

                            if local:
                                try :
                                    d=skbio.alignment.local_pairwise_align_protein(seq,seq_ref)
                                except :
                                    print("Problem with sequences ")
                                    print(seq)
                                    print(seq_ref)

                            else:
                                try :
                                    d=skbio.alignment.global_pairwise_align_protein(seq,seq_ref)
                                except :
                                    print("Problem with sequences ")
                                    print(seq)
                                    print(seq_ref)
                            norm=get_self_alignement(seq_ref,AA_type,AA_scores)

                            # Here the zero index is just to have the same dimensions as the FD and IDRs
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]={}
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['all']=[{}]
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['all'][0]['norm']=str(d[1]/norm)
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['all'][0]['raw']=str(d[1])

                            # Getting the relevant IDR information
                            temp_ref=IDRs[Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq_id']]
                            temp=IDRs[Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']]

                            ################### Folded domain #########################
                            # Geting the folded domains bounds
                            bounds_folded_ref=[[0]]

                            N_folded_ref=0
                            for i in range(len(temp_ref['dis_doms']['bounds'])):
                                bounds_folded_ref[N_folded_ref]+=[temp_ref['dis_doms']['bounds'][N_folded_ref][0]]
                                bounds_folded_ref+=[[temp_ref['dis_doms']['bounds'][N_folded_ref][1]]]
                                N_folded_ref+=1
                            bounds_folded_ref[-1]+=[len(seq_ref_raw)]

                            bounds_folded=[[0]]
                            N_folded=0
                            for i in range(len(temp['dis_doms']['bounds'])):
                                bounds_folded[N_folded]+=[temp['dis_doms']['bounds'][N_folded][0]]
                                bounds_folded+=[[temp['dis_doms']['bounds'][N_folded][1]]]
                                N_folded+=1
                            bounds_folded[-1]+=[len(seq_raw)]

                            bounds_dis=temp['dis_doms']['bounds']
                            bounds_dis_ref=temp_ref['dis_doms']['bounds']

                            #### Okay this is the hard bit : I will try to offset the sequences based on the difference of length of IDR

                            bounds_folded=np.array(bounds_folded)
                            bounds_folded_ref=np.array(bounds_folded_ref)
                            # I will need to save the effective bounds for comparison

                            bounds_IDR_eff=[]
                            #first dimension is for ref. second dimension is for seq

                            offset=0

                            folded_ind_ref=[]
                            folded_ind=[]
                            disodered_ind_ref=[]
                            disodered_ind=[]

                            for i in range(min(len(bounds_folded),len(bounds_folded_ref))):
                                #this is for the folded domains
                                dom_len_ref=bounds_folded_ref[i][1]-bounds_folded_ref[i][0]
                                dom_len=bounds_folded[i][1]-bounds_folded[i][0]
                                offset+=dom_len_ref-dom_len

                                #This is because of the "zero filling" of folded domains
                                if i>0:
                                    #Now for the disordered domains
                                    dom_len_ref=bounds_folded_ref[i][0]-bounds_folded_ref[i-1][1]
                                    dom_len=bounds_folded[i][0]-bounds_folded[i-1][1]
                                    offset+=dom_len_ref-dom_len

                                    # I should replace this with a percent overlap

                                    if test_overlap(bounds_dis[i-1]+offset,bounds_dis_ref[i-1],min_overlap_fraction):
                                        disodered_ind_ref+=[i-1]
                                        disodered_ind+=[i-1]

                                if test_overlap(bounds_folded[i]+offset,bounds_folded_ref[i],min_overlap_fraction):
                                    folded_ind_ref+=[i]
                                    folded_ind+=[i]


                            N_folded_eff=min(N_folded_ref,N_folded)

                            if N_folded_eff==0:
                                continue

                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['FDs']=[]
                            if not cross_check :
                                for FD in folded_ind_ref:
                                    bound1=bounds_folded[FD][0]
                                    bound2=bounds_folded[FD][1]

                                    bound1_ref=bounds_folded_ref[FD][0]
                                    bound2_ref=bounds_folded_ref[FD][1]

                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['FDs']+=[{}]


                                    FD_seq_ref_raw=Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq'][bound1_ref:bound2_ref]
                                    FD_seq_raw=Orthology_all[orths][orga][N_orth]['version'][vers]['seq'][bound1:bound2]
                                    if clean_sequence:
                                        FD_seq_ref_raw = clean_seq(FD_seq_ref_raw)
                                        FD_seq_raw = clean_seq(FD_seq_raw)

                                    FD_seq_ref=skbio.sequence.Protein(FD_seq_ref_raw)
                                    FD_seq = skbio.sequence.Protein(FD_seq_raw)

                                    if local :
                                        d=skbio.alignment.local_pairwise_align_protein(FD_seq,FD_seq_ref)
                                    else:
                                        d=skbio.alignment.global_pairwise_align_protein(FD_seq,FD_seq_ref)

                                    norm=get_self_alignement(FD_seq_ref,AA_type,AA_scores)

                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['FDs'][-1]['norm']=str(d[1]/norm)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['FDs'][-1]['raw']=str(d[1])

                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['FDs'][-1]['region_ref']=[str(bound1_ref),str(bound2_ref)]
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['FDs'][-1]['region']=[str(bound1),str(bound2)]

                            ################### IDRs #########################
                            # This is to make a single loop for both the refrerence and the comapred sequnece
                            N_dis=min(len(temp_ref['dis_doms']['bounds']),len(temp['dis_doms']['bounds']))
                            if N_dis==0:
                                continue

                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs']=[]
                            if not cross_check :
                                for DR in disodered_ind_ref:
                                    bound1=bounds_dis[DR][0]
                                    bound2=bounds_dis[DR][1]

                                    bound1_ref=bounds_dis_ref[DR][0]
                                    bound2_ref=bounds_dis_ref[DR][1]

                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs']+=[{}]
                                    IDR_seq_ref_raw=Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq'][bound1_ref:bound2_ref]
                                    IDR_seq_raw = Orthology_all[orths][orga][N_orth]['version'][vers]['seq'][bound1:bound2]

                                    if clean_sequence:
                                        IDR_seq_ref_raw = clean_seq(IDR_seq_ref_raw)
                                        IDR_seq_raw = clean_seq(IDR_seq_raw)

                                    IDR_seq_ref=skbio.sequence.Protein(IDR_seq_ref_raw)
                                    IDR_seq = skbio.sequence.Protein(IDR_seq_raw)

                                    if local :
                                        d=skbio.alignment.local_pairwise_align_protein(IDR_seq,IDR_seq_ref)
                                    else:
                                        d=skbio.alignment.global_pairwise_align_protein(IDR_seq,IDR_seq_ref)

                                    norm=get_self_alignement(IDR_seq_ref,AA_type,AA_scores)

                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['norm']=str(d[1]/norm)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['raw']=str(d[1])

                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['region_ref']=[str(bound1_ref),str(bound2_ref)]
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['region']=[str(bound1),str(bound2)]

                                    RG_ref=float(IDRs[Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq_id']]['dis_doms']['rgs'][DR])
                                    RG_oth=float(IDRs[Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']]['dis_doms']['rgs'][DR])

                                    RE_ref=float(IDRs[Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq_id']]['dis_doms']['res'][DR])
                                    RE_oth=float(IDRs[Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']]['dis_doms']['res'][DR])

                                    AS_ref=float(IDRs[Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq_id']]['dis_doms']['asphs'][DR])
                                    AS_oth=float(IDRs[Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']]['dis_doms']['asphs'][DR])

                                    NU_ref=float(IDRs[Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq_id']]['dis_doms']['asphs'][DR])
                                    NU_oth=float(IDRs[Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']]['dis_doms']['asphs'][DR])

                                    PR_ref=float(IDRs[Orthology_all[orths][ref_org][N_orth_ref]['version'][vers_ref]['seq_id']]['dis_doms']['prefs'][DR])
                                    PR_oth=float(IDRs[Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']]['dis_doms']['prefs'][DR])

                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['Rg_ratio']=str(RG_oth/RG_ref)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['Re_ratio']=str(RE_oth/RE_ref)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['As_ratio']=str(AS_oth/AS_ref)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['Nu_ratio']=str(NU_oth/NU_ref)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]['IDRs'][-1]['Pref_ratio']=str(PR_oth/PR_ref)

    # for i in Sequence_homology_all:
    #     for j in Sequence_homology_all[i]:
    #         for k in Sequence_homology_all[i][j]:
    #             for l in Sequence_homology_all[i][j][k]:
    #                 for m in Sequence_homology_all[i][j][k][l]:
    #                 with open('test.txt','w') as f:
    #                     json.dump(Sequence_homology_all[i][j][k][l],f)

    with open(output_file,'w') as f:
        json.dump(Sequence_homology_all,f)
