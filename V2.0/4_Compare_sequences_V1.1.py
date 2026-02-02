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
import Orthology_utils as OU

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
    parser.add_argument("--file_Seq_Prop","-fi",help='File name where the file is a json file containing created by the Get_IDR program, contaning IDR domain boundaries')
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

    if args.file_Seq_Prop :
        filename_Seq_Prop=args.file_Seq_Prop
    else :
        filename_Seq_Prop='Sequence_properties.json'

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
    SeqProp=OU.Seq_Prop_Manager(filename_Seq_Prop)
    Sequence_homology_all={}

    Sequence_homology_all[ref_org]={oth:{} for oth in all_orga}
    val_types_IDRs=['rg','re','asph','nu','pref']
    val_types_charge=['pI']
    val_types_seq=['NCPR', 'FCR','f-','f+','Kappa']
    length_ratio=0.8
    # labels=["all","IDRs","FDs"]

    for orths in Orthology_all:
        # I should make an orthology score matrix here
        # I need a reference organism, which will carry the scores for the others
        # They will be in order
        # I need an all organisms AND orthologs AND all versions AND all Seq_Prop comparison (Yeah that's a lot)
        if not ref_org in (Orthology_all[orths].keys()):
            continue
        for orga in Orthology_all[orths]:
            #  Dimension for each orga
            for orth_id in Orthology_all[orths][orga]:
                if orth_id=='N/A':
                    continue
                # dimension for each ortholog
                for orth_ref_id in Orthology_all[orths][ref_org]:
                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id]={}
                    if orth_id!='ENSMUSG00000032718':
                        continue
                     # These are version of the saem gene (i.e. alleles)
                    for seq_id in Orthology_all[orths][orga][orth_id]:
                        if seq_id=='name':
                            continue
                        seq_raw=Orthology_all[orths][orga][orth_id][seq_id]
                        if clean_sequence:
                            seq_raw=clean_seq(seq_raw)
                        for ref_id in Orthology_all[orths][ref_org][orth_ref_id]:
                            if ref_id=='name':
                                continue
                            seq_ref_raw=Orthology_all[orths][ref_org][orth_ref_id][ref_id]
                            if clean_sequence:
                                seq_ref_raw=clean_seq(seq_ref_raw)

                            if min(len(seq_ref_raw),len(seq_raw))/max(len(seq_ref_raw),len(seq_raw))<length_ratio :
                                continue

                            score,norm=OU.get_homology_score(seq_raw,seq_ref_raw,local=local)
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]={}
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["all"]={}
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]={}
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]['Homology']=score
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]['Homology_ratio']=norm
                            for val_type in val_types_seq:
                                val_ref=SeqProp.get_values(ref_id,"all",str(0)+'_'+str(len(seq_ref_raw)),val_type)
                                val_oth=SeqProp.get_values(seq_id,"all",str(0)+'_'+str(len(seq_raw)),val_type)
                                Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))][val_type]=str(val_oth-val_ref)

                            bounds,bounds_labels=SeqProp.get_bounds(seq_id,"FDs")
                            bounds_ref,bounds_labels_ref=SeqProp.get_bounds(ref_id,"FDs")

                            bound_match,match_score,bounds_not_folded,bounds_not_folded_ref=OU.find_matching_folded_domains(bounds,bounds_ref,seq_raw,seq_ref_raw)
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["FDs"]={}
                            for m in range(len(bound_match)):
                                Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["FDs"][bound_match[m]]={}
                                Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["FDs"][bound_match[m]]['Homology']=str(match_score[m][0])
                                Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["FDs"][bound_match[m]]['Homology_ratio']=str(match_score[m][1])
                                for val_type in val_types_seq:
                                    val_ref=SeqProp.get_values(ref_id,"FDs",bound_match[m].split('_&_')[0],val_type)
                                    val_oth=SeqProp.get_values(seq_id,"FDs",bound_match[m].split('_&_')[1],val_type)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["FDs"][bound_match[m]][val_type]=str(val_oth-val_ref)

                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["IDRs"]={}
                            Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["NFDs"]={}

                            bounds_dis_ref,bounds_labels_dis_ref=SeqProp.get_bounds(ref_id,"IDRs")
                            bounds_dis,bounds_labels_dis=SeqProp.get_bounds(seq_id,"IDRs")

                            for m in range(len(bounds_not_folded)):

                                seq_ref=seq_ref_raw[bounds_not_folded_ref[m,0]:bounds_not_folded_ref[m,1]]
                                seq_oth=seq_raw[bounds_not_folded[m,0]:bounds_not_folded[m,1]]

                                score,norm=OU.get_homology_score(seq_ref,seq_oth,local)
                                bounds_label=str(bounds_not_folded_ref[m,0])+'_'+str(bounds_not_folded_ref[m,1])+'_&_'+str(bounds_not_folded[m,0])+'_'+str(bounds_not_folded[m,1])

                                Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["NFDs"][bounds_label]={}
                                Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["NFDs"][bounds_label]['Homology']=str(score)
                                Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id][ "NFDs"][bounds_label]['Homology_ratio']=str(norm)
                                if bounds_label.split('_&_')[0] in bounds_labels_dis_ref and bounds_label.split('_&_')[1] in bounds_labels_dis:
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["IDRs"][bounds_label]={}
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["IDRs"][bounds_label]['Homology']=str(score)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["IDRs"][bounds_label]['Homology_ratio']=str(norm)
                                    for val_type in val_types_seq:
                                        val_ref=SeqProp.get_values(ref_id,"IDRs",bounds_label.split('_&_')[0],val_type)
                                        val_oth=SeqProp.get_values(seq_id,"IDRs",bounds_label.split('_&_')[1],val_type)
                                        Sequence_homology_all[ref_org][orga][orth_ref_id+'_'+orth_id][ref_id+'_'+seq_id]["IDRs"][bounds_label][val_type]=str(val_oth-val_ref)
                            # Here I should only compare IDRS that correspond to NFDs
                            # This means for every couple of NFDs, check that they each correspond to an IDR
                            # If they do, save the ratio fo the ensemble predictions into the homology json
                            # From here on, we will recombine the IDRs
    with open(output_file,'w') as f:
        json.dump(Sequence_homology_all,f)
