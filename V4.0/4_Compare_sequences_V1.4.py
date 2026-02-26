import matplotlib
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
import json
import Orthology_utils as OU
import os
import numpy as np
import math

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Compares sequences of orthologs between a reference species and a number of other species. 
    Comparisons includes comparing sequence features, sequence ensemble prediction and charge sequence features, as well as calcualting the homology between sequences.""")
    parser.add_argument("--orthology_file","-of",help='File name where the file is a json file containing created by the Ortholog program')
    parser.add_argument("--sequences_file", "-sf",help='Name of the sequences database json output file. Default is Sequences.json')
    parser.add_argument("--properties_file", "-pf",help='Name of the sequence properties file. Default is Sequence_properties.json')
    parser.add_argument("--homologies_file","-hf",help='Name of the homology database file. Default is Gene_homologies.json')
    parser.add_argument("--min_len_fraction","-mlf",help='Minimum length fraction. Best to keep 0, build the database for everything, and discard later, however, if the algorithm runs too slowly, this is a simple way of speeding it up.')
    parser.add_argument("--reference_specie","-rs",help="Name of the reference specie to which sequences are compared")
    parser.add_argument("--additional_species","-as",help="Name of additional species which are compared to the reference specie.",default=[],nargs='+')
    parser.add_argument("--do_homo","-dh",help="Whether to calculate the homology score. 1 is on, 0 is off. Default is 1 (on)")
    parser.add_argument("--align_type","-at",help="Alignment type for the homology. Can be NW (Needleman–Wunsch, i.e. global) or SW (Smith–Waterman, i.e. local). Default is SW.")

    ############### NEW
    parser.add_argument("--seq_labels","-sl",help="Labels for the sequence annotation. Must correspond to name given in other programs. Default is \"all\", which is the whole sequence. Can be multiple arguments.")
    parser.add_argument("--iso_compare","-ic",help="How to compare isoforms. Option are : MLF (Minimum Length Fraction) Compares all protein isoform that have a pair length ratio above the value specified; MLO (Maximum Length Only) : only compare the longest isoforms between the two species. This is what is done in Ensembl. Default is MLO.")
    args=parser.parse_args()
    # Last thing to do :
    # Make sure that the choosing of longest isoform for both species is possible at the plotting point.
    # Implement comparison based on species gene list
    # Transfer name of species from gProfiler to Ensembl id


    if args.iso_compare :
        if args.iso_compare=='MLO':
            use_MLF=False
        elif args.iso_compare=='MLF':
            use_MLF=True
        else :
            print("Invalid argument for align_type")
            quit()
    else :
        use_MLF=False

    if args.align_type :
        if args.align_type=='SW':
            local=True
        elif args.align_type=='NW':
            local=False
        else :
            print("Invalid argument for align_type")
            quit()
    else :
        local=True

    if args.do_homo :
        try :
            if int(args.do_homo)==1 :
                do_homo=True
            elif int(args.do_homo)==0 :
                do_homo=False
            else :
                print("Wrong value for do_homo. Can only be 1 or 0")
        except :
            print("Wrong format for do_homo")
            quit()
    else :
        do_homo=True

    if args.min_len_fraction :
        try :
            length_ratio=float(args.min_len_fraction)
        except :
            print("Wrong format for mlf")
        if length_ratio<0 or length_ratio>1.:
            print("Wrong value for mlf. Please use a number between 0 and 1")
    else :
        length_ratio=0.0
        print("No argument for mlf, defaulting to "+str(length_ratio))

    if args.orthology_file :
        orthology_file=args.orthology_file
    else :
        orthology_file='Orthology.json'

    if args.properties_file :
        properties_file=args.properties_file
    else :
        properties_file='Sequence_properties.json'

    if args.sequences_file :
        sequences_file=args.sequences_file
    else :
        sequences_file='Sequences.json'

    if args.homologies_file :
        homologies_file=args.homologies_file
    else :
        homologies_file='Homology.json'

    if args.reference_specie :
        ref_org=args.reference_specie
    else :
        print("You must specify a reference specie")
        quit(1)

    if args.additional_species :
        other_organism=args.additional_species
    else :

        print("You need at least two ogranism to compare")
        quit()

    if args.seq_labels :
        seq_labels=args.seq_labels
    else :
        seq_labels=["all"]

    clean_sequence=True

    all_orga=[ref_org]+other_organism

    print('Loading orthology file')
    f=open(orthology_file)
    Orthology_all=json.load(f)

    print('Loading sequences file')
    f=open(sequences_file)
    Sequences_all=json.load(f)

    print('Loading sequence properties file')
    SeqProp=OU.Seq_Prop_Manager(properties_file)
    print(SeqProp)
    print('Loading homology file')
    if os.path.exists(homologies_file):
        # I should load the existing file in case one wants to do a difference comparison
        f=open(homologies_file)
        Sequence_homology_all=json.load(f)
    else :
        Sequence_homology_all={}

    if not ref_org in Sequence_homology_all :
        Sequence_homology_all[ref_org]={}

    val_types_IDRs=['rg','re','asph','nu','pref']
    val_types_charge=['pI']
    val_types_seq=['NCPR', 'FCR','f-','f+','Kappa']

    for orths in Orthology_all:
        # I should make an orthology score matrix here
        # I need a reference organism, which will carry the scores for the others
        # They will be in order
        # I need an all organisms AND orthologs AND all versions AND all Seq_Prop comparison (Yeah that's a lot)
        if not ref_org in (Orthology_all[orths].keys()):
            continue
        for orga in Orthology_all[orths]:
            if not orga in Sequence_homology_all[ref_org]:
                Sequence_homology_all[ref_org][orga]={}
            #  Dimension for each orga
            for orth_ref_id in Orthology_all[orths][ref_org]:
                if orth_ref_id=='N/A':
                    continue
                if not orth_ref_id in Sequence_homology_all[ref_org][orga].keys():
                    Sequence_homology_all[ref_org][orga][orth_ref_id]={}
                # dimension for each ortholog
                for orth_id in Orthology_all[orths][orga]:
                    if orth_id=='N/A':
                        continue
                    if not orth_id in Sequence_homology_all[ref_org][orga][orth_ref_id].keys():
                        Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id]={}
                    # Change of strategy : we will now list all of the sequences, then go over construct a pairs list
                    seq_list_ref=[]
                    seq_list_ref_len=[]

                    for ref_id in Orthology_all[orths][ref_org][orth_ref_id]:

                        seq_ref_raw=Sequences_all[ref_id]

                        if clean_sequence:
                            seq_ref_raw=OU.clean_seq(seq_ref_raw)

                        seq_list_ref+=[ref_id]
                        seq_list_ref_len+=[len(seq_ref_raw)]

                    seq_list_ref=np.array(seq_list_ref)
                    seq_list_ref_len=np.array(seq_list_ref_len)

                    seq_list_oth=[]
                    seq_list_oth_len=[]
                    for seq_id in Orthology_all[orths][orga][orth_id]:
                        seq_raw=Sequences_all[seq_id]
                        if clean_sequence:
                            seq_raw=OU.clean_seq(seq_raw)
                        seq_list_oth+=[seq_id]
                        seq_list_oth_len+=[len(seq_raw)]

                    seq_list_oth=np.array(seq_list_oth)
                    seq_list_oth_len=np.array(seq_list_oth_len)

                    if use_MLF :
                        list_pairs=[[],[]]
                        for s1 in range(len(seq_list_ref)):
                            for s2 in range(len(seq_list_oth)):
                                if  min(seq_list_ref_len[s1],seq_list_oth_len[s2])/max(seq_list_ref_len[s1],seq_list_oth_len[s2])<length_ratio:
                                    continue
                                list_pairs[0]+=[seq_list_ref[s1]]
                                list_pairs[1]+=[seq_list_oth[s2]]

                    else :
                        s1=np.argmax(seq_list_ref_len)
                        s2=np.argmax(seq_list_oth_len)
                        list_pairs=[[seq_list_ref[s1]],[seq_list_oth[s2]]]

                    for sp in range(len(list_pairs[0])):
                        ref_id=list_pairs[0][sp]
                        seq_id=list_pairs[1][sp]

                        seq_ref_raw=Sequences_all[ref_id]
                        seq_raw=Sequences_all[seq_id]

                        if "all" in seq_labels :
                            if not ref_id in Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id]:
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id]={}
                            if not seq_id in Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id]:
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]={}
                            if not 'all' in Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]:
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"]={}


                            bounds_label=str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))

                            if not bounds_label in Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"]:
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][bounds_label]={}

                            if do_homo:
                                if not 'Homology' in Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]:
                                    score,norm=OU.get_homology_score(seq_raw,seq_ref_raw,local=local)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]['Homology']=str(score)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]['Homology_ratio']=str(norm)

                            for val_type in val_types_seq:
                                val_ref=SeqProp.get_values(ref_id,"all",str(0)+'_'+str(len(seq_ref_raw)),val_type)
                                val_oth=SeqProp.get_values(seq_id,"all",str(0)+'_'+str(len(seq_raw)),val_type)
                                if not math.isnan(val_ref) and not math.isnan(val_oth):
                                    Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))][val_type]=str(val_oth-val_ref)

                        if 'FDs' in seq_labels:
                            bounds,bounds_labels=SeqProp.get_bounds(seq_id,"FDs")
                            bounds_ref,bounds_labels_ref=SeqProp.get_bounds(ref_id,"FDs")

                            bound_match,match_score,bounds_not_folded,bounds_not_folded_ref=OU.find_matching_folded_domains(bounds,bounds_ref,seq_raw,seq_ref_raw)
                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["FDs"]={}
                            for m in range(len(bound_match)):
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["FDs"][bound_match[m]]={}
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["FDs"][bound_match[m]]['Homology']=str(match_score[m][0])
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["FDs"][bound_match[m]]['Homology_ratio']=str(match_score[m][1])
                                for val_type in val_types_seq:
                                    val_ref=SeqProp.get_values(ref_id,"FDs",bound_match[m].split('_&_')[0],val_type)
                                    val_oth=SeqProp.get_values(seq_id,"FDs",bound_match[m].split('_&_')[1],val_type)
                                    if not math.isnan(val_ref) and not math.isnan(val_oth):
                                        Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["FDs"][bound_match[m]][val_type]=str(val_oth-val_ref)

                        if "IDRs" in seq_labels or "NFDs" in seq_labels :
                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["IDRs"]={}
                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["NFDs"]={}

                            bounds_dis_ref,bounds_labels_dis_ref=SeqProp.get_bounds(ref_id,"IDRs")
                            bounds_dis,bounds_labels_dis=SeqProp.get_bounds(seq_id,"IDRs")

                            for m in range(len(bounds_not_folded)):
                                seq_ref=seq_ref_raw[bounds_not_folded_ref[m,0]:bounds_not_folded_ref[m,1]]
                                seq_oth=seq_raw[bounds_not_folded[m,0]:bounds_not_folded[m,1]]

                                score,norm=OU.get_homology_score(seq_ref,seq_oth,local)
                                bounds_label=str(bounds_not_folded_ref[m,0])+'_'+str(bounds_not_folded_ref[m,1])+'_&_'+str(bounds_not_folded[m,0])+'_'+str(bounds_not_folded[m,1])

                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["NFDs"][bounds_label]={}
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["NFDs"][bounds_label]['Homology']=str(score)
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id][ "NFDs"][bounds_label]['Homology_ratio']=str(norm)
                                if bounds_label.split('_&_')[0] in bounds_labels_dis_ref and bounds_label.split('_&_')[1] in bounds_labels_dis:
                                    Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["IDRs"][bounds_label]={}
                                    Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["IDRs"][bounds_label]['Homology']=str(score)
                                    Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["IDRs"][bounds_label]['Homology_ratio']=str(norm)
                                    for val_type in val_types_seq:
                                        val_ref=SeqProp.get_values(ref_id,"IDRs",bounds_label.split('_&_')[0],val_type)
                                        val_oth=SeqProp.get_values(seq_id,"IDRs",bounds_label.split('_&_')[1],val_type)
                                        if not math.isnan(val_ref) and not math.isnan(val_oth):
                                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["IDRs"][bounds_label][val_type]=str(val_oth-val_ref)

    with open(homologies_file,'w') as f:
        json.dump(Sequence_homology_all,f)
