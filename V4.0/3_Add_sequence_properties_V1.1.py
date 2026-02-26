import json
import argparse
import numpy as np
import Orthology_utils as OU

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Adds sequence properties to the property file. Sequence properties can be : (1) sequence features 
    (fraction of charge residue (FCR), fraction of positive and negative residues, Net Charge Per Residue from primary sequence (Does not account for pH)) 
    (2) Sequence ensemble properties prediction from Sparrow (only strictly valid for IDRs), which include end to end distance, radius of gyration, asphericity, 
    ect... or (3) pH specific charge properties prediction based on MEDOC (Isoelectric point and Net Charge, although mostly valid for disordered regions, is always done for everything).""")
    parser.add_argument("--orthology_file","-of",help='Name of the orthology database json output file. Default is Orthology.json')
    parser.add_argument("--sequences_file","-sf",help='Name of the sequences database json output file. Default is Sequences.json')
    parser.add_argument("--properties_file","-pf",help='Name of the sequence properties file. Default is Sequence_properties.json')
    parser.add_argument("--do_seq_feat","-dsf",help='Whether to analyze sequence features. 0 if off, 1 is on. Default is 1 (on)')
    parser.add_argument("--do_seq_ens","-dse",help='Whether to predict sequence ensemble properties. 0 if off, 1 is on. Default is 1 (on)')
    parser.add_argument("--do_seq_pH","-dsp",help='Whether to perform pH calculations. 0 if off, 1 is on. Default is 1 (on)')
    parser.add_argument("--force_pred","-fp",help='Whether to force calculations on every ensemble types. By default, the ensemble prediction are limited to all non-folded domains. 0 if off, 1 is on. Default is 0 (off)')
    parser.add_argument("--pH_val","-pHv",help='Values at which to calculate the effective charge. Can have multiple arguments',default=[],nargs='+')
    args = parser.parse_args()

    if args.do_seq_pH:
            if args.do_seq_pH=='0':
                do_seq_pH=False
            elif args.do_seq_pH=='1':
                do_seq_pH=True
            else :
                print("Wrong value for do_seq_pH. Please use 0 and 1.")
                quit()
    else:
        do_seq_pH=True

    if args.force_pred:
        if args.force_pred=='0':
            force_pred=False
        elif args.force_pred=='1':
            force_pred=True
        else:
            print("Wrong value for force_pred. Please use 0 and 1.")
            quit()
    else:
        force_pred=True

    if args.pH_val:
        try :
            pH_ref=[]
            for i in range(len(args.pH_val)):
                pH_ref+=[float(args.pH_val[i])]
        except :
            print("Wrong value(s) for pH_val.")

        if do_seq_pH and  len(pH_ref)!=0 :
            print("Incompatible arguments : you must have do_seq_pH turned on if you want to compute the charge.")
            quit()

    else:
        pH_ref=[]

    if args.do_seq_feat:
            if args.do_seq_feat=='0':
                do_seq_feat=False
            elif args.do_seq_feat=='1':
                do_seq_feat=True
            else :
                print("Wrong value for do_seq_feat. Please use 0 and 1.")
                quit()
    else:
        do_seq_feat=True

    if args.do_seq_ens:
            if args.do_seq_ens=='0':
                do_seq_ens=False
            elif args.do_seq_ens=='1':
                do_seq_ens=True
            else :
                print("Wrong value for do_seq_ens. Please use 0 and 1.")
                quit()
    else:
        do_seq_ens=True

    if args.orthology_file :
        filename=args.orthology_file
    else :
        filename='Orthology.json'

    if args.properties_file :
        properties_file=args.properties_file
    else :
        properties_file='Sequence_properties.json'

    if args.sequences_file :
        sequences_file=args.sequences_file
    else :
        sequences_file='Sequences.json'

    clean_sequence=True

    f=open(filename)
    Orthology_all=json.load(f)

    SeqProp=OU.Seq_Prop_Manager(properties_file)

    f=open(sequences_file)
    Sequences_all=json.load(f)

    if do_seq_pH:
        # MEDOC part
        pH_range=[0,14]
        pH=np.arange(pH_range[0],pH_range[1]+0.01,0.01,dtype=np.longdouble)
        OU.initialize_additive_DF_array_MEDOC_public()

    force_pred=False
    labels=["all","IDRs","FDs","NFDs"]

    for orths in Orthology_all:
        # I should make an orthology score matrix here
        # I need a reference organism, which will carry the scores for the others
        # They will be in order
        # I need an all organisms AND orthologs AND all versions AND all IDRs comparison (Yeah that's a lot)
        # For the adding sequence properties, I do not need the loop on the ref and compa, I will just do that at the
        # same time as the homology
        for orga in Orthology_all[orths]:
            #  dimension for each orga
            for orth_id in Orthology_all[orths][orga]:
                # These are version of the saem gene (i.e. isoforms)
                for seq_id in Orthology_all[orths][orga][orth_id]:
                    seq_raw=Sequences_all[seq_id]
                    if clean_sequence :
                        seq_raw=OU.clean_seq(seq_raw)
                    for label in labels:
                        bounds,bounds_labels=SeqProp.get_bounds(seq_id,label)
                        if bounds is None:
                            continue
                        for i in range(len(bounds)):
                            if bounds[i,0]==bounds[i,1]:
                                continue
                            if do_seq_pH:
                                SeqProp.add_pH_information(seq_id,seq_raw,pH_ref,pH,label,bounds_labels[i])

                            if do_seq_feat:
                                SeqProp.add_sequence_properties(seq_id,seq_raw,label,bounds_labels[i])

                            if do_seq_ens:
                                if force_pred :
                                    SeqProp.add_ensemble_properties(seq_id,seq_raw,label,bounds_labels[i])
                                elif label!="IDRs"  and label!="NFDs" :
                                    SeqProp.add_ensemble_properties(seq_id,seq_raw,label,bounds_labels[i])

                    bounds_fd=np.array(SeqProp.get_bounds(seq_id,"FDs")[0])
                    SeqProp.add_percent_folded(seq_id,seq_raw,bounds_fd)
    SeqProp.save()
