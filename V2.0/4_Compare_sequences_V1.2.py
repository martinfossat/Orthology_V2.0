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
import Orthology_utils as OU

AA_type,AA_scores=FU.get_self_homology_score()
if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Checks the orthologs of a gene name list.\n
    The user needs to specify an original species name, corresponding to the gene name list, and a reference specie, which may be used for 
    homology comparison, in other steps. Additional species may be given, but those two first are required.\n
    Species name must follow the id in gProfiler (https://biit.cs.ut.ee/gprofiler/page/organism-list).\n
    The output is a file containing the gene name and gene IDs for all species, and that is required to use subsequent programs.""")
    parser.add_argument("--orthology_file","-of",help='File name where the file is a json file containing created by the Ortholog program')
    parser.add_argument("--sequences_file", "-sf",help='Name of the sequences database json output file. Default is Sequences.json')
    parser.add_argument("--properties_file", "-pf",help='Name of the sequence properties file. Default is Sequence_properties.json')
    parser.add_argument("--homologies_file","-hf",help='Name of the homology database file. Default is Gene_homologies.json')
    parser.add_argument("--min_len_fraction","-mif",help='Minimum length fraction. Best to keep 0, build the database for everything, and discard later')
    parser.add_argument("--reference_specie","-rs",help="Reference specie")
    parser.add_argument("--additional_species","-as",help="Additional species",default=[],nargs='+')
    #parser.add_argument("--max_size_factor","-msf",help="Factor for the length of the ortholog array compared to the input gene list size. Must be integer, bigger number is slower, but if many ortholog exists, may be necessary")
    #parser.add_argument("--delete_cross_refs","-dcr",help="Whether gene name that share an ortholog should be delete, so they don't appear twice.  1 is delete, 0 is keep.")
    args = parser.parse_args()

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
        homologies_file='Gene_homologies.json'

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

    f=open(orthology_file)
    Orthology_all=json.load(f)

    f=open(sequences_file)
    Sequences_all=json.load(f)
    local=True
    SeqProp=OU.Seq_Prop_Manager(properties_file)
    Sequence_homology_all={}

    Sequence_homology_all[ref_org]={oth:{} for oth in all_orga}
    val_types_IDRs=['rg','re','asph','nu','pref']
    val_types_charge=['pI']
    val_types_seq=['NCPR', 'FCR','f-','f+','Kappa']

    # labels=["all","IDRs","FDs"]

    for orths in Orthology_all:
        # I should make an orthology score matrix here
        # I need a reference organism, which will carry the scores for the others
        # They will be in order
        # I need an all organisms AND orthologs AND all versions AND all Seq_Prop comparison (Yeah that's a lot)
        if not ref_org in (Orthology_all[orths].keys()):
            continue

        for orga in Orthology_all[orths]:
            print(orga)
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
                     # These are version of the same gene (i.e. isoform)
                    for ref_id in Orthology_all[orths][ref_org][orth_ref_id]:
                        if ref_id=='name':
                            continue
                        if not ref_id in Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id].keys():
                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id]={}
                    
                        seq_ref_raw=Sequences_all[ref_id]
                        if clean_sequence:
                            seq_ref_raw=OU.clean_seq(seq_ref_raw)
                        
                        for seq_id in Orthology_all[orths][orga][orth_id]:
                            if seq_id=='name':
                                continue
                            if not seq_id in Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id].keys():
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]={}
                            seq_raw=Sequences_all[seq_id]
                            if clean_sequence:
                                seq_raw=OU.clean_seq(seq_raw)

                            if min(len(seq_ref_raw),len(seq_raw))/max(len(seq_ref_raw),len(seq_raw))<length_ratio :
                                continue

                            score,norm=OU.get_homology_score(seq_raw,seq_ref_raw,local=local)
                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"]={}
                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]={}
                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]['Homology']=str(score)
                            Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))]['Homology_ratio']=str(norm)

                            for val_type in val_types_seq:
                                val_ref=SeqProp.get_values(ref_id,"all",str(0)+'_'+str(len(seq_ref_raw)),val_type)
                                val_oth=SeqProp.get_values(seq_id,"all",str(0)+'_'+str(len(seq_raw)),val_type)
                                Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["all"][str(0)+'_'+str(len(seq_ref_raw))+'_&_'+str(0)+'_'+str(len(seq_raw))][val_type]=str(val_oth-val_ref)

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
                                    Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["FDs"][bound_match[m]][val_type]=str(val_oth-val_ref)

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
                                        Sequence_homology_all[ref_org][orga][orth_ref_id][orth_id][ref_id][seq_id]["IDRs"][bounds_label][val_type]=str(val_oth-val_ref)
                            # Here I should only compare IDRs that correspond to NFDs
                            # This means for every couple of NFDs, check that they each correspond to an IDR
                            # If they do, save the ratio fo the ensemble predictions into the homology json
                            # From here on, we will recombine the IDRs

    with open(homologies_file,'w') as f:
        json.dump(Sequence_homology_all,f)
