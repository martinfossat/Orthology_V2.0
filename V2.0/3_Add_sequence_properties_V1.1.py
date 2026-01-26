
from sparrow.predictors import batch_predict
import json
import matplotlib
import Fossat_utils as FU
import qcan_utils as QC
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

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Checks the orthologs of a gene name list.\n
    The user needs to specify an original species name, corresponding to the gene name list, and a reference specie, which may be used for 
    homology comparison, in other steps. Additional species may be given, but those two first are required.\n
    Species name must follow the id in gProfiler (https://biit.cs.ut.ee/gprofiler/page/organism-list).\n
    The output is a file containing the gene name and gene IDs for all species, and that is required to use subsequent programs.""")
    parser.add_argument("--file_orthologs","-fo",help='File name where the file is a json file containing created by the Ortholog program')
    parser.add_argument("--file_IDRs","-fi",help='File name where the file is a json file containing created by the Get_IDR program, contaning IDR domain boundaries')
    parser.add_argument("--max_size_factor","-msf",help="Factor for the length of the ortholog array compared to the input gene list size. Must be integer, bigger number is slower, but if many orhtolog exists, may be necessary")
    parser.add_argument("--delete_cross_refs","-dcr",help="Whether gene name that share an ortholog should be delete, so they don't appear twice.  1 is delete, 0 is keep.")
    args = parser.parse_args()

    if args.file_orthologs :
        filename=args.file_orthologs
    else :
        filename='Gene_orthology.json'

    if args.file_IDRs :
        filename_IDRs=args.file_IDRs
    else :
        filename_IDRs='Sequence_properties.json'

    file_path=os.path.realpath(__file__)
    with open(os.path.dirname(file_path)+"/g_Profiler_Organisms_names_dic.json", "r") as fp:
        dictionary_organisms_gprofiler = json.load(fp)
    clean_sequence=True
    dictionary_organisms_gprofiler_inv=dict(zip(dictionary_organisms_gprofiler.values(), dictionary_organisms_gprofiler.keys()))

    f=open(filename)
    Orthology_all=json.load(f)

    f=open(filename_IDRs)
    SeqProp=OU.Seq_Prop_Manager(filename_IDRs)

    # MEDOC part
    pH_range=[0,14]
    pH=np.arange(pH_range[0],pH_range[1]+0.01,0.01,dtype=np.longdouble)
    QC.initialize_additive_DF_array_MEDOC_public()
    pH_ref=[7.4]

    do_seq_prop=True
    do_seq_pH=True
    do_ensemble_prop=True
    force_prediction=True

    local=True


    labels=["all","IDRs","FDs"]

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
                # These are version of the saem gene (i.e. alleles)
                for seq_id in Orthology_all[orths][orga][orth_id]:
                    if seq_id=='name':
                        continue
                    seq_raw=Orthology_all[orths][orga][orth_id][seq_id]
                    if clean_sequence :
                        seq_raw=clean_seq(seq_raw)
                    for label in labels:
                        bounds,bounds_labels=SeqProp.get_bounds(seq_id,label)
                        for i in range(len(bounds)):
                            if do_seq_pH:
                                if force_prediction and label!="IDRs":
                                    SeqProp.add_pH_information(seq_id,seq_raw,pH_ref,pH,label,bounds_labels[i])
                                else :
                                    SeqProp.add_pH_information(seq_id,seq_raw,pH_ref,pH,label,bounds_labels[i])
                            if do_seq_prop:
                                SeqProp.add_sequence_properties(seq_id,seq_raw,label,bounds_labels[i])
                            if do_ensemble_prop:
                                if force_prediction and label!="IDRs":
                                    SeqProp.add_ensemble_properties(seq_id,seq_raw,label,bounds_labels[i])
                                else :
                                    SeqProp.add_ensemble_properties(seq_id,seq_raw,label,bounds_labels[i])
                    bounds_fd=np.array(SeqProp.get_bounds(seq_id,"FDs")[0])
                    SeqProp.add_percent_folded(seq_id,seq_raw,bounds_fd)
    SeqProp.save()
