import skbio
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
        # 6. Sequence Physics Properties (FU)
        fn,fp,NCPR,FCR,kappa=FU.get_seq_properties(seq,window_size=5)

        self.SeqProp[seq_id][label][bounds_label]['NCPR']=str(NCPR)
        self.SeqProp[seq_id][label][bounds_label]['FCR']=str(FCR)
        self.SeqProp[seq_id][label][bounds_label]['f+']=str(fp)
        self.SeqProp[seq_id][label][bounds_label]['f-']=str(fn)
        self.SeqProp[seq_id][label][bounds_label]['Kappa']=str(kappa)

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

    # parser.add_argument("--original_specie","-os",help="Original specie")
    parser.add_argument("--reference_specie","-rs",help="Reference specie")
    parser.add_argument("--additional_species","-as",help="Additional species",default=[],nargs='+')
    parser.add_argument("--max_size_factor","-msf",help="Factor for the length of the otrholog array compared to the input gene list size. Must be integer, bigger number is slower, but if many orhtolog exists, may be necessary")
    parser.add_argument("--delete_cross_refs","-dcr",help="Whether gene name that share an ortholog should be delete, so they don't appear twice.  1 is delete, 0 is keep.")
    args = parser.parse_args()

    if args.file_orthologs :
        filename=args.file_orthologs
    else :
        filename='Gene_orthology.json'

    if args.file_IDRs :
        filename_IDRs=args.file_IDRs
    else :
        filename_IDRs='Proteins_IDR_properties.json'


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

    f=open(filename_IDRs)
    #IDRs=json.load(f)
    SeqProp=Seq_Prop_Manager(filename_IDRs)

    # MEDOC part
    pH_range=[0,14]
    pH=np.arange(pH_range[0],pH_range[1]+0.01,0.01,dtype=np.longdouble)
    QC.initialize_additive_DF_array_MEDOC_public()
    pH_ref=[7.4]



    do_seq_prop=True
    do_seq_pH=True
    local=True
    N1=0
    cross_check=False

    labels=["all","IDRs","FDs"]

    for orths in Orthology_all:
        print(100*N1/len(Orthology_all))
        N1+=1
        # I should make an orthology score matrix here
        # I need a reference organism, which will carry the scores for the others
        # They will be in order
        # I need an all organisms AND orthologs AND all versions AND all IDRs comparison (Yeah that's a lot)
        # For the adding sequence properties, I do not need the loop on the ref and compa, I will just do that at the
        # same time as the homology
        for orga in Orthology_all[orths]:
            #  dimension for each orga
            for N_orth in Orthology_all[orths][orga]:
                # dimension for each ortolog
                orth_id=Orthology_all[orths][orga][N_orth]['gene_id']
                # These are version of the saem gene (i.e. alleles)
                for vers in Orthology_all[orths][orga][N_orth]['version']:
                    seq_id=Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']
                    seq_raw=Orthology_all[orths][orga][N_orth]['version'][vers]['seq']
                    if clean_sequence :
                        seq_raw=clean_seq(seq_raw)
                    for label in labels:
                        bounds,bounds_labels=SeqProp.get_bounds(seq_id,label)
                        for i in range(len(bounds)):
                            if do_seq_pH:
                                SeqProp.add_pH_information(seq_id,seq_raw,pH_ref,pH,label,bounds_labels[i])
                            if do_seq_prop:
                                SeqProp.add_sequence_properties(seq_id,seq_raw,label,bounds_labels[i])

    SeqProp.save()

