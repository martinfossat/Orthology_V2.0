import json
import metapredict as MP
import argparse
import os
import Orthology_utils as OU

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Gets the predicted IDR regions, and predict their ensemble properties""")
    parser.add_argument("--orthology_file","-of",help='Name of the orthology file')
    parser.add_argument("--properties_file","-pf",help='Name of the sequence properties file. Default is Sequence_properties.json')
    parser.add_argument("--sequences_file", "-sf",help='Name of the sequences database json output file. Default is Sequences.json')
    # parser.add_argument("--species","-as",help="Additional species",default=[],nargs='+')
    args=parser.parse_args()

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

    file_path=os.path.realpath(__file__)
    with open(os.path.dirname(file_path)+"/g_Profiler_Organisms_names_dic.json", "r") as fp:
        dictionary_organisms_gprofiler = json.load(fp)
    dictionary_organisms_gprofiler_inv=dict(zip(dictionary_organisms_gprofiler.values(), dictionary_organisms_gprofiler.keys()))

    # Opening JSON file
    f=open(orthology_file)
    Orthology_all=json.load(f)

    f=open(sequences_file)
    Sequences_all=json.load(f)

    Seq_Prop={}
    N1=0
    for orths in Orthology_all:
        N1+=1
        for orga in Orthology_all[orths]:
            for gene_id in Orthology_all[orths][orga]:
                for seq_id in Orthology_all[orths][orga][gene_id]:
                    seq=OU.clean_seq(Sequences_all[seq_id])

                    Seq_Prop[seq_id]={}
                    Seq_Prop[seq_id]['all']={}
                    Seq_Prop[seq_id]['all'][str(0)+'_'+str(len(seq))]={}
                    Seq_Prop[seq_id]['IDRs']={}

                    dis_bounds=MP.predict_disorder_domains(seq).disordered_domain_boundaries
                    fd_bounds=OU.get_bounds_inverted(dis_bounds,seq)
                    for bounds in dis_bounds:
                        label_bounds=str(bounds[0])+'_'+str(bounds[1])
                        Seq_Prop[seq_id]['IDRs'][label_bounds]={}

                    Seq_Prop[seq_id]['FDs']={}
                    for bounds in fd_bounds:
                        label_bounds=str(bounds[0])+'_'+str(bounds[1])
                        Seq_Prop[seq_id]['FDs'][label_bounds]={}

        with open(properties_file,'w') as f:
            json.dump(Seq_Prop,f)