import json
import metapredict as MP
import argparse
import Orthology_utils as OU

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Gets the predicted IDR regions, and predict their ensemble properties""")
    parser.add_argument("--properties_file","-pf",help='Name of the sequence properties file. Default is Sequence_properties.json')
    parser.add_argument("--sequences_file", "-sf",help='Name of the sequences database json output file. Default is Sequences.json')
    args=parser.parse_args()

    if args.properties_file :
        properties_file=args.properties_file
    else :
        properties_file='Sequence_properties.json'

    if args.sequences_file :
        sequences_file=args.sequences_file
    else :
        sequences_file='Sequences.json'

    f=open(sequences_file)
    Sequences_all=json.load(f)

    Seq_Prop={}
    N1=0

    for seq_id in Sequences_all:
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