import json
import metapredict as MP
from sparrow.predictors import batch_predict
import argparse
import os
import Orthology_utils as OU


if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Gets the predicted IDR regions, and predict their ensemble properties""")
    parser.add_argument("--file_name","-f",help='File name containing the name of the genes in the "original" species in the first column, other columns are ignored.')
    parser.add_argument("--output_file_name","-of",help='Name of the outputfile. Default is Gene_names.txt')
    parser.add_argument("--original_specie","-os",help="Original specie")
    parser.add_argument("--additional_species","-as",help="Additional species",default=[],nargs='+')
    parser.add_argument("--delete_x","-dx",help="Replace symbol X in proteins sequences by an emtpy character")
    parser.add_argument("--delete_u","-du",help="Replace symbol U in proteins sequences by an emtpy character")
    parser.add_argument("--delete_any","-da",help="Replace symbol * in proteins sequences by an emtpy character")
    parser.add_argument("--max_size_factor","-msf",help="Factor for the length of the otrholog array compared to the input gene list size. Must be integer, bigger number is slower, but if many orhtolog exists, may be necessary")
    args = parser.parse_args()

    if args.file_name :
        filename=args.file_name
    else :
        print("You must specify an input file name")
        quit(1)

    if args.output_file_name :
        output_file=args.output_file_name
    else :
        output_file='Sequence_properties.json'

    if args.delete_x :
        try :
            X_is_None=bool(args.delete_x)
        except :
            print("Invalid option for replacing X")
            X_is_None=True
    else :
        X_is_None=True

    if args.delete_u :
        try :
            U_is_None=bool(args.delete_u)
        except :
            print("Invalid option for replacing X")
            U_is_None=True
    else :
        U_is_None=True

    if args.delete_any :
        try :
            Any_is_None=bool(args.delete_any)
        except :
            print("Invalid option for replacing X")
            Any_is_None=True
    else :
        Any_is_None=True

    # X_is_None=True
    # U_is_None=True
    # Any_is_None=True
    # Whether to delete the X (undefined AAs)

    file_path=os.path.realpath(__file__)
    with open(os.path.dirname(file_path)+"/g_Profiler_Organisms_names_dic.json", "r") as fp:
        dictionary_organisms_gprofiler = json.load(fp)
    dictionary_organisms_gprofiler_inv=dict(zip(dictionary_organisms_gprofiler.values(), dictionary_organisms_gprofiler.keys()))

    # Opening JSON file
    f=open(filename)
    Orthology_all=json.load(f)

    Seq_Prop={}
    N1=0
    for orths in Orthology_all:
        print(100*N1/len(Orthology_all))
        N1+=1
        for orga in Orthology_all[orths]:
            for gene_id in Orthology_all[orths][orga]:
                for seq_id in Orthology_all[orths][orga][gene_id]:
                    if seq_id=='name':
                        continue
                    seq=Orthology_all[orths][orga][gene_id][seq_id]
                    if X_is_None :
                        seq=seq.replace('X','')
                    if U_is_None :
                        seq=seq.replace('U','')
                    if Any_is_None :
                        seq=seq.replace('*','')

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

        with open(output_file,'w') as f:
            json.dump(Seq_Prop,f)