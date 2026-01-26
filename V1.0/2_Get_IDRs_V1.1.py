import json
import metapredict as MP
from sparrow.predictors import batch_predict
import argparse
import os
def get_bounds_folded(dis_bounds,seq):
    if len(dis_bounds)==0:
        return [[0,len(seq)]]
    if dis_bounds[0][0]==0:
        bounds_folded=[]
    else :
        bounds_folded=[[0]]
    for b in range(len(dis_bounds)):
        if dis_bounds[b][0]!=0:
            bounds_folded[-1]+=[dis_bounds[b][0]]
        if dis_bounds[b][1]!=len(seq):
            bounds_folded+=[[dis_bounds[b][1]]]
    if dis_bounds[b][1]!=len(seq):
        bounds_folded[-1]+=[len(seq)]
    return bounds_folded

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
        output_file='Proteins_IDR_properties.json'

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

    IDRs={}
    N1=0
    for orths in Orthology_all:
        print(100*N1/len(Orthology_all))
        N1+=1
        for orga in Orthology_all[orths]:
            for N_orth in Orthology_all[orths][orga]:
                for vers in Orthology_all[orths][orga][N_orth]['version']:
                    Orthology_all[orths][orga][N_orth]['version'][vers]['dis_doms']={}
                    seq=Orthology_all[orths][orga][N_orth]['version'][vers]['seq']
                    if X_is_None :
                        seq=seq.replace('X','')
                    if U_is_None :
                        seq=seq.replace('U','')
                    if Any_is_None :
                        seq=seq.replace('*','')

                    seq_id=Orthology_all[orths][orga][N_orth]['version'][vers]['seq_id']
                    IDRs[seq_id]={}
                    IDRs[seq_id]['all']={}
                    IDRs[seq_id]['all'][str(0)+'_'+str(len(seq))]={}

                    IDRs[seq_id]['IDRs']={}
                    dis_bounds=MP.predict_disorder_domains(seq).disordered_domain_boundaries
                    fd_bounds=get_bounds_folded(dis_bounds,seq)
                    # IDRs[seq_id]['IDRs']['bounds']=dis_bounds
                    # IDRs[seq_id]['IDRs']['rgs']=[]
                    # IDRs[seq_id]['IDRs']['res']=[]
                    # IDRs[seq_id]['IDRs']['asphs']=[]
                    # IDRs[seq_id]['IDRs']['nus']=[]
                    # IDRs[seq_id]['IDRs']['prefs']=[]
                    temp={}
                    for bounds in dis_bounds:
                        label_bounds=str(bounds[0])+'_'+str(bounds[1])
                        temp[label_bounds]=seq[bounds[0]:bounds[1]]
                        IDRs[seq_id]['IDRs'][label_bounds]={}
                        IDRs[seq_id]['IDRs'][label_bounds]['rg']=str(batch_predict.batch_predict(temp,network='scaled_rg')[label_bounds][1])
                        IDRs[seq_id]['IDRs'][label_bounds]['re']=str(batch_predict.batch_predict(temp,network='scaled_re')[label_bounds][1])
                        IDRs[seq_id]['IDRs'][label_bounds]['asph']=str(batch_predict.batch_predict(temp,network='asphericity')[label_bounds][1])
                        IDRs[seq_id]['IDRs'][label_bounds]['nu']=str(batch_predict.batch_predict(temp,network='scaling_exponent')[label_bounds][1])
                        IDRs[seq_id]['IDRs'][label_bounds]['pref']=str(batch_predict.batch_predict(temp,network='prefactor')[label_bounds][1])

                    IDRs[seq_id]['FDs']={}

                    for bounds in fd_bounds:
                        label_bounds=str(bounds[0])+'_'+str(bounds[1])
                        temp[label_bounds]=seq[bounds[0]:bounds[1]]
                        IDRs[seq_id]['FDs'][label_bounds]={}

                        # IDRs[seq_id]['IDRs'][label_bounds]['rgs']=batch_predict.batch_predict(temp,network='scaled_rg')[label_bounds]
                        # IDRs[seq_id]['IDRs'][label_bounds]['res']=batch_predict.batch_predict(temp,network='scaled_re')[label_bounds]
                        # IDRs[seq_id]['IDRs'][label_bounds]['asphs']=batch_predict.batch_predict(temp,network='asphericity')[label_bounds]
                        # IDRs[seq_id]['IDRs'][label_bounds]['nus']=batch_predict.batch_predict(temp,network='scaling_exponent')[label_bounds]
                        # IDRs[seq_id]['IDRs'][label_bounds]['prefs']=batch_predict.batch_predict(temp,network='prefactor')[label_bounds]

                    # temp={}
                    # for i in range(len(IDRs[seq_id]['dis_doms']['bounds'])):
                    #     inds=IDRs[seq_id]['dis_doms']['bounds'][i]
                    #     temp[str(inds[0])+'_'+str(inds[1])]=seq[inds[0]:inds[1]]
                    #
                    #     rg=batch_predict.batch_predict(temp,network='scaled_rg')
                    #     re=batch_predict.batch_predict(temp,network='scaled_re')
                    #     asph=batch_predict.batch_predict(temp,network='asphericity')
                    #     nu=batch_predict.batch_predict(temp,network='scaling_exponent')
                    #     pref=batch_predict.batch_predict(temp,network='prefactor')
                    #
                    #     IDRs[seq_id]['dis_doms']['rgs']+=[str(rg[str(inds[0])+'_'+str(inds[1])][1])]
                    #     IDRs[seq_id]['dis_doms']['res']+=[str(re[str(inds[0])+'_'+str(inds[1])][1])]
                    #     IDRs[seq_id]['dis_doms']['asphs']+=[str(asph[str(inds[0])+'_'+str(inds[1])][1])]
                    #     IDRs[seq_id]['dis_doms']['nus']+=[str(nu[str(inds[0])+'_'+str(inds[1])][1])]
                    #     IDRs[seq_id]['dis_doms']['prefs']+=[str(pref[str(inds[0])+'_'+str(inds[1])][1])]

        with open(output_file,'w') as f:
            json.dump(IDRs,f)