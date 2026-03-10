import json
import argparse
import os
import Orthology_utils as OU


label_dic={"rg":"Predicted IDR radius of gyration",
           're':"Predicted IDR end to end distance",
           'as':"Predicted IDR asphericity normalized",
           'nu':"Predicted IDR scaling exponent normalized",
           'pref':"Predicted IDR scaling prefactor normalized",
           'n':"IDR length",
           'pI':"Predicted isoelectric point"}

label_dic_short={"rg":"Radius of gyration",
                 're':"End to end distance",
                 'asph':"Asphericity",
                 'nu':"Scaling exponent",
                 'pref':"Scaling prefactor ",
                 'n':"IDR length",
                 'pI':"Isoelectric point"}

# First is start, second is end, third is width
# None is get min or max
bin_properties={"rg":[None,None,5],
                're':[None,None,1],
                'asph':[0,1,0.01],
                'nu':[None,None,0.01],
                'pref':[None,None,0.1],
                'n':[None,None,0.1],
                'pI':[0,14,0.1],
                'Kappa':[0,1,0.01],
                'FCR':[0,1,0.01],
                'NCPR':[0,1,0.01],
                'pct_folded':[0,100,1],
                'NCPR':[0,1,0.05],
                'f-':[0,1,0.01],
                'f+':[0,1,0.01],
                'n':[None,None,1]}

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""This program reads the Orthology and Sequence properties files and plot all relevant properties individually for each species.""")
    parser.add_argument("--orthology_file","-of",help='Name of the orthology database json output file. Default is Orthology.json')
    parser.add_argument("--properties_file", "-pf",help='Name of the sequence properties file. Default is Sequence_properties.json')

    ############### NEW
    parser.add_argument("--seq_labels","-sl",help="Labels for the sequence annotation. Must correspond to name given in other programs. Default is \"all\", which is the whole sequence. Can be multiple arguments.",default=[],nargs='+')
    parser.add_argument("--annotation_file","-af",help='Name of the gene annotation json file. Default is Gene_annotations.json')
    parser.add_argument("--annotation_type","-at",help='Name of the annotation type (PC, BP, MF, SML, SAL). Type -at list_types to get details',default=[],nargs='+')
    parser.add_argument("--annotation_names","-an",help='Annotation names for each corresponding annotation type.',default=[],nargs='+')
    parser.add_argument("--plot_list", "-pl", help='List genes to be plotted separately. Genes must be in columns, each first element of the column being the label for the gene list.')

    args = parser.parse_args()

    OU.check_and_create_rep('Plots')
    if args.plot_list :
        plot_list_file=args.plot_list
        data=OU.load_file(plot_list_file)
        plot_list=[[]for j in range(len(data[0].split('\t')))]
        plot_list_name=[data[0].split('\t')[j].replace('\n','') for j in range(len(data[0].split('\t')))]
        for i in range(1,len(data)):
            temp=data[i].split('\t')
            for j in range(len(temp)):
                to_add=temp[j].replace('\n','')
                if to_add!='':
                    plot_list[j]+=[temp[j].replace('\n','')]
    else :
        plot_list_name=[]

    if args.orthology_file :
        orthology_file=args.orthology_file
    else :
        orthology_file='Orthology.json'

    if args.properties_file :
        properties_file=args.properties_file
    else :
        properties_file='Sequence_properties.json'


    if args.seq_labels :
        seq_labels=args.seq_labels
    else :
        seq_labels=["all"]

    if args.annotation_file:
        annotation_file=args.annotation_file
    else:
        annotation_file='Gene_annotations.json'

    if args.annotation_type :
        annotation_type=args.annotation_type
        if annotation_type[0]=="list_types":
            dict=OU.get_gene_annotation_label_dic()
            print("Annotation types : ")
            for i in dict:
                print(dict[i]+" : "+i)
            quit()
    else :
        annotation_type=[]

    if args.annotation_names :
        annotation_names=args.annotation_names
    else :
        annotation_names=[]

    if len(annotation_names)!=len(annotation_type):
        print("You must specify as many annotation types as annotation names")
        quit()
    annotation={}
    for i in range(len(annotation_type)):
        if not annotation_type[i] in annotation:
            annotation[annotation_type[i]]=[]
        annotation[annotation_type[i]]+=[annotation_names[i]]
    del annotation_type,annotation_names

    if os.path.exists(annotation_file):
        print('Loading gene annotations file')
        f=open(annotation_file)
        Gene_annotation=json.load(f)
    else :
        Gene_annotation={}
    print('Loading orthology file')
    f=open(orthology_file)
    orthology=json.load(f)

    print('Loading sequence properties file')
    # This contains the precomputed sequence feature, ensemble properties and charge properties
    f=open(properties_file)
    properties=json.load(f)

    region_types=seq_labels
    save_all_prop={}

    species=['hsapiens','drerio','mmusculus']

    props_pH=['pI']
    prop_feat=['NCPR',"FCR","Kappa",]
    prop_ens=["asph","rg","re"]
    prop_ext=['q_at_pH_7.4']
    props=props_pH+prop_feat+prop_ens+prop_ext

    for region_type in region_types:
        masks_all={}
        for at in annotation:
            masks_all[at]={}
            for an in annotation[at]:
                masks_all[at][an]={}
                for orga in species:
                    masks_all[at][an][orga]=[]

        print("Retrieving data "+region_type)
        for orth in orthology:
            for orga in orthology[orth]:
                for orth_id in orthology[orth][orga]:
                    # Retrieving data
                    for prot_id in orthology[orth][orga][orth_id]:
                        if not region_type in properties[prot_id]:
                            continue
                        if not region_type in save_all_prop:
                            save_all_prop[region_type]={}

                        for region in properties[prot_id][region_type]:
                            for prop in props:
                                if not prop in save_all_prop[region_type]:
                                    save_all_prop[region_type][prop]={}
                                if not orga in save_all_prop[region_type][prop]:
                                    save_all_prop[region_type][prop][orga]=[]
                                if prop in properties[prot_id][region_type][region]:

                                    save_all_prop[region_type][prop][orga]+=[float(properties[prot_id][region_type][region][prop])]
                                else :
                                    save_all_prop[region_type][prop][orga]+=[float('NaN')]
                            # Making the masks
                            for at in annotation:
                                for an in range(len(annotation[at])):
                                    if orth_id in Gene_annotation:
                                        if Gene_annotation[orth_id][at] is None:
                                            masks_all[at][annotation[at][an]][orga]+=[False]
                                            continue
                                    else:
                                        masks_all[at][annotation[at][an]][orga]+=[False]
                                        continue

                                    if annotation[at][an] in Gene_annotation[orth_id][at]:
                                        masks_all[at][annotation[at][an]][orga]+=[True]
                                    else:
                                        masks_all[at][annotation[at][an]][orga]+=[False]

        print("Plotting "+region_type)
        for at in annotation:
            for an in range(len(annotation[at])):
                if not region_type in save_all_prop:
                    continue
                print("Plotting "+at+'_'+annotation[at][an])
                masked_database=OU.apply_mask(save_all_prop[region_type],masks_all[at][annotation[at][an]])
                OU.plot_hist(masked_database,region_type,species,label_dic,bin_properties)

                ensembles=['FCR','f+','f-','Kappa','asph','pct_folded']
                name=at+"_"+annotation[at][an]+'_Sequence_features'
                OU.plot_bar_charts(masked_database,region_type,species,ensembles,name,label_dic_short)

                ensembles=['NCPR']
                name=at+"_"+annotation[at][an]+'_NCPR'
                OU.plot_bar_charts(masked_database,region_type,species,ensembles,name,label_dic_short)

                ensembles=['q_at_pH_7.4']
                name=at+"_"+annotation[at][an]+'_Charge_at_pH'
                OU.plot_bar_charts(masked_database,region_type,species,ensembles,name,label_dic_short)

                ensembles=['pI']
                name=at+"_"+annotation[at][an]+'_Isoelectric_point'
                OU.plot_bar_charts(masked_database,region_type,species,ensembles,name,label_dic_short)

                ensembles=['rg','re','n']
                name=at+"_"+annotation[at][an]+'_Predicted_Ensemble_properties'
                OU.plot_bar_charts(masked_database,region_type,species,ensembles,name,label_dic_short)

        OU.plot_hist(save_all_prop[region_type],region_type,species,label_dic,bin_properties)

        ensembles=['NCPR','FCR','f+','f-','Kappa','pct_folded','asph']
        name='Sequence_features'
        OU.plot_bar_charts(save_all_prop[region_type],region_type,species,ensembles,name,label_dic_short)

        ensembles=['q_at_pH_7.4']
        name='Charge_at_pH'
        OU.plot_bar_charts(save_all_prop[region_type],region_type,species,ensembles,name,label_dic_short)

        ensembles=['pI']
        name='Isoelectric_point'
        OU.plot_bar_charts(save_all_prop[region_type],region_type,species,ensembles,name,label_dic_short)

        ensembles=['rg','re','n']
        name='Predicted_IDR_properties'
        OU.plot_bar_charts(save_all_prop[region_type],region_type,species,ensembles,name,label_dic_short)

