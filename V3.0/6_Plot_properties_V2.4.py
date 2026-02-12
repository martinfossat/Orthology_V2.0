import json
import Fossat_utils as FU
import argparse
import Orthology_utils as OU

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Checks the orthologs of a gene name list.\n
    The user needs to specify an original species name, corresponding to the gene name list, and a reference specie, which may be used for 
    homology comparison, in other steps. Additional species may be given, but those two first are required.\n
    Species name must follow the id in gProfiler (https://biit.cs.ut.ee/gprofiler/page/organism-list).\n
    The output is a file containing the gene name and gene IDs for all species, and that is required to use subsequent programs.""")
    parser.add_argument("--orthology_file","-of",help='Name of the orthology database json output file. Default is Orthology.json')
    parser.add_argument("--properties_file", "-pf",help='Name of the sequence properties file. Default is Sequence_properties.json')
    args = parser.parse_args()

    FU.check_and_create_rep('Plots')

    if args.orthology_file :
        orthology_file=args.orthology_file
    else :
        orthology_file='Orthology.json'

    if args.properties_file :
        properties_file=args.properties_file
    else :
        properties_file='Sequence_properties.json'

    f=open(orthology_file)
    orthology=json.load(f)
    # This contains the precomputed sequence feature, ensemble properties and charge properties
    f=open(properties_file)
    properties=json.load(f)


    region_types=['all','IDRs','NFDs','FDs']
    save_all_prop={}
    species=[]
    for orth in orthology:
        for orga in orthology[orth]:
            if not orga in species:
                species+=[orga]
            for orth_id in orthology[orth][orga]:
                for prot_id in orthology[orth][orga][orth_id]:
                    for region_type in properties[prot_id]:
                        for region in properties[prot_id][region_type]:
                            for prop in properties[prot_id][region_type][region]:
                                if not prop in save_all_prop:
                                    save_all_prop[prop]={}
                                if not region_type in save_all_prop[prop]:
                                    save_all_prop[prop][region_type]={}
                                if not orga in save_all_prop[prop][region_type]:
                                    save_all_prop[prop][region_type][orga]=[]
                                save_all_prop[prop][region_type][orga]+=[float(properties[prot_id][region_type][region][prop])]

    # colors=[plt.cm.gist_rainbow(i/(len(species))) for i in range(len(species))]
    # colors_ens=['g','r','orange']

    # type_ensemble=[]

    label_dic={"rg" : "Predicted IDR radius of gyration",
               're': "Predicted IDR end to end distance",
               'as': "Predicted IDR asphericity normalized",
               'nu': "Predicted IDR scaling exponent normalized",
               'pref':"Predicted IDR scaling prefactor normalized",
               'n' :"IDR length",
               'pI' : "Predicted isoelectric point"}

    label_dic_short={"rg" : "Radius of gyration",
               're': "End to end distance",
               'asph': "Asphericity",
               'nu': "Scaling exponent",
               'pref':"Scaling prefactor ",
               'n':"IDR length",
               'pI' : "Isoelectric point"}
    # First is start, second is end, third is width
    # None is get min or max
    bin_properties={"rg" : [None,None,5],
               're': [None,None,1],
               'asph': [0,1,0.01],
               'nu': [None,None,0.01],
               'pref':[None,None,0.1],
               'n':[None,None,0.1],
               'pI' : [0,14,0.1],
               'Kappa':[0,1,0.01],
                'FCR':[0,1,0.01],
                'NCPR':[0,1,0.01],
                'pct_folded':[0,100,1],
                'NCPR':[0,1,0.05],
                'f-':[0,1,0.01],
                'f+':[0,1,0.01],
                'n': [None,None,1]}

    OU.plot_hist(save_all_prop,region_types,species,label_dic,bin_properties)

    ensembles=['NCPR','FCR','f+','f-','Kappa','pct_folded']
    name='Sequence_features'
    OU.plot_bar_charts(save_all_prop,region_types,species,ensembles,name,label_dic_short)

    ensembles=['pI','q_at_pH_7.4']
    name='pH properties'
    OU.plot_bar_charts(save_all_prop,region_types,species,ensembles,name,label_dic_short)
    #OU.plot_hist(save_all_prop,region_types,bin_properties)

    ensembles=['rg','re','asph','n']
    name='Predicted IDR properties'
    OU.plot_bar_charts(save_all_prop,region_types,species,ensembles,name,label_dic_short)


