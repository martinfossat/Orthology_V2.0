import json
import matplotlib
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
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
#File was renamed from 4_Get_homology_scores_V1.5.py
def plot_homology_scores(save_homo,prefix):
    plot_mean=np.zeros((len(region_types),len(species)))
    plot_std=np.zeros((len(region_types),len(species)))
    reg_n=0
    max_len=0
    for region in region_types:
        i=0
        plt.figure(region)
        for specie in save_homo[region]:
            temp_save_all_ens=np.array([])
            temp_save_all_all=np.array([])
            #This will contain all of the
            if not specie in species :
                continue
            for name in save_homo[region][specie]:
                # Each of these are an orthologs
                if len(save_homo[region][specie][name])==0:
                    continue

                temp_ens_mean=np.array(save_homo[region][specie][name])
                temp_all_mean=np.array(save_homo['all'][specie][name])

                # Getting rid of the nan, in both arrays

                nans_bool=np.invert(np.isnan(temp_ens_mean))
                temp_ens_mean=temp_ens_mean[nans_bool]
                temp_all_mean=temp_all_mean[nans_bool]

                # We could take the top pct here
                # This should be based on the all score only
                top_ind=int(np.ceil(len(temp_all_mean)*(1-pct_top_orth)))
                ind_sort=np.array(np.argsort(temp_all_mean)[min(top_ind,len(temp_all_mean)-1):])
                # top_ind=int(np.ceil(len(temp_alleles)*(1-pct_top_alle)))
                # save_homo[region_type][orga][orth_ref]+=[np.mean(np.sort(np.array(temp_alleles))[min(top_ind,len(temp_alleles)-1):])]
                temp_ens_mean=temp_ens_mean[ind_sort]
                temp_all_mean=temp_all_mean[ind_sort]

                temp_save_all_ens=np.append(temp_save_all_ens,temp_ens_mean)
                temp_save_all_all=np.append(temp_save_all_all,temp_all_mean)
                max_len=max(max_len,len(temp_ens_mean))

            if len(temp_save_all_ens)>=1 :
                bins=np.arange(0,np.amax(temp_save_all_ens),binwidth)
                plt.hist(temp_save_all_ens,bins=bins,label=specie,histtype='step',linewidth=2)
                #This is the mean over orthologs
                plot_mean[reg_n,i]=np.mean(temp_save_all_ens)
                plot_std[reg_n,i]=np.std(temp_save_all_ens)
            i+=1

        if max_len!=1:
            plt.legend()
            plt.savefig('Hist_homologies'+region_types[reg_n]+'.pdf')
        plt.close()
        reg_n+=1

    x=[i for i in range(len(save_homo))]
    N_offset=1
    for i in range(len(species)):
        if i!=0:
            pre='_'
        else :
            pre=''
        for j in range(len(region_types)):
            pos=i-0.5+(j+1)/(len(region_types)+N_offset)
            plt.bar(pos,plot_mean[j][i],color=colors_ens[j],width=1/(len(region_types)+N_offset),label=pre+region_types[j],yerr=plot_std[j][i],capsize=5)

    xticks=[i for i in range(len(species))]
    plt.title('Average relative homology to hsapiens')
    plt.ylabel('Mean relative homology score')
    plt.xticks(xticks,species)
    plt.legend()
    plt.savefig(prefix+'Homologies_bar_plot.pdf')
    plt.close()

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Checks the orthologs of a gene name list.\n
    The user needs to specify an original species name, corresponding to the gene name list, and a reference specie, which may be used for 
    homology comparison, in other steps. Additional species may be given, but those two first are required.\n
    Species name must follow the id in gProfiler (https://biit.cs.ut.ee/gprofiler/page/organism-list).\n
    The output is a file containing the gene name and gene IDs for all species, and that is required to use subsequent programs.""")
    parser.add_argument("--orthology_file","-of",help='Name of the orthology database json output file. Default is Orthology.json')
    parser.add_argument("--properties_file", "-pf",help='Name of the sequence properties file. Default is Sequence_properties.json')
    # parser.add_argument("--reference_specie","-rs",help="Reference specie")
    # parser.add_argument("--top_allele_fraction","-taf",help="Top fraction of alleles that are kept using overall homology as a metric. 0 is only most homologous, 1 is all. ")
    # parser.add_argument("--top_ortholog_fraction","-tof",help="Top fraction of orthologs that are kept using overall homology as a metric. 0 is only most homologous, 1 is all. ")
    # parser.add_argument("--additional_species","-as",help="Additional species",default=[],nargs='+')
    args = parser.parse_args()

    # if args.top_allele_fraction :
    #     try :
    #         top_allele_fraction=float(args.top_allele_fraction)
    #         if top_allele_fraction<0 or top_allele_fraction>1.:
    #             print("Wrong value for top_allele_fraction. Please use a number between 0 and 1")
    #             quit()
    #     except :
    #         print("Wrong format for mof")
    # else :
    #     top_allele_fraction=0.0

    # if args.top_ortholog_fraction :
    #     try :
    #         top_ortholog_fraction=float(args.top_ortholog_fraction)
    #         if top_ortholog_fraction<0 or top_ortholog_fraction>1.:
    #             print("Wrong value for top_ortholog_fraction. Please use a number between 0 and 1")
    #             quit()
    #     except :
    #         print("Wrong format for mof")
    # else :
    #     top_ortholog_fraction=1.0


    if args.orthology_file :
        orthology_file=args.orthology_file
    else :
        orthology_file='Orthology.json'

    if args.properties_file :
        properties_file=args.properties_file
    else :
        properties_file='Sequence_properties.json'

    #
    # if args.reference_specie :
    #     ref_org=args.reference_specie
    # else :
    #     print("You must specify a reference specie")
    #     quit(1)

    # if args.additional_species :
    #     other_organism=args.additional_species
    # else :
    #     other_organism=[]

    file_path=os.path.realpath(__file__)

    with open(os.path.dirname(file_path)+"/g_Profiler_Organisms_names_dic.json", "r") as fp:
        dictionary_organisms_gprofiler = json.load(fp)

    dictionary_organisms_gprofiler_inv=dict(zip(dictionary_organisms_gprofiler.values(), dictionary_organisms_gprofiler.keys()))
    #
    # all_orga=[ref_org]+other_organism
    # species=all_orga
    #
    # pct_top_orth=top_ortholog_fraction
    # pct_top_alle=top_allele_fraction

    f=open(orthology_file)
    orthology=json.load(f)

    f=open(properties_file)
    properties=json.load(f)

    # This contains the precomputed IDR ensemble properties from ALBATROSS
    # ref_orga=ref_org
    region_types=['all','IDRs','NFDs','FDs']
    #save_all_prop={region_types[j]:{species[i]:{} for i in range(len(species))} for j in range(len(region_types))}
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


    binwidth=0.05

    colors=[plt.cm.gist_rainbow(i/(len(species))) for i in range(len(species))]
    colors_ens=['g','r','orange']

    type_ensemble=['Rg_ratio','Re_ratio','As_ratio','Nu_ratio','Pref_ratio','N_ratio']

    save_ensemble={species[j]:{type_ensemble[i]:[] for i in range(len(type_ensemble))} for j in range(len(species))}

    label_dic={"Rg_ratio" : "Predicted IDR radius of gyration normalized by "+ref_orga,
               'Re_ratio': "Predicted IDR end to end distance normalized by "+ref_orga,
               'As_ratio': "Predicted IDR asphericity normalized by "+ref_orga,
                'Nu_ratio': "Predicted IDR scaling exponent normalized by "+ref_orga,
               'Pref_ratio':"Predicted IDR scaling prefactor normalized by "+ref_orga,
               'N_ratio':"IDR length normalized by "+ref_orga}

    label_dic_short={"Rg_ratio" : "Radius\nof gyration",
               'Re_ratio': "End to end\ndistance",
               'As_ratio': "Asphericity",
               'Nu_ratio': "Scaling\nexponent",
               'Pref_ratio':"Scaling\nprefactor ",
               'N_ratio':"IDR\nlength"}

    for region_type in region_types:
        i=0
        plt.close()
        N_offset=1
        for prop in save_all_prop:
            j=0
            for orga in save_all_prop[prop][region_type]:
                if i!=0:
                    pre='_'
                else :
                    pre=''

                pos=i-0.5+(j+1)/(len(species)+N_offset)
                plt.figure('barchart')
                plt.bar(pos,np.mean(save_all_prop[prop][region_type][orga]),color=colors[j],width=1/(len(species)+N_offset),label=pre+orga,yerr=np.std(save_all_prop[prop][region_type][orga]),capsize=3)
                plt.figure(prop)
                bins=np.arange(0,np.amax(save_all_prop[prop][region_type][orga]),binwidth)

                bins_center=(bins[1:]+bins[:-1])/2.
                hist,bins_temp=np.histogram(save_all_prop[prop][region_type][orga],bins=bins)
                hist=hist/np.amax(hist)
                plt.step(bins_center,hist,linewidth=1,color=colors[j])
                j+=1

            plt.xlabel(label_dic[prop])
            i+=1
            plt.legend()
            plt.savefig("Hist_ensemble_"+prop+'.pdf')

        plt.close()
        plt.figure('barchart')
        plt.title('Predicted IDR ensemble properties of orthologs relative to '+str(ref_orga))
        plt.legend()
        xticks=[i for i in range(len(type_ensemble))]
        plt.xticks(xticks,[label_dic_short[type_ensemble[i]] for i in range(len(type_ensemble))])
        plt.savefig("Ensemble_properties_homologies"+region_type+".pdf")
