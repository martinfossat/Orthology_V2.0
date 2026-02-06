import json
import Fossat_utils as FU
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
import Orthology_utils as OU

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Checks the orthologs of a gene name list.\n
    The user needs to specify an original species name, corresponding to the gene name list, and a reference specie, which may be used for 
    homology comparison, in other steps. Additional species may be given, but those two first are required.\n
    Species name must follow the id in gProfiler (https://biit.cs.ut.ee/gprofiler/page/organism-list).\n
    The output is a file containing the gene name and gene IDs for all species, and that is required to use subsequent programs.""")
    parser.add_argument("--homologies_file","-hf",help='Name of the homology database file. Default is Gene_homologies.json')
    parser.add_argument("--name_file","-nf",help='Name of the names database json output file. Default is Names.json')
    parser.add_argument("--reference_specie","-rs",help="Reference specie")
    parser.add_argument("--top_specie","-ts",help="Top specie for the comparison (top_specie/norm_specie)")
    parser.add_argument("--norm_specie","-ns",help="Norm specie for the comparison (top_specie/norm_specie)")
    parser.add_argument("--top_iso_fraction","-tif",help="Top fraction of isoforms that are kept using overall homology as a metric. 0 is only most homologous, 1 is all.")
    parser.add_argument("--top_ortholog_fraction","-tof",help="Top fraction of orthologs that are kept using overall homology as a metric. 0 is only most homologous, 1 is all. ")
    parser.add_argument("--min_len_ratio","-mlr",help="Mimimun length ratio between the species and the reference species for orthologs to be plotted")
    parser.add_argument("--bin_width","-bw",help="Width of the bins for histograms. Default is 0.2")

    args = parser.parse_args()

    if args.min_len_ratio :
        try :
            min_len_ratio=float(args.min_len_ratio)
            if min_len_ratio<0 or min_len_ratio>1.:
                print("Wrong value for min_len_ratio. Please use a number between 0 and 1")
                quit()
        except :
            print("Wrong format for min_len_ratio")
            quit()
    else :
        min_len_ratio=0.0
    if args.bin_width :
        try :
            binwidth=float(args.bin_width)

        except :
            print("Wrong format for bin bin_width")
            quit()
    else :
        binwidth=0.02

    if args.top_iso_fraction :
        try :
            top_iso_fraction=float(args.top_allele_fraction)
            if top_iso_fraction<0 or top_iso_fraction>1.:
                print("Wrong value for top_allele_fraction. Please use a number between 0 and 1")
                quit()
        except :
            print("Wrong format for mof")
    else :
        top_iso_fraction=0.0

    if args.top_ortholog_fraction :
        try :
            top_ortholog_fraction=float(args.top_ortholog_fraction)
            if top_ortholog_fraction<0 or top_ortholog_fraction>1.:
                print("Wrong value for top_ortholog_fraction. Please use a number between 0 and 1")
                quit()
        except :
            print("Wrong format for mof")
            quit()
    else :
        top_ortholog_fraction=0.0

    if args.homologies_file :
        homologies_file=args.homologies_file
    else :
        homologies_file='Homology.json'

    if args.name_file :
        name_file=args.name_file

    else :
        name_file='Names.json'

    if args.reference_specie :
        ref_orga=args.reference_specie
    else :
        print("You must specify a reference specie")
        quit(1)

    if args.top_specie :
        specie_top=args.top_specie
    else :
        print("You must specify a top specie")
        quit(1)

    if args.norm_specie :
        specie_norm=args.norm_specie
    else :
        print("You must specify a top specie")
        quit(1)

    file_path=os.path.realpath(__file__)
    with open(os.path.dirname(file_path)+"/g_Profiler_Organisms_names_dic.json", "r") as fp:
        dictionary_organisms_gprofiler = json.load(fp)

    dictionary_organisms_gprofiler_inv=dict(zip(dictionary_organisms_gprofiler.values(), dictionary_organisms_gprofiler.keys()))

    f=open(homologies_file)
    save_homo=json.load(f)

    region_types=['all']#,'IDRs','FDs']

    save_all_compare=[]
    save_all_ids_ref=[]
    save_all_ids_top = []
    save_all_ids_norm = []
    save_all_compare_norm=[]
    save_all_compare_top = []

    for label in region_types:
        for orth_ref in save_homo[ref_orga][specie_top]:
            if not orth_ref in save_homo[ref_orga][specie_norm].keys():
                # Skip if an orhtolog exist in only one of the species that we are comparing
                continue
            save_temp_top,save_temp_ids_top=OU.get_all_homologies(ref_orga,orth_ref,label,specie_top,save_homo,min_len_ratio,top_iso_fraction)
            save_temp_norm,save_temp_ids_norm=OU.get_all_homologies(ref_orga,orth_ref,label,specie_norm,save_homo,min_len_ratio,top_iso_fraction)

            save_temp_top,best_top_id=OU.get_top_x_pct(save_temp_top,top_ortholog_fraction,names=save_temp_ids_top)
            save_temp_norm,best_norm_id=OU.get_top_x_pct(save_temp_norm,top_ortholog_fraction,names=save_temp_ids_norm)

            if not ((best_top_id is None ) or (best_norm_id is None)):
                # There are cases where the the ortholog does not exist in one organism, so don't use
                save_all_ids_top+=[best_top_id]
                save_all_ids_norm+=[best_norm_id]
                save_all_ids_ref+=[orth_ref]
                save_all_compare+=[np.mean(save_temp_top)/np.mean(save_temp_norm)]
                save_all_compare_top+=[np.mean(save_temp_top)]
                save_all_compare_norm+=[np.mean(save_temp_norm)]

    # List into array conversion
    save_all_ids_ref=np.array(save_all_ids_ref)
    save_all_ids_top=np.array(save_all_ids_top)
    save_all_ids_norm=np.array(save_all_ids_norm)
    save_all_compare=np.array(save_all_compare)
    save_all_compare_top= np.array(save_all_compare_top)
    save_all_compare_norm= np.array(save_all_compare_norm)

    ################ Plotting 2D histograms
    tick_sep=0.2
    x_ticks=np.round(np.arange(0,1+tick_sep, tick_sep),3)
    y_ticks =np.round(np.arange(0,1+tick_sep, tick_sep),3)
    ylabel='Normalized homology ratio '+specie_norm+' to '+ref_orga
    xlabel='Normalized homology ratio '+specie_top+' to '+ref_orga
    name='2D_hist_2_species.pdf'
    OU.plot_2d_hist(save_all_compare_norm,save_all_compare_top,xlabel,ylabel,x_ticks,y_ticks,name,binwidth)

    x_ticks=np.round(np.arange(0, np.nanmax(save_all_compare_top)+tick_sep, tick_sep),3)
    y_ticks =np.round(np.arange(0, np.nanmax(save_all_compare)+tick_sep, tick_sep),3)
    xlabel='Normalized homology ratio '+specie_top+' to '+ref_orga
    ylabel='Normalized homology ratio ('+specie_top+'/'+specie_norm+')'
    name='2D_hist_top.pdf'
    OU.plot_2d_hist(save_all_compare,save_all_compare_top,xlabel,ylabel,x_ticks,y_ticks,name,binwidth)

    x_ticks=np.round(np.arange(0, np.nanmax(save_all_compare_norm)+tick_sep, tick_sep),3)
    y_ticks =np.round(np.arange(0, np.nanmax(save_all_compare)+tick_sep, tick_sep),3)
    xlabel='Normalized homology ratio '+specie_norm+' to '+ref_orga
    ylabel='Normalized homology ratio ('+specie_top+'/'+specie_norm+')'
    name='2D_hist_norm.pdf'
    OU.plot_2d_hist(save_all_compare,save_all_compare_norm,xlabel,ylabel,x_ticks,y_ticks,name,binwidth)

    # Plotting 1D histograms
    bins=np.arange(0,1+binwidth,binwidth)
    plt.title('Normalized homology to hsapiens')
    plt.xlim(0,1)
    plt.hist(save_all_compare_norm,bins=bins,histtype='step',color='orange',label=specie_norm)
    plt.hist(save_all_compare_top,bins=bins,histtype='step',color='g',label=specie_top)
    plt.legend()
    plt.ylabel("Number of occurences")
    plt.xlabel("Normalized homology to "+ref_orga)
    plt.savefig('Hist_both_species.pdf')
    plt.close()

    bins = np.arange(0, np.nanmax(save_all_compare)+binwidth, binwidth)
    plt.xlim(0,np.nanmax(save_all_compare)+binwidth)
    plt.hist(save_all_compare,bins=bins,histtype='step',color='k')
    plt.ylabel("Number of occurences")
    plt.xlabel('Normalized homology to '+ref_orga+' ratio ('+specie_top+'/'+specie_norm+')')
    plt.savefig('Hist_ratio_'+specie_top+'_'+specie_norm+'.pdf')
    plt.close()

    # Sorting for the text file
    inds_sorted=np.argsort(save_all_compare)
    save_all_ids_ref=save_all_ids_ref[inds_sorted]
    save_all_ids_top=save_all_ids_top[inds_sorted]
    save_all_ids_norm=save_all_ids_norm[inds_sorted]
    save_all_compare=save_all_compare[inds_sorted]
    save_all_compare_top=save_all_compare_top[inds_sorted]
    save_all_compare_norm=save_all_compare_norm[inds_sorted]
    # Getting the names of the genes from the file
    f=open(name_file)
    save_names=json.load(f)

    # Wrtting the sorted homology comaprison
    W=ref_orga+'_id\t'+ref_orga+'_name\t'+specie_top+'_id\t'+specie_top+'_name\t'+specie_top+'_score\t'+specie_norm+'_id\t'+specie_norm+'_name\t'+specie_norm+'_score\t'+specie_top+'/'+specie_norm+'\n'
    for i in range(len(save_all_ids_ref)):
        W+=str(save_all_ids_ref[i])+'\t'+str(save_names[save_all_ids_ref[i]])+'\t'+str(save_all_ids_top[i])+'\t'+str(save_names[save_all_ids_top[i]])+'\t'+str(save_all_compare_top[i])+'\t'+str(save_all_ids_norm[i])+'\t'+str(save_names[save_all_ids_norm[i]])+'\t'+str(save_all_compare_norm[i])+'\t'+str(save_all_compare[i])+'\n'

    FU.write_file('Sorted_ids.txt',W)

