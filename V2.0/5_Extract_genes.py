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
    parser.add_argument("--homologies_file_name","-fh",help='Name of the input homologies json file. Default is Homologies.txt')
    parser.add_argument("--orthology_file_name","-fo",help='Name of the Orthology json file. Default is Homologies.txt')
    parser.add_argument("--reference_specie","-rs",help="Reference specie")
    parser.add_argument("--top_iso_fraction","-taf",help="Top fraction of isoforms that are kept using overall homology as a metric. 0 is only most homologous, 1 is all.")
    parser.add_argument("--top_ortholog_fraction","-tof",help="Top fraction of orthologs that are kept using overall homology as a metric. 0 is only most homologous, 1 is all. ")
    # Must add the top and norm specie
    parser.add_argument("--min_len_ratio","-mlr",help="")
    parser.add_argument("--additional_species","-as",help="Additional species",default=[],nargs='+')
    args = parser.parse_args()
    name_file='Names.json'
    if args.top_iso_fraction :
        try :
            top_iso_fraction=float(args.top_allele_fraction)
            if top_iso_fraction<0 or top_iso_fraction>1.:
                print("Wrong value for top_allele_fraction. Please use a number between 0 and 1")
                quit()
        except :
            print("Wrong format for mof")
    else :
        top_iso_fraction=1.0

    if args.top_ortholog_fraction :
        try :
            top_ortholog_fraction=float(args.top_ortholog_fraction)
            if top_ortholog_fraction<0 or top_ortholog_fraction>1.:
                print("Wrong value for top_ortholog_fraction. Please use a number between 0 and 1")
                quit()
        except :
            print("Wrong format for mof")
    else :
        top_ortholog_fraction=0.0

    if args.homologies_file_name :
        filename=args.homologies_file_name
    else :
        filename='Gene_homologies.json'

    if args.orthology_file_name :
        filename_ortho=args.homologies_file_name
    else :
        filename_ortho='Gene_orthology.json'


    if args.reference_specie :
        ref_org=args.reference_specie
    else :
        print("You must specify a reference specie")
        quit(1)

    if args.additional_species :
        other_organism=args.additional_species
    else :
        other_organism=[]
    ref_orga = ref_org
    file_path=os.path.realpath(__file__)

    with open(os.path.dirname(file_path)+"/g_Profiler_Organisms_names_dic.json", "r") as fp:
        dictionary_organisms_gprofiler = json.load(fp)

    dictionary_organisms_gprofiler_inv=dict(zip(dictionary_organisms_gprofiler.values(), dictionary_organisms_gprofiler.keys()))

    all_orga=[ref_org]+other_organism
    species=all_orga

    pct_top_orth=top_ortholog_fraction
    #pct_top_alle=top_allele_fraction

    f=open(filename)
    save_homo=json.load(f)
    # This contains the precomputed IDR ensemble properties from ALBATROSS
    region_types=['all']#,'IDRs','FDs']
    # I need to save the top of each score as a function of the top of all other scores
    # What I can do is keep the scores in lists with the same index.
    specie_norm='mmusculus'
    specie_top = 'drerio'

    plt.close()
    save_all_compare=[]
    save_all_ids_ref=[]
    save_all_ids_top = []
    save_all_ids_norm = []

    save_all_compare_norm=[]
    save_all_compare_top = []
    if pct_top_orth !=0. :
        print("This is not currently possible")
        exit()

    for label in region_types:
        for orth_ref in save_homo[ref_orga][specie_top]:
            if not orth_ref in save_homo[ref_orga][specie_norm].keys():
                # Skip if an orhtolog exist in only one of the species that we are comparing
                continue
            save_temp_top=[]
            save_temp_norm=[]
            #Norm
            # For each orthologs, find all the scores for the top and norm species.
            save_temp_ids_top=[]
            save_temp_ids_norm=[]
            for orth in save_homo[ref_orga][specie_top][orth_ref]:
                temp_top=[]
                for prot_ref in save_homo[ref_orga][specie_top][orth_ref][orth]:
                    for prot in save_homo[ref_orga][specie_top][orth_ref][orth][prot_ref]:
                        temp_region=[]
                        names_iso_temp=[]
                        for region in save_homo[ref_orga][specie_top][orth_ref][orth][prot_ref][prot][label]:
                            temp_region+=[float(save_homo[ref_orga][specie_top][orth_ref][orth][prot_ref][prot][label][region]['Homology_ratio'])]
                        if len(temp_region)!=0:
                            temp_top+=[np.mean(temp_region)]
                save_temp_top+=[np.mean(OU.get_top_x_pct(temp_top,top_iso_fraction)[0])]
                save_temp_ids_top+=[orth]

            # For each orthologs, find all the scores for the top and norm species.
            for orth in save_homo[ref_orga][specie_norm][orth_ref]:
                temp_top=[]
                for prot_ref in save_homo[ref_orga][specie_norm][orth_ref][orth]:
                    for prot in save_homo[ref_orga][specie_norm][orth_ref][orth][prot_ref]:
                        temp_region=[]
                        for region in save_homo[ref_orga][specie_norm][orth_ref][orth][prot_ref][prot][label]:
                            temp_region+=[float(save_homo[ref_orga][specie_norm][orth_ref][orth][prot_ref][prot][label][region]['Homology_ratio'])]
                        if len(temp_region)!=0:
                            temp_top+=[np.mean(temp_region)]

                save_temp_norm+=[np.mean(OU.get_top_x_pct(temp_top,top_iso_fraction)[0])]
                save_temp_ids_norm+=[orth]

            save_temp_top,best_top_id=OU.get_top_x_pct(save_temp_top,top_ortholog_fraction,names=save_temp_ids_top)
            save_temp_norm,best_norm_id=OU.get_top_x_pct(save_temp_norm,top_ortholog_fraction,names=save_temp_ids_norm)

            # Now I need to get the IDs of the best allele, best orthologs, for each ortholog

            save_all_ids_top+=[best_top_id]
            save_all_ids_norm+=[best_norm_id]
            save_all_ids_ref+=[orth_ref]
            save_all_compare+=[np.mean(save_temp_top)/np.mean(save_temp_norm)]
            save_all_compare_top+=[np.mean(save_temp_top)]
            save_all_compare_norm+=[np.mean(save_temp_norm)]

    save_all_ids_ref=np.array(save_all_ids_ref)
    save_all_ids_top=np.array(save_all_ids_top)
    save_all_ids_norm=np.array(save_all_ids_norm)
    save_all_compare=np.array(save_all_compare)
    save_all_compare_top= np.array(save_all_compare_top)
    save_all_compare_norm= np.array(save_all_compare_norm)


    plt.figure('2D_hist_norm',figsize=(4,4))
    binwidth = 0.01

    bins_x=np.arange(0, np.nanmax(save_all_compare), binwidth)
    bins_y=np.arange(0, 1, binwidth)
    yticks_lab=[i for i in range(len(bins_x))]
    tick_sep=1
    x_ticks=np.arange(0, np.nanmax(save_all_compare), tick_sep)
    y_ticks = np.arange(0, np.nanmax(save_all_compare), tick_sep)
    hist=np.histogram2d(save_all_compare,save_all_compare_norm,bins=[bins_x,bins_y])[0]

    plt.yticks(x_ticks,x_ticks)
    plt.xticks(y_ticks,y_ticks)

    plt.ylabel('Normalized homology ratio ('+specie_top+'/'+specie_norm+')')
    plt.xlabel('Normalized homology ratio '+specie_norm+' to '+ref_org)

    plt.imshow(hist,aspect='auto',origin='lower',extent=(0.,np.nanmax(save_all_compare_norm),0.,np.nanmax(save_all_compare)),interpolation='None',cmap='Spectral')
    plt.ylim(0,3)
    plt.savefig('2D_hist_norm.pdf')
    plt.close()

    plt.figure('2D_hist_top',figsize=(4,4))
    hist=np.histogram2d(save_all_compare,save_all_compare_top,bins=[bins_x,bins_y])[0]
    plt.yticks(x_ticks,x_ticks)
    plt.xticks(y_ticks, y_ticks)
    plt.ylabel('Normalized homology ratio ('+specie_top+'/'+specie_norm+')')
    plt.xlabel('Normalized homology ratio '+specie_top+' to '+ref_org)
    plt.imshow(hist,aspect='auto',origin='lower',extent=(0.,np.nanmax(save_all_compare_top),0.,np.nanmax(save_all_compare)),interpolation='None',cmap='Spectral')
    plt.ylim(0,3)
    plt.savefig('2D_hist_top.pdf')
    plt.close()

    plt.figure('2D_hist_norm',figsize=(4,4))

    bins_x=np.arange(0, 1, binwidth)
    bins_y=np.arange(0, 1, binwidth)
    hist=np.histogram2d(save_all_compare_norm,save_all_compare_top,bins=[bins_x,bins_y])[0]

    plt.yticks(x_ticks,x_ticks)
    plt.xticks(y_ticks, y_ticks)
    plt.ylabel('Normalized homology ratio '+specie_norm+' to '+ref_org)
    plt.xlabel('Normalized homology ratio '+specie_top+' to '+ref_org)
    plt.plot([0,1],[0,1],linestyle='--',color='k')
    plt.imshow(hist,aspect='auto',origin='lower',interpolation='None',cmap='Spectral',extent=(0.,1,0.,1))

    plt.savefig('2D_hist_2_species.pdf')
    plt.close()


    f=open(name_file)
    save_names=json.load(f)

    inds_sorted=np.argsort(save_all_compare)
    save_all_ids_ref=save_all_ids_ref[inds_sorted]
    save_all_ids_top=save_all_ids_top[inds_sorted]
    save_all_ids_norm=save_all_ids_norm[inds_sorted]
    save_all_compare=save_all_compare[inds_sorted]
    save_all_compare_top=save_all_compare_top[inds_sorted]
    save_all_compare_norm=save_all_compare_norm[inds_sorted]




    W=ref_orga+'_id\t'+ref_orga+'_name\t'+specie_top+'_id\t'+specie_top+'_name\t'+specie_top+'_score\t'+specie_norm+'_id\t'+specie_norm+'_name\t'+specie_norm+'_score\t'+specie_top+'/'+specie_norm+'\n'
    for i in range(len(save_all_ids_ref)):
        W+=str(save_all_ids_ref[i])+'\t'+str(save_names[save_all_ids_ref[i]])+'\t'+str(save_all_ids_top[i])+'\t'+str(save_names[save_all_ids_top[i]])+'\t'+str(save_all_compare_top[i])+'\t'+str(save_all_ids_norm[i])+'\t'+str(save_names[save_all_ids_norm[i]])+'\t'+str(save_all_compare_norm[i])+'\t'+str(save_all_compare[i])+'\n'

    FU.write_file('Sorted_ids.txt',W)

    bins = np.arange(0, np.nanmax(save_all_compare_norm), binwidth)
    plt.title('Normalized homology to hsapiens ')
    plt.xlim(0,3)
    plt.hist(save_all_compare_norm,bins=bins)
    plt.savefig(specie_norm+'_hist.pdf')
    plt.close()

    bins = np.arange(0, np.nanmax(save_all_compare_top), binwidth)
    plt.title('Normalized homology to hsapiens '+specie_top)
    plt.xlim(0,3)
    plt.hist(save_all_compare_top,bins=bins)
    plt.savefig(specie_top+'_hist.pdf')
    plt.close()

    bins = np.arange(0, np.nanmax(save_all_compare), binwidth)
    plt.title('Normalized homology to hsapiens '+specie_norm)
    plt.xlim(0,3)
    plt.hist(save_all_compare,bins=bins)
    plt.savefig('Compare_species.pdf')
    plt.close()
