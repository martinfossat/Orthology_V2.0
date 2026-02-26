import json
import matplotlib
matplotlib.use("pgf")
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
pgf_with_latex = {
    "text.usetex": True,
    "pgf.preamble":
        r'\usepackage{color}',
    "font.family": "Arial" }
matplotlib.rcParams.update(pgf_with_latex)
import matplotlib.pyplot as plt
import numpy as np
import argparse
import Orthology_utils as OU
import os
global annotation_type,annotation_names,Gene_annotation

if __name__=="__main__":
    global annotation_type,annotation_names,Gene_annotation
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Plot the sequence comparison between two species with regard to a third reference specie.""")
    parser.add_argument("--homologies_file","-hf",help='Name of the homology database file. Default is Gene_homologies.json')
    parser.add_argument("--name_file","-nf",help='Name of the names database json output file. Default is Names.json')
    parser.add_argument("--reference_specie","-rs",help="Reference specie")
    parser.add_argument("--top_specie","-ts",help="Top specie for the comparison (top_specie/norm_specie)")
    parser.add_argument("--norm_specie","-ns",help="Norm specie for the comparison (top_specie/norm_specie)")
    parser.add_argument("--top_iso_fraction","-tif",help="Top fraction of isoforms that are kept using overall homology as a metric. 0 is only most homologous, 1 is all.")
    parser.add_argument("--top_orth_fraction","-tof",help="Top fraction of orthologs that are kept using overall homology as a metric. 0 is only most homologous, 1 is all. ")
    parser.add_argument("--min_len_ratio","-mlr",help="Mimimun length ratio between the species and the reference species for orthologs to be plotted")
    parser.add_argument("--bin_width","-bw",help="Width of the bins for histograms. Default is 0.2")
    parser.add_argument("--factor_len_ratio","-flr",help="Whether to multiply the homology by the length ratio. Can be 1 (on) or 0 (off). Default is 0")

    ############### NEW
    parser.add_argument("--seq_labels","-sl",help="Labels for the sequence annotation. Must correspond to name given in other programs. Default is \"all\", which is the whole sequence. Can be multiple arguments.",default=[],nargs='+')
    parser.add_argument("--annotation_file","-af",help='Name of the gene annotation json file. Default is Gene_annotations.json')
    parser.add_argument("--annotation_type","-at",help='Name of the annotation type (PC, BP, MF, SML, SAL). Type -at list_types to get details',default=[],nargs='+')
    parser.add_argument("--annotation_names","-an",help='Annotation names for each corresponding annotation type.',default=[],nargs='+')
    parser.add_argument("--iso_compare","-ic",help="How to compare isoforms. Option are : MLF (Minimum Length Fraction) Compares all protein isoform that have a pair length ratio above the value specified; MLO (Maximum Length Only) : only compare the longest isoforms between the two species. This is what is done in Ensembl. Default is MLO.")
    parser.add_argument("--plot_list", "-pl", help='List genes to be plotted separately. Genes must be in columns, each first element of the column being the label for the gene list.')
    args = parser.parse_args()

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

    if args.iso_compare :
        if args.iso_compare=='MLO':
            use_MLF=False
        elif args.iso_compare=='MLF':
            use_MLF=True
        else :
            print("Invalid argument for align_type")
            quit()
    else :
        use_MLF=False

    if args.factor_len_ratio:
            if args.factor_len_ratio=='0':
                factor_len_ratio=False
            elif args.factor_len_ratio=='1':
                factor_len_ratio=True
            else :
                print("Wrong value for factor_length_ratio. Please use 0 and 1")
                quit()
    else:
        factor_len_ratio=False

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
                print("Wrong value for top_iso_fraction. Please use a number between 0 and 1")
                quit()
        except :
            print("Wrong format for tif")
            quit()
    else :
        top_iso_fraction=0.0

    if args.top_orth_fraction :
        try :
            top_ortholog_fraction=float(args.top_orth_fraction)
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
    OU.initialize_gene_annotation_names(annotation)

    OU.check_and_create_rep('Plots')
    
    print('Loading homology file')
    f=open(homologies_file)
    save_homo=json.load(f)
    if os.path.exists(annotation_file):
        print('Loading gene annotations file')
        f=open(annotation_file)
        Gene_annotation=json.load(f)
    else :
        Gene_annotation={}
    OU.initialize_gene_annotation(Gene_annotation)

    #Now that the clas is created, I need to create a object for each specific class, that is, the combination of seq_labels, gene label type and gene label names
    AllSeqLab={}
    save_means={}
    # okay big change in V2.0, the gene annotation changes everything
    for label in seq_labels:
        OU.check_and_create_rep('Plots/'+label)
        OU.check_and_create_rep('Plots/'+label+'/Homology/')
        AllSeqLab[label]=OU.Homology_Compare_Manager()

        print("Retrieving data "+label)
        for orth_ref in save_homo[ref_orga][specie_top]:
            if not orth_ref in save_homo[ref_orga][specie_norm].keys():
                # Skip if an ortholog exist in only one of the species that we are comparing
                continue

            save_temp_top,temp_len_ratio_top,save_temp_ids_top=OU.get_all_homologies(ref_orga,orth_ref,label,specie_top,save_homo,min_len_ratio,top_iso_fraction,factor_length_ratio=factor_len_ratio,MLO_only=not use_MLF)
            save_temp_norm,temp_len_ratio_norm,save_temp_ids_norm=OU.get_all_homologies(ref_orga,orth_ref,label,specie_norm,save_homo,min_len_ratio,top_iso_fraction,factor_length_ratio=factor_len_ratio,MLO_only=not use_MLF)

            save_temp_top,temp_len_ratio_top,best_top_id=OU.get_top_x_pct(save_temp_top,temp_len_ratio_top,top_ortholog_fraction,names=save_temp_ids_top)
            save_temp_norm,temp_len_ratio_norm,best_norm_id=OU.get_top_x_pct(save_temp_norm,temp_len_ratio_norm,top_ortholog_fraction,names=save_temp_ids_norm)

            if not ((best_top_id is None ) or (best_norm_id is None)):
                # There are cases where the ortholog does not exist in one organism, so don't use
                AllSeqLab[label].add_entry(best_top_id,best_norm_id,orth_ref,save_temp_top,save_temp_norm,temp_len_ratio_top,temp_len_ratio_norm)

        print("Plotting "+label)
        AllSeqLab[label].sort_all()

        save_means[label]={}
        for at in annotation:
            save_means[label][at]={}
            for an in annotation[at]:
                compare,compare_top,compare_norm=AllSeqLab[label].get_compare_masked(at,an)
                OU.plot_2D_hists_wrapper(compare,compare_top,compare_norm,specie_top,specie_norm,ref_orga,binwidth,label,at+'_'+an)
                save_means[label][at][an]={}
                save_means[label][at][an][specie_top]=[np.mean(compare_top),np.std(compare_top)]
                save_means[label][at][an][specie_norm]=[np.mean(compare_norm),np.std(compare_norm)]

        compare,compare_top,compare_norm=AllSeqLab[label].get_compare()
        OU.plot_2D_hists_wrapper(compare,compare_top,compare_norm,specie_top,specie_norm,ref_orga,binwidth,label,"")


        i=0
        for at in annotation:
            len_tmp=len(annotation[at])+1
            len_tmp2=len_tmp+2
            width=1./(len_tmp2)
            offset=width*(len_tmp-1)/2.
            colors=[plt .cm.nipy_spectral(i/(len_tmp)) for i in range(len_tmp)]
            j=0
            for an in annotation[at]:

                plt.bar((j)/len_tmp2-offset,save_means[label][at][an][specie_top][0],yerr=save_means[label][at][an][specie_top][1],label=an,color=colors[j],width=width)
                plt.bar(1+(j)/len_tmp2-offset,save_means[label][at][an][specie_norm][0],yerr=save_means[label][at][an][specie_norm][1],color=colors[j],width=width)
                j+=1
            i+=1
            plt.bar((j)/len_tmp2-offset,np.mean(compare_top),yerr=np.std(compare_top),label='all',color=colors[j],width=width)
            plt.bar(1+(j)/len_tmp2-offset,np.mean(compare_norm),yerr=np.std(compare_norm),color=colors[j],width=width)
            plt.legend()
            xlabs=[specie_top,specie_norm]
            xticks=[0,1]
            plt.xticks(xticks,xlabs)
            plt.ylabel('Relative homology to '+ref_orga)
            OU.check_and_create_rep('Plots/'+label+'/Homology/')
            plt.savefig('Plots/'+label+'/Homology/'+at+'.pdf')



        # Getting the names of the genes from the file
        f=open(name_file)
        save_names=json.load(f)
        # Writing the sorted homology comparison
        W_tmp_ref=ref_orga+'_id\t'+ref_orga+'_name'
        W_tmp_top=specie_top+'_id\t'+specie_top+'_name\t'+'Isoform_pair\t'+specie_top+'_score'
        W_tmp_norm=specie_norm+'_id\t'+specie_norm+'_name\t'+'Isoform_pair\t'+specie_norm+'_score'

        W=W_tmp_ref+'\t'+W_tmp_top+'\t'+W_tmp_norm+'\t'+specie_top+'/'+specie_norm+'\n'

        # There is one line per reference gene, implying the reference gene is the same for each other organism, but
        # However, the reference isoform used may not be the same in the two organism, so we put it twice
        save_all_compare,save_all_ids_top,save_all_ids_norm,save_all_ids_ref,save_all_compare_top,save_all_compare_norm,save_all_len_ratio=AllSeqLab[label].get_all()
        for i in range(len(save_all_ids_ref)):
            W_tmp_ref=str(save_all_ids_ref[i])+'\t'+str(save_names[save_all_ids_ref[i]])

            split_tmp=save_all_ids_top[i].split('__')
            gene_id=split_tmp[0]
            W_tmp_top=str(gene_id)+'\t'+str(save_names[gene_id])+'\t'+split_tmp[1]+'\t'+str(save_all_compare_top[i])


            split_tmp=save_all_ids_norm[i].split('__')
            gene_id=split_tmp[0]
            W_tmp_norm=str(gene_id)+'\t'+str(save_names[gene_id])+'\t'+split_tmp[1]+'\t'+str(save_all_compare_norm[i])

            W+=W_tmp_ref+'\t'+W_tmp_top+'\t'+W_tmp_norm+'\t'+str(save_all_compare[i])+'\n'
        OU.write_file('Sorted_ids_'+label+'.txt',W)




