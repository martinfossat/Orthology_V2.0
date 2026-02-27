import matplotlib

matplotlib.use("pgf")
from matplotlib.backends.backend_pgf import FigureCanvasPgf

matplotlib.backend_bases.register_backend('pdf',FigureCanvasPgf)
pgf_with_latex={
    "text.usetex":True,
    "pgf.preamble":
        r'\usepackage{color}',
    "font.family":"Arial"}
matplotlib.rcParams.update(pgf_with_latex)
import matplotlib.pyplot as plt
import numpy as np
from .File_tools import *

def plot_2d_hist(x,y,xlabel,ylabel,x_ticks,y_ticks,name,binwidth):
    if len(x)==0 or len(y)==0:
        return None
    plt.figure(figsize=(4,4))
    plt.subplot2grid((8,8),(1,0),colspan=7,rowspan=7)

    bins_x=np.arange(0,np.nanmax(x)+binwidth,binwidth)
    bins_y=np.arange(0,np.nanmax(y)+binwidth,binwidth)

    hist=np.histogram2d(x,y,bins=[bins_x,bins_y])[0]

    plt.xticks(x_ticks,x_ticks)
    plt.yticks(y_ticks,y_ticks)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.plot([0,1],[0,1],linestyle='--',color='k')
    plt.imshow(hist,aspect='auto',origin='lower',interpolation='None',cmap='Spectral',extent=(0.,np.nanmax(y),0.,np.nanmax(x)))

    plt.subplot2grid((8,8),(1,7),colspan=1,rowspan=7)

    label_empty=['' for i in range(len(x_ticks))]

    plt.yticks(x_ticks,label_empty)
    plt.xticks([],[])

    plt.hist(x,bins=bins_x,orientation='horizontal',histtype='step',color='g')
    plt.ylim(0,np.nanmax(x))

    plt.subplot2grid((8,8),(0,0),colspan=7,rowspan=1)

    plt.yticks([],[])

    label_empty=['' for i in range(len(y_ticks))]

    plt.xticks(y_ticks,label_empty)

    plt.hist(y,bins=bins_y,histtype='step',color='g')
    plt.xlim(0,np.nanmax(y))

    plt.savefig(name)
    plt.close()


def plot_hist(save_all_prop,region_type,species,label_dic,bin_properties):
    colors=[plt.cm.gist_rainbow(i/(len(species))) for i in range(len(species))]
    type_ensemble=[]
    check_and_create_rep('Plots/'+region_type)
    for prop in save_all_prop:
        if not prop in bin_properties:
            continue
        if not prop in type_ensemble :
            type_ensemble+=[prop]
        j=0
        if not region_type in save_all_prop[prop]:
            continue
        for orga in save_all_prop[prop]:
            # if prop in bin_properties:
            temp_tuple=()
            #Need to remove and do explicitly, this nis not necessary
            for i in range(2):
                if bin_properties[prop][i] is None :
                    if i==0:
                        modulo=np.amin(save_all_prop[prop][orga])%bin_properties[prop][2]
                        temp_tuple+=(np.amin(save_all_prop[prop][orga])-modulo-bin_properties[prop][2],)

                    else :
                        modulo=np.amax(save_all_prop[prop][orga])%bin_properties[prop][2]
                        temp_tuple+=(np.amax(save_all_prop[prop][orga])+(1-modulo)+bin_properties[prop][2],)
                else :
                    if i==0:
                        temp_tuple+=(bin_properties[prop][i]-bin_properties[prop][2],)
                    else :
                        temp_tuple+=(bin_properties[prop][i]+bin_properties[prop][2],)
            temp_tuple+=(bin_properties[prop][2],)
            bins=np.arange(temp_tuple[0],temp_tuple[1],temp_tuple[2])

            bins_center=(bins[1:]+bins[:-1])/2.
            hist,bins_temp=np.histogram(save_all_prop[prop][orga],bins=bins)
            if len(hist)>0:
                hist=hist/np.amax(hist)
                plt.step(bins_center,hist,linewidth=1,color=colors[j],label=orga)
                j+=1
        if prop in label_dic:
            xlabel=label_dic[prop]
        else :
            xlabel=prop
        plt.xlabel(xlabel)
        plt.xlim(bin_properties[prop][0],bin_properties[prop][1])
        plt.legend()
        plt.savefig('Plots/'+region_type+"/Hist_ensemble_"+prop+'_'+region_type+'.pdf')
        plt.close()

def plot_bar_charts(save_all_prop,region_type,species,ensembles,name,label_dic_short):
    import matplotlib
    matplotlib.use("pgf")
    from matplotlib.backends.backend_pgf import FigureCanvasPgf
    matplotlib.backend_bases.register_backend('pdf',FigureCanvasPgf)
    pgf_with_latex={
        "text.usetex":True,
        "pgf.preamble":
            r'\usepackage{color}',
        "font.family":"Arial"}
    matplotlib.rcParams.update(pgf_with_latex)
    import matplotlib.pyplot as plt

    colors=[plt.cm.gist_rainbow(i/(len(species))) for i in range(len(species))]

    type_ensemble=[]

    i=0
    print('Plots/'+region_type)
    check_and_create_rep('Plots/'+region_type)
    N_offset=1
    for prop in save_all_prop:
        if not prop in ensembles:
            continue
        if not prop in type_ensemble:
            type_ensemble+=[prop]
        j=0
        for orga in save_all_prop[prop]:

            if i!=0:
                pre='_'
            else:
                pre=''
            pos=i-0.5+(j+1)/(len(species)+N_offset)
            if prop=='pct_folded':
                norm=100.
            else :
                norm=1.
            value=np.mean(np.array(save_all_prop[prop][orga])/norm)
            val_std=np.std(np.array(save_all_prop[prop][orga])/norm)
            plt.bar(pos,value,color=colors[j],width=1/(len(species)+N_offset),label=pre+orga,yerr=val_std,capsize=3)
            j+=1
        i+=1
    plt.title(name.replace('_',' '))
    plt.legend()
    xticks=[i for i in range(len(type_ensemble))]
    x_labels=[]
    for t in range(len(type_ensemble)):
        if type_ensemble[t] in label_dic_short:
            x_labels+=[label_dic_short[type_ensemble[t]]]
        else :
            x_labels+=[type_ensemble[t].replace("_"," ")]
    plt.xticks(xticks,x_labels,rotation=45,ha='right',va='top')
    plt.tight_layout()
    print('Plots/'+region_type+"/Sequence_properties")
    check_and_create_rep('Plots/'+region_type+"/Sequence_properties")
    print('Plots/'+region_type+"/Sequence_properties/"+name+".pdf")
    plt.savefig('Plots/'+region_type+"/Sequence_properties/"+name+".pdf")
    plt.close()

def plot_2D_hists_wrapper(compare,compare_top,compare_norm,specie_top,specie_norm,ref_orga,binwidth,seq_label,gene_label):
    check_and_create_rep('Plots/')
    check_and_create_rep('Plots/'+seq_label)
    check_and_create_rep('Plots/'+seq_label+'/Homology/')
    check_and_create_rep('Plots/'+seq_label+'/Homology/'+gene_label)
    ################ Plotting 2D histograms
    tick_sep=0.2
    x_ticks=np.round(np.arange(0,1+tick_sep,tick_sep),3)
    y_ticks=np.round(np.arange(0,1+tick_sep,tick_sep),3)
    ylabel='Normalized homology ratio '+specie_norm+' to '+ref_orga
    xlabel='Normalized homology ratio '+specie_top+' to '+ref_orga
    name='Plots/'+seq_label+'/Homology/'+gene_label+'/2D_hist_2_species.pdf'
    plot_2d_hist(compare_norm,compare_top,xlabel,ylabel,x_ticks,y_ticks,name,binwidth)

    x_ticks=np.round(np.arange(0,np.nanmax(compare_top)+tick_sep,tick_sep),3)
    y_ticks=np.round(np.arange(0,np.nanmax(compare)+tick_sep,tick_sep),3)
    xlabel='Normalized homology ratio '+specie_top+' to '+ref_orga
    ylabel='Normalized homology ratio ('+specie_top+'/'+specie_norm+')'
    name='Plots/'+seq_label+'/Homology/'+gene_label+'/2D_hist_top.pdf'
    plot_2d_hist(compare,compare_top,xlabel,ylabel,x_ticks,y_ticks,name,binwidth)

    x_ticks=np.round(np.arange(0,np.nanmax(compare_norm)+tick_sep,tick_sep),3)
    y_ticks=np.round(np.arange(0,np.nanmax(compare)+tick_sep,tick_sep),3)
    xlabel='Normalized homology ratio '+specie_norm+' to '+ref_orga
    ylabel='Normalized homology ratio ('+specie_top+'/'+specie_norm+')'
    name='Plots/'+seq_label+'/Homology/'+gene_label+'/2D_hist_norm.pdf'
    plot_2d_hist(compare,compare_norm,xlabel,ylabel,x_ticks,y_ticks,name,binwidth)

    # Plotting 1D histograms
    bins=np.arange(0,1+binwidth,binwidth)
    plt.title('Normalized homology to hsapiens')
    plt.xlim(0,1)
    plt.hist(compare_norm,bins=bins,histtype='step',color='orange',label=specie_norm)
    plt.hist(compare_top,bins=bins,histtype='step',color='g',label=specie_top)
    plt.legend()
    plt.ylabel("Number of occurences")
    plt.xlabel("Normalized homology to "+ref_orga)
    plt.savefig('Plots/'+seq_label+'/Homology/'+gene_label+'/Hist_both_species.pdf')
    plt.close()

    bins=np.arange(0,np.nanmax(compare)+binwidth,binwidth)
    plt.xlim(0,np.nanmax(compare)+binwidth)
    plt.xlim(0,2)
    plt.hist(compare,bins=bins,histtype='step',color='k')
    plt.ylabel("Number of occurences")
    plt.xlabel('Normalized homology to '+ref_orga+' ratio ('+specie_top+'/'+specie_norm+')')
    plt.savefig('Plots/'+seq_label+'/Homology/'+gene_label+'/Hist_ratio_'+specie_top+'_'+specie_norm+'.pdf')
    plt.close()
