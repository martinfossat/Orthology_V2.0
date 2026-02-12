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


def plot_hist(save_all_prop,region_types,species,label_dic,bin_properties):
    colors=[plt.cm.gist_rainbow(i/(len(species))) for i in range(len(species))]
    type_ensemble=[]
    for region_type in region_types:
        check_and_create_rep('Plots/'+region_type)
        for prop in save_all_prop:
            if not prop in bin_properties:
                continue
            if not prop in type_ensemble :
                type_ensemble+=[prop]
            j=0
            if not region_type in save_all_prop[prop]:
                continue
            for orga in save_all_prop[prop][region_type]:
                # if prop in bin_properties:
                temp_tuple=()
                #Need to remove and do explicitly, this nis not necessary
                for i in range(2):
                    if bin_properties[prop][i] is None :

                        if i==0:
                            modulo=np.amin(save_all_prop[prop][region_type][orga])%bin_properties[prop][2]
                            temp_tuple+=(np.amin(save_all_prop[prop][region_type][orga])-modulo-bin_properties[prop][2],)

                        else :
                            modulo=np.amax(save_all_prop[prop][region_type][orga])%bin_properties[prop][2]
                            temp_tuple+=(np.amax(save_all_prop[prop][region_type][orga])+(1-modulo)+bin_properties[prop][2],)
                    else :
                        if i==0:
                            temp_tuple+=(bin_properties[prop][i]-bin_properties[prop][2],)
                        else :
                            temp_tuple+=(bin_properties[prop][i]+bin_properties[prop][2],)
                temp_tuple+=(bin_properties[prop][2],)
                bins=np.arange(temp_tuple[0],temp_tuple[1],temp_tuple[2])

                bins_center=(bins[1:]+bins[:-1])/2.
                hist,bins_temp=np.histogram(save_all_prop[prop][region_type][orga],bins=bins)
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

def plot_bar_charts(save_all_prop,region_types,species,ensembles,name,label_dic_short):
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
    for region_type in region_types:
        i=0
        check_and_create_rep('Plots/'+region_type)
        N_offset=1
        for prop in save_all_prop:
            if not prop in ensembles:
                continue
            if not prop in type_ensemble:
                type_ensemble+=[prop]
            j=0
            if not region_type in save_all_prop[prop]:
                continue
            for orga in save_all_prop[prop][region_type]:

                if i!=0:
                    pre='_'
                else:
                    pre=''
                pos=i-0.5+(j+1)/(len(species)+N_offset)
                if prop=='pct_folded':
                    norm=100.
                else :
                    norm=1.
                value=np.mean(np.array(save_all_prop[prop][region_type][orga])/norm)
                val_std=np.std(np.array(save_all_prop[prop][region_type][orga])/norm)
                plt.bar(pos,value,color=colors[j],width=1/(len(species)+N_offset),label=pre+orga,yerr=val_std,capsize=3)
                j+=1
            i+=1
        plt.title('Predicted ensemble properties')
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
        plt.savefig('Plots/'+region_type+"/Ensemble_properties_"+name+".pdf")
        plt.close()

