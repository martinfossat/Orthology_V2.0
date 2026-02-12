import numpy as np
import requests
import json
import os
import argparse
import Orthology_utils as OU

server = "https://rest.ensembl.org/"

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Checks the orthologs for a file containing a list of gene names. The user needs to specify an original species name, corresponding to the gene name list, and a number of additional species may be used for 
    homology comparison, in subsequent algorithm.\nSpecies name must follow the id in gProfiler (https://biit.cs.ut.ee/gprofiler/page/organism-list). This program output an orthology file, a gene name file and a sequences files, all in json format.""")
    parser.add_argument("--gene_file","-f",help='Name of test file containing the name of the genes in the "original" species in the first column, other columns are ignored.')
    parser.add_argument("--original_specie","-os",help="Name of the species corresponding to the gene names provided in the input file. Must respect the naming convention in gProfiler (i.e. the id on https://biit.cs.ut.ee/gprofiler/page/organism-list ).This option can take multiple arguments.")
    parser.add_argument("--additional_species","-as",help="Additional species. Can take multiple arguments.",default=[],nargs='+')
    parser.add_argument("--orthology_file","-of",help='Name of the orthology database json output file. Default is Orthology.json')
    parser.add_argument("--sequences_file","-sf",help='Name of the sequences database json output file. Default is Sequences.json')
    parser.add_argument("--names_file","-nf",help='Name of the names database json output file. Default is Names.json')
    args = parser.parse_args()

    if args.gene_file :
        filename=args.gene_file
    else :
        print("You must specify an input file name")
        quit(1)

    if args.orthology_file :
        orthology_file=args.orthology_file

    else :
        orthology_file='Orthology.json'

    if args.sequences_file :
        sequences_file=args.sequences_file

    else :
        sequences_file='Sequences.json'

    if args.names_file :
        name_file=args.names_file

    else :
        name_file='Names.json'

    if args.original_specie :
        original_organism=args.original_specie
    else :
        print("You must specify an original specie")
        quit(1)

    if args.additional_species :
        other_organism=args.additional_species
    else :
        other_organism=[]

    dictionary_organisms_gprofiler=OU.get_orga_dic()
    dictionary_organisms_gprofiler_inv=dict(zip(dictionary_organisms_gprofiler.values(),dictionary_organisms_gprofiler.keys()))
    organisms_types=[]
    organisms_all=[original_organism]+other_organism
    ind_org=0

    for i in range(len(organisms_all)):
        if not dictionary_organisms_gprofiler_inv[organisms_all[i]] :
            print('One of the organisms does not match the name used in gProfiler, please only use gProfiler ID names.')
            quit(1)



    server="https://rest.ensembl.org/"

    raw_data=OU.read_file(filename)
    data=np.array([raw_data[i][0] for i in range(len(raw_data))])

    save_gene_name=np.empty((len(data)*4,len(organisms_all)),dtype="<U30")
    save_gene_id =np.empty((len(data)*4,len(organisms_all)),dtype="<U30")
    N_ids=np.zeros((len(data)*4),dtype=int)
    max_query_size=1000000000 # Just a large number, I don't know that there is actually a limit
    all_gene_names=np.ndarray.tolist(data[1:])
    MAX_same=0
    MAX_temp=0
    last_line_change=0
    save_prev_name=''
    save_prev_id=''

    Names_all={}
    Orthology_all={}

    print('Getting ortholog names')
    for t in range(1,len(organisms_all)):
        line=0
        N_q=0
        last_line_change=0
        print("Sending query to gProfiler for "+dictionary_organisms_gprofiler_inv[organisms_all[t]]+" ("+str(t)+"/"+str(len(organisms_all)-1)+")")

        #Gets the names from the input file
        gene_names=np.ndarray.tolist(data[N_q*max_query_size:min((N_q+1)*max_query_size,len(data[:]))])
        #get the ortholog names from the profiler website
        out=OU.query_profiler(organisms_all[ind_org],organisms_all[t],gene_names)
        print("Query received (size="+str(len(out['result']))+")")

        for i in range(len(out['result'])):
            name_org=out['result'][i]['incoming'].upper()
            id_org=out['result'][i]['converted'].upper()
            if out['result'][i]['name'] is not None:
                name_new=out['result'][i]['name'].upper()
            id_new=out['result'][i]['ortholog_ensg'].upper()
            N=out['result'][i]['n_result']

            # If that gene does not exist yet
            if name_org not in Orthology_all :
                ind=ind_org
                Orthology_all[name_org]={}
                Orthology_all[name_org][organisms_all[ind]]={}
                Orthology_all[name_org][organisms_all[ind]][id_org]=[]
                #Orthology_all[name_org][organisms_all[ind]][id_org]['name']=name_org
                Names_all[id_org]=name_org
            # If that organism in that gene does not exist yet
            if organisms_all[t] not in Orthology_all[name_org] :
                Orthology_all[name_org][organisms_all[t]]={}
            Orthology_all[name_org][organisms_all[t]][id_new]=[]
            #Orthology_all[name_org][organisms_all[t]][id_new]['name']=name_new
            Names_all[id_new]=name_new

    # Ok now we have all the orthologs, and we can check the sequences
    # For that we first get all gene ids, and we then ask the server, then we assign the protein in the json
    # This is to have only one query,and thus limit server usage

    queries_gene_id_list=[]
    queries_seq_list=[]
    queries_seq_id_list=[]

    all_ids=OU.get_all_ids(Orthology_all)
    LEN_MAX=50 # This si the maxx rom the ensemble rest API
    N=0
    #Turns out, the server would rather have lots of connection than big requests....
    while N*LEN_MAX<len(all_ids):
        print("Requesting sequences "+str(N*100*LEN_MAX/len(all_ids)))
        ext='sequence/id/'+'?multiple_sequences=1;type=protein'
        headers={"Content-Type":"application/json","Accept":"application/json"}
        r=requests.post(server+ext,headers=headers,data='{ "ids" : '+str(all_ids[N*LEN_MAX:(N+1)*LEN_MAX]).replace("'",'"')+' }')
        seqs_ref=r.json()
        N+=1
        # This is to avoid the no results coming back from the query (is a dict instead of a list)
        if type(seqs_ref)!=type([]):
            continue

        # Now to find efficiently the location of each of the output sequence
        # For this I can construct an array of tuples, where the tuple is the location of the
        for i in range(len(seqs_ref)):
            queries_gene_id_list += [seqs_ref[i]['query']]
            queries_seq_list += [seqs_ref[i]['seq']]
            queries_seq_id_list+=[seqs_ref[i]['id']]

    queries_gene_id_list=np.array(queries_gene_id_list)
    queries_seq_list=np.array(queries_seq_list)
    queries_seq_id_list=np.array(queries_seq_id_list)

    Sequences_all={}

    # Loop across gene names
    for orths in Orthology_all:
        # Loop across organisms
        for orga in Orthology_all[orths]:
            for gene_id in Orthology_all[orths][orga]:
                out=np.argwhere(gene_id==queries_gene_id_list)
                if len(out)>0:
                    inds=out[:,0]
                    for vers in inds :
                        seq_id=queries_seq_id_list[vers]
                        Sequences_all[seq_id]=queries_seq_list[vers]
                        Orthology_all[orths][orga][gene_id]+=[seq_id]


    with open(orthology_file,'w') as f:
        json.dump(Orthology_all,f)
    with open(name_file,'w') as f:
        json.dump(Names_all,f)
    with open(sequences_file,'w') as f:
        json.dump(Sequences_all,f)

