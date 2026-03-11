import numpy as np
import requests
import json
import argparse
import Orthology_utils as OU
import time

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
    parser.add_argument("--division","-d",help='Division. Default is Vertebrates. Can be Fungi, Plants, Protists, Metazoa, Bacteria or Vertebrates')
    parser.add_argument("--list_species","-ls",help='Whether to print the list of species in the specified division.')
    args = parser.parse_args()

    if args.list_species :
        if args.list_species=='1':
            print_list=True
        else:
            print_list=False
    else :
        print_list=False

    if args.division :
        division=args.division
    else :
        division='Vertebrates'

    dictionary_organisms=OU.get_orga_dic_ensembl(division)
    dictionary_organisms_inv=dict(zip(dictionary_organisms.values(),dictionary_organisms.keys()))

    if print_list :
        print('List of species')
        for sp in dictionary_organisms:
            print(sp,'\t\t\t\t',dictionary_organisms[sp])
        quit()

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

    organisms_types=[]
    organisms_all=[original_organism]+other_organism

    for i in range(len(organisms_all)):
        if not organisms_all[i] in dictionary_organisms_inv :
            print('One of the organisms ('+organisms_all[i]+') does not match the name used in Ensembl, please only use Ensembl ID names, or specify a division (-d option).')
            quit(1)

    # Standard headers for Ensembl
    headers={"Content-Type":"application/json"}

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
    write_partial=True
    W_does_not_exit=''

    try :
        f=open(orthology_file)
        Orthology_all=json.load(f)
    except :
        Orthology_all={}

    Sequences_all={}


    headers={"Content-Type":"application/json"}
    server="https://rest.ensembl.org/"
    print('Getting ortholog names')
    for t in range(1,len(organisms_all)):
        for name_org in data:
            if name_org not in Orthology_all:
                Orthology_all[name_org]={}
                Orthology_all[name_org][original_organism]={}

            if organisms_all[t] not in Orthology_all[name_org] :
                Orthology_all[name_org][organisms_all[t]]={}

            # If not means we have already done that
            if len(Orthology_all[name_org][organisms_all[t]])!=0:
                print("Already have")
                print(name_org)
                # This is if there is already an entry : that means we have done this gene succesfully
                continue

            ext="/homology/symbol/"+original_organism+"/"+name_org+"?"+"target_species="+organisms_all[t]+";type=orthologues"+";compara="+division

            r1=requests.get(server+ext,headers=headers)

            N_wait=2
            while r1.status_code==429:
                sleep_time=float(r1.headers.get("Retry-After", 1))
                time.sleep(sleep_time+N_wait)
                r1=requests.get(server+ext, headers=headers)  # Retry once
                print("Max connection wait (Orthologs)")
                N_wait+=1
            if not r1.ok:
                if r1.status_code==400 :
                    W_does_not_exit+=name_org+'\n'
                else :
                    print(r1.status_code)
                continue

            decoded=r1.json()
            print(decoded)
            if len(decoded["data"])==0:
                continue

            # Here we are assuming the first gene is the best gene.
            # This may not be the case for a very small subset of proteins
            id_org=decoded["data"][0]['id']
            Orthology_all[name_org][original_organism][id_org]=[]
            for num in range(len(decoded["data"][0]['homologies'])):
                pct_id=decoded["data"][0]['homologies'][num]['target']['perc_id']
                prot_id_org=decoded["data"][0]['homologies'][num]['source']['protein_id']
                id_new=decoded["data"][0]['homologies'][num]['target']['id']
                prot_id_new=decoded["data"][0]['homologies'][num]['target']['protein_id']

                Orthology_all[name_org][organisms_all[t]][id_new]=[]

            if write_partial:
                with open(orthology_file, 'w') as f:
                    json.dump(Orthology_all, f)



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
        N_wait=2
        while r.status_code==429:
            sleep_time=float(r.headers.get("Retry-After", 1))
            time.sleep(sleep_time+N_wait)
            r=requests.get(server+ext, headers=headers)  # Retry once
            print("Max connection wait (Sequences)")
            N_wait+=1
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
                        if not seq_id in Orthology_all[orths][orga][gene_id]:
                            Orthology_all[orths][orga][gene_id]+=[seq_id]




    with open(orthology_file,'w') as f:
        json.dump(Orthology_all,f)

    with open(sequences_file,'w') as f:
        json.dump(Sequences_all,f)
