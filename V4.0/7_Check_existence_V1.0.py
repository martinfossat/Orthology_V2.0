import json
import argparse
import Orthology_utils as OU
import numpy as np


if __name__=="__main__":
    ################## Parser declaration ######################

    parser = argparse.ArgumentParser(description="""This program takes the orhtology database and checks which species have an orthologs""")
    parser.add_argument("--species","-s",help='Specie to check the existence of orthologs for',nargs='+',default=[])
    args = parser.parse_args()

    if args.species :
        species=args.species
    else :
        print("You must specify at least one specie")
        quit()

    f=open('Orthology.json')
    Orthology_all=json.load(f)

    Exists=np.zeros((len(Orthology_all),len(species)),dtype=int)
    all_any=np.zeros((len(Orthology_all),2),dtype=int)
    n=0
    W='Name\t'
    for specie in species :
        W+=specie+'\t'
    W+='all\tany\n'


    for orths in Orthology_all:
        tmp_list=[]
        W+=orths+'\t'
        for s in range(len(species)):
            orga=species[s]
            if orga in Orthology_all[orths]:
                if len(Orthology_all[orths][orga])==1 and 'N/A' in Orthology_all[orths][orga]:
                    # Exists[n][s]=0
                    W+=str(0)+'\t'

                else :
                    W+=str(1)+'\t'
                    Exists[n][s]=1
            else :
                W+=str(0)+'\t'
                Exists[n][s]=0
        all_any[n,0]=int(np.all(Exists[n,:]==1))
        all_any[n,1]=int(np.any(Exists[n,:]==1))
        W+=str(all_any[n,0])+'\t'+str(all_any[n,1])+'\n'
        n+=1
    W+='sum\t'
    for s in range(len(species)):
        W+=str(np.sum(Exists[:,s]))+'\t'
    W+=str(np.sum(all_any[:,0]))+'\t'+str(np.sum(all_any[:,1]))

    OU.write_file('Exists.txt',W)

    #print(Exists)



