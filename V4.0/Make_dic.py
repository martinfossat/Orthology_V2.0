import json


def get_orga_dic_ensembl():
    import requests
    server="https://rest.ensembl.org"
    ext="/info/species?"
    r=requests.get(server+ext,headers={"Content-Type":"application/json"})

    species_list=[s['name'] for s in r.json()['species']]
    name_list=[n['display_name'] for n in r.json()['species']]

    dic_out={}
    for i,j in zip(species_list,name_list):
        dic_out[j]=i

    with open('Species_Ensembl.json','w') as f:
        json.dump(dic_out,f)


    return dic_out
get_orga_dic_ensembl()

