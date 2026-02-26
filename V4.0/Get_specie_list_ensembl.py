import requests
server = "https://rest.ensembl.org"
ext = "/info/species?"
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

species_list = [s['name'] for s in r.json()['species']]
name_list=[n['display_name'] for n in r.json()['species']]
W='{'
for i,j in zip(species_list,name_list):
    W+='\t"'+i+'" : "'+j+'",\n'
W=W[:-2]+'}'
print(W)