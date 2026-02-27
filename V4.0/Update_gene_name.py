import requests
import Orthology_utils as OU



def get_current_symbol(species, old_symbol):
    server = "https://rest.ensembl.org"
    # This endpoint finds all Ensembl objects linked to the symbol
    ext = f"/xrefs/symbol/{species}/{old_symbol}?"

    r = requests.get(server + ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        return None

    decoded = r.json()
    print(decoded)
    # Usually, the first result for 'gene' type is what you need
    for entry in decoded:
        if entry['type'] == 'gene':
            return entry['id'] # This is the ENSG ID (The "Truth")

data=OU.read_file('Gene_names.txt')

for i in range(len(data)):
    old=data[i][0]
    new=get_current_symbol("danio_rerio",old)
    print(old,new)
    input()

