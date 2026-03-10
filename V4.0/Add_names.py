import json
import requests
import time
f=open("Orthology.json","r")
Orthology_all=json.load(f)

try :
    f=open("Names.json","r")
    Names_all=json.load(f)
except :
    Names_all={}

headers={"Content-Type":"application/json"}
server="https://rest.ensembl.org/"
N=0
for orths in Orthology_all:
    print(100*N/len(Orthology_all))
    N+=1
    for orga in Orthology_all[orths]:
        for orth_id in Orthology_all[orths][orga]:
            if not orth_id in Names_all:
                # Looking up names of orthologs
                ext="lookup/id/"+orth_id+"?expand=1"
                r1=requests.get(server+ext, headers=headers)

                N_wait=2
                while r1.status_code==429:
                    sleep_time=float(r1.headers.get("Retry-After", 1))
                    time.sleep(sleep_time+N_wait)
                    r1=requests.get(server+ext, headers=headers)  # Retry once
                try :
                    name_new=r1.json()['display_name']
                    Names_all[orth_id]=name_new
                except :
                    Names_all[orth_id]="N/A"

with open('Names.json','w') as f:
    json.dump(Names_all, f)