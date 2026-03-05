import requests
import Orthology_utils as OU
import json

import numpy as np
import time
absent=OU.read_file('Not_in_danio_rerio.txt')
server="https://rest.ensembl.org/"
headers={"Content-Type":"application/json"}
url = 'https://zfin.org/downloads/aliases.txt'
dic={}

response=requests.get(url)
response.raise_for_status()
data=response.text.splitlines()
for abs in range(len(absent)):
    print(100*abs/len(absent))
    for i in range(len(data)):
        found1=False
        found2=False
        tmp=np.array(data[i].split('\t'))
        for j in range(len(tmp)):
            if absent[abs][0].lower()==tmp[j].lower():
                found1=True
                break

        if found1 :
            for j in range(2,len(tmp)):
                ext="/homology/symbol/danio_rerio/"+tmp[j]+"?target_species=homo_sapiens;type=orthologues"
                r=requests.get(server+ext,headers=headers)
                N_wait=2
                while r.status_code==429:
                    sleep_time=float(r.headers.get("Retry-After", 1))
                    time.sleep(sleep_time+N_wait)
                    r=requests.get(server+ext, headers=headers)
                    N_wait+=1
                    print("Waiting for access")
                if not r.ok:
                    continue
                decoded=r.json()
                try :
                    if len(decoded['data'][0]['homologies'])>0:
                        found2=True
                        if absent[abs][0]!= tmp[j]:
                            dic[absent[abs][0]]=tmp[j]

                except :
                    print(r)
                    print(decoded)
                    print()
                    continue
    print(dic)
with open('Name_dic.json','w') as f:
    json.dump(dic,f)