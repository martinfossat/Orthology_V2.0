import requests
import Orthology_utils as OU
import json
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
        tmp=data[i].split('\t')
        for j in range(len(tmp)):
            if absent[abs][0].lower()==tmp[j].lower():
                found1=True
                break

        if found1 :
            for j in range(2,len(tmp)):
                # try :
                ext="/homology/symbol/danio_rerio/"+tmp[j]+"?target_species=homo_sapiens;type=orthologues"
                #ext="/homology/symbol/danio_rerio/"+tmp[j]+"?target_species=danio_rerio;type=paralogues"
                r=requests.get(server+ext,headers=headers)

                decoded=r.json()

                if not r.ok:
                    #r.raise_for_status()
                    continue

                # ext="lookup/id/"+id_new+"?expand=1"
                # r2=requests.get(server+ext,headers=headers)
                #
                # try :
                #     r2.raise_for_status()
                #     name_new=r2.json()['display_name']
                # except ConnectionError:
                #     name_new='N/A_(CE)'
                # except Timeout :
                #     name_new='N/A_(TO)'
                # except HTTPError :
                #     name_new='N/A_('+str(HTTPError.r2.status_code)+')'
                # except RequestException:
                #     name_new='N/A_('+str(RequestException)+')'
                # print(dic)
                # print()
                # print(tmp[j],absent[abs][0])
                # print(decoded['data'][0]['id'])
                # input(decoded['data'][0])

                if len(decoded['data'][0]['homologies'])>0:

                    found2=True
                    dic[absent[abs][0]]=tmp[j]
                    # print(dic)
                    # input(decoded['data'][0]['homologies'])
                # except :
                #     continue
    print(dic)
with open('Name_dic.json','w') as f:
    json.dump(dic,f)