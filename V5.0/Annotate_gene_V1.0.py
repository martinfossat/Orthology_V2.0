import json
import argparse
import Orthology_utils as OU
import requests
from requests.exceptions import HTTPError, ConnectionError, Timeout, RequestException

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""""")
    parser.add_argument("--orthology_file","-of",help='Name of the orthology database json input file. Default is Orthology.json')
    parser.add_argument("--annotation_file","-af",help='Name of the gene annotation json file. Default is Gene_annotations.json')
    parser.add_argument("--verbose","-v",help='Verbose (print percent progress), can be 1 (On) or 0 (Off). Default is 1.')
    args = parser.parse_args()

    if args.orthology_file:
        orthology_file=args.orthology_file
    else:
        orthology_file='Orthology.json'

    if args.annotation_file:
        annotation_file=args.annotation_file
    else:
        annotation_file='Gene_annotations.json'

    if args.verbose :
        try :
            if int(args.verbose)==1 :
                verbose=True
            elif int(args.verbose)==0 :
                verbose=False
            else :
                print("Wrong value for verbose. Can only be 1 or 0")
        except :
            print("Wrong format for verbose")

    else :
        verbose=True

    orga_ref="hsapiens"
    Gene_annotation={}
    dict=OU.get_gene_annotation_label_dic()
    print('Loading file')
    f=open(orthology_file)
    orthology=json.load(f)

    Nkeys=len(orthology.keys())
    i=0
    for name in orthology:
        print("{:.2f}".format(100*float(i)/Nkeys)+" % done")
        i+=1
        orga_ref="hsapiens"
        if not orga_ref in orthology[name]:
            continue

        for orth in orthology[name][orga_ref]:
            if orth=="N/A":
                continue
            if orth in Gene_annotation :
                continue
            url='https://www.proteinatlas.org/'+orth+'.json'
            try:
                r2=requests.post(url=url,timeout=10)
                r2.raise_for_status()

            except HTTPError as http_err:
                print(f"HTTP error occurred: {http_err}")  # e.g., 404 Not Found
                continue
            except ConnectionError as conn_err:
                print(f"Error connecting: {conn_err}")  # e.g., DNS failure or refused connection
                continue
            except Timeout as timeout_err:
                print(f"Timeout error: {timeout_err}")  # The server took too long to respond
                continue
            except RequestException as req_err:
                print(f"An ambiguous error occurred: {req_err}")
                continue
            out=r2.json()
            Gene_annotation[orth]={}
            for lab in dict:
                Gene_annotation[orth][dict[lab]]=out[lab]

    print('Writing file')
    with open(annotation_file,'w') as f:
        json.dump(Gene_annotation,f)
