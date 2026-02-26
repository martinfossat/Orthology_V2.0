import json
import argparse
import Orthology_utils as OU

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

    print('Loading orthology file')
    f=open(orthology_file)
    orthology=json.load(f)

    print('Loading annotation file')
    f=open(annotation_file)
    Gene_annotation=json.load(f)

    lateral_transfer=True
    Nkeys=len(orthology.keys())
    i=0
    j=0
    for name in orthology:
        print("{:.2f}".format(100*float(i)/Nkeys)+" % done")
        i+=1
        orga_ref="hsapiens"
        if not orga_ref in orthology[name]:
            continue
        for orga in orthology[name]:
            #lateral transfer means if there is two orthologs in the reference species, do you transfer the annotations to one another
            if (not lateral_transfer) and orga==orga_ref:
                continue
            for orth_ref in orthology[name][orga_ref]:
                if orth_ref=="N/A":
                    continue

                if not orth_ref in Gene_annotation:
                    continue

                for orth in orthology[name][orga]:
                    if orth=="N/A":
                        continue

                    if not orth in Gene_annotation :
                        Gene_annotation[orth]={}

                    for lab in dict :
                        if not dict[lab] in Gene_annotation[orth]:
                            Gene_annotation[orth][dict[lab]]=[]
                        print(Gene_annotation[orth_ref][dict[lab]])
                        if Gene_annotation[orth_ref][dict[lab]] is None:
                            continue
                        for an in range(len(Gene_annotation[orth_ref][dict[lab]])):

                            if Gene_annotation[orth][dict[lab]] is None :
                                Gene_annotation[orth][dict[lab]]=[Gene_annotation[orth_ref][dict[lab]][an]]
                            elif not Gene_annotation[orth_ref][dict[lab]][an] in Gene_annotation[orth][dict[lab]]:
                                Gene_annotation[orth][dict[lab]]+=[Gene_annotation[orth_ref][dict[lab]][an]]
    print(j)
    print('Writing file')
    with open(annotation_file,'w') as f:
        json.dump(Gene_annotation,f)
