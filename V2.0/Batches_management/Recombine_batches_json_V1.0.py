import json
import argparse
Homologies_rec={}
Names_rec={}
Seq_Prop_rec={}
Ortho_rec={}
if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Gets the predicted IDR regions, and predict their ensemble properties""")
    # parser.add_argument("--name_prefix","-np",help='Prefix of the file names ')
    # parser.add_argument("--name_suffix","-ns",help='Suffix of the file names ')
    parser.add_argument("--number_of_files", "-n",help='Number of files to recombine')
    # parser.add_argument("--outfile_name", "-on",help='Name of the output file ')
    args = parser.parse_args()
    # if args.name_prefix:
    #     prefix=args.name_prefix
    # else :
    #     prefix='./'
    #
    # if args.name_suffix:
    #     suffix=args.name_suffix
    # else :
    #     suffix='.json'
    #
    # if args.outfile_name:
    #     outfile_name=args.outfile_name
    # else:
    #     outfile_name='./Homology.json'

    if args.number_of_files:
        try :
            MAX_file=int(args.number_of_files)

        except :
            print("Wrong value for argument number_of_files")
            exit(1)
    else :
        print("Must specify number of files to recombine")
        exit(0)

    for i in range(0,MAX_file):
        # Homologies
        file_name='./Batches/Homology/N_'+str(i)+'.json'
        try :
            f=open(file_name)
            homologies=json.load(f)
        except :
            print('Could not open '+file_name)
            continue

        for ref_orga in homologies:
            if not ref_orga in Homologies_rec.keys():
                Homologies_rec[ref_orga]={}
            for orga in homologies[ref_orga]:
                if not orga in Homologies_rec[ref_orga].keys():
                    Homologies_rec[ref_orga][orga]={}
                for orth_homo_id in homologies[ref_orga][orga]:
                    Homologies_rec[ref_orga][orga][orth_homo_id]=homologies[ref_orga][orga][orth_homo_id]

        # Gene names
        file_name='./Batches/Names/N_'+str(i)+'.json'
        try :
            f=open(file_name)
            names=json.load(f)
        except :
            print('Could not open '+file_name)
            continue
        for n in names:
            Names_rec[n]=names[n]
        # Sequence properties
        file_name='./Batches/Seq_prop/N_'+str(i)+'.json'

        try :
            f=open(file_name)
            seq_prop=json.load(f)
        except :
            print('Could not open '+file_name)
            continue
        for n in seq_prop:
            Seq_Prop_rec[n]=seq_prop[n]

        # Orthology
        file_name='./Batches/Ortho/N_'+str(i)+'.json'
        try :
            f=open(file_name)
            ortho=json.load(f)
        except :
            print('Could not open '+file_name)
            continue
        for n in ortho:
            Ortho_rec[n]=ortho[n]

    with open('Homology.json','w') as f:
        json.dump(Homologies_rec,f)
        
    with open('Names.json','w') as f:
        json.dump(Names_rec,f)

    with open('Sequence_properties.json','w') as f:
        json.dump(Seq_Prop_rec,f)

    with open('Orthology.json','w') as f:
        json.dump(Ortho_rec,f)