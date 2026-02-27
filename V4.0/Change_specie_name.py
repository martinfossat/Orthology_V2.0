import json
import argparse
import Orthology_utils as OU

orga_dic={"hsapiens": "homo_sapiens",
            "homo_sapiens":"homo_sapiens",
          "drerio": "danio_rerio",
          "danio_rerio" : "danio_rerio",
          "mmusculus": "mus_musculus",
          "mus_musculus":"mus_musculus"}

if __name__=="__main__":
    ################## Parser declaration ######################
    parser=argparse.ArgumentParser(
        description="""Gets the predicted IDR regions, and predict their ensemble properties""")
    parser.add_argument("--number_of_files", "-n", help='Number of files to recombine')
    args=parser.parse_args()

    if args.number_of_files:
        try:
            MAX_file=int(args.number_of_files)

        except:
            print("Wrong value for argument number_of_files")
            exit(1)
    else:
        print("Must specify number of files to recombine")
        exit(0)
    OU.check_and_create_rep('./Batches_new')
    OU.check_and_create_rep('./Batches_new/Homology')
    OU.check_and_create_rep('./Batches_new/Ortho')
    for i in range(MAX_file):
        Homologies_rec={}
        Ortho_rec={}
        file_name='./Batches/Homology/N_'+str(i)+'.json'

        try:
            f=open(file_name)
            homologies=json.load(f)
        except:
            print('Could not open '+file_name)
            continue

        for ref_orga in homologies:
            if not orga_dic[ref_orga] in Homologies_rec.keys():
                Homologies_rec[orga_dic[ref_orga]]={}
            for orga in homologies[ref_orga]:
                if not orga_dic[orga] in Homologies_rec[orga_dic[ref_orga]].keys():
                    Homologies_rec[orga_dic[ref_orga]][orga_dic[orga]]={}
                for orth_homo_id in homologies[ref_orga][orga]:
                    Homologies_rec[orga_dic[ref_orga]][orga_dic[orga]][orth_homo_id]=homologies[ref_orga][orga][orth_homo_id]

        with open('./Batches_new/Homology/N_'+str(i)+'.json', 'w') as f:
            json.dump(Homologies_rec, f)

        # Orthology
        file_name='./Batches/Ortho/N_'+str(i)+'.json'
        try:
            f=open(file_name)
            ortho=json.load(f)
        except:
            print('Could not open '+file_name)
            continue

        for n in ortho:
            Ortho_rec[n]={}
            for orga1 in ortho[n]:
                Ortho_rec[n][orga_dic[orga1]]=ortho[n][orga1]

        with open('./Batches_new/Ortho/N_'+str(i)+'.json', 'w') as f:
            json.dump(Ortho_rec, f)