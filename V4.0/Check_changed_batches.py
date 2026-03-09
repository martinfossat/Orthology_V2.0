import json
import argparse
import Orthology_utils as OU
Homologies_rec={}
Names_rec={}
Seq_Prop_rec={}
Ortho_rec={}

orga_dic={"hsapiens": "homo_sapiens",
            "homo_sapiens":"homo_sapiens",
          "drerio": "danio_rerio",
          "danio_rerio" : "danio_rerio",
          "mmusculus": "mus_musculus",
          "mus_musculus":"mus_musculus"}

if __name__=="__main__":
    ################## Parser declaration ######################
    parser = argparse.ArgumentParser(description="""Gets the predicted IDR regions, and predict their ensemble properties""")
    parser.add_argument("--number_of_files", "-n",help='Number of files to recombine')
    args = parser.parse_args()

    f=open('Name_dic.json')
    dic_changed=json.load(f)

    pre='./Batches/Input_names/batch_'
    post='.txt'
    start=0
    end=2500
    W=''
    N=0
    for i in range(start,end):
        data=OU.read_file(pre+str(i)+post)
        for j in range(len(data)):
            print(data[j][0])
            if data[j][0] in dic_changed :
                W+=str(i)+'\n'
                N+=1
                break

    OU.write_file('To_redo.txt',W)
    print(len(dic_changed))
    print(N)
    exit()