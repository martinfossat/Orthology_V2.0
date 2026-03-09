import json
import Orthology_utils as OU
import numpy as np

f=open('Name_dic.json')
Name_dic=json.load(f)

filename='Gene_names.txt'
raw_data=OU.read_file(filename)
data=np.array([raw_data[i][0] for i in range(len(raw_data))])
W_new=''
for i in range(len(data)):

    if data[i] in Name_dic :
        W_new+=Name_dic[data[i]]
    else :
        W_new+=data[i]

    W_new+='\n'
OU.write_file('Gene_names_new.txt',W_new)
print(W_new)