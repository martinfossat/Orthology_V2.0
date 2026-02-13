import os


def check_and_create_rep(directory,silent=True):
    try :
        os.makedirs(directory)
    except :
        if os.path.exists(directory) and not silent:
            print('Could not create directory ('+directory+') : already exists')
        if not os.path.exists(directory):
            print('Error creating directory '+directory)
def write_file(file_name,W,silent=True):
    try :
        with open(file_name,'w') as f :
            f.writelines(W)
            f.close()
        if not silent :
            print(file_name+" written.")
    except :
        print("Could not open "+file_name)


def load_file(file_name,silent=False):
    try :
        with open(file_name,'r') as f :
            data=f.readlines()
            f.close()
    except :
        if silent==False:
            print("Could not open "+file_name)
    return data
def read_file(file_name,silent=False,split='') :
    data=load_file(file_name,silent=silent)
    len1=len(data)
    arr=[[] for i in range(len1)]
    for i in range(len1):
        if split!='':
            temp=data[i].split(split)
        else :
            temp=data[i].split()
        arr[i]+=temp
    return arr
