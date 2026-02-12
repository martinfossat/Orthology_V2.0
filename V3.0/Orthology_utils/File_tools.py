
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
