# -*- coding: utf-8 -*-
import os

def listdir(path, list_name):  #传入存储的list
    for file in os.listdir(path):
        file_path = os.path.join(path, file)
        if os.path.isdir(file_path):
            listdir(file_path, list_name)
        else:
            file_path=file_path.replace("\\","/")
            file_path=file_path.replace("../r","http://www.hbpding.com/CLIPA/r")
            list_name.append(file_path)
            
res=[]
listdir('..\\r\\data',res)

f = open("./allRdatas.txt","w",encoding="utf-8")
f.write("\n".join(res))
f.flush()
f.close()