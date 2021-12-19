import os
import shutil
import copy
import json
from PIL import Image

def parsehelp(folder,catalog):
    help_tmp={}
    files=os.listdir(folder)
    id=0
    for f in files:
        filepath=folder+"/"+f
        if not os.path.isfile(filepath):
            print("intodir"+filepath)
            catalog_new=copy.deepcopy(catalog)
            catalog_new.append(f)
            help_tmp[f]=parsehelp(filepath,catalog_new)
        elif "_title.txt" in f:
            print("intotitle:"+filepath)
            help_tmp["00_title"]=readFile(filepath)[0]
        elif "_plain.txt" in f:
            print("intoplain:"+filepath)
            help_tmp[f.split(".")[0]]=readFile(filepath)
        elif "_italic.txt" in f:
            print("intoitalic:"+filepath)
            help_tmp[f.split(".")[0]]=readFile(filepath)
        elif "_bold.txt" in f:
            print("intobold:"+filepath)
            help_tmp[f.split(".")[0]]=readFile(filepath)
        elif "_table.txt" in f:
            print("intobold:"+filepath)
            help_tmp[f.split(".")[0]]=readFile(filepath)
        elif ".png" in f:
            print("intopng:"+filepath)
            
            figname="-".join(catalog)+"_"+str(Image.open(filepath).size[0])+"_"+str(id)+".png"
            id+=1
            shutil.copy(filepath,"../help/"+figname)
            help_tmp[f.split(".")[0]+"_png"]=figname
        
    return(help_tmp)



def readFile(filename):
    with open(filename, 'r') as f:
        return(f.readlines())
oldfiles=os.listdir('../help/')
for of in oldfiles:
    os.remove("../help/"+of)

catalog=[]
help=parsehelp("../myScripts/help",catalog)

f = open("../help/index.json","w",encoding="utf-8")
f.write(json.dumps(help))
f.flush()
f.close()

print("执行完毕，检查一下有没有报错，按回车退出")
input()