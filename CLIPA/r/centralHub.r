
requestGUID<-commandArgs(T)[1]

result<-list()
start<-Sys.time()
fit<-try({
	setwd("../r")
	library(rjson)
	#write.table(userParams_raw,file="userParams_raw.txt")
	source("W_functions.R")
	source("L_functions.R")
	userParams<-fromJSON(file=paste0("../plotCache/userParams_",requestGUID,".json"))
	#requestGUID<-userParams$requestGUID
	jsonpath<-userParams$path
	result$Model<-userParams$Model
	result$requestGUID<-userParams$requestGUID
	#write(toJSON(userParams),file="../plotCache/userParams_debug.json")

	
	drugs<-userParams$foldList_drugs
	#sldkjls<-errortest
	gene<-userParams$singleGeneInput
	
	

	if(userParams$use=="jsShow"){
		res<-processData(userParams)
		save(res,file=paste0("../plotCache/data_jsShow_",userParams$plotGUID,".rdata"))
		result$res<-res
		result$code=0
	}else if(userParams$use=="download"){
		filetype<-userParams$filetype
		if(userParams$Model!="DepMap,Fusions"){
			if(filetype==".pdf"){
				FigFunc<-"pdf"
			}else if(filetype==".eps"){
				FigFunc<-"postscript"
			}else if(filetype==".tiff"){
				FigFunc<-"tiff"
			}else if(filetype==".svg"){
				FigFunc<-"svg"
			}
			Figpath<-paste0("../plotCache/",userParams$plotGUID,filetype)
			urlpath<-paste0("./plotCache/",userParams$plotGUID,filetype)
			data_jsShow<-paste0("../plotCache/data_jsShow_",userParams$plotGUID,".rdata")
			if(!file.exists(data_jsShow)){
				res<-processData(userParams)
			}else{
				load(data_jsShow)
			}
			FigOutput(res,FigFunc,Figpath)
		}else{
			cellLines<-userParams$foldList_cellLines
			urlpath<-paste0("./r/data/DepMap/Fusion/",cellLines,"_Fusion",filetype)
		}
		result$Figpath<-urlpath
		result$code=0
		
	}else if(userParams$use=="checkGenes"){
		result$res<-checkGenes(userParams$genes)
	}

	
})

result$debug$Delay<-as.numeric(difftime(Sys.time(),start,unit="sec"))
if('try-error' %in% class(fit)){
	co<-attr(fit,"condition")
	#ver<-shell("Rscript --version",intern=TRUE)
	result$error<-paste0(co[1],"<br>----------------<br>Detailed error:<br>",co[2],"<br>--------------<br>")
	
}



write(toJSON(result),file=jsonpath)
print("done")
