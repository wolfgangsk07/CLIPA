library(rjson)


get_userParams<-function(userParams_raw){
	userParams_raw<-gsub("[+]","%20",userParams_raw)
	#write(userParams_raw,file="debug.txt")
	userParams<-fromJSON(URLdecode(userParams_raw))
	#write(URLdecode(userParams_raw),file="../userParams.json")
	#userParams<-fromJSON(file="../userParams.json")
}


get_Drugs<-function(userParams){
	return(names( userParams$foldList_drugs)[ unlist(userParams$foldList_drugs)])
}

get_cellLines<-function(userParams){
	return(names( userParams$foldList_cellines)[ unlist(userParams$foldList_cellines)])
}


processData<-function(userParams){
	Dataset<-strsplit(userParams$Model,",")[[1]][1]
	Model<-strsplit(userParams$Model,",")[[1]][2]
	if(userParams$Model %in% c("GDSC,Drug-Omics","NCI60,Drug-Omics","DepMap,Drug-Omics")){
		drugs<-userParams$foldList_drugs
		gene<-userParams$singleGeneInput
		cor_method<-userParams$cor_method
		mutation_method<-userParams$mutation_method
		MHT<-userParams$MHT
		pcutoff<-userParams$pcutoff
		cellLines<-userParams$foldList_cellLines
		#if(userParams$omicsSelector %in% c("mRNA","methylation_promoter","methylation_body","methylation_site","CNV","mutation")){
			res<-Drug_Omics_response(Dataset=Dataset,Omics=userParams$omicsSelector,gene=gene,drug_input=drugs,cell_lines=cellLines,
				statistic_method=MHT,cut_off=pcutoff,cor_method=cor_method,mutation_method=mutation_method)
		#}else{
		#	res<-"error"
		#}
	}else if(userParams$Model %in% c("GDSC,Drug-Pathway","NCI60,Drug-Pathway","DepMap,Drug-Pathway")){
		drugs<-userParams$foldList_drugs
		cellLines<-userParams$foldList_cellLines
		cor_method<-userParams$cor_method
		mutation_method<-userParams$mutation_method
		MHT<-userParams$MHT
		pcutoff<-userParams$pcutoff
		if(userParams$foldList_geneSets$type=="MSigDB"){
			MSigDB_Name<-userParams$foldList_geneSets$MSigDB_Name
			res<-Drug_Omics_response(Dataset=Dataset,Omics="GSVA",drug_input=drugs,cell_lines=cellLines,geneset=MSigDB_Name,statistic_method=MHT,cut_off=pcutoff,cor_method=cor_method,mutation_method=mutation_method,preset=TRUE)
		}else if(userParams$foldList_geneSets$type=="customGeneList"){
			gene<-userParams$foldList_geneSets$genes
			res<-Drug_Omics_response(Dataset=Dataset,Omics="GSVA",drug_input=drugs,cell_lines=cellLines,gene=gene,statistic_method=MHT,cut_off=pcutoff,cor_method=cor_method,mutation_method=mutation_method,preset=FALSE)
		}else{
			res<-list()
			res$error<-"unknown gsva type"
		}
		
		
	}else if(userParams$Model %in% c("GDSC,Omics-Omics (cis-regulation)","NCI60,Omics-Omics (cis-regulation)")){
		gene<-userParams$singleGeneInput
		cellLines<-userParams$foldList_cellLines
		cor_method<-userParams$cor_method
		mutation_method<-userParams$mutation_method
		MHT<-userParams$MHT
		pcutoff<-userParams$pcutoff
		res<-MultiOmics_calculation(Dataset=Dataset,gene=gene,cell_lines=cellLines,
			statistic_method=MHT,cut_off=pcutoff,cor_method=cor_method,mutation_method=mutation_method)
	}else if(userParams$Model %in% c("DepMap,Omics-Omics (cis-regulation)")){
		gene<-userParams$singleGeneInput
		cellLines<-userParams$foldList_cellLines
		cor_method<-userParams$cor_method
		mutation_method<-userParams$mutation_method
		MHT<-userParams$MHT
		pcutoff<-userParams$pcutoff
		res<-DepMap_MultiOmics_calculation(Dataset=Dataset,gene=gene,cell_lines=cellLines,
		statistic_method=MHT,cut_off=pcutoff,cor_method=cor_method,mutation_method=mutation_method)
	}else if(userParams$Model %in% c("GDSC,Cancer Subtypes","NCI60,Cancer Subtypes","DepMap,Cancer Subtypes")){
		genes<-userParams$genes
		omics<-userParams$omics
		cellLines<-userParams$foldList_cellLines
		clusterAlgorithm<-userParams$clusterAlgorithm
		groupNum<-as.numeric(userParams$groupNum)
		res<-Subcluster_calculation(Dataset=Dataset,Omics=omics,genesets=genes,cell_lines=cellLines,Num_clust=groupNum,cluster_algorithm=clusterAlgorithm,txtPath=paste0("../plotCache/clust_result_",userParams$plotGUID,".csv"))
	}else if(userParams$Model %in% c("GDSC,Cancer Types Summary","NCI60,Cancer Types Summary","DepMap,Cancer Types Summary")){
		if(userParams$omicsSelector %in% c("mRNA","methylation_promoter","methylation_body","methylation_site","protein","microRNA","CNV_continuous","Chromatin","Metabolomics","Metastatic","CRISPR","RNAi_combine")){
			gene<-userParams$singleGeneInput
			res<-Boxplot_calculation(Dataset=Dataset,gene=gene,Omics=userParams$omicsSelector)
		}else if(userParams$omicsSelector %in% c("CNV","mutation")){
			gene<-userParams$singleGeneInput
			res<-SummBarplot_calculation(Dataset=Dataset,gene=gene,Omics=userParams$omicsSelector)
		}
	}else if(userParams$Model %in% c("GDSC,Omics-Omics (trans-regulation)","NCI60,Omics-Omics (trans-regulation)","DepMap,Omics-Omics (trans-regulation)")){
			Omics<-userParams$multiOmicsSelector$Omics
			Genes<-userParams$multiOmicsSelector$Genes
			Rank<-as.numeric(userParams$multiOmicsSelector$Rank)
			Omics<-Omics[order(Rank,decreasing=FALSE)]
			Genes<-Genes[order(Rank,decreasing=FALSE)]
			cellLines<-userParams$foldList_cellLines
			reference<-c(Genes[1],Omics[1])
			cor_method<-userParams$cor_method
			mutation_method<-userParams$mutation_method
			MHT<-userParams$MHT
			pcutoff<-userParams$pcutoff
			res<-Trans_Omics_calculation(Dataset=Dataset,genes_sets=Genes,Omics=Omics,cell_lines=cellLines,reference=reference,
				statistic_method=MHT,cut_off=pcutoff,cor_method=cor_method,mutation_method=mutation_method)

	}else if(userParams$Model %in% c("GDSC,Machine Learning in Predicting Drug Response")){
			Omics<-userParams$omics
			Genes<-userParams$genes
			cellLines<-userParams$foldList_cellLines
			drugs<-userParams$foldList_drugs
			Algorithm<-userParams$Algorithm
			
			res<-Drug_response_model_calculation(Dataset=Dataset,genesets=Genes,Omics=Omics,cell_lines=cellLines,drug_input=drugs,Algorithm=Algorithm)
	
	}else if(userParams$Model %in% c("DepMap,Fusions")){

			cellLines<-userParams$foldList_cellLines

     
			res<-list()
			res$commomParams$rightMargin<-0
     
			res$commomParams$leftMargin<-0
			res$commomParams$width<-7*300
			res$commomParams$DocHeight<-7*300
			res$plotElements$svg$plotType<-"drawRawSVG"
			header='<svg id="Fusionsvg" width="100%" height="100%" xmlns="http://www.w3.org/2000/svg" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:svgjs="http://svgjs.com/svgjs" text-rendering="optimizeLegibility" viewBox="0 0 504 504" style="background:white">'
			svgRaw<-readLines(paste0("./data/DepMap/Fusion/",cellLines,"_Fusion.svg"))[-1]
			svgRaw[1]<-header
			res$plotElements$svg$data<-paste(svgRaw,collapse="")
	}else{
		res<-list()
		res$error<-"这个模块还没写好"
	}
	return(res)
}

checkGenes<-function(genes_raw){
	
	res<-list()
	if(FALSE){
	genes_raw<-gsub(",",", ",genes_raw)
	genes_raw<-gsub(",","\n",genes_raw)
	genes_raw<-gsub(",","\t",genes_raw)
	genes<-strsplit(genes_raw,",")[[1]]
	load("./data/checkGenes.RData")
	
	if(length(genes)==1){
		if(genes %in% checkGenes$ourGenes){
			res$text<-"The submitted gene has passed the inspection!"
			res$color<-"green"
		}else if(genes %in% checkGenes$officalGenes){
			res$text<-"The submitted gene out of our range!"
			res$color<-"red"
		}else{
			res$text<-"The input gene is invalid!"
			res$color<-"red"
		}
	}else if(length(genes)>1){
		invalidGenes<-c()
		outRangeGenes<-c()
		for(i in genes){
			if(i %in% checkGenes$ourGenes){
				next
			}else if(i %in% checkGenes$officalGenes){
				outRangeGenes<-c(outRangeGenes,i)
			}else{
				invalidGenes<-c(invalidGenes,i)
			}
		}
		if(length(outRangeGenes)==0&length(invalidGenes)==0){
			res$text<-"All the submitted genes have passed the inspection!"
			res$color<-"green"
		}else if(length(outRangeGenes)==0&length(invalidGenes)!=0){
		
			res$text<-paste("The following gene(s) are invalid:\n",paste(invalidGenes,collapse=", "))
			res$color<-"red"
		}else if(length(outRangeGenes)!=0&length(invalidGenes)==0){
			res$text<-paste("The following gene(s) are out of our range:\n",paste(invalidGenes,collapse=", "))
			res$color<-"red"
		}else if(length(outRangeGenes)!=0&length(invalidGenes)!=0){
			res$text<-paste("The following gene(s) are invalid:\n",paste(invalidGenes,collapse=", "),"\nThe following gene(s) are out of our range:\n",paste(invalidGenes,collapse=", "))
			res$color<-"red"
		}
	}
	
	}
	res$text<-"The submitted gene has passed the inspection!"
	res$color<-"green"
	return(res)
}

#结果转图形参数(小弟负责）=================

heatmap_translate<-function(calcData=NULL){
	#计算矩阵数据
	if(is.null(calcData)){
		colorPool<-colorRampPalette(c("red","black","green"))(200)
		calcData<-list()
		calcData$data[["class1"]]$data<-sapply(1:30,function(x)sample(colorPool,20))
		calcData$data[["class1"]]$name<-"class1"
		calcData$data[["class1"]]$color<-"yellow"
		calcData$data[["class1"]]$cellNames<-paste0(sample(names(precip),30),"_",round(rnorm(1,100,100),0))
		calcData$data[["class2"]]$data<-sapply(1:5,function(x)sample(colorPool,20))
		calcData$data[["class2"]]$name<-"class2"
		calcData$data[["class2"]]$color<-"blue"
		calcData$data[["class2"]]$cellNames<-paste0(sample(names(precip),5),"_",round(rnorm(1,100,100),0))
		calcData$data[["class3"]]$data<-sapply(1:10,function(x)sample(colorPool,20))
		calcData$data[["class3"]]$name<-"class3"
		calcData$data[["class3"]]$color<-"grey"
		calcData$data[["class3"]]$cellNames<-paste0(sample(names(precip),10),"_",round(rnorm(1,100,100),0))
		calcData$geneNames<-sample(names(precip),20)
		clust_result<-data.frame(clust_result=c(rep(calcData$data[["class1"]]$name,ncol(calcData$data[["class1"]]$data)),rep(calcData$data[["class2"]]$name,ncol(calcData$data[["class2"]]$data)),rep(calcData$data[["class3"]]$name,ncol(calcData$data[["class3"]]$data))))
		rownames(clust_result)<-c(calcData$data[["class1"]]$cellNames,calcData$data[["class2"]]$cellNames,calcData$data[["class3"]]$cellNames)
		calcData$table<-clust_result
		calcData$txtPath<-"../plotCache/clust_result.txt"
		calcData$color_range<-c("#1C577E","white","#B51F23")
	}
	
	dpi<-300
	#save(calcData,file="calcData.rda")
	#可调节图像参数
	xlength<-sum(sapply(calcData$data,function(x)ncol(x$data)))
	if(xlength>50){
		rectWidth<-2000/xlength
	}else{
		rectWidth<-40
	}
	
	ylength<-nrow(calcData$data[[1]]$data)
	if(ylength>20){
		rectHeigth<-800/ylength
	}else{
		rectHeigth<-40
	}
	if(rectWidth<20){
		showcellNames<-FALSE
	}else{
		showcellNames<-TRUE
	}
	
	if(rectHeigth<25){
		showgeneNames<-FALSE
	}else{
		showgeneNames<-TRUE
	}
	
	classInterval<-10

	classTitleHeight<-30
	classTitleBottomInterval<-5
	margin<-5
	colBarWidth<-30
	colBarHeight<-100
	colBarMargin<-30
	#定义公共画图参数======================================================================================
	if(showgeneNames){
		geneNames_width<-max(sapply(calcData$geneNames,function(x)getStrLen(x,0.5)))*dpi
	}else{
		geneNames_width<-0
	}
	if(showcellNames){
		cellNames_height<-max(sapply(unlist(sapply(calcData$data,function(c)c$cellNames)),function(x)getStrLen(x,0.5)))*dpi
	}else{
		cellNames_height<-0
	}
     res<-list()
	 res$commomParams$rightMargin<-margin
     
	 res$commomParams$leftMargin<-margin
     res$commomParams$width<-sum(sapply(calcData$data,function(x)ncol(x$data)))*rectWidth+classInterval*(length(calcData$data)-1)+margin+geneNames_width+colBarWidth+colBarMargin*2
	 plotHeight<-max(sapply(calcData$data,function(x)nrow(x$data)))*rectHeigth+classTitleBottomInterval+classTitleHeight
	 res$commomParams$DocHeight<-plotHeight+margin*2+cellNames_height

	#获得每个方块参数
	classLeft<-0
	for(i in 1:length(calcData$data)){
		
		for(j in 1:nrow(calcData$data[[i]]$data)){
			for(k in 1:length(calcData$data[[i]]$data[j,])){
				res$plotElements[[paste(i,j,k,sep="_")]]$plotType="drawRect"
				res$plotElements[[paste(i,j,k,sep="_")]]$color<-calcData$data[[i]]$data[j,k]
				res$plotElements[[paste(i,j,k,sep="_")]]$left<-classLeft+(k-1)*rectWidth+geneNames_width
				res$plotElements[[paste(i,j,k,sep="_")]]$bottom<-res$commomParams$DocHeight-j*rectHeigth-classTitleBottomInterval-classTitleHeight-margin
				res$plotElements[[paste(i,j,k,sep="_")]]$width<-rectWidth
				res$plotElements[[paste(i,j,k,sep="_")]]$height<-rectHeigth
				res$plotElements[[paste(i,j,k,sep="_")]]$className<-"justHover"
				
					
				if(k==1&showgeneNames){ #定义基因名
					res$plotElements[[paste("gene",i,j,k)]]$plotType="drawText"
					res$plotElements[[paste("gene",i,j,k)]]$text=calcData$geneNames[j]
					res$plotElements[[paste("gene",i,j,k)]]$x=geneNames_width-margin
					res$plotElements[[paste("gene",i,j,k)]]$y=res$commomParams$DocHeight-j*rectHeigth-classTitleBottomInterval-classTitleHeight-margin+rectHeigth/2
					res$plotElements[[paste("gene",i,j,k)]]$cex=0.5
					res$plotElements[[paste("gene",i,j,k)]]$adj=c(1,0.5)
					res$plotElements[[paste("gene",i,j,k)]]$color="black"
					res$plotElements[[paste("gene",i,j,k)]]$srt=0
				}
				
				if(j==1&showcellNames){ #定义细胞名
					res$plotElements[[paste("cell",i,j,k)]]$plotType="drawText"
					res$plotElements[[paste("cell",i,j,k)]]$text=calcData$data[[i]]$cellNames[k]
					res$plotElements[[paste("cell",i,j,k)]]$x=classLeft+(k-1)*rectWidth+rectWidth/2+geneNames_width
					res$plotElements[[paste("cell",i,j,k)]]$y=cellNames_height
					res$plotElements[[paste("cell",i,j,k)]]$cex=0.5
					res$plotElements[[paste("cell",i,j,k)]]$adj=c(1,0.5)
					res$plotElements[[paste("cell",i,j,k)]]$color="black"
					res$plotElements[[paste("cell",i,j,k)]]$srt=90
				}
			}
		
		
		}
		classLeft<-classLeft+ncol(calcData$data[[i]]$data)*rectWidth+classInterval
	}
	#定义抬头============
	classLeft<-0
	for(i in 1:length(calcData$data)){
		res$plotElements[[paste0(calcData$data[[i]]$name,i)]]$plotType="drawRect"
		res$plotElements[[paste0(calcData$data[[i]]$name,i)]]$color<-calcData$data[[i]]$color
		res$plotElements[[paste0(calcData$data[[i]]$name,i)]]$left<-classLeft+geneNames_width
		res$plotElements[[paste0(calcData$data[[i]]$name,i)]]$bottom<-res$commomParams$DocHeight-classTitleHeight-margin
		res$plotElements[[paste0(calcData$data[[i]]$name,i)]]$width<-ncol(calcData$data[[i]]$data)*rectWidth
		res$plotElements[[paste0(calcData$data[[i]]$name,i)]]$height<-classTitleHeight
		
		
		res$plotElements[[paste0(calcData$data[[i]]$name,"text",i)]]$plotType="drawText"
		res$plotElements[[paste0(calcData$data[[i]]$name,"text",i)]]$text=calcData$data[[i]]$name
		res$plotElements[[paste0(calcData$data[[i]]$name,"text",i)]]$x=classLeft+ncol(calcData$data[[i]]$data)*rectWidth/2+geneNames_width
		res$plotElements[[paste0(calcData$data[[i]]$name,"text",i)]]$y=res$commomParams$DocHeight-classTitleHeight/2-margin
		res$plotElements[[paste0(calcData$data[[i]]$name,"text",i)]]$cex=0.5
		res$plotElements[[paste0(calcData$data[[i]]$name,"text",i)]]$adj=c(0.5,0.5)
		res$plotElements[[paste0(calcData$data[[i]]$name,"text",i)]]$color="white"
		res$plotElements[[paste0(calcData$data[[i]]$name,"text",i)]]$srt=0
		
		
		classLeft<-classLeft+ncol(calcData$data[[i]]$data)*rectWidth+classInterval
	}
	#画colbar
		
		res$plotElements$colBar$plotType="drawColBar"
		res$plotElements$colBar$text=c("High","Low")
		res$plotElements$colBar$color_range<-calcData$color_range
		res$plotElements$colBar$left<-classLeft+colBarMargin+geneNames_width
		res$plotElements$colBar$bottom<-res$commomParams$DocHeight-colBarHeight-colBarMargin*2
		res$plotElements$colBar$width<-colBarWidth
		res$plotElements$colBar$height<-colBarHeight
	#定义表格========
	classLeft<-0
	cellNames_list<-list()
	lineHeight<-35
	tableInnerMargin<-10
	downloadBarHeight<-60
	titleWidth<-150
	res$table$downloadBarHeight<-40
	res$table$downloadBarWidth<-330
	res$table$downloadBarLeft<-res$commomParams$leftMargin
	res$table$downloadBarBaseline<--downloadBarHeight+5
	res$table$url<-calcData$txtPath
	contentWidth<-res$commomParams$width-margin-titleWidth
	for(i in 1:length(calcData$data)){
		
		cellNames<-calcData$data[[i]]$cellNames
	
		
		index<-1
		cellNames_collapsed<-c()
		while(length(cellNames)!=0){
			cellNames_tmp<-c()
			for(c in 1:length(cellNames)){
				cellNames_tmp<-c(cellNames_tmp,cellNames[c])
				if(getStrLen(paste(cellNames_tmp,collapse=", "),0.5)*dpi>=contentWidth-tableInnerMargin*2){
					if(c==1){
						cellNames<-cellNames[-1]
					}else{
						cellNames_tmp<-cellNames_tmp[1:(length(cellNames_tmp)-1)]
						cellNames<-cellNames[c:length(cellNames)]
					}
					
					cellNames_collapsed<-c(cellNames_collapsed,paste(cellNames_tmp,collapse=", "))
					index<-index+1
					break
				}
				if(c==length(cellNames)){
					
					cellNames<-c()
					cellNames_collapsed<-c(cellNames_collapsed,paste(cellNames_tmp,collapse=", "))
					index<-index+1
				}
			}
		}
		cellNames_collapsed[-length(cellNames_collapsed)]<-paste0(cellNames_collapsed[-length(cellNames_collapsed)],",")
		cellNames_list[[i]]<-cellNames_collapsed
	}
	if(!showcellNames){
		res$table$height<-sum(sapply(cellNames_list,length)*lineHeight+tableInnerMargin*2)+margin+downloadBarHeight
		
		res$table$plotElements[["tableOutline"]]$plotType="drawRect"
		#res$table$plotElements[["tableOutline"]]$color<-"white"
		res$table$plotElements[["tableOutline"]]$border<-"black"
		res$table$plotElements[["tableOutline"]]$left<-0
		res$table$plotElements[["tableOutline"]]$bottom<-margin-res$table$height
		res$table$plotElements[["tableOutline"]]$width<-res$commomParams$width-margin
		res$table$plotElements[["tableOutline"]]$height<-res$table$height-margin-downloadBarHeight
		
		res$table$plotElements[["vline"]]$ys<-c(margin-res$table$height,-downloadBarHeight)
		res$table$plotElements[["vline"]]$xs<-rep(titleWidth,2)
		res$table$plotElements[["vline"]]$plotType<-"drawPolyline"
		res$table$plotElements[["vline"]]$color<-"black"
		res$table$plotElements[["vline"]]$lwd<-0.5
		
		innerBase<--downloadBarHeight
		for(i in 1:length(cellNames_list)){
			if(i!=length(cellNames_list)){
				res$table$plotElements[[paste("innerline",i)]]$ys<-rep(innerBase-length(cellNames_list[[i]])*lineHeight-tableInnerMargin*2,2)
				res$table$plotElements[[paste("innerline",i)]]$xs<-c(0,res$commomParams$width-margin)
				res$table$plotElements[[paste("innerline",i)]]$plotType<-"drawPolyline"
				res$table$plotElements[[paste("innerline",i)]]$color<-"black"
				res$table$plotElements[[paste("innerline",i)]]$lwd<-0.5
			}
			res$table$plotElements[[paste("className",i,j,sep="_")]]$plotType="drawText"
			res$table$plotElements[[paste("className",i,j,sep="_")]]$text=names(calcData$data)[i]
			res$table$plotElements[[paste("className",i,j,sep="_")]]$x=titleWidth/2
			res$table$plotElements[[paste("className",i,j,sep="_")]]$y=innerBase-length(cellNames_list[[i]])*lineHeight/2-tableInnerMargin
			res$table$plotElements[[paste("className",i,j,sep="_")]]$cex=0.5
			res$table$plotElements[[paste("className",i,j,sep="_")]]$adj=c(0.5,0.5)
			res$table$plotElements[[paste("className",i,j,sep="_")]]$color="black"
			res$table$plotElements[[paste("className",i,j,sep="_")]]$srt=0
			for(j in 1:length(cellNames_list[[i]])){
				res$table$plotElements[[paste("table",i,j,sep="_")]]$plotType="drawText"
				res$table$plotElements[[paste("table",i,j,sep="_")]]$text=cellNames_list[[i]][j]
				res$table$plotElements[[paste("table",i,j,sep="_")]]$x=titleWidth+tableInnerMargin
				res$table$plotElements[[paste("table",i,j,sep="_")]]$y=innerBase-lineHeight*j-tableInnerMargin
				res$table$plotElements[[paste("table",i,j,sep="_")]]$cex=0.5
				res$table$plotElements[[paste("table",i,j,sep="_")]]$adj=c(0,1)
				res$table$plotElements[[paste("table",i,j,sep="_")]]$color="black"
				res$table$plotElements[[paste("table",i,j,sep="_")]]$srt=0
			}
			innerBase<-innerBase-length(cellNames_list[[i]])*lineHeight-tableInnerMargin*2
		}
	}else{
		res$table$height<-margin+downloadBarHeight
	}
	#save(res,file="res")
	return(res)
	
}
plot_translate<-function(calcData=NULL){
	if(is.null(calcData)){
		calcData<-list()
		calcData$type$a$x<--10:10*25
		calcData$type$a$y<-rep(50,21)
		calcData$type$a$shape<-"circle"
		calcData$type$a$color<-"green"
		calcData$type$a$size<-10
		
		calcData$type$d$x<-rep(50,21)
		calcData$type$d$y<--10:10*25
		calcData$type$d$shape<-"rect"
		calcData$type$d$color<-"red"
		calcData$type$d$size<-10
		
		calcData$type$c$x<-rnorm(200,0,1)*100
		calcData$type$c$y<-rnorm(200,0,1)*100
		calcData$type$c$shape<-"circle"
		calcData$type$c$color<-"blue"
		calcData$type$c$size<-10
		
		calcData$xlab="xlab"
		calcData$ylab="ylab"
	}

	 #定义公共画图参数======================================================================================
	marginTop<-10
	marginBotoom<-100
	plotHeight<-1000
	plotWidth<-1000
    res<-list()
	res$commomParams$rightMargin<-20
	res$commomParams$leftMargin<-100
    res$commomParams$width<-plotWidth+res$commomParams$rightMargin
	res$commomParams$DocHeight<-plotHeight+marginBotoom+marginTop
	#定义PLOT参数
	xmin<-0
	xmax<-0
	ymax<-0
	ymin<-0
	
	for(i in 1:length(calcData$type)){
		xmin<-min(c(xmin,calcData$type[[i]]$x))
		xmax<-max(c(xmax,calcData$type[[i]]$x))
		ymin<-min(c(ymin,calcData$type[[i]]$y))
		ymax<-max(c(ymax,calcData$type[[i]]$y))
	}
	print(paste(xmin,xmax,ymin,ymax))
	yseg<-ceiling((ymax-ymin)/100)*100/5
	new_ymin<-floor(ymin/yseg)*yseg
	new_ymax<-ceiling(ymax/yseg)*yseg
	xseg<-ceiling((xmax-xmin)/100)*100/5
	new_xmin<-floor(xmin/xseg)*xseg
	new_xmax<-ceiling(xmax/xseg)*xseg
	print(paste(new_xmin,new_xmax,new_ymin,new_ymax))
	xscale<-plotWidth/(new_xmax-new_xmin)
	yscale<-plotHeight/(new_ymax-new_ymin)
	for(i in 1:length(calcData$type)){
		res$plotElements[[paste("plot",i)]]$plotType<-"drawPlot"
		res$plotElements[[paste("plot",i)]]$x<-(calcData$type[[i]]$x-new_xmin)*xscale
		res$plotElements[[paste("plot",i)]]$y<-(calcData$type[[i]]$y-new_ymin)*yscale+marginBotoom
		res$plotElements[[paste("plot",i)]]$color<-calcData$type[[i]]$color
		res$plotElements[[paste("plot",i)]]$shape<-calcData$type[[i]]$shape
		res$plotElements[[paste("plot",i)]]$size<-calcData$type[[i]]$size
	}
	#定义坐标轴
	
	
	res$plotElements$axis$plotType<-"drawAxis"
	res$plotElements$axis$x$labs<-(new_xmin/xseg+1):(new_xmax/xseg+1)*xseg
	res$plotElements$axis$x$ats<-(1:6)*plotWidth/6
	res$plotElements$axis$x$startEnd<-c(0,plotWidth)
	res$plotElements$axis$x$baseline<-marginBotoom
	res$plotElements$axis$x$srt<-0
	
	res$plotElements$axis$y$labs<-(new_ymin/yseg+1):(new_ymax/yseg+1)*yseg
	res$plotElements$axis$y$ats<-1:6*plotHeight/6+marginBotoom
	res$plotElements$axis$y$startEnd<-c(marginBotoom,plotHeight+marginBotoom)
	res$plotElements$axis$y$x<-0
	res$plotElements$axis$y$srt<-0
	return(res)
}

SummBarplot_translate<-function(BarData=NULL){
	#计算数据=======================
	if(is.null(BarData)){
		BarData<-list()
		
		BarData$Type<-sample(c(LETTERS,letters),28)
		for(i in 1:length(BarData$Type)){
			BarData$data[[BarData$Type[i]]]<-abs(rnorm(4,0,0.1))
		}
		BarData$colors<-rainbow(4)
		BarData$anno<-c("NOR","LOSS","GAIN","AMP")
		BarData$Y_text<-"test"
	}
	#save(BarData,file="BarData.rda")
	Type<-BarData$Type
	maxValue<-max(sapply(BarData$data,function(x)sum(sapply(x,as.numeric))))
	fitedY<-autoFitAxis(0,maxValue*100)
	Scale<-100/tail(fitedY,1)

	 #定义公共画图参数======================================================================================
	 dpi<-300
	 topMargin<-20
	 if(!is.null(BarData$legend)){
		legendHeight<-20
	 }else{
		legendHeight<-0
	 }
     res<-list()
	 res$commomParams$rightMargin<-10
     res$commomParams$DocHeight<-1000
     res$commomParams$width<-2000
	 
	 baseline<-max(sapply(Type,function(x)getStrLen(x,0.5)))*dpi/2^0.5+20
     plotHeight<-res$commomParams$DocHeight-baseline-topMargin-legendHeight
	 
	 barWidthVSboxInterval<-0.5
     barWidth<-res$commomParams$width/length(Type)*barWidthVSboxInterval
	 res$commomParams$leftMargin<-max(c(max(sapply(Type,function(x)getStrLen(x,0.5))*dpi/2^0.5-(1:length(Type)-1)*barWidth/barWidthVSboxInterval),getStrLen("100%",0.5)*dpi+60))
	#定义坐标轴===========================================
	res$plotElements$axis$plotType<-"drawAxis"
	res$plotElements$axis$x$labs<-Type
	res$plotElements$axis$x$ats<-c(rep(0,length(Type)))
	res$plotElements$axis$x$startEnd<-c(0,res$commomParams$width)
	res$plotElements$axis$x$baseline<-baseline
	res$plotElements$axis$x$srt<-45
	
	res$plotElements$axis$y$labs<-paste0(fitedY,"%")
	res$plotElements$axis$y$ats<-fitedY/tail(fitedY,1)*plotHeight+baseline
	res$plotElements$axis$y$startEnd<-c(baseline,plotHeight+baseline)
	res$plotElements$axis$y$x<-0
	res$plotElements$axis$y$srt<-0
	res$plotElements$axis$y$ylab<-BarData$Y_text
	if(getStrLen(BarData$Y_text,0.75)*dpi/2>plotHeight/2){
		res$plotElements$axis$y$ylabRescale<-(plotHeight/2)/(getStrLen(BarData$Y_text,0.75)*dpi/2)
	}
	#处理每个柱子
	
	for (i in Type) {
		res$plotElements[[i]]$popupText$cell_lines_text<-list()
		for(t in names(BarData$data[[i]])){
			res$plotElements[[i]]$data[[t]]=BarData$data[[i]][[t]]*plotHeight*Scale
			res$plotElements[[i]]$color[[t]]<-BarData$colors[[i]][[t]]
			res$plotElements[[i]]$popupText$cell_lines_text[[t]] <- BarData$cell_lines_text[[i]][[t]]
			if(!is.null(BarData$Sub_Text[[i]])){
				res$plotElements[[i]]$popupText$Sub_Text[[t]] <- BarData$Sub_Text[[i]][[t]]
			}
		}
		
		res$plotElements[[i]]$barWidth=barWidth
		res$plotElements[[i]]$x=which(Type==i)*barWidth/barWidthVSboxInterval-barWidth/barWidthVSboxInterval/2
		res$plotElements$axis$x$ats[which(Type==i)]<-res$plotElements[[i]]$x
		res$plotElements[[i]]$plotType = "SummBarplot"
		res$plotElements[[i]]$baseline=baseline
		
		res$plotElements[[i]]$popupText$main<-BarData$Main_Text[[i]]
		
		
	}
	if(!is.null(BarData$legend)){
		res$plotElements$legend<-legend_translate(unlist(BarData$legend),res$commomParams$DocHeight-legendHeight-topMargin,res$commomParams$width-100)
	}#save(BarData,file="../plotCache/BarData.rda")
	#res$warnings<-c("I don't known what to warn for!","This is also a warning test!")
	 return(res)
}
legend_translate<-function(legendData=NULL,baseline,right){
	if(is.null(legendData)){
		legendData<-rainbow(4)
		names(legendData)<-c("Deletion","Loss","Gain","Amplification")
	}
	dpi<-300
	#legend默认右对齐
	rectSize<-20
	InnerMargin<-5
	rect_textMargin<-2
	totalLength<-sum(sapply(names(legendData),function(x)getStrLen(x,0.5)*dpi))+length(legendData)*(rectSize+InnerMargin+rect_textMargin)
	leftStart<-right-totalLength
	legend<-list()
	offset<-0
	for(i in 1:length(legendData)){
		legend$legend[[names(legendData)[i]]]$baseline<-baseline
		legend$legend[[names(legendData)[i]]]$size<-rectSize
		legend$legend[[names(legendData)[i]]]$color<- as.character(legendData[i])
		legend$legend[[names(legendData)[i]]]$left<-leftStart+offset
		legend$legend[[names(legendData)[i]]]$text<-names(legendData)[i]
		offset<-offset+rectSize+InnerMargin+rect_textMargin+getStrLen(names(legendData)[i],0.5)*dpi
	}
	legend$plotType<-"drawLegend"
	return(legend)
	


}

modelRoc_translate<-function(modelRocData=NULL){
	if(is.null(modelRocData)){
		load("../myScripts/modelrocDEMO.RData")
		modelRocData<-list()
		modelRocData$specificities<-modelroc$specificities
		modelRocData$sensitivities<-modelroc$sensitivities
		modelRocData$ci<-as.numeric(modelroc$ci)
		modelRocData$model<-"NNet"
	}
	#定义公共画图参数======================================================================================
	dpi<-300
	topMargin<-12
	botomMargin<-130
	leftMargin<-130
	rightMargin<-12
    res<-list()
	res$commomParams$rightMargin<-rightMargin
    res$commomParams$DocHeight<-600+botomMargin
    res$commomParams$width<-600
	res$commomParams$leftMargin<-leftMargin
	baseline<-botomMargin
	#定义坐标轴===========================================
	res$plotElements$axis$plotType<-"drawAxis"
	res$plotElements$axis$x$labs<-format(0:5/5,nsmall=1)
	res$plotElements$axis$x$ats<-0:5/5*(res$commomParams$width-rightMargin)
	res$plotElements$axis$x$startEnd<-c(0,res$commomParams$width-rightMargin)
	res$plotElements$axis$x$baseline<-baseline
	res$plotElements$axis$x$srt<-0
	res$plotElements$axis$x$xlab<-"1-Specificity (FPR)"
	
	res$plotElements$axis$y$labs<-format(0:5/5,nsmall=1)
	res$plotElements$axis$y$ats<-0:5/5*(res$commomParams$DocHeight-topMargin-botomMargin)+baseline
	res$plotElements$axis$y$startEnd<-c(baseline,res$commomParams$DocHeight-topMargin)
	res$plotElements$axis$y$x<-0
	res$plotElements$axis$y$srt<-0
	res$plotElements$axis$y$ylab<-"Sensitivity (TPR)"
	res$plotElements$axis$box<-TRUE
	#开始作图
	res$plotElements$polyline$xs<-round((1-modelRocData$specificities)*(res$commomParams$width-rightMargin),1)
	res$plotElements$polyline$ys<-round(modelRocData$sensitivities*(res$commomParams$DocHeight-topMargin-botomMargin)+baseline,1)
	res$plotElements$polyline$plotType<-"drawPolyline"
	res$plotElements$polyline$color<-"black"
	res$plotElements$polyline$lwd<-1.5
	res$plotElements$polyline$className<-"justHover"
	
	res$plotElements$crossline$ys<-c(0,1)*(res$commomParams$DocHeight-topMargin-botomMargin)+baseline
	res$plotElements$crossline$xs<-c(0,1)*(res$commomParams$width-rightMargin)
	res$plotElements$crossline$plotType<-"drawPolyline"
	res$plotElements$crossline$color<-"black"
	res$plotElements$crossline$lwd<-0.5
	res$plotElements$crossline$lty<-"dashed"
	
	if(as.numeric(modelRocData$ci[2])>=0.5){
		textX<-0.6
		textY<-0.2
	}else{
		textX<-0.4
		textY<-0.8
	}
	
	res$plotElements[["text"]]$plotType="drawText"
	res$plotElements[["text"]]$text=paste(modelRocData$model,"AUC:", sprintf("%0.3f",modelRocData$ci[2])," 95%CI:", sprintf("%0.3f",modelRocData$ci[1]),"-",sprintf("%0.3f",modelRocData$ci[3]))
	res$plotElements[["text"]]$x=textX*(res$commomParams$width-rightMargin)
	res$plotElements[["text"]]$y=textY*(res$commomParams$DocHeight-topMargin-botomMargin)+baseline
	res$plotElements[["text"]]$cex=0.4
	res$plotElements[["text"]]$adj=c(0.5,0.5)
	res$plotElements[["text"]]$color="black"
	res$plotElements[["text"]]$srt=0
					
					
	return(res)

}
boxplot_translate<-function(calcData=NULL){
	if(is.null(calcData)){
		load("../myScripts/Single_Omics.rda")
		calcData<-list()
		calcData$data<-Single_Omics 
		calcData$ylab<-"llsjdkljfsldjfl"
	}
	res<-list()
  dpi<-300
  Single_Omics<-calcData$data
  if((getQuantile(mean(Single_Omics[,2]),Single_Omics[,2])>0.7|getQuantile(mean(Single_Omics[,2]),Single_Omics[,2])<0.3)&min(Single_Omics[,2])>0){
	calcData$ylab<-paste(calcData$ylab,"(log2)")
	Single_Omics[,2]<-log(Single_Omics[,2],2)
  
  }
  Type <- unique(Single_Omics$cancer_type)
  calc<-list()
  min_expr<-c()
  max_expr<-c()
  colorPool<-c('#E41A1C', '#A73C52', '#6B5F88', '#3780B3', '#3F918C', '#47A266',
    '#53A651', '#6D8470', '#87638F', '#A5548D', '#C96555', '#ED761C',
    '#FF9508', '#FFC11A', '#FFEE2C', '#EBDA30', '#CC9F2C', '#AD6428',
    '#BB614F', '#D77083', '#F37FB8', '#DA88B3', '#B990A6', '#999999')
  #colors<-rainbow(length(Type))
  colors<-colorPool[sapply(1:length(Type)-1,function(x)(x %% length(colorPool)))+1]
  Q2_forOder<-c()
  for (i in Type) {
    Single_Omics_tmp <- Single_Omics[Single_Omics$cancer_type == i,]
    if(length(grep("-",names(Single_Omics_tmp)[2]))==1|
       length(grep("\\/",names(Single_Omics_tmp)[2]))==1|
       length(grep("\\(",names(Single_Omics_tmp)[2]))==1|
       length(grep("\\:",names(Single_Omics_tmp)[2]))==1){
      fmla <- as.formula(paste0("`",names(Single_Omics_tmp)[2],"`","~","cancer_type"))
    }else{
      fmla <- as.formula(paste(names(Single_Omics_tmp)[2],"~","cancer_type"))
    }
    boxplot_result <- boxplot(fmla,Single_Omics_tmp,plot=FALSE)
	calc[[i]]$orgData$data<-as.numeric(Single_Omics_tmp[,2])
	calc[[i]]$orgData$cellNames<-rownames(Single_Omics_tmp)
    calc[[i]]$abnormal <- as.numeric(boxplot_result$out)
    boxplot_result <- as.numeric(boxplot_result$stats)
    calc[[i]]$Q3 <- boxplot_result[4]
    calc[[i]]$Q2 <- boxplot_result[3]
    Q2_forOder<-c(Q2_forOder,calc[[i]]$Q2)
    calc[[i]]$Q1 <- boxplot_result[2]
    calc[[i]]$min_observation <- boxplot_result[1]
    calc[[i]]$max_observation <- boxplot_result[5]
    min_expr<-min(c(min_expr,calc[[i]]$min_observation))
    max_expr<-max(c(max_expr,calc[[i]]$max_observation))
  }
  fitedY<-autoFitAxis(min_expr,max_expr)
  Type<-Type[order(Q2_forOder)]
  calc<-calc[order(Q2_forOder)]
  #定义公共画图参数======================================================================================
  innerMargin<-20
  
  res$commomParams$rightMargin<-10
  res$commomParams$DocHeight<-1000
  #res$commomParams$width<-2000
  boxWidth<-60-length(Type)
  baseline<-max(sapply(Type,function(x)getStrLen(x,0.5)))*dpi/2^0.5+20
  plotHeight<-res$commomParams$DocHeight-baseline-innerMargin
  
  boxWidthVSboxInterval<-0.6
  #boxWidth<-(res$commomParams$width-res$commomParams$rightMargin)/length(Type)*boxWidthVSboxInterval
  res$commomParams$width<-boxWidth/boxWidthVSboxInterval*length(Type)+res$commomParams$rightMargin
  res$commomParams$leftMargin<-max(c(130,max(sapply(Type,function(x)getStrLen(x,0.5))*dpi/2^0.5-(1:length(Type)-1)*boxWidth/boxWidthVSboxInterval)))
  #定义坐标轴===========================================
  res$plotElements$axis$plotType<-"drawAxis"
  res$plotElements$axis$x$labs<-Type
  res$plotElements$axis$x$ats<-c(rep(0,length(Type)))
  res$plotElements$axis$x$startEnd<-c(0,res$commomParams$width-res$commomParams$rightMargin)
  res$plotElements$axis$x$baseline<-baseline
  res$plotElements$axis$x$srt<-45
  res$plotElements$axis$y$labs<-fitedY
  res$plotElements$axis$y$ats<-0:(length(fitedY)-1)/(length(fitedY)-1)*plotHeight+baseline
  res$plotElements$axis$y$startEnd<-c(baseline,plotHeight+baseline)
  res$plotElements$axis$y$x<-0
  res$plotElements$axis$y$srt<-0
  res$plotElements$axis$y$ylab<-calcData$ylab
  if(getStrLen(calcData$ylab,0.75)*dpi/2>plotHeight/2){
		res$plotElements$axis$y$ylabRescale<-(plotHeight/2)/(getStrLen(calcData$ylab,0.75)*dpi/2)
	}
  #处理每个箱子的参数
  for (i in Type) {
    Qs<-c(calc[[i]]$min_observation,calc[[i]]$Q1,calc[[i]]$Q2,calc[[i]]$Q3,calc[[i]]$max_observation)
    Qs<-(Qs-fitedY[1])/(tail(fitedY,1)-fitedY[1])*plotHeight
	
    res$plotElements[[i]]$Qs=Qs
    #res$plotElements[[i]]$abnormal= (calc[[i]]$abnormal-min_expr)/(max_expr-min_expr)*plotHeight
    res$plotElements[[i]]$boxWidth=boxWidth
    res$plotElements[[i]]$x=which(Type==i)*boxWidth/boxWidthVSboxInterval-boxWidth/boxWidthVSboxInterval/2
    res$plotElements$axis$x$ats[which(Type==i)]<-res$plotElements[[i]]$x
    res$plotElements[[i]]$plotType = "Boxplot"
    res$plotElements[[i]]$baseline=baseline
    res$plotElements[[i]]$color<-colors[which(Type==i)]
	res$plotElements[[i]]$plotData$ys=round((calc[[i]]$orgData$data-fitedY[1])/(tail(fitedY,1)-fitedY[1])*plotHeight)
	res$plotElements[[i]]$plotData$xs=round(sample(1:100,length(calc[[i]]$orgData$data),replace=TRUE)/100*boxWidth+res$plotElements[[i]]$x-boxWidth/2)
	res$plotElements[[i]]$plotData$cellNames=calc[[i]]$orgData$cellNames
	res$plotElements[[i]]$plotData$r=boxWidth/3.5
	res$plotElements[[i]]$plotData$width=boxWidth/boxWidthVSboxInterval
  }
	return(res)
}
#画图子函数(小弟负责）===================
drawPolyline<-function(plotElement){
	if(!is.null(plotElement$lty)){
		lty<-plotElement$lty
	}else{
		lty<-"solid"
	}
	points(plotElement$xs,plotElement$ys,type="l",lwd=plotElement$lwd,col=plotElement$color,lty=lty)
}
drawPlot<-function(plotElement){
	if(plotElement$shape=="circle"){
		require(plotrix)
		for(i in 1:length(plotElement$x)){
			draw.circle(plotElement$x[i],plotElement$y[i],plotElement$size/2,10,col=plotElement$color,border=NA)
		}
	
	}else if(plotElement$shape=="rect"){
		for(i in 1:length(plotElement$x)){
			rect(plotElement$x[i],plotElement$y[i],plotElement$x[i]+plotElement$size,plotElement$y[i]+plotElement$size,col=plotElement$color,border=NA)
		}
	}
}

SummBarplot<-function(plotElement){
	subBaseline<-0
	for(i in 1:length(plotElement$data)){
		rect(plotElement$x-plotElement$barWidth/2,plotElement$baseline+subBaseline,plotElement$x+plotElement$barWidth/2,plotElement$baseline+subBaseline+plotElement$data[[i]],col=plotElement$color[[i]],border=NA)
		subBaseline<-subBaseline+plotElement$data[[i]]
	}

}

numericalBarplot<-function(plotElement){
  if(is.null(plotElement$resize)){
    
    data<-normalize_base0(plotElement$data)
  }else{
    data<-plotElement$data
  }
  data[which(is.na(data))]<-0
  a <- as.numeric(data)*plotElement$buckHeight
  aa<-rep(0,length(a)*4)
  aa[(1:length(a))*4-3]<-0
  aa[(1:length(a))*4-2]<-a
  aa[(1:length(a))*4-1]<-a
  aa[(1:length(a))*4-0]<-0
  aa<- aa+plotElement$baseline
  xx<-rep(0,length(a)*4)
  xx[(1:length(a))*4-3]<-c(1:length(a))*(res$commomParams$barWidth+res$commomParams$barInterval)-res$commomParams$barWidth/2
  xx[(1:length(a))*4-2]<-c(1:length(a))*(res$commomParams$barWidth+res$commomParams$barInterval)-res$commomParams$barWidth/2
  xx[(1:length(a))*4-1]<-c(1:length(a))*(res$commomParams$barWidth+res$commomParams$barInterval)+res$commomParams$barWidth/2
  xx[(1:length(a))*4-0]<-c(1:length(a))*(res$commomParams$barWidth+res$commomParams$barInterval)+res$commomParams$barWidth/2
  polygon(x=xx,y=aa,col=plotElement$pillarColor,border=NA)
  text(plotElement$leftText,x=-2,y=plotElement$baseline,cex = 0.4,adj=c(1,0),font=2)
  if(!is.null(plotElement$rightText)){
    text(plotElement$rightText,x=length(a)*(res$commomParams$barWidth+res$commomParams$barInterval)+5,y=plotElement$baseline,col = plotElement$rightTextColor,cex = 0.4,adj=c(0,0),font=2)
  }
}
categoryBarplot<-function(plotElement){
  categorys<-names(plotElement$color)
  numOfCategorys<-c()
  colOfCategorys<-c()
  for(i in categorys){
    numOfCategorys<-c(numOfCategorys,sum(as.character(plotElement$data) == i))
    colOfCategorys<-c(colOfCategorys,plotElement$color[[i]])
  }
  position <-c(0,cumsum(numOfCategorys))
  for(i in 2:length(position)){
    
    endPoint<-position[i]
    startPoint<-position[i-1]+1
    if(endPoint<startPoint)next
    data1 <- rep(0.8,times=endPoint-startPoint+1)
    
    a <- as.numeric(data1)*plotElement$buckHeight
    aa<-rep(0,length(a)*4)
    aa[(1:length(a))*4-3]<-0
    aa[(1:length(a))*4-2]<-a
    aa[(1:length(a))*4-1]<-a
    aa[(1:length(a))*4-0]<-0
    aa<- aa+plotElement$baseline
    
    xx<-rep(0,length(a)*4)
    xx[(1:length(a))*4-3]<-c(startPoint:endPoint)*(res$commomParams$barWidth+res$commomParams$barInterval)-res$commomParams$barWidth/2
    xx[(1:length(a))*4-2]<-c(startPoint:endPoint)*(res$commomParams$barWidth+res$commomParams$barInterval)-res$commomParams$barWidth/2
    xx[(1:length(a))*4-1]<-c(startPoint:endPoint)*(res$commomParams$barWidth+res$commomParams$barInterval)+res$commomParams$barWidth/2
    xx[(1:length(a))*4-0]<-c(startPoint:endPoint)*(res$commomParams$barWidth+res$commomParams$barInterval)+res$commomParams$barWidth/2
    polygon(x=xx,y=aa,col=colOfCategorys[i-1],border=NA)
  }		
  text(plotElement$leftText,x=-2,y=plotElement$baseline,cex = 0.4,adj=c(1,0),font=2)
  if(!is.null(plotElement$rightText)){
    text(plotElement$rightText,x=length(a)*(res$commomParams$barWidth+res$commomParams$barInterval)+5,y=plotElement$baseline,col = plotElement$rightTextColor,cex = 0.4,adj=c(0,0),font=2)
}
}
drawLegend<-function(legend){
  if(!is.null(legend$comment)){
    text(legend$comment$text,x=legend$left,y=legend$comment$baseline,cex = 0.4,adj=c(0,0),font=2)
  }
  
  if(!is.null(legend$title)){
    text(legend$title$text,x=legend$title$left,y=legend$title$baseline,cex=0.4,adj=c(0,0),font=2)
  }
  
  if(!is.null(legend$legend)){
    for(i in legend$legend){
      rect(i$left,i$baseline,i$left+i$size,i$baseline+i$size,col = i$color,bty="n",border =NA)
      text(i$text,x=i$left+i$size+i$size/4,y=i$baseline,cex=0.4,adj=c(0,0),font=2)
    }
    
  }
}

Boxplot<-function(plotElement){

	#plotElements$type="Boxplot"
	#plotElements$baseline=100
	 #这里五个长度的向量代表下边缘，Q1,Q2,Q3，上边缘的Y坐标
	#plotElements$Qs=c(1,2,3,4,5)  
	#plotElements$abnormal=c(9,10,11)
	#plotElements$boxWidth=10
	#plotElements$x=2
	

	baseline<-plotElement$baseline
	rect(plotElement$x-plotElement$boxWidth/2,plotElement$Qs[2]+baseline,plotElement$x+plotElement$boxWidth/2,plotElement$Qs[4]+baseline,lwd=0.5,col=plotElement$color)
	points(c(plotElement$x-plotElement$boxWidth/2,plotElement$x+plotElement$boxWidth/2),c(plotElement$Qs[3],plotElement$Qs[3])+baseline,type="l",lwd=0.5)
	points(c(plotElement$x,plotElement$x),c(plotElement$Qs[4],plotElement$Qs[5])+baseline,type="l",lwd=0.5)
	points(c(plotElement$x,plotElement$x),c(plotElement$Qs[1],plotElement$Qs[2])+baseline,type="l",lwd=0.5)
	points(c(plotElement$x,plotElement$x),c(plotElement$Qs[1],plotElement$Qs[2])+baseline,type="l",lwd=0.5)
	points(c(plotElement$x-plotElement$boxWidth/2,plotElement$x+plotElement$boxWidth/2),c(plotElement$Qs[1],plotElement$Qs[1])+baseline,type="l",lwd=0.5)
	points(c(plotElement$x-plotElement$boxWidth/2,plotElement$x+plotElement$boxWidth/2),c(plotElement$Qs[5],plotElement$Qs[5])+baseline,type="l",lwd=0.5)
	if(!is.null(plotElement$abnormal)){
		for(i in plotElement$abnormal){
			points(plotElement$x,i+baseline)
		}
	}
}

drawAxis<-function(plotElement){
	if(!is.null(plotElement$x)){
		points(plotElement$x$startEnd,c(plotElement$x$baseline,plotElement$x$baseline),type="l",lwd=0.5)
		if(plotElement$x$srt!=0){
			adj<-c(1,0)
		}else{
			adj<-c(0.5,1)
		}
		for(i in 1:length(plotElement$x$ats)){
			points(rep(plotElement$x$ats[i],2),c(plotElement$x$baseline,plotElement$x$baseline-10),type="l",lwd=0.5)
			text(plotElement$x$labs[i],x=plotElement$x$ats[i],y=plotElement$x$baseline-20,srt=plotElement$x$srt,adj=adj,cex=0.5)
		}
		if(!is.null(plotElement$x$xlab)){
			text(plotElement$x$xlab,x=mean(plotElement$x$startEnd),y=plotElement$x$baseline-75,srt=0,adj=c(0.5,0.5),cex=0.75,font=2)
		}
	}

	if(!is.null(plotElement$y)){
		points(rep(plotElement$y$x,2),plotElement$y$startEnd,type="l",lwd=0.5)
		for(i in 1:length(plotElement$y$ats)){
			points(c(plotElement$y$x,plotElement$y$x-10),rep(plotElement$y$ats[i],2),type="l",lwd=0.5)
			text(plotElement$y$labs[i],y=plotElement$y$ats[i],x=plotElement$y$x-20,srt=plotElement$y$srt,adj=c(1,0.5),cex=0.5,lwd=0.5)
		}
		if(!is.null(plotElement$y$ylab)){
			if(!is.null(plotElement$y$ylabRescale)){
				ylabRescale<-plotElement$y$ylabRescale
			}else{
				ylabRescale<-1
			}
			text(plotElement$y$ylab,y=mean(plotElement$y$startEnd),x=plotElement$y$x-95,srt=90,adj=c(0.5,0.5),cex=0.75*ylabRescale,font=2)
		}
	}
	if(!is.null(plotElement$box)){
		if(plotElement$box){
			points(c(plotElement$y$x,plotElement$x$startEnd[2],plotElement$x$startEnd[2]),c(plotElement$y$startEnd[2],plotElement$y$startEnd[2],plotElement$y$startEnd[1]),type="l",lwd=0.5)
		}
	}
}
drawRect<-function(plotElement){
	if(is.null(plotElement$border)){
		border<-NA
	}else{
		border<-plotElement$border
	}
	if(is.null(plotElement$color)){
		color<-NA
	}else{
		color<-plotElement$color
	}
	rect(plotElement$left,plotElement$bottom,plotElement$left+plotElement$width,plotElement$bottom+plotElement$height,col=color,border=border)
}
#figType="pdf"
#filename="1.pdf"
drawText<-function(plotElement){
	text(plotElement$text,x=plotElement$x,y=plotElement$y,cex=plotElement$cex,col=plotElement$color,adj=plotElement$adj,srt=plotElement$srt)

}
drawColBar<-function(plotElement){
	gradientCol<-colorRampPalette(c(plotElement$color_range[1],plotElement$color_range[2],plotElement$color_range[3]))(100)
	rect(rep(plotElement$left,length(gradientCol)),
		(1:length(gradientCol)-1)*plotElement$height/length(gradientCol)+plotElement$bottom,
		rep(plotElement$left+plotElement$width,length(gradientCol)),
		1:length(gradientCol)*plotElement$height/length(gradientCol)+plotElement$bottom,
		col=gradientCol,border=NA)
	text(plotElement$text[1],x=plotElement$left+plotElement$width/2,y=plotElement$bottom+plotElement$height+10,cex = 0.4,adj=c(0.5,0))
	text(plotElement$text[2],x=plotElement$left+plotElement$width/2,y=plotElement$bottom-10,cex = 0.4,adj=c(0.5,1))
}

#其他支持函数(小弟负责）==============



getStrLen<-function(str,cex=1){
  
  lettersSize<-c(11.1,11.1,12,12,11.1,10.2,13,12,4.6,8.4,
                 11.1,9.3,13.9,12,13,11.1,13,12,11.1,10.2,
                 12,11.1,15.7,11.1,11.1,10.2,
                 9.3,9.1,8.4,9.3,9.3,4.6,9.3,9.3,3.7,3.7,8.4,3.7,13.9,9.3,9.3,9.3,9.3,5.5,8.4,4.6,9.3,8.4,12,8.3,8.3,8.3,
                 9.7,9.7,9.3,4.7,4.7,14.9,5.6,5.6,9.3,9.3,9.3,9.3,9.3,9.3,9.3,9.3,9.3,9.3,4.7,9.7)/100
  symbols<-c("-","+","_",",",".","%","(",")",0:9," ","?")
  text<-c(LETTERS,letters,symbols)
  names(lettersSize)<-text
  fullLengh<-0
  for(i in 1:nchar(str)){
    if(substring(str,i,i) %in% names(lettersSize)){
      fullLengh<-lettersSize[substring(str,i,i)]+fullLengh
    }else{
      fullLengh<-mean(lettersSize)+fullLengh
      print(paste(substring(str,i,i)," is unknown, replaced with average length!"))
    }
  }
  return(as.numeric(fullLengh*cex))
  
}

normalize_base0<-function(x){
  x_clean<-unlist(na.omit(x))
  x<-x-min(x_clean)
  return(x/(max(x_clean)-min(x_clean)))
}

normalize_central0<-function(x){
  x_clean<-na.omit(x)
  x<-x-min(x_clean)
  return(x/(max(x_clean)-min(x_clean))-0.5)
}

autoFitAxis<-function(bottom,top){

	range<-top-bottom
	scale<-10^floor(log(range,10))
	rangeAdj<-range/scale
	if(rangeAdj<1.5){
		seg<-0.25
	}else if(rangeAdj<4){
		seg<-0.5
	}else if(rangeAdj<7){
		seg<-1
	}else{
		seg<-2
	}
	return(seq(floor(bottom/(seg*scale))*(seg*scale),ceiling(top/(seg*scale))*(seg*scale),seg*scale))
}
getQuantile<-function(num,v){
	through<-TRUE
	v<-v[order(v)]
	for(i in 1:length(v)){
		if(num<v[i]){
			through<-FALSE
			break
		}
	}
	if(through){
		return(1)
	}else{
		return((i-1)/length(v))
	}
	

}