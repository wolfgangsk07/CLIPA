Drug_Omics_response <- function (Dataset="GDSC",Omics = "mRNA",gene= "TP53",
                                 genesets = c("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF",
                                              "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", "TUBB6" , "VEGFA"),
                                 drug_input	= c("AZD5363","PD0325901"),
                                 cell_lines=c("A253", "BB49-HNC","Detroit562") ,
                                 Omics_pillar_color = "#34495e",Drug_pillar_color="#928a97",Mutation_color="#000091",No_mutation_color ="grey",
                                 DEL_color = "#000091",LOSS_color = "#4343b5",NOR_color = "#8e8e8e",GAIN_color = "#ce5555",AMP_color = "#ba0000",
                                 resistant_color="#ba0000",sensitive_color="#000091",cut_off=0.2,cor_method="spearman",mutation_method="wilcoxon",
                                 statistic_method="FDR",preset = TRUE){
  options(stringsAsFactors = F)
  
  if(Omics == "GSVA"){
	if(preset == TRUE){
    preset_GSVA_path <- paste0("./data/",Dataset,"/",Dataset,"_gsva_es.RData")
    load(preset_GSVA_path)
    gsva_es <- get(ls()[ls()==paste0(Dataset,"_gsva_es")])
	gene <- genesets
    Final_Omics <- gsva_es[cell_lines,c(gene,"cell_lines")]
  }else{
    library(GSVA)
    File_path <- paste0("./data/",Dataset,"/",Dataset,"_mRNA.RData")
    load(File_path)
    genes <-list()
    gene <- "GSVA_score"
    genes$GSVA_score <- genesets
    Final_Omics <- get(ls()[grep(paste0(Dataset,"_Omics"),ls())])
    gsva_es <- gsva(as.matrix(Final_Omics), genes)
    gsva_es <- as.data.frame(t(gsva_es))
    gsva_es$cell_line <- row.names(gsva_es)
    Final_Omics <- gsva_es[cell_lines,]
  }}else{
    name_end <- ifelse(Omics =="methylation_site",6,2)
    Omics_file_name <-paste0("./data/",Dataset,"/",Omics,"/",substr(gene,1,name_end),".RData")
    if (file.exists(Omics_file_name)) {
      load(Omics_file_name) 
    }else{
      res <- list()
      res$error <- paste0("This Omics file could not be found")
      return(res)
    }
    if( gene %in% names(list)){
      Final_Omics <- as.data.frame(list[c("cell_lines",gene)])
      row.names(Final_Omics) <- Final_Omics$cell_lines
      Final_Omics <- Final_Omics[cell_lines,]
      #Final_Omics[,gene] <- as.numeric(Final_Omics[,gene])
    }else{
      res <- list()
      res$error <- paste0("This selected gene is not found in ",Dataset, " ", Omics," dataset")
      return(res)
    }
  }
  if(length(grep("-",gene))==1){
    names(Final_Omics)[-grep("cell_lines",names(Final_Omics))] <- 
      gsub("\\.","-",names(Final_Omics)[-grep("cell_lines",names(Final_Omics))])
  }
  Final_Omics_tmp <- Final_Omics[order(Final_Omics[,gene],na.last = F),]
  drug_response_path <- paste0("./data/",Dataset,"/",Dataset,"_drug_response.RData")
  load(drug_response_path)
  
   Final_drug_response <- get(ls()[which(ls()==paste0(Dataset,"_drug_response"))])
  Final_drug_response$cell_line_drug <- row.names(Final_drug_response)
  Final_drug_tmp <- Final_drug_response[row.names(Final_Omics_tmp),c(drug_input,"cell_line_drug")]
  final <- cbind(Final_drug_tmp,Final_Omics_tmp)
  final <- final[,-grep("cell_line",names(final))]
  names(final) <- gsub("\\.","-",names(final))
  if(Omics=="mutation"){
    Mutation_cell <- row.names(final)[final[,grep(gene,names(final))]=="Mutation"]
    None_Mutation_cell <- row.names(final)[final[,grep(gene,names(final))]=="None"]
    if( length(Mutation_cell) == 0 ){
      res <- list()
      res$error <- "There is no mutation of this gene in your selected cell lines"
      return(res)
    }else{
      pro <- matrix(NA,length(drug_input),2)
      row.names(pro) <- drug_input
      colnames(pro) <- c("FDR","P")
      test_method <- ifelse(mutation_method =="wilcoxon","wilcox.test","t.test")
      pro <- as.data.frame(pro)
      for (i in drug_input) {
        tmp <- get(test_method)(final[Mutation_cell,i],
                                final[None_Mutation_cell,i])
        pro[i,2]<- tmp$p.value
      }
      
      pro_new <- pro[order(pro$P),]
      pro_new$FDR <- p.adjust(pro_new$P)
      
      
      if (statistic_method =="P") {
        pro_new$biaoji <- ifelse(pro_new$P < 0.0001,"****",ifelse(pro_new$P < 0.001,"***",
                                                                  ifelse(pro_new$P < 0.01,"**",
                                                                         ifelse(pro_new$P < 0.05,"*",""))))
        pro_new$col <- ifelse(pro_new$P < 0.05,"black","#928a97")
        pro_new$P_new  <- sprintf("%0.3f",  pro_new$P)
        pro_new$wenben <- paste0("P = ",pro_new$P_new," ",pro_new$biaoji)
      }else{
        pro_new$biaoji <- ifelse(pro_new$FDR < 0.0001,"****",ifelse(pro_new$FDR < 0.001,"***",
                                                                  ifelse(pro_new$FDR < 0.01,"**",
                                                                         ifelse(pro_new$FDR < 0.05,"*",""))))
        pro_new$col <- ifelse(pro_new$FDR < 0.05,"black","#928a97")
        pro_new$FDR_new  <- sprintf("%0.3f",  pro_new$FDR)
        pro_new$wenben <- paste0("FDR = ",pro_new$FDR_new," ",pro_new$biaoji)
      }
      pro_new$wenben <- gsub("= 0.000","< 0.001",pro_new$wenben)
      }}else{
      pro <- matrix(NA,length(drug_input),2)
      row.names(pro) <- drug_input
      colnames(pro) <- c("cor","P")
      pro <- as.data.frame(pro)
      for (i in drug_input) {
        tmp <- cor.test(final[,i],final[,gene],method = cor_method)
        pro[i,1] <-tmp$estimate
        pro[i,2]<- tmp$p.value
      }
      pro$cor1 <- sprintf("%0.3f", pro$cor)
      pro_new <- pro[order(pro$P),]
      pro_new$FDR <- p.adjust(pro_new$P)
      
      
      if (statistic_method =="P") {
        pro_new$biaoji <- ifelse(pro_new$P < 0.0001,"****",ifelse(pro_new$P < 0.001,"***",
                                                                  ifelse(pro_new$P < 0.01,"**",ifelse(pro_new$P < 0.05,"*",""))))
      }else{
        pro_new$biaoji <- ifelse(pro_new$FDR < 0.0001,"****",ifelse(pro_new$FDR < 0.001,"***",
                                                                  ifelse(pro_new$FDR < 0.01,"**",ifelse(pro_new$FDR < 0.05,"*",""))))
      }
      pro_new$col <- ifelse(abs(pro_new$cor) < cut_off,"#928a97","black")
      pro_new$qianzui <- ifelse(cor_method=="spearman","Rho = ","r = ")
      pro_new$wenben <- paste0(pro_new$qianzui,pro_new$cor1," ",pro_new$biaoji)
      if(Omics %in% c("mRNA","methylation_promoter","methylation_body","methylation_site","GSVA","protein", "microRNA")) {
        final[!is.na(final[,gene]),gene] <- seq(0.3,0.7,length.out =length(final[!is.na(final[,gene]),gene] ))
        final[is.na(final[,gene]),gene] <- 0
      }else if(Omics=="CNV"){
        final[,gene] <- final[,gene]
      }
    }
  
  final <- final[,c(drug_input,gene)]
  if(Dataset=="GDSC"){
    load("./data/GDSC/GDSC_drug_response_Binary.RData")
	GDSC_drug_response_Binary$cell_lines <- row.names(GDSC_drug_response_Binary)
    Drug_response_Binary <- GDSC_drug_response_Binary[row.names(final),c("cell_lines",drug_input)]
  }
  pro_new <- pro_new[drug_input,]
  pro_new$P <- round(pro_new$P,3)
  pro_new$FDR <- round(pro_new$FDR,3)
  
  if(Omics=="mutation"){
    final[,-dim(final)[2]] <- round(final[,-dim(final)[2]],3)
  }else{
    final <- round(final,3)
    pro_new$cor <- round(pro_new$cor,3)
  }
  res<-list()
  res$Dataset <- Dataset
  res$Omics <- Omics
  
  
  #开始组织画图参数======================================================================================
  buckHeight <- 25
  buckInterv <- buckHeight/10
  
  if(length(final[,gene])>600){
	barWidth = 2.25*600/length(final[,gene])
	barInterval = 0.75*600/length(final[,gene])
  }else{
	barWidth = 2.25
	barInterval = 0.75
  }

  #leftMargin=280
  
  rightMargin=180
  bottomMargin=10
  topMargin=10
  baseline=0+bottomMargin
  dpi <- 300
  
  if(Omics %in% c("methylation_site","GSVA","microRNA","Metabolomics","Chromatin")){
    OmicsTag<-""
  }else if (Omics=="mRNA"){
    OmicsTag<-" mRNA"
  }else if (Omics=="mutation"){
    OmicsTag<-" mutation"
  }else if (Omics=="methylation_body"){
    OmicsTag<-" body"
  }else if(Omics=="methylation_promoter"){
    OmicsTag<-" promoter"
  }else if(Omics=="CNV"){
    OmicsTag<-" CNV"
  }else if(Omics=="protein"){
    OmicsTag<-" protein"
  }else if(Omics=="CRISPR"){
    OmicsTag<-" CRISPR"
  }else if(Omics=="RNAi_combine"){
    OmicsTag<-" RNAi"
  }else if(Omics=="CNV_continuous"){
    OmicsTag<-" CNV"
  }else if(Omics=="Metastatic"){
    OmicsTag<-" metastatisis"
  }
  
  
  leftMargin<-max(sapply(c(drug_input,paste0(gene,OmicsTag)),function(x)getStrLen(x,0.5)))*dpi+5
  
  
  #当画幅宽度小于700，则不再继续减小画幅宽度
  minWidth<-round(getStrLen("* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001",0.5)*dpi)
  
  minWidth<-700
	currWidth<-length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin
  if(currWidth<minWidth){
    leftMargin<-leftMargin+(minWidth-currWidth)/2
    rightMargin<-rightMargin+(minWidth-currWidth)/2
  }
  
  
  #定义公共画图参数======================================================================================
  res$commomParams$barWidth<-barWidth
  res$commomParams$barInterval<-barInterval
  res$commomParams$leftMargin<-leftMargin
  res$commomParams$rightMargin<-rightMargin
  res$commomParams$width<-(barWidth+barInterval)*length(final[,gene])
  
  #定义legend参数======================================================================================
  if(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin==minWidth){
    legendLeft=-leftMargin+30
  }else if(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin<870){
    legendLeft=(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin-minWidth)/2-leftMargin+30
  }else{
    legendLeft=0
  }

  if(Omics=="CNV"){
    
    res$plotElements$legend$plotType<-"drawLegend"
    
    res$plotElements$legend$left<-legendLeft
    res$plotElements$legend$comment$baseline<-baseline
    res$plotElements$legend$comment$text<-"* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001"
    baseline<-baseline+buckHeight+buckInterv
    
    res$plotElements$legend$title$text<-"CNV type:"
    res$plotElements$legend$title$baseline<-baseline
    res$plotElements$legend$title$left<-legendLeft
    lenOfTitle<-110
    legendOffset<-0
    res$plotElements$legend$legend[["-2"]]$baseline<-baseline
    res$plotElements$legend$legend[["-2"]]$size<-15
    res$plotElements$legend$legend[["-2"]]$color<-DEL_color
    res$plotElements$legend$legend[["-2"]]$left<-legendLeft+lenOfTitle
    res$plotElements$legend$legend[["-2"]]$text<-"DEL"
    
	legendOffset<-legendOffset+getStrLen("DEL",0.5)*300+15
    res$plotElements$legend$legend[["-1"]]$baseline<-baseline
    res$plotElements$legend$legend[["-1"]]$size<-15
    res$plotElements$legend$legend[["-1"]]$color<-LOSS_color
    res$plotElements$legend$legend[["-1"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend$legend[["-1"]]$text<-"LOSS"
    
	legendOffset<-legendOffset+getStrLen("LOSS",0.5)*300+15
    res$plotElements$legend$legend[["0"]]$baseline<-baseline
    res$plotElements$legend$legend[["0"]]$size<-15
    res$plotElements$legend$legend[["0"]]$color<-NOR_color
    res$plotElements$legend$legend[["0"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend$legend[["0"]]$text<-"NOR"
    
    
    legendOffset<-legendOffset+getStrLen("NOR",0.5)*300+15
    res$plotElements$legend$legend[["1"]]$baseline<-baseline
    res$plotElements$legend$legend[["1"]]$size<-15
    res$plotElements$legend$legend[["1"]]$color<-GAIN_color
    res$plotElements$legend$legend[["1"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend$legend[["1"]]$text<-"GAIN"
    
	legendOffset<-legendOffset+getStrLen("GAIN",0.5)*300+15
    res$plotElements$legend$legend[["2"]]$baseline<-baseline
    res$plotElements$legend$legend[["2"]]$size<-15
    res$plotElements$legend$legend[["2"]]$color<-AMP_color
    res$plotElements$legend$legend[["2"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend$legend[["2"]]$text<-"AMP"
    
    baseline<-baseline+buckHeight*1.5+buckInterv
  }else if(Omics=="mutation") {
    res$plotElements$legend$plotType<-"drawLegend"
    res$plotElements$legend$left<-legendLeft
    res$plotElements$legend$comment$baseline<-baseline
    res$plotElements$legend$comment$text<-"* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001"
    baseline<-baseline+buckHeight+buckInterv
    
    res$plotElements$legend$title$text<-"Type:"
    res$plotElements$legend$title$baseline<-baseline
    res$plotElements$legend$title$left<-legendLeft
    lenOfTitle<-60
    
    res$plotElements$legend$legend[["Mutation"]]$baseline<-baseline
    res$plotElements$legend$legend[["Mutation"]]$size<-15
    res$plotElements$legend$legend[["Mutation"]]$color<- Mutation_color
    res$plotElements$legend$legend[["Mutation"]]$left<-legendLeft+lenOfTitle
    res$plotElements$legend$legend[["Mutation"]]$text<-"Mutation"
    
    res$plotElements$legend$legend[["None"]]$baseline<-baseline
    res$plotElements$legend$legend[["None"]]$size<-15
    res$plotElements$legend$legend[["None"]]$color<- No_mutation_color
    res$plotElements$legend$legend[["None"]]$left<-legendLeft+lenOfTitle+110
    res$plotElements$legend$legend[["None"]]$text<-"No alteration"
    baseline<-baseline+buckHeight*1.5+buckInterv
  }else{
    res$plotElements$legend$plotType<-"drawLegend"
    res$plotElements$legend$left<-legendLeft
    res$plotElements$legend$comment$text<-"* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001"
    res$plotElements$legend$comment$baseline<-baseline
    baseline<-baseline+buckHeight*1.5+buckInterv
  }
  
  
  #定义正文参数（组学）======================================================================================
  if(Omics=="CNV"){
    lab<-c("-2","-1","0","1","2")
    color_tmp<-c(DEL_color,LOSS_color,NOR_color,GAIN_color,AMP_color)
    for(i in 1:length(lab)){
      res$plotElements[[paste0("gene_",i)]]$plotType<-"numericalBarplot"
      res$plotElements[[paste0("gene_",i)]]$pillarColor<-color_tmp[i]
      res$plotElements[[paste0("gene_",i)]]$baseline<-baseline
      CNVdata_tmp<-rep(0,length(final[,gene]))
      CNVdata_tmp[which(final[,gene]==lab[i])]<-0.8
      res$plotElements[[paste0("gene_",i)]]$data <- CNVdata_tmp
      res$plotElements[[paste0("gene_",i)]]$buckHeight<-buckHeight
      res$plotElements[[paste0("gene_",i)]]$pillarClassName<-"CellAndDrugResistTags"
      res$plotElements[[paste0("gene_",i)]]$resize<-"false"
      if(i==1){
        res$plotElements[[paste0("gene_",i)]]$leftText <- paste0(gene,OmicsTag)
      }
    }
  }else if(Omics=="mutation"){
    lab<-c("Mutation","None")
    color_tmp<-c(Mutation_color,No_mutation_color)
    for(i in 1:length(lab)){
      res$plotElements[[paste0("gene_",i)]]$plotType<-"numericalBarplot"
      res$plotElements[[paste0("gene_",i)]]$pillarColor<-color_tmp[i]
      res$plotElements[[paste0("gene_",i)]]$baseline<-baseline
      CNVdata_tmp<-rep(0,length(final[,gene]))
      CNVdata_tmp[which(final[,gene]==lab[i])]<-0.8
      res$plotElements[[paste0("gene_",i)]]$data <- CNVdata_tmp
      res$plotElements[[paste0("gene_",i)]]$buckHeight<-buckHeight
      res$plotElements[[paste0("gene_",i)]]$pillarClassName<-"CellAndDrugResistTags"
      res$plotElements[[paste0("gene_",i)]]$resize<-"false"
      if(i==1){
        res$plotElements[[paste0("gene_",i)]]$leftText <- paste0(gene,OmicsTag)
      }
    }
  }else{
    res$plotElements$gene$plotType<-"numericalBarplot"
    res$plotElements$gene$pillarColor<-Omics_pillar_color
    res$plotElements$gene$baseline<-baseline
    res$plotElements$gene$leftText <- gsub("_"," ",paste0(gene,OmicsTag))
    res$plotElements$gene$data <- final[,gene]
    res$plotElements$gene$buckHeight<-buckHeight
    res$plotElements$gene$pillarClassName<-"CellAndDrugResistTags"
    res$plotElements$gene$resize<-"false"
  }
  
  if (Omics=="methylation_site") {
    methylation_site_annotation <- paste0("./data/GDSC/methylation_site_annotation/",substr(gene,1,name_end),"_ann.RData")
    load(methylation_site_annotation)
    final_annotation <- list[[gene]]
    names(final_annotation)<-list[["item"]]
    final_annotation[2] <- gsub(" ","",final_annotation[2] )
    res$plotElements$gene$leftTextAdditionInfos <- toJSON(final_annotation)
    res$plotElements$gene$leftTextClassName <- "MethSiteTags"
  }
  baseline<-baseline+buckInterv+buckHeight
  library(rjson)
  
  #row.names(final)==row.names(Drug_response_Binary)
  #定义正文参数（药物）======================================================================================
  for(i in drug_input ){
    res$plotElements[[i]]$plotType<-"numericalBarplot"
    res$plotElements[[i]]$data <- final[,i]
    if(Dataset=="GDSC"){
      res$plotElements[[i]]$pillarAdditionInfos <- sapply(Drug_response_Binary[,i],function(x)toJSON(c(i,x)),USE.NAMES=FALSE)
    }
    res$plotElements[[i]]$leftText <- i
    res$plotElements[[i]]$rightText <- pro_new$wenben[which(row.names(pro_new)==i)]
    res$plotElements[[i]]$rightTextColor <- pro_new$col[which(row.names(pro_new)==i)]
    res$plotElements[[i]]$pillarColor<-Drug_pillar_color
    res$plotElements[[i]]$baseline<-baseline
    res$plotElements[[i]]$buckHeight<-buckHeight
    res$plotElements[[i]]$leftTextClassName<-"drugTags"
    res$plotElements[[i]]$pillarClassName<-"CellAndDrugResistTags"
    baseline<-baseline+buckInterv+buckHeight
  }
  
  #收尾获得文档高度，并返回======================================================================================
  res$commomParams$DocHeight<-baseline+topMargin
  res$cell_line <- row.names(final)
  #save(res,file="r:/debug.rda")
  return(res)
}

######

FigOutput <- function(res,figType="pdf",filename="1.pdf"){
  #res为上一步产生的计算结果
  #figType为需要输出的图像类型，如"pdf","svg","postscript","png","jpeg","bitmap"等
  #filename为保存的文件名
  #Num_drug <- length(res$gene$data)
  
  dpi <- 300
  
  
  
  if(figType=="tiff"){
    tiff(file=filename,width =(res$commomParams$leftMargin+res$commomParams$width+res$commomParams$rightMargin)/dpi,
         height = res$commomParams$DocHeight/dpi,compression="lzw",units="in",res=1200)
  }else{
    get(figType)(file=filename,width =(res$commomParams$leftMargin+res$commomParams$width+res$commomParams$rightMargin)/dpi,
                 height = res$commomParams$DocHeight/dpi)
  }
  
  par(mar=c(0,0,0,0),omi=c(0,0,0,0))
  plot(-99999,-99999,xlim=c(-res$commomParams$leftMargin,res$commomParams$width+res$commomParams$rightMargin),
       ylim=c(0,res$commomParams$DocHeight),xaxs="i",yaxs="i",bty="n",xaxt="n",yaxt="n",xlab="",ylab=""
  )
  
  
  
  for(i in res$plotElements){
    get(i$plotType)(i)
    
  }
  
  dev.off()
}

#计算函数（大哥负责）=============
MultiOmics_calculation <- function (Dataset="GDSC", gene= "TP53",cell_lines=c("A253", "BB49-HNC","Detroit562") ,
                                    mRNA_pillar_color = "#34495e",methylation_site_pillar_color="#928a97",
                                    TMB_pillar_color="#34495e",ploidy_pillar_color="#928a97",
                                    Mutation_color="#000091",No_mutation_color ="grey",
                                    MSI_color="#000091",MSS_color ="grey",DEL_color = "#000091",
                                    LOSS_color = "#4343b5",NOR_color = "#8e8e8e",GAIN_color = "#ce5555",
                                    AMP_color = "#ba0000",protein_pillar_color = "#928a97",
                                    cut_off=0.2,cor_method="spearman",mutation_method="wilcoxon",statistic_method="FDR"){
  options(stringsAsFactors = F)
	#-----mRNA
    mRNA_file_name <- paste0("./data/",Dataset,"/mRNA/",substr(gene,1,2),".RData")
    load(mRNA_file_name)
    mRNA_Omics <- as.data.frame(list[c("cell_lines",gene)])
	names(mRNA_Omics) <- gsub("\\.","-",names(mRNA_Omics))
    row.names(mRNA_Omics) <- mRNA_Omics$cell_lines
    names(mRNA_Omics)[names(mRNA_Omics)==gene] <- paste0(gene,"_mRNA")
    #-----CNV
    CNV_file_name <- paste0("./data/",Dataset,"/CNV/",substr(gene,1,2),".RData")
    load(CNV_file_name)
    CNV_Omics <- as.data.frame(list[c("cell_lines",gene)])
	names(CNV_Omics) <- gsub("\\.","-",names(CNV_Omics))
    row.names(CNV_Omics) <- CNV_Omics$cell_lines
    names(CNV_Omics)[names(CNV_Omics)==gene] <- paste0(gene,"_CNV")
    #-------mutation
    mutation_file_name <- paste0("./data/",Dataset,"/mutation/",substr(gene,1,2),".RData")
    load(mutation_file_name)
    mutation_Omics <- as.data.frame(list[c("cell_lines",gene)])
	names(mutation_Omics) <- gsub("\\.","-",names(mutation_Omics))
    row.names(mutation_Omics) <- mutation_Omics$cell_lines
    names(mutation_Omics)[names(mutation_Omics)==gene] <- paste0(gene,"_mutation")
    #--------methylation site annotation
    methylation_site_annotation <- paste0("./data/GDSC/methylation_gene_annotation/",substr(gene,1,2),"_methylation.RData")
    load(methylation_site_annotation)
    CpG <- list[[gene]]$methylation_site
    #--------methylation_site
    methylation_site_file_name <- paste0("./data/",Dataset,"/methylation_site/",substr(CpG[1],1,6),".RData")
    load(methylation_site_file_name)
    methylation_site_Omics <- as.data.frame(list[c("cell_lines",CpG[1])])
    row.names(methylation_site_Omics) <- methylation_site_Omics$cell_lines
    names(methylation_site_Omics)[names(methylation_site_Omics)=="cell_lines"] <- 
      paste0("cell_lines_",CpG[1])
    if (length(CpG) >1 ) {
      for (i in 2:length(CpG)) {
        methylation_site_file_name <- paste0("./data/",Dataset,"/methylation_site/",substr(CpG[i],1,6),".RData")
        load(methylation_site_file_name)
        methylation_site_Omics_tmp <- as.data.frame(list[c("cell_lines",CpG[i])])
        row.names(methylation_site_Omics_tmp) <- methylation_site_Omics_tmp$cell_lines
        names(methylation_site_Omics_tmp)[names(methylation_site_Omics_tmp)=="cell_lines"] <- 
          paste0("cell_lines_",CpG[i])
        methylation_site_Omics <- cbind(methylation_site_Omics,methylation_site_Omics_tmp)
      }}else{
        methylation_site_Omics <- methylation_site_Omics
      }
  methylation_site_Omics <- methylation_site_Omics[,-grep("cell_line",names(methylation_site_Omics))]
  methylation_site_Omics$cell_lines <- row.names(methylation_site_Omics)
  #----添加蛋白组学
  if(Dataset== "NCI60"){
    gene_information_path <- paste0("./data/",Dataset,"/",Dataset,"_gene_information.RData")
    load(gene_information_path)
    gene_information <- get(ls()[which(ls()==paste0(Dataset,"_gene_information"))])
    if(gene %in% gene_information$protein$Gene){
      protein_file_name <- paste0("./data/",Dataset,"/protein/",substr(gene,1,2),".RData")
      load(protein_file_name)
      protein_Omics <- as.data.frame(list[c("cell_lines",gene)])
	  names(protein_Omics) <- gsub("\\.","-",names(protein_Omics))
      row.names(protein_Omics) <- protein_Omics$cell_lines
      names(protein_Omics)[names(protein_Omics)==gene] <- paste0(gene,"_protein")
    }else{
      res$error_proteomics <- paste0("This selected gene is not found in NCI60 proteomics dataset")
    }}
  #--------mutational_burden
  if(Dataset=="GDSC"){
    load("data/GDSC/GDSC_cell_annotation.RData")
    TMB <- GDSC_cell_annotation[,c("Cell_Name","mutational_burden")]
    row.names(TMB) <- TMB$Cell_Name
    names(TMB)[names(TMB)=="Cell_Name"] <- "cell_lines"
    names(TMB)[names(TMB) == "mutational_burden"] <- "TMB"
  }else if(Dataset== "NCI60"){
    load("data/NCI60/NCI60_cell_annotation.RData")
    TMB <- NCI60_cell_annotation[,c("Cell_Name","mutational_burden")]
    row.names(TMB) <- TMB$Cell_Name
    names(TMB)[names(TMB)=="Cell_Name"] <- "cell_lines"
    names(TMB)[names(TMB) == "mutational_burden"] <- "TMB"
  }
  #------ploidy
  if(Dataset=="GDSC"){
    load("data/GDSC/GDSC_cell_annotation.RData")
    ploidy <- GDSC_cell_annotation[,c("Cell_Name","ploidy")]
    row.names(ploidy) <- ploidy$Cell_Name
    names(ploidy)[names(ploidy)=="Cell_Name"] <- "cell_lines"
  }else if(Dataset== "NCI60"){
    load("data/NCI60/NCI60_cell_annotation.RData")
    ploidy <- NCI60_cell_annotation[,c("Cell_Name","ploidy")]
    row.names(ploidy) <- ploidy$Cell_Name
    names(ploidy)[names(ploidy)=="Cell_Name"] <- "cell_lines"
  }
  #------MSI_Status
  if(Dataset=="GDSC"){
    load("data/GDSC/GDSC_cell_annotation.RData")
    MSI <- GDSC_cell_annotation[,c("Cell_Name","MSI_Status")]
    row.names(MSI) <- MSI$Cell_Name
    names(MSI)[names(MSI)=="Cell_Name"] <- "cell_lines"
  }else if(Dataset== "NCI60"){
    load("data/NCI60/NCI60_cell_annotation.RData")
    MSI <- NCI60_cell_annotation[,c("Cell_Name","MSI_Status")]
    row.names(MSI) <- MSI$Cell_Name
    names(MSI)[names(MSI)=="Cell_Name"] <- "cell_lines"
  }
  #----合并数据
  if(length(grep("protein_Omics",ls()))==1){
    protein_Omics[60,] <- c("MDA-N",NA)
    row.names(protein_Omics)[60] <- "MDA-N"
    protein_Omics <- protein_Omics[cell_lines,]
    protein_Omics[,paste0(gene,"_protein")] <- as.numeric(protein_Omics[,paste0(gene,"_protein")])
  }
  mutation_Omics <- mutation_Omics[cell_lines,]
  library(dplyr)
  if(length(grep("protein_Omics",ls()))==1){
    multiomics <- left_join(left_join(left_join(left_join(left_join(left_join(left_join(mutation_Omics,CNV_Omics),mRNA_Omics),
                                                methylation_site_Omics),TMB),ploidy),MSI),protein_Omics)}else{
                                                  multiomics <- left_join(left_join(left_join(left_join(left_join(left_join(mutation_Omics,CNV_Omics),mRNA_Omics),
                                                                                        methylation_site_Omics),TMB),ploidy),MSI)  
                                                }
  
  row.names(multiomics) <- multiomics$cell_lines
  multiomics <- multiomics[,-grep("cell_line",names(multiomics))]
  #------------计算相关系数
  pro <- matrix(NA,dim(multiomics)[2]-1,2)
  row.names(pro) <- names(multiomics)[-grep("mRNA",names(multiomics))]
  colnames(pro) <- c("cor","P")
  pro <- as.data.frame(pro)
  corrlation_variable <- row.names(pro)[-c(grep("mutation",row.names(pro)),grep("MSI_Status",row.names(pro)))]
  for (i in corrlation_variable) {
    tmp <- cor.test(as.numeric(multiomics[,i]),
                    as.numeric(multiomics[,grep("mRNA",names(multiomics))]),
                    method = cor_method)
    pro[i,1] <-tmp$estimate
    pro[i,2]<- tmp$p.value
  }
  Mutation_cell <- row.names(multiomics)[multiomics[,grep("mutation",names(multiomics))]=="Mutation"]
  None_Mutation_cell <- row.names(multiomics)[multiomics[,grep("mutation",names(multiomics))]=="None"]
  Mutation_cell <- Mutation_cell[!is.na(Mutation_cell)]
  None_Mutation_cell <- None_Mutation_cell[!is.na(None_Mutation_cell)]
  if( length(Mutation_cell) == 0 ){
    res$error_mutation <- "There is no mutation of this gene in your selected cell lines"
    pro[grep("mutation",row.names(pro)),2] <- NA
  }else{
    test_method <- ifelse(mutation_method =="wilcoxon","wilcox.test","t.test")
    tmp <- get(test_method)(multiomics[Mutation_cell,grep("mRNA",names(multiomics))],
                            multiomics[None_Mutation_cell,grep("mRNA",names(multiomics))])
    pro[grep("mutation",row.names(pro)),2] <- tmp$p.value
  }
  table(multiomics$MSI_Status)
  #MSI MSS 
  MSS_cell <- row.names(multiomics)[multiomics[,grep("MSI",names(multiomics))]=="MSS"]
  MSI_cell <- row.names(multiomics)[multiomics[,grep("MSI",names(multiomics))]=="MSI"]
  MSI_cell <- MSI_cell[!is.na(MSI_cell)]
  MSS_cell <- MSS_cell[!is.na(MSS_cell)]
  if(length(MSI_cell) == 0 ){
  res <- list()
    res$error <- "Your selected cell lines are all in MSS status"
    pro[grep("MSI",row.names(pro)),2] <- NA
  }else{
    test_method <- ifelse(mutation_method =="wilcoxon","wilcox.test","t.test")
    tmp <- get(test_method)(multiomics[MSS_cell,grep("mRNA",names(multiomics))],
                            multiomics[MSI_cell,grep("mRNA",names(multiomics))])
    pro[grep("MSI",row.names(pro)),2] <- tmp$p.value
  }
  pro$cor1 <- sprintf("%0.3f", pro$cor)
  pro_new <- pro[order(pro$P),]
  pro_new$FDR <- p.adjust(pro_new$P)
  pro_new$P_new  <- sprintf("%0.3f",  pro_new$P)
  pro_new$FDR_new  <- sprintf("%0.3f",  pro_new$FDR)
  pro_new$col <- NA
  pro_new$col <- ifelse(abs(pro_new$cor) < cut_off,"#928a97","black")
  if(statistic_method=="P"){
  pro_new$biaoji <- ifelse(pro_new$P < 0.0001,"****",ifelse(pro_new$P < 0.001,"***",
                                                            ifelse(pro_new$P < 0.01,"**",ifelse(pro_new$P < 0.05,"*",""))))
		pro_new[grep("mutation",row.names(pro_new)),"col"] <-  ifelse(pro_new[grep("mutation",row.names(pro_new)),"P"] < 0.05,
                                                                "black","#928a97")
  pro_new[grep("MSI",row.names(pro_new)),"col"] <-  ifelse(pro_new[grep("MSI",row.names(pro_new)),"P"] < 0.05,
                                                           "black","#928a97")													
															
															}else if(statistic_method=="FDR"){
		pro_new$biaoji <- ifelse(pro_new$FDR < 0.0001,"****",ifelse(pro_new$FDR < 0.001,"***",
                                                            ifelse(pro_new$FDR < 0.01,"**",ifelse(pro_new$FDR < 0.05,"*",""))))
			pro_new[grep("mutation",row.names(pro_new)),"col"] <-  ifelse(pro_new[grep("mutation",row.names(pro_new)),"FDR"] < 0.05,
                                                                "black","#928a97")
  pro_new[grep("MSI",row.names(pro_new)),"col"] <-  ifelse(pro_new[grep("MSI",row.names(pro_new)),"FDR"] < 0.05,
                                                           "black","#928a97")}												
	mutaion_MSI_location <- c(grep("mutation",row.names(pro_new)),grep("MSI",row.names(pro_new)))
  
  
  pro_new$qianzui <- NA
  pro_new$wenben <- NA
  
  
  
    pro_new[-mutaion_MSI_location,"qianzui"] <- ifelse(cor_method=="spearman","Rho = ","r = ")
    if(statistic_method=="P"){
	pro_new[mutaion_MSI_location,"qianzui"] <- "P = "
	pro_new[mutaion_MSI_location,"wenben"] <- 
      paste0("P = ",pro_new[mutaion_MSI_location,"P_new"],
             pro_new[mutaion_MSI_location,"biaoji"])
	}else if(statistic_method=="FDR"){
	pro_new[mutaion_MSI_location,"qianzui"] <- "FDR = "
	pro_new[mutaion_MSI_location,"wenben"] <- 
      paste0("FDR = ",pro_new[mutaion_MSI_location,"FDR_new"],
             pro_new[mutaion_MSI_location,"biaoji"])}
    
	
	pro_new[-mutaion_MSI_location,"wenben"] <- 
      paste0(pro_new[-mutaion_MSI_location,"qianzui"],pro_new[-mutaion_MSI_location,"cor1"],
             pro_new[-mutaion_MSI_location,"biaoji"])
    
    pro_new$wenben <- gsub("= 0.000","< 0.001",pro_new$wenben)
  if(is.na(pro_new[grep("MSI",row.names(pro_new)),"P"])){
    pro_new[grep("MSI",row.names(pro_new)),"wenben"] <- ""
  }
  if(is.na(pro_new[grep("mutation",row.names(pro_new)),"P"])){
    pro_new[grep("mutation",row.names(pro_new)),"wenben"] <- ""
  }
  #----------mRNA 赋值
  multiomics <- multiomics[order(multiomics[,paste0(gene,"_mRNA")],na.last = F),]
  multiomics[!is.na(multiomics[,grep("mRNA",names(multiomics))]),grep("mRNA",names(multiomics))] <- 
    seq(0.3,0.7,length.out =length(multiomics[!is.na(multiomics[,grep("mRNA",names(multiomics))]),
                                              grep("mRNA",names(multiomics))] ))
  multiomics[is.na(multiomics[,grep("mRNA",names(multiomics))]),grep("mRNA",names(multiomics))] <- 0
  multiomics[,-c(grep("mutation",names(multiomics)),grep("CNV",names(multiomics)),
                 grep("MSI",names(multiomics)))] <- 
    round(multiomics[,-c(grep("mutation",names(multiomics)),grep("CNV",names(multiomics)),
                         grep("MSI",names(multiomics)))],3)
  pro_new <- pro_new[names(multiomics)[-grep("mRNA",names(multiomics))],]
  #---存储计算结果
  res<-list()
  buckHeight <- 25
  buckInterv <- buckHeight/10
   if(length(cell_lines)>600){
	barWidth = 2.25*600/length(cell_lines)
	barInterval = 0.75*600/length(cell_lines)
  }else{
	barWidth = 2.25
	barInterval = 0.75
  }
  rightMargin=180
  bottomMargin=10
  topMargin=25
  baseline=0+bottomMargin
  dpi <- 300
  leftMargin<-max(sapply(names(multiomics),function(x)getStrLen(x,0.5)))*dpi+5
  #当画幅宽度小于570，则不再继续减小画幅宽度
   minWidth<-round(getStrLen("* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001",0.5)*dpi)
  
  minWidth<-700
	currWidth<-length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin
  if(currWidth<minWidth){
    leftMargin<-leftMargin+(minWidth-currWidth)/2
    rightMargin<-rightMargin+(minWidth-currWidth)/2
  }
  #定义公共画图参数======================================================================================
  res$commomParams$barWidth<-barWidth
  res$commomParams$barInterval<-barInterval
  res$commomParams$leftMargin<-leftMargin
  res$commomParams$rightMargin<-rightMargin
  res$commomParams$width<-(barWidth+barInterval)*length(cell_lines)
  #定义legend参数======================================================================================
  if(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin==minWidth){
    legendLeft=-leftMargin+30
  }else if(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin<870){
    legendLeft=(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin-minWidth)/2-leftMargin+30
  }else{
    legendLeft=0
  }

  
  
  ######
  
  res$plotElements$legend_text$plotType <-"drawLegend"
  
  res$plotElements$legend_text$left<-legendLeft
  res$plotElements$legend_text$comment$baseline<- baseline
  res$plotElements$legend_text$comment$text<-"* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001"
  
  baseline<-baseline+buckHeight+buckInterv
  
  res$plotElements$legend_CNV$plotType <-"drawLegend"
  res$plotElements$legend_CNV$title$text<-"CNV type:"
  res$plotElements$legend_CNV$title$baseline<- baseline
  res$plotElements$legend_CNV$title$left<- legendLeft
  lenOfTitle<-110
  legendOffset<-0
  res$plotElements$legend_CNV$legend[["-2"]]$baseline<- baseline
  res$plotElements$legend_CNV$legend[["-2"]]$size <- 15
  res$plotElements$legend_CNV$legend[["-2"]]$color<-DEL_color
  res$plotElements$legend_CNV$legend[["-2"]]$left<-legendLeft+lenOfTitle
  res$plotElements$legend_CNV$legend[["-2"]]$text<-"DEL"
  
  legendOffset<-legendOffset+getStrLen("DEL",0.5)*300+15
  res$plotElements$legend_CNV$legend[["-1"]]$baseline<-baseline
  res$plotElements$legend_CNV$legend[["-1"]]$size<-15
  res$plotElements$legend_CNV$legend[["-1"]]$color<-LOSS_color
  res$plotElements$legend_CNV$legend[["-1"]]$left<-legendLeft+lenOfTitle+legendOffset
  res$plotElements$legend_CNV$legend[["-1"]]$text<-"LOSS"
  
  legendOffset<-legendOffset+getStrLen("LOSS",0.5)*300+15
  res$plotElements$legend_CNV$legend[["0"]]$baseline<-baseline
  res$plotElements$legend_CNV$legend[["0"]]$size<-15
  res$plotElements$legend_CNV$legend[["0"]]$color<-NOR_color
  res$plotElements$legend_CNV$legend[["0"]]$left<-legendLeft+lenOfTitle+legendOffset
  res$plotElements$legend_CNV$legend[["0"]]$text<-"NOR"
  
  
  legendOffset<-legendOffset+getStrLen("NOR",0.5)*300+15
  res$plotElements$legend_CNV$legend[["1"]]$baseline<-baseline
  res$plotElements$legend_CNV$legend[["1"]]$size<-15
  res$plotElements$legend_CNV$legend[["1"]]$color<-GAIN_color
  res$plotElements$legend_CNV$legend[["1"]]$left<-legendLeft+lenOfTitle+legendOffset
  res$plotElements$legend_CNV$legend[["1"]]$text<-"GAIN"
  
  legendOffset<-legendOffset+getStrLen("GAIN",0.5)*300+15
  res$plotElements$legend_CNV$legend[["2"]]$baseline<-baseline
  res$plotElements$legend_CNV$legend[["2"]]$size<-15
  res$plotElements$legend_CNV$legend[["2"]]$color<-AMP_color
  res$plotElements$legend_CNV$legend[["2"]]$left<-legendLeft+lenOfTitle+legendOffset
  res$plotElements$legend_CNV$legend[["2"]]$text<-"AMP"
  
  baseline<-baseline+buckHeight+buckInterv
  
  res$plotElements$legend_mutation$plotType <-"drawLegend"
  res$plotElements$legend_mutation$title$text<-"Mutation type:"
  res$plotElements$legend_mutation$title$baseline<-baseline
  res$plotElements$legend_mutation$title$left<-legendLeft
  lenOfTitle<-150
  
  legendOffset<-0
  res$plotElements$legend_mutation$legend[["Mutation"]]$baseline<-baseline
  res$plotElements$legend_mutation$legend[["Mutation"]]$size<-15
  res$plotElements$legend_mutation$legend[["Mutation"]]$color<- Mutation_color
  res$plotElements$legend_mutation$legend[["Mutation"]]$left<-legendLeft+lenOfTitle
  res$plotElements$legend_mutation$legend[["Mutation"]]$text<-"Mutation"
  
  legendOffset<-legendOffset+getStrLen("Mutation",0.5)*300+15
  res$plotElements$legend_mutation$legend[["None"]]$baseline<-baseline
  res$plotElements$legend_mutation$legend[["None"]]$size<-15
  res$plotElements$legend_mutation$legend[["None"]]$color<- No_mutation_color
  res$plotElements$legend_mutation$legend[["None"]]$left<-legendLeft+lenOfTitle+legendOffset
  res$plotElements$legend_mutation$legend[["None"]]$text<-"No alteration"
  
  
  baseline<-baseline+buckHeight+buckInterv
  
  res$plotElements$legend_MSI$plotType <-"drawLegend"
  res$plotElements$legend_MSI$title$text<-"MSI Status:"
  res$plotElements$legend_MSI$title$baseline<-baseline
  res$plotElements$legend_MSI$title$left<-legendLeft
  lenOfTitle<-120
  
  legendOffset<-0
  res$plotElements$legend_MSI$legend[["MSI"]]$baseline<-baseline
  res$plotElements$legend_MSI$legend[["MSI"]]$size<-15
  res$plotElements$legend_MSI$legend[["MSI"]]$color<- MSI_color
  res$plotElements$legend_MSI$legend[["MSI"]]$left<-legendLeft+lenOfTitle
  res$plotElements$legend_MSI$legend[["MSI"]]$text<-"MSI"
  
  legendOffset<-legendOffset+getStrLen("MSI",0.5)*300+15
  res$plotElements$legend_MSI$legend[["MSS"]]$baseline<-baseline
  res$plotElements$legend_MSI$legend[["MSS"]]$size<-15
  res$plotElements$legend_MSI$legend[["MSS"]]$color<- MSS_color
  res$plotElements$legend_MSI$legend[["MSS"]]$left<-legendLeft+lenOfTitle+legendOffset
  res$plotElements$legend_MSI$legend[["MSS"]]$text<-"MSS"
  
  baseline<-baseline+buckHeight*1.5+buckInterv
  
  #-----------绘制正文的参数
  
  #mRNA
  res$plotElements$mRNA$plotType<-"numericalBarplot"
  res$plotElements$mRNA$pillarColor<- mRNA_pillar_color
  res$plotElements$mRNA$baseline<-baseline
  res$plotElements$mRNA$leftText <- paste0(gene," mRNA")
  res$plotElements$mRNA$data <- multiomics[,grep("mRNA",names(multiomics))]
  res$plotElements$mRNA$buckHeight<-buckHeight
  res$plotElements$mRNA$rightText <- ""
  res$plotElements$mRNA$pillarClassName<-"CellAndDrugResistTags"
  res$plotElements$mRNA$resize<-"false"
  
  baseline<-baseline+buckInterv+buckHeight
  
  #CNV
  lab<-c("-2","-1","0","1","2")
  color_tmp<-c(DEL_color,LOSS_color,NOR_color,GAIN_color,AMP_color)
  for(i in 1:length(lab)){
    
    res$plotElements[[paste0("CNV_",i)]]$plotType<-"numericalBarplot"
    res$plotElements[[paste0("CNV_",i)]]$pillarColor<-color_tmp[i]
    res$plotElements[[paste0("CNV_",i)]]$baseline<-baseline
    CNVdata_tmp<-rep(0,length(multiomics[,grep("CNV",names(multiomics))]))
    CNVdata_tmp[which(multiomics[,grep("CNV",names(multiomics))]==lab[i])]<-0.8
    res$plotElements[[paste0("CNV_",i)]]$data <- CNVdata_tmp
    res$plotElements[[paste0("CNV_",i)]]$buckHeight<-buckHeight
    res$plotElements[[paste0("CNV_",i)]]$pillarClassName<-"CellAndDrugResistTags"
    res$plotElements[[paste0("CNV_",i)]]$resize<-"false"
    if(i==1){
      res$plotElements[[paste0("CNV_",i)]]$leftText <- paste0(gene," CNV")
      res$plotElements[[paste0("CNV_",i)]]$rightText <- pro_new$wenben[grep("CNV",row.names(pro_new))]
      res$plotElements[[paste0("CNV_",i)]]$rightTextColor <- pro_new$col[grep("CNV",row.names(pro_new))]
      
    }
  }
  
  
  
  baseline<-baseline+buckInterv+buckHeight
  
  #Mutation
  lab<-c("Mutation","None")
  color_tmp<-c(Mutation_color,No_mutation_color)
  for(i in 1:length(lab)){
    res$plotElements[[paste0("mutation_",i)]]$plotType<-"numericalBarplot"
    res$plotElements[[paste0("mutation_",i)]]$pillarColor<-color_tmp[i]
    res$plotElements[[paste0("mutation_",i)]]$baseline<-baseline
    mutationdata_tmp<-rep(0,length(multiomics[,grep("mutation",names(multiomics))]))
    mutationdata_tmp[which(multiomics[,grep("mutation",names(multiomics))]==lab[i])]<-0.8
    res$plotElements[[paste0("mutation_",i)]]$data <- mutationdata_tmp
    res$plotElements[[paste0("mutation_",i)]]$buckHeight<-buckHeight
    res$plotElements[[paste0("mutation_",i)]]$pillarClassName<-"CellAndDrugResistTags"
    res$plotElements[[paste0("mutation_",i)]]$resize<-"false"
    if(i==1){
      res$plotElements[[paste0("mutation_",i)]]$leftText <- paste0(gene," Mutation")
      res$plotElements[[paste0("mutation_",i)]]$rightText <- pro_new$wenben[grep("mutation",row.names(pro_new))]
      res$plotElements[[paste0("mutation_",i)]]$rightTextColor <- pro_new$col[grep("mutation",row.names(pro_new))]
    }
  }
  baseline<-baseline+buckInterv+buckHeight
  
  #--protien
  if(length(grep("protein_Omics",ls()))==1){
    res$plotElements$protein$plotType<-"numericalBarplot"
    res$plotElements$protein$pillarColor<- protein_pillar_color
    res$plotElements$protein$baseline<-baseline
    res$plotElements$protein$leftText <- paste0(gene," protein")
    res$plotElements$protein$data <- multiomics[,grep("protein",names(multiomics))]
    res$plotElements$protein$buckHeight<- buckHeight
    res$plotElements$protein$rightText <- pro_new$wenben[grep("protein",row.names(pro_new))]
    res$plotElements$protein$rightTextColor <-pro_new$col[grep("protein",row.names(pro_new))]
    res$plotElements$protein$pillarClassName<-"CellAndDrugResistTags"
    #res$plotElements$protein$resize<-Ture
    baseline<-baseline+buckInterv+buckHeight
  }
  
  #methylation site
  library(rjson)
  #定义正文参数（甲基化位点）======================================================================================
  for(i in CpG ){
    res$plotElements[[i]]$plotType<-"numericalBarplot"
    res$plotElements[[i]]$data <- multiomics[,i]
    res$plotElements[[i]]$leftText <- i
    res$plotElements[[i]]$rightText <- pro_new$wenben[which(row.names(pro_new)==i)]
    res$plotElements[[i]]$rightTextColor <- pro_new$col[which(row.names(pro_new)==i)]
    res$plotElements[[i]]$pillarColor<- methylation_site_pillar_color
    res$plotElements[[i]]$baseline<- baseline
    res$plotElements[[i]]$buckHeight<- buckHeight
    res$plotElements[[i]]$leftTextClassName <- "MethSiteTags"
    res$plotElements[[i]]$pillarClassName<-"CellAndDrugResistTags"
    res$plotElements[[i]]$resize<-"false"
    methylation_site_annotation_file_name <- paste0("./data/GDSC/methylation_site_annotation/",substr(i,1,6),"_ann.RData")
    load(methylation_site_annotation_file_name)
    methylation_site_annotation <- list[[i]]
    names(methylation_site_annotation)<-list[["item"]]
    methylation_site_annotation[2] <- gsub(" ","",methylation_site_annotation[2])
    res$plotElements[[i]]$leftTextAdditionInfos <- toJSON(methylation_site_annotation)
    baseline<-baseline+buckInterv+buckHeight
  }
  #临床数据
  #TMB
  res$plotElements$TMB$plotType<-"numericalBarplot"
  res$plotElements$TMB$pillarColor<- TMB_pillar_color
  res$plotElements$TMB$baseline<- baseline
  res$plotElements$TMB$leftText <- "TMB"
  res$plotElements$TMB$rightText <- pro_new$wenben[grep("TMB",row.names(pro_new))]
  res$plotElements$TMB$rightTextColor <- pro_new$col[grep("TMB",row.names(pro_new))]
  #试运行
  res$plotElements$TMB$data <- log2(multiomics[,grep("TMB",names(multiomics))])
  res$plotElements$TMB$buckHeight<-buckHeight
  res$plotElements$TMB$pillarClassName<-"CellAndDrugResistTags"
  #res$plotElements$clinical$TMB$resize<-"false"
  
  baseline<-baseline+buckInterv+buckHeight
  #-----------ploidy
  res$plotElements$ploidy$plotType<-"numericalBarplot"
  res$plotElements$ploidy$pillarColor<- ploidy_pillar_color
  res$plotElements$ploidy$baseline<- baseline
  res$plotElements$ploidy$leftText <- "Ploidy"
  res$plotElements$ploidy$rightText <- pro_new$wenben[grep("ploidy",row.names(pro_new))]
  res$plotElements$ploidy$rightTextColor <- pro_new$col[grep("ploidy",row.names(pro_new))]
  res$plotElements$ploidy$data <- multiomics[,grep("ploidy",names(multiomics))]
  res$plotElements$ploidy$buckHeight<-buckHeight
  res$plotElements$ploidy$pillarClassName<-"CellAndDrugResistTags"
  #res$plotElements$clinical$ploidy$resize<-"false"
  
  baseline<-baseline+buckInterv+buckHeight
  
  #MSI
  lab<-c("MSI","MSS")
  color_tmp<-c(MSI_color,MSS_color)
  for(i in 1:length(lab)){
    res$plotElements[[paste0("MSI_",i)]]$plotType<-"numericalBarplot"
    res$plotElements[[paste0("MSI_",i)]]$pillarColor<-color_tmp[i]
    res$plotElements[[paste0("MSI_",i)]]$baseline<- baseline
    MSIdata_tmp<-rep(0,length(multiomics[,grep("MSI",names(multiomics))]))
    MSIdata_tmp[which(multiomics[,grep("MSI",names(multiomics))]==lab[i])]<-0.8
    res$plotElements[[paste0("MSI_",i)]]$data <- MSIdata_tmp
    res$plotElements[[paste0("MSI_",i)]]$buckHeight<-buckHeight
    res$plotElements[[paste0("MSI_",i)]]$pillarClassName<-"CellAndDrugResistTags"
    res$plotElements[[paste0("MSI_",i)]]$resize<-"false"
    if(i==1){
      res$plotElements[[paste0("MSI_",i)]]$leftText <- "MSI Status"
      res$plotElements[[paste0("MSI_",i)]]$rightText <- pro_new$wenben[grep("MSI",row.names(pro_new))]
      res$plotElements[[paste0("MSI_",i)]]$rightTextColor <- pro_new$col[grep("MSI",row.names(pro_new))]
    }
  }
  #收尾获得文档高度，并返回======================================================================================
  res$commomParams$DocHeight<- baseline+topMargin
  res$cell_line <- row.names(multiomics)
  return(res)
  
}  



##########
Boxplot_calculation <- function (Dataset="GDSC",Omics = "mRNA",gene= "TP53"){
  options(stringsAsFactors = F)
  cell_annotation_path <-  paste0("./data/",Dataset,"/",Dataset,"_cell_annotation.RData")
  load(cell_annotation_path)
  cell_annotation <- get(ls()[grep(paste0(Dataset,"_cell_annotation"),ls())])
  cell_annotation <- cell_annotation[,c("Cell_Name","cancer_type")]
  row.names(cell_annotation) <- cell_annotation$Cell_Name
  if (Omics  %in% c("Penetrance","Metastatic","Chromatin","Metabolomics","Metastatic")) {
    Omics_file_name <- paste0("./data/",Dataset,"/",Dataset,"_",Omics,".RData")
  }else{
    name_end <- ifelse(Omics =="methylation_site",6,2)
    Omics_file_name <- paste0("./data/",Dataset,"/",Omics,"/",substr(gene,1,name_end),".RData")
  }
  if (file.exists(Omics_file_name)) {
    load(Omics_file_name)
  }else{
    res <- list()
    res$error <- paste0("This ", Omics ," file could not be found")
    return(res)
  }
  if (Omics %in% c("Penetrance","Metastatic","Chromatin","Metabolomics")) {
    Single_Omics <- get(ls()[ls()==paste0(Dataset,"_Omics")])
    Single_Omics$cell_lines <- row.names(Single_Omics)
    Single_Omics <- Single_Omics[,c("cell_lines",gene)]
    names(Single_Omics)[names(Single_Omics)==gene] <- paste0(gene,"_",Omics)
    names(Single_Omics)[names(Single_Omics)=="cell_lines"] <- paste0("cell_lines","_",Omics)
  }else{
    if( gene %in% names(list)){
      Single_Omics <- as.data.frame(list[c("cell_lines",gene)])
      row.names(Single_Omics) <- Single_Omics$cell_lines
      if(length(grep("-",gene))==1){
        names(Single_Omics)[-grep("cell_lines",names(Single_Omics))] <- 
          gsub("\\.","-",names(Single_Omics)[-grep("cell_lines",names(Single_Omics))])
      }
      names(Single_Omics)[names(Single_Omics)==gene] <- paste0(gene,"_",Omics)
      names(Single_Omics)[names(Single_Omics)=="cell_lines"] <- paste0("cell_lines","_",Omics)
    }else{
      res <- list()
      res$error <- paste0("This selected gene is not found in ",Dataset, " ", Omics," dataset")
      return(res)
    }}
  cross <- intersect(row.names(Single_Omics),row.names(cell_annotation))
  Single_Omics <- Single_Omics[cross,]
  cell_annotation <- cell_annotation[cross,]
  Single_Omics <- cbind(cell_annotation,Single_Omics)
  Single_Omics <- Single_Omics[,-grep("Cell",names(Single_Omics),ignore.case = T)]
  Single_Omics <- Single_Omics[!is.na(Single_Omics[,grep(gene,names(Single_Omics))]),]
  if(Omics == "microRNA"){
    names(Single_Omics)[names(Single_Omics)== paste0(gene,"_",Omics)] <- 
      gsub("-miR-","-mir-",names(Single_Omics)[names(Single_Omics)== paste0(gene,"_",Omics)] )
  }
  ####
  Freq <- as.data.frame.array(table(Single_Omics$cancer_type))
  names(Freq) <- "Freq"
  naru <- row.names(Freq)[Freq$Freq>3]
  Single_Omics <- Single_Omics[Single_Omics$cancer_type%in% naru,]
  ###
  
  if(Omics == "mRNA"){
  ylab <- paste0("Expression level of ", gene)
}else if(Omics == "methylation_promoter"){
  ylab <- paste0("promoter DNA methylation level of ", gene)
}else if(Omics == "methylation_body"){
  ylab <- paste0("Gene body DNA methylation level of ", gene)
}else if(Omics == "methylation_site"){
  ylab <- paste0("DNA methylation level of ", gene)
}else if(Omics == "protein"){
  ylab <- paste0("Protein level of ", gene)
}else if(Omics == "microRNA"){
  ylab <- paste0("Expression level of ", gene)
}else if(Omics == "Penetrance"){
  ylab <- "The penetrance of cell lines"
}else if(Omics == "Metastatic"){
  if(gene != "all"){
  ylab <- paste0("The ",gene," metastatic potential of cell lines")
  }else{
  ylab <- paste0("The total metastatic potential of cell lines")
  }
}else if(Omics == "RNAi_combine"){
  ylab <- paste0("RNAi dependency of ",gene)
}else if(Omics == "CRISPR"){
  ylab <- paste0("CRISPR-Cas9 dependency of ",gene)
}else if(Omics == "Chromatin"){
  ylab <- paste0("The level of ",gene ," modification")
}else if(Omics == "Metabolomics"){
  ylab <- paste0("The Metabolomic level of ",gene)
}else if(Omics == "CNV_continuous"){
    ylab <- paste0("The CNV level of ",gene)
  }
  calcData<-list()
  calcData$data<-Single_Omics
  calcData$ylab<-ylab
  
  return(boxplot_translate(calcData))
}
#######

SummBarplot_calculation<-function(Dataset="GDSC",gene="TP53",Omics="CNV",
                                  Mutation_color="#008000",
                                  DEL_color = "#000091",LOSS_color = "#4343b5",
                                  GAIN_color = "#ce5555",AMP_color = "#ba0000"){
  options(stringsAsFactors = F)
  cell_annotation_path <-  paste0("./data/",Dataset,"/",Dataset,"_cell_annotation.RData")
  load(cell_annotation_path)
  cell_annotation <- get(ls()[grep(paste0(Dataset,"_cell_annotation"),ls())])
  cell_annotation <- cell_annotation[,c("Cell_Name","cancer_type")]
  row.names(cell_annotation) <- cell_annotation$Cell_Name
  name_end <- ifelse(Omics =="methylation_site",6,2)
  Omics_file_name <- paste0("./data/",Dataset,"/",Omics,"/",substr(gene,1,name_end),".RData")
  if (file.exists(Omics_file_name)) {
    load(Omics_file_name)
  }else{
    res <- list()
    res$error <- paste0("This ", Omics ," file could not be found")
    return(res)
  }
  if( gene %in% names(list)){
    Single_Omics <- as.data.frame(list[c("cell_lines",gene)])
    row.names(Single_Omics) <- Single_Omics$cell_lines
    names(Single_Omics)[-grep("cell_lines",names(Single_Omics))] <- 
      gsub("\\.","-",names(Single_Omics)[-grep("cell_lines",names(Single_Omics))])
    names(Single_Omics)[names(Single_Omics)==gene] <- paste0(gene,"_",Omics)
    names(Single_Omics)[names(Single_Omics)=="cell_lines"] <- paste0("cell_lines","_",Omics)
  }else{
    res <- list()
    res$error <- paste0("This selected gene is not found in ",Dataset, " ", Omics," dataset")
    return(res)
  }
  cross <- intersect(row.names(Single_Omics),row.names(cell_annotation))
  Single_Omics <- Single_Omics[cross,]
  cell_annotation <- cell_annotation[cross,]
  Single_Omics <- cbind(cell_annotation,Single_Omics)
  Single_Omics <- Single_Omics[,-grep("Cell",names(Single_Omics),ignore.case = T)]
  BarData<-list()
  if(Omics=="CNV"){
    index <- which(names(Single_Omics)==paste0(gene,"_",Omics))
    Single_Omics[,index][Single_Omics[,index]=="-1" ] <- "Loss"
    Single_Omics[,index][Single_Omics[,index]=="-2" ] <- "Deletion"
    Single_Omics[,index][Single_Omics[,index]=="1" ] <- "Gain"
    Single_Omics[,index][Single_Omics[,index]=="0" ] <- "Normal"
    Single_Omics[,index][Single_Omics[,index]=="2" ] <- "Amplification"
    Single_Omics_paixu <- Single_Omics
    Single_Omics_paixu[,index] <- ifelse(Single_Omics_paixu[,index]  =="Normal","None","Altered")
    Single_Omics_paixu <- table(Single_Omics_paixu$cancer_type,Single_Omics_paixu[,index])
    Single_Omics_paixu<- as.data.frame.array(Single_Omics_paixu)
    Single_Omics_paixu$Total <- rowSums(Single_Omics_paixu)
    Single_Omics_paixu$frequency <- Single_Omics_paixu$Altered/Single_Omics_paixu$Total
    Single_Omics_paixu <- Single_Omics_paixu[order(Single_Omics_paixu$frequency,decreasing = T),]
    BarData$Type<-row.names(Single_Omics_paixu)
    Total_change_Number <- sum(Single_Omics[,paste0(gene,"_",Omics)] !="Normal")
    if(Total_change_Number>0){
      Total_change_Number <- Total_change_Number
    }else{
      res <- list()
      res$error <- paste0("There is no CNV alteration of this selected gene in ",Dataset," dataset")
      return(res)
    }
    for(i in 1:length(BarData$Type)){
      Single_Omics_tmp <- Single_Omics[Single_Omics$cancer_type == BarData$Type[i],]
      tmp <- as.data.frame(table(Single_Omics_tmp[, paste0(gene,"_",Omics)]))
      tmp$Var1 <- as.character(tmp$Var1 )
      tmp$frequency <- round(tmp$Freq/dim(Single_Omics_tmp)[1],3)
      tmp <- tmp[tmp$Var1 !="Normal",]
      queshi <- setdiff(c("Amplification","Gain","Loss","Deletion"),tmp$Var1)
      if(length(queshi !=0)){
        pro <- as.data.frame(matrix(NA,length(queshi),3))
        names(pro) <- names(tmp)
        pro$Var1 <- queshi
        pro$Freq <- rep(0,length(queshi))
        pro$frequency <- rep(0,length(queshi))
        tmp <- rbind(tmp,pro)
      }else{
        tmp <- tmp
      }
      row.names(tmp) <- tmp$Var1
      tmp <- tmp[c("Deletion","Loss","Gain","Amplification"),]
      data_tmp<-lapply(tmp$frequency,function(x)x)
      names(data_tmp)<-c("Deletion","Loss","Gain","Amplification")
      BarData$data[[BarData$Type[i]]] <- data_tmp
	  color_tmp<-lapply(c(DEL_color,LOSS_color,GAIN_color,AMP_color),function(x)x)
      names(color_tmp)<-c("Deletion","Loss","Gain","Amplification")
      BarData$colors[[BarData$Type[i]]] <-color_tmp
      Total_change_percent_type <- paste0(round(sum(tmp$Freq)/dim(Single_Omics_tmp)[1]*100,2),"%")
      BarData$Main_Text[[BarData$Type[i]]] <- paste0(gene," altered in ", sum(tmp$Freq),  " of ", dim(Single_Omics_tmp)[1]," cases (",Total_change_percent_type,")")
      Index <- which(names(Single_Omics_tmp)==paste0(gene,"_",Omics))
      Sub_Text<-list()
      cell_lines_text_tmp<-list()
      for (j in c("Deletion","Loss","Gain","Amplification")) {
        cell_lines_text_tmp[[j]] <- row.names(Single_Omics_tmp)[Single_Omics_tmp[Index]==j]
        change_percent_type <- paste0(round(tmp$Freq[tmp$Var1==j]/dim(Single_Omics_tmp)[1]*100,2),"%")
        change_Num_type <- tmp$Freq[tmp$Var1==j]
        case <- ifelse(change_Num_type >1,"cases","case")
        if(j == "Deletion"){
          Wenben <- "Deletion in "
        }else if(j == "Loss"){
          Wenben <- "Loss in "
        }else if(j == "Gain"){
          Wenben <- "Gain in "
        }else if(j == "Amplification"){
          Wenben <- "Amplification in "
        }
        Sub_Text[[j]] <- paste0(Wenben,change_percent_type," (",tmp$Freq[tmp$Var1==j]," ",case,")")
      }
      BarData$Sub_Text[[BarData$Type[i]]]<-Sub_Text
      BarData$cell_lines_text[[BarData$Type[i]]]<-cell_lines_text_tmp
	  BarData$legend<-color_tmp
    }
  }else{
    index <- which(names(Single_Omics)==paste0(gene,"_",Omics))
    Single_Omics_paixu <- table(Single_Omics$cancer_type,Single_Omics[,index])
    Single_Omics_paixu<- as.data.frame.array(Single_Omics_paixu)
    Single_Omics_paixu$Total <- rowSums(Single_Omics_paixu)
    Single_Omics_paixu$frequency <- Single_Omics_paixu$Mutation/Single_Omics_paixu$Total
    Single_Omics_paixu <- Single_Omics_paixu[order(Single_Omics_paixu$frequency,decreasing = T),]
    BarData$Type<-row.names(Single_Omics_paixu)
    Total_Mutation_number <- length(row.names(Single_Omics)[Single_Omics[,grep(gene,names(Single_Omics))]=="Mutation"])
    if(Total_Mutation_number>0){
      Total_Mutation_number <- Total_Mutation_number
    }else{
      res <- list()
      res$error <- paste0("There is no mutation alteration of this selected gene in ",Dataset," dataset")
      return(res)
    }
    for(i in 1:length(BarData$Type)){
      Single_Omics_tmp <- Single_Omics[Single_Omics$cancer_type == BarData$Type[i],]
      Mutation_number_type <- length(row.names(Single_Omics_tmp)[Single_Omics_tmp[,grep(gene,names(Single_Omics_tmp))]=="Mutation"])
      Total_number_type <- dim(Single_Omics_tmp)[1]
      BarData$data[[BarData$Type[i]]]$mut <-Mutation_number_type/Total_number_type
      BarData$colors[[BarData$Type[i]]]$mut <-Mutation_color
      BarData$Main_Text[[BarData$Type[i]]] <- paste0(gene," altered in ", Mutation_number_type,  " of ", Total_number_type," cases")
      Index <- which(names(Single_Omics_tmp)==paste0(gene,"_",Omics))
      BarData$cell_lines_text[[BarData$Type[i]]]$mut <- row.names(Single_Omics_tmp)[Single_Omics_tmp[Index]=="Mutation"]
    }
  }
  BarData$Y_text <- "Alteration Frequency"
  BarData$Omics <- Omics
  #将结果转化成图形
  res<-SummBarplot_translate(BarData)
  
  return(res)
}
############
##########
Subcluster_calculation <- function(Dataset="GDSC",Omics = "mRNA",
                                   genesets=  c("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF",
                                                "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", 
                                                "TUBB6" , "VEGFA"),cluster_algorithm="km",
                                   Min_color = "#1C577E",Mid_color = "white",
                                   Max_color = "#B51F23",Num_clust=3,
                                   cell_lines =c("A253","BB30-HNC","BB49-HNC","BHY","BICR10",
                                                 "BICR22","BICR31","BICR78","Ca9-22","CAL-27",
                                                 "CAL-33","Detroit562","DOK","FADU","NCI-H3118",
                                                 "HN","HO-1-N-1","HO-1-u-1","HSC-2","HSC-3",
                                                 "HSC-4","JHU-011","JHU-022","KON","KOSC-2",
                                                 "LB771-HNC","OSC-19","OSC-20","PCI-15A","PCI-30",
                                                 "PCI-38","PCI-4B","PCI-6A","PE-CA-PJ15","RPMI-2650",
                                                 "SAS","SAT","SCC-15","SCC-25","SCC-4",
                                                 "SCC-9","SKN-3","COLO-680N","EC-GI-10","ESO26",
                                                 "ESO51","FLO-1","HCE-4","KYAE-1","KYSE-140",
                                                 "KYSE-150","KYSE-180","KYSE-220","KYSE-270","KYSE-410",
                                                 "KYSE-450","KYSE-50","KYSE-510","KYSE-520","KYSE-70",
                                                 "OACM5-1","OACp4C","OE19","OE21","OE33",
                                                 "SK-GT-4","TE-1","TE-10","TE-11","TE-12",
                                                 "TE-15","TE-4","TE-5","TE-6","TE-8",
                                                 "TE-9","T-T","CESS","CTV-1","GDM-1",
                                                 "HEL","HL-60","KASUMI-1","KG-1","KMOE-2",
                                                 "ME-1","ML-2","MOLM-13","MONO-MAC-6","NKM-1",
                                                 "NOMO-1","OCI-AML2","OCI-AML3","OCI-AML5","OCI-M1",
                                                 "P31-FUJ","PL-21","QIMR-WIL","SIG-M5","THP-1"),
								txtPath="../plotCache/clust_result.txt"){
  options(stringsAsFactors = F)
  if(Omics =="methylation_site"){
    gene_information_path <- paste0("data/",Dataset,"/",Dataset,"_methylation_site_infomation.RData")
    load(gene_information_path)
    gene_information <- get(ls()[grep(paste0(Dataset,"_methylation_site_infomation"),ls())])
    chaji <- setdiff(genesets,gene_information) 
    genesets <- intersect(genesets,gene_information)
  }else{
    gene_information_path <- paste0("data/",Dataset,"/",Dataset,"_gene_information.RData")
    load(gene_information_path)
    gene_information <- get(ls()[grep(paste0(Dataset,"_gene_information"),ls())])
    chaji <- setdiff(genesets,gene_information[[Omics]]$Gene) 
    genesets <- intersect(genesets,gene_information[[Omics]]$Gene)
  }
  name_end <- ifelse(Omics =="methylation_site",6,2)
  Omics_file_name <- paste0("./data/",Dataset,"/",Omics,"/",substr(genesets[1],1,name_end),".RData")
  load(Omics_file_name)
  Final_Omics <- as.data.frame(list[c("cell_lines",genesets[1])])
  row.names(Final_Omics) <- Final_Omics$cell_lines
  names(Final_Omics)[names(Final_Omics)=="cell_lines"] <- 
    paste0("cell_lines_",genesets[1])
  if (length(genesets) >1 ) {
    for (i in 2:length(genesets)) {
      Omics_file_name <- paste0("./data/",Dataset,"/",Omics,"/",substr(genesets[i],1,name_end),".RData")
      load(Omics_file_name)
      Final_Omics_tmp <- as.data.frame(list[c("cell_lines",genesets[i])])
      row.names(Final_Omics_tmp) <- Final_Omics_tmp$cell_lines
      names(Final_Omics_tmp)[names(Final_Omics_tmp)=="cell_lines"] <- 
        paste0("cell_lines_",genesets[i])
      Final_Omics <- cbind(Final_Omics,Final_Omics_tmp)
    }}else{
      Final_Omics <- Final_Omics
    }
  Final_Omics <- Final_Omics[cell_lines,-grep("cell_lines",names(Final_Omics))]
  names(Final_Omics) <- gsub("\\.","-",names(Final_Omics))
  Final_Omics <- Final_Omics[rowMeans(is.na(Final_Omics))< 0.15,]
  Final_Omics_scale <- as.data.frame(scale(Final_Omics))
  
  if(sum(is.na(Final_Omics_scale))>0){
    library("impute")
    mat=as.matrix(Final_Omics_scale)
    mat=impute.knn(mat)
    Final_Omics_scale= as.data.frame(mat$data)
  }
  if(cluster_algorithm =="km"){
    kc <- kmeans(Final_Omics_scale, Num_clust)
    clust_result <- kc$cluster
  }else if(cluster_algorithm =="hc"){
    dist.r = dist(Final_Omics_scale, method="euclidean")
    hc.r = hclust(dist.r, method = "ward.D") 
    clust_result <- cutree(hc.r, k = Num_clust)
  }else if(cluster_algorithm =="pam"){
    library(cluster)
    fit_pam=pam(Final_Omics_scale,Num_clust)
    clust_result <- fit_pam$clustering
  }
  clust_result <- clust_result[order(clust_result)]
  clust_result <- as.data.frame(clust_result)
  clust_result_tmp <- as.data.frame(table(clust_result))
  Num_while <- 1
  while (sum(clust_result_tmp$Freq<=1)>0) {
    if(cluster_algorithm =="km"){
      kc <- kmeans(Final_Omics_scale, Num_clust)
      clust_result <- kc$cluster
    }else if(cluster_algorithm =="hc"){
      dist.r = dist(Final_Omics_scale, method="euclidean")
      hc.r = hclust(dist.r, method = "ward.D") 
      clust_result <- cutree(hc.r, k = Num_clust)
    }else if(cluster_algorithm =="pam"){
      library(cluster)
      fit_pam=pam(Final_Omics_scale,Num_clust)
      clust_result <- fit_pam$clustering
    }
    clust_result_tmp <- as.data.frame(table(clust_result))
    Num_while <- Num_while+1
    if(Num_while>3){
      break
    }
  }
  res <- list()
  if(sum(clust_result_tmp$Freq<=1)>0){
    res$error <- paste0("There is only one cell line in a sub-group of the clustering result. Please try again or select other cell lines.")
    return(res)
  }
  clust_result <- as.data.frame(clust_result)
  clust_result_table <- clust_result
  clust_result_table$cell <- row.names(clust_result_table)
  clust_result_table$clust_result <- as.roman(clust_result_table$clust_result)
  clust_result_table <- clust_result_table[,c(2,1)]
  write.csv(clust_result_table,file=txtPath, row.names =F,quote=F) 
  
  library(reshape2)
  library(tidyr)
  Final_Omics_scale$cell_lines <- row.names(Final_Omics_scale)
  Final_Omics_scale <- melt(Final_Omics_scale)
  names(Final_Omics_scale)[names(Final_Omics_scale)=="variable"] <- "Gene"
  names(Final_Omics_scale)[names(Final_Omics_scale)=="value"] <- "Scale_value"
  # 分割基因表达值，步长为0.01（如果希望颜色更加细腻可以步长缩短，但是没有必要）
  Final_Omics_scale$range <- cut(Final_Omics_scale$Scale_value,
                                 breaks = c(min(Final_Omics_scale$Scale_value)-0.01,seq(-1,
                                                                                   1,0.01),max(Final_Omics_scale$Scale_value)))  
  rangeMat1 <- levels(Final_Omics_scale$range)
  rbPal1 <- colorRampPalette(colors = c(Min_color,Mid_color,Max_color)) 
  col.vec1 <- rbPal1(length(rangeMat1))
  names(col.vec1) <- rangeMat1 # 产生配对的颜色向量
  Final_Omics_scale$color <- col.vec1[as.character(Final_Omics_scale$range)] # 匹配每个区间对应的颜色
  head(Final_Omics_scale)
  Final_Omics_scale <- Final_Omics_scale[,c("cell_lines","Gene","color")]
  Final_Omics_scale <- spread(Final_Omics_scale,cell_lines,color)
  row.names(Final_Omics_scale) <- Final_Omics_scale$Gene
  Final_Omics_scale <- Final_Omics_scale[,-1]
  color <- c("#1C577E", "#198888", "#B51F23", "#9C06B7")
  calcData<-list()
  for (i in 1:Num_clust) {
    j <- as.character(as.roman(i))
    calcData$data[[j]]$data <- Final_Omics_scale[, row.names(clust_result)[clust_result$clust_result==i]]
    calcData$data[[j]]$name <- j
    calcData$data[[j]]$color <-color[i]
	calcData$data[[j]]$cellNames<-row.names(clust_result)[clust_result$clust_result==i]
  }
  calcData$geneNames<- row.names(Final_Omics_scale)
  calcData$table<-clust_result
  calcData$txtPath<-txtPath
  calcData$color_range <- c(Min_color,Mid_color,Max_color)
  #save(clust_result,file="../plotCache/clust_result.rda")
  #将结果转化成图形
  res<-heatmap_translate(calcData)
  if(length(chaji)>0){
    #res <- list()
    res$error <- paste0(chaji, " is not found in ",Dataset," ",Omics ," dataset")
  }
  for (i in 1:Num_clust) {
    j <- as.character(as.roman(i))
    res[[j]]$cell_lines <- row.names(clust_result)[clust_result$clust_result==i]
  }
  res
  return(res)
}
###


Trans_Omics_calculation <- function (Dataset="GDSC",Omics = c("mRNA","mRNA","mRNA","methylation_promoter",
                                                              "methylation_promoter","methylation_promoter",
                                                              "methylation_body","methylation_body","methylation_body"
                                                              ,"CNV","CNV","CNV","mutation","mutation",
                                                              "mutation","methylation_site",
                                                              "methylation_site","methylation_site"),
                                     genes_sets = c("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF",
                                                    "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", 
                                                    "TUBB6" , "VEGFA","cg03641722","cg09436562","cg01316516"),
                                     reference=c("ACOT7","mRNA"),
                                     cell_lines =c("A253","BB30-HNC","BB49-HNC","BHY","BICR10",
                                                   "BICR22","BICR31","BICR78","Ca9-22","CAL-27",
                                                   "CAL-33","Detroit562","DOK","FADU","NCI-H3118",
                                                   "HN","HO-1-N-1","HO-1-u-1","HSC-2","HSC-3",
                                                   "HSC-4","JHU-011","JHU-022","KON","KOSC-2",
                                                   "LB771-HNC","OSC-19","OSC-20","PCI-15A","PCI-30",
                                                   "PCI-38","PCI-4B","PCI-6A","PE-CA-PJ15","RPMI-2650",
                                                   "SAS","SAT","SCC-15","SCC-25","SCC-4",
                                                   "SCC-9","SKN-3","COLO-680N","EC-GI-10","ESO26",
                                                   "ESO51","FLO-1","HCE-4","KYAE-1","KYSE-140",
                                                   "KYSE-150","KYSE-180","KYSE-220","KYSE-270","KYSE-410",
                                                   "KYSE-450","KYSE-50","KYSE-510","KYSE-520","KYSE-70",
                                                   "OACM5-1","OACp4C","OE19","OE21","OE33",
                                                   "SK-GT-4","TE-1","TE-10","TE-11","TE-12",
                                                   "TE-15","TE-4","TE-5","TE-6","TE-8",
                                                   "TE-9","T-T","CESS","CTV-1","GDM-1",
                                                   "HEL","HL-60","KASUMI-1","KG-1","KMOE-2",
                                                   "ME-1","ML-2","MOLM-13","MONO-MAC-6","NKM-1",
                                                   "NOMO-1","OCI-AML2","OCI-AML3","OCI-AML5","OCI-M1",
                                                   "P31-FUJ","PL-21","QIMR-WIL","SIG-M5","THP-1"),
                                     Omics_pillar_color = "#34495e",reference_pillar_color="#34495e",
                                     Mutation_color="#000091",No_mutation_color ="grey",
                                     DEL_color = "#000091",LOSS_color = "#4343b5",
                                     NOR_color = "#8e8e8e",GAIN_color = "#ce5555",AMP_color = "#ba0000",
                                     resistant_color="#ba0000",sensitive_color="#000091",
                                     cut_off=0.2,cor_method="spearman",
                                     statistic_method="FDR",mutation_method ="wilcoxon"){
  options(stringsAsFactors = F)
  names(genes_sets) <- Omics
  res<-list()
  gene_information_path <- paste0("data/",Dataset,"/",Dataset,"_gene_information.RData")
  load(gene_information_path)
  gene_information <- get(ls()[grep(paste0(Dataset,"_gene_information"),ls())])
  
  variable <- intersect(c("mRNA","methylation_promoter","RNAi_combine","CNV","methylation_body",
                          "methylation_site",
                          "CNV_continuous","mutation","CRSPR","protein","microRNA","Chromatin","Metabolomics",
                          "Metastatic"),Omics) 
  Final <- list()
  for (i in variable) {
    name_end <- ifelse(i =="methylation_site",6,2)
    if(i=="methylation_site"){
      methylation_site_infomation_path <- paste0("data/",Dataset,"/",Dataset,"_methylation_site_infomation.RData")
      load(methylation_site_infomation_path)
      methylation_site_infomation <- get(ls()[grep(paste0(Dataset,"_methylation_site_infomation"),ls())])
      present_genes <- genes_sets[names(genes_sets)==i]
      chaji <- setdiff(present_genes,methylation_site_infomation) 
      genesets_selected <- intersect(present_genes,methylation_site_infomation)
      file_name <- paste0("./data/",Dataset,"/", i,"/",substr(genesets_selected[1],1,name_end),".RData")
    }else{
      present_genes <- genes_sets[names(genes_sets)==i]
      chaji <- setdiff(present_genes,gene_information[[i]]$Gene) 
      genesets_selected <- intersect(present_genes,gene_information[[i]]$Gene)
      file_name <- paste0("./data/",Dataset,"/", i,"/",substr(genesets_selected[1],1,name_end),".RData")
    }
    if(length(genesets_selected)==0){
      res$warning <-   paste0(present_genes," is not found in GDSC ",i," Dataset")
      next
    }
    if (length(chaji)>0) {
      res[[paste0("error_",i)]] <-   paste0(chaji," is not found in GDSC ",i," Dataset")
    }
    load(file_name)
    Final_Omics <- as.data.frame(list[c("cell_lines",genesets_selected[1])])
    names(Final_Omics)[names(Final_Omics)!="cell_lines"] <- genesets_selected[1]
    row.names(Final_Omics) <- Final_Omics$cell_lines
    names(Final_Omics)[names(Final_Omics)=="cell_lines"] <- 
      paste0("cell_lines_",genesets_selected[1])
    if (length(genesets_selected) >1 ) {
      for (j in 2:length(genesets_selected)) {
        file_name <- paste0("./data/",Dataset,"/", i,"/",substr(genesets_selected[j],1,name_end),".RData")
        load(file_name)
        Final_Omics_tmp <- as.data.frame(list[c("cell_lines",genesets_selected[j])])
        names(Final_Omics_tmp)[names(Final_Omics_tmp)!="cell_lines"] <- genesets_selected[j]
        row.names(Final_Omics_tmp) <- Final_Omics_tmp$cell_lines
        names(Final_Omics_tmp)[names(Final_Omics_tmp)=="cell_lines"] <- 
          paste0("cell_lines_",genesets_selected[j])
        #print(sum(row.names(Final_Omics_tmp)==row.names(Final_Omics)))
        Final_Omics <- cbind(Final_Omics,Final_Omics_tmp)
      }}else{
        Final_Omics <- Final_Omics
      }
    #Final_Omics <- Final_Omics[,-grep("cell_lines",names(Final_Omics))]
    names(Final_Omics) <- gsub("\\.","-",names(Final_Omics))
    names(Final_Omics)[-grep("cell_lines",names(Final_Omics))] <- paste0(names(Final_Omics)[-grep("cell_lines",names(Final_Omics))],"&",i)
    names(Final_Omics)[grep("cell_lines",names(Final_Omics))] <- "cell_lines"
    Final_Omics <- Final_Omics[,!duplicated(names(Final_Omics))]
    Final[[i]] <- Final_Omics
  }
  if(Dataset %in% c("GDSC","NCI60")){
    TP53_mutation_file_name <- paste0("./data/",Dataset,"/mutation/TP.RData")
    load(TP53_mutation_file_name)
    Final_variable <- as.data.frame(list[c("cell_lines","TP53")])
    for (i in 1:length(Final)) {
      Final_variable <- dplyr::left_join(Final_variable,Final[[i]])
    }
    row.names(Final_variable) <- Final_variable[,which(names(Final_variable)=="cell_lines")]
    Final_variable <- Final_variable[,-c(which(names(Final_variable)=="cell_lines"),
                                         which(names(Final_variable)=="TP53"))]
  }else{
    load("data/DepMap/DepMap_demo_data.RData")
    Final_variable <- DepMap_demo_data
    for (i in 1:length(Final)) {
      Final_variable <- dplyr::left_join(Final_variable,Final[[i]])
    }
    row.names(Final_variable) <- Final_variable[,which(names(Final_variable)=="cell_lines")]
    Final_variable <- Final_variable[,-c(which(names(Final_variable)=="cell_lines"),
                                         which(names(Final_variable)=="demo_data"))]
  }
  multiomics <- Final_variable[cell_lines,]
  multiomics <- multiomics[rowMeans(is.na(multiomics))<1,]
  cell_lines <- row.names(multiomics)
  head(multiomics)
  reference_gene <- paste0(reference[1],"&",reference[2])
  
  #------------计算相关系数
  pro <- matrix(NA,dim(multiomics)[2]-1,2)
  row.names(pro) <- names(multiomics)[-grep(reference_gene,names(multiomics))]
  colnames(pro) <- c("cor","P")
  pro <- as.data.frame(pro)
  if(length(grep("mutation",row.names(pro)))>0){
    corrlation_variable <- row.names(pro)[-c(grep("mutation",row.names(pro)))]
  }else{
    corrlation_variable <- row.names(pro)
  }
  
  if(reference[2]=="mutation"){
    Mutation_cell <- row.names(multiomics)[multiomics[,reference_gene]=="Mutation"]
    None_Mutation_cell <-row.names(multiomics)[multiomics[,reference_gene]=="None"]
    if(length(Mutation_cell)==0){
      res <- list()
      res <-  paste0("There is no mutation of ",
                     reference[1]," in your selected cell lines")
      return(res)
    }else{
      for (i in corrlation_variable) {
        test_method <- ifelse(mutation_method =="wilcoxon","wilcox.test","t.test")
        tmp <- get(test_method)(multiomics[Mutation_cell,i],
                                multiomics[None_Mutation_cell,i])
        pro[i,2] <- tmp$p.value
      }
    }
    mutation_index <- setdiff(grep("mutation",names(multiomics)),grep(reference_gene,names(multiomics)))
    if(length(mutation_index)>0){
      for (j in mutation_index) {
        tmp <- chisq.test(multiomics[,reference_gene],multiomics[,j])
        pro[names(multiomics)[j],2] <- tmp$p.value
      }}
    pro_new <- pro[order(pro$P),]
    pro_new$FDR <- p.adjust(pro_new$P)
    pro_new$P_new  <- sprintf("%0.3f",  pro_new$P)
    na_index <- which(is.na(pro_new[,"P"]))
    if(statistic_method=="FDR"){
      pro_new$biaoji <- ifelse(pro_new$FDR < 0.0001,"****",ifelse(pro_new$FDR < 0.001,"***",
                                                                  ifelse(pro_new$FDR < 0.01,"**",ifelse(pro_new$FDR < 0.05,"*",""))))
      pro_new$col <-  ifelse(pro_new$FDR < 0.05,"black","#928a97")
      pro_new$FDR_new <-  sprintf("%0.3f",  pro_new$FDR)
      pro_new$wenben <- paste0("FDR = ", pro_new$FDR_new ,pro_new$biaoji)
    }else{
      pro_new$biaoji <- ifelse(pro_new$P < 0.0001,"****",ifelse(pro_new$P < 0.001,"***",
                                                                ifelse(pro_new$P < 0.01,"**",ifelse(pro_new$P < 0.05,"*",""))))
      pro_new$col <-  ifelse(pro_new$P < 0.05,"black","#928a97")
      pro_new$wenben <- paste0("P = ", pro_new$P_new ,pro_new$biaoji)
    }
    pro_new$wenben[na_index] <- ""
  }else{
    if(length(corrlation_variable)>0){
      for (i in corrlation_variable) {
        tmp <- cor.test(as.numeric(multiomics[,i]),
                        as.numeric(multiomics[,grep(reference_gene,names(multiomics))]),
                        method = cor_method)
        pro[i,1] <-tmp$estimate
        pro[i,2]<- tmp$p.value
      }
    }
    
    mutation_index <- grep("mutation",names(multiomics))
    if(length(mutation_index)>0){
      for (i in mutation_index) {
        Mutation_cell <- row.names(multiomics)[multiomics[,i]=="Mutation"]
        None_Mutation_cell <-row.names(multiomics)[multiomics[,i]=="None"]
        if( length(Mutation_cell) == 0 ){
          res$error_mutation <- paste0("There is no mutation of ",
                                       gsub("_mutation","",names(multiomics)[i])," in your selected cell lines")
          pro[names(multiomics)[i],2] <- NA
        }else{
          test_method <- ifelse(mutation_method =="wilcoxon","wilcox.test","t.test")
          tmp <- get(test_method)(multiomics[Mutation_cell,reference_gene],
                                  multiomics[None_Mutation_cell,reference_gene])
          pro[names(multiomics)[i],2] <- tmp$p.value
        }
      }}
    pro$cor1 <- sprintf("%0.3f", pro$cor)
    pro_new <- pro[order(pro$P),]
    pro_new$FDR <- p.adjust(pro_new$P)
    pro_new$P_new  <- sprintf("%0.3f",  pro_new$P)
    pro_new$FDR_new  <- sprintf("%0.3f",  pro_new$FDR)
    pro_new$col <- ifelse(abs(pro_new$cor) < cut_off,"#928a97","black")
    if(statistic_method=="P"){
      pro_new$biaoji <- ifelse(pro_new$P < 0.0001,"****",ifelse(pro_new$P < 0.001,"***",
                                                                ifelse(pro_new$P < 0.01,"**",ifelse(pro_new$P < 0.05,"*",""))))
      pro_new[grep("mutation",row.names(pro_new)),"col"] <-  ifelse(pro_new[grep("mutation",row.names(pro_new)),"P"] < 0.05,
                                                                    "black","#928a97")
      
    }else{pro_new$biaoji <- ifelse(pro_new$FDR < 0.0001,"****",ifelse(pro_new$FDR < 0.001,"***",
                                                                      ifelse(pro_new$FDR < 0.01,"**",ifelse(pro_new$FDR < 0.05,"*",""))))
    pro_new[grep("mutation",row.names(pro_new)),"col"] <-  ifelse(pro_new[grep("mutation",row.names(pro_new)),"FDR"] < 0.05,
                                                                  "black","#928a97")												  
    }
    mutaion_location <- grep("mutation",row.names(pro_new))
    
    
    pro_new$qianzui <- NA
    pro_new$wenben <- NA
    na_index <- which(is.na(pro_new[,"P"]))
    if (length(mutaion_location)>0) {
      pro_new[-mutaion_location,"qianzui"] <- ifelse(cor_method=="spearman","Rho = ","r = ")
      if(statistic_method=="P"){pro_new[mutaion_location,"qianzui"] <- "P = "
      pro_new[mutaion_location,"wenben"] <- 
        paste0("P = ",pro_new[mutaion_location,"P_new"],pro_new[mutaion_location,"biaoji"])
      }else{pro_new[mutaion_location,"qianzui"] <- "FDR = "
      pro_new[mutaion_location,"wenben"] <- 
        paste0("FDR = ",pro_new[mutaion_location,"FDR_new"],pro_new[mutaion_location,"biaoji"])}
      pro_new[-mutaion_location,"wenben"] <- 
        paste0(pro_new[-mutaion_location,"qianzui"],pro_new[-mutaion_location,"cor1"],
               pro_new[-mutaion_location,"biaoji"])
      
    }else{
      pro_new[,"qianzui"] <- ifelse(cor_method=="spearman","Rho = ","r = ")
      pro_new[,"wenben"] <- 
        paste0(pro_new[,"qianzui"],pro_new[,"cor1"],
               pro_new[,"biaoji"])
    }
    
   
    
    
    pro_new$wenben <- gsub("= 0.000","< 0.001",pro_new$wenben)
    pro_new$wenben[na_index] <- ""
  }
  
  multiomics <- multiomics[order(multiomics[,reference_gene],na.last = F),]
  
  if(reference[2] %in% c("mRNA","methylation_promoter","RNAi_combine","methylation_body",
                         "methylation_site",
                         "CNV_continuous","CRSPR","protein","microRNA","Chromatin","Metabolomics",
                         "Metastatic")){
    multiomics[!is.na(multiomics[,reference_gene]),reference_gene] <- 
      seq(0.3,0.7,length.out =length(  multiomics[!is.na(multiomics[,reference_gene]),reference_gene]))
  }
  if(reference[2] %in% c("CNV","mutation")){
    multiomics[is.na(multiomics[,reference_gene]),reference_gene] <- "No_alteration"
  }else{
    multiomics[is.na(multiomics[,reference_gene]),reference_gene] <- 0
  }
  
  multiomics[,-c(grep("mutation",names(multiomics)),grep("CNV",names(multiomics)))] <- 
    round(multiomics[,-c(grep("mutation",names(multiomics)),grep("CNV",names(multiomics)))],3)
  pro_new <- pro_new[names(multiomics)[-grep(reference_gene,names(multiomics))],]
  
  #---存储计算结果
  res<-list()
  buckHeight <- 25
  buckInterv <- buckHeight/10
  if(length(cell_lines)>600){
    barWidth = 2.25*600/length(cell_lines)
    barInterval = 0.75*600/length(cell_lines)
  }else{
    barWidth = 2.25
    barInterval = 0.75
  }
  rightMargin=180
  bottomMargin=10
  topMargin=25
  baseline=0+bottomMargin
  dpi <- 300
  ########
  if(Omics %in% c("methylation_site","microRNA","Metabolomics","Chromatin")){
    OmicsTag<-""
  }else if (Omics=="mRNA"){
    OmicsTag<-" mRNA"
  }else if (Omics=="mutation"){
    OmicsTag<-" mutation"
  }else if (Omics=="methylation_body"){
    OmicsTag<-" body"
  }else if(Omics=="methylation_promoter"){
    OmicsTag<-" promoter"
  }else if(Omics=="CNV"){
    OmicsTag<-" CNV"
  }else if(Omics=="protein"){
    OmicsTag<-" protein"
  }else if(Omics=="CRISPR"){
    OmicsTag<-" CRISPR"
  }else if(Omics=="RNAi_combine"){
    OmicsTag<-" RNAi"
  }else if(Omics=="CNV_continuous"){
    OmicsTag<-" CNV"
  }else if(Omics=="Metastatic"){
    OmicsTag<-" metastatisis"
  }
  
  
  leftMargin<-max(sapply(c(paste0(gene,OmicsTag)),function(x)getStrLen(x,0.5)))*dpi+5
  
  
  
  #########
  leftMargin<-max(sapply(names(multiomics),function(x)getStrLen(x,0.5)))*dpi+5
  #当画幅宽度小于570，则不再继续减小画幅宽度
  minWidth<-round(getStrLen("* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001",0.5)*dpi)
  
  minWidth<-700
  currWidth<-length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin
  if(currWidth<minWidth){
    leftMargin<-leftMargin+(minWidth-currWidth)/2
    rightMargin<-rightMargin+(minWidth-currWidth)/2
  }
  #定义公共画图参数======================================================================================
  res$commomParams$barWidth<-barWidth
  res$commomParams$barInterval<-barInterval
  res$commomParams$leftMargin<-leftMargin
  res$commomParams$rightMargin<-rightMargin
  res$commomParams$width<-(barWidth+barInterval)*length(cell_lines)
  #定义legend参数======================================================================================
  if(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin==minWidth){
    legendLeft=-leftMargin+30
  }else if(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin<870){
    legendLeft=(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin-minWidth)/2-leftMargin+30
  }else{
    legendLeft=0
  }
  
  #####
  
  res$plotElements$legend_text$plotType <-"drawLegend"
  
  res$plotElements$legend_text$left<-legendLeft
  res$plotElements$legend_text$comment$baseline<- baseline
  res$plotElements$legend_text$comment$text<-"* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001"
  
  baseline<-baseline+buckHeight+buckInterv
  if("CNV" %in% Omics){
    res$plotElements$legend_CNV$plotType <-"drawLegend"
    res$plotElements$legend_CNV$title$text<-"CNV type:"
    res$plotElements$legend_CNV$title$baseline<- baseline
    res$plotElements$legend_CNV$title$left<- legendLeft
    lenOfTitle<-110
    
	legendOffset<-0
    res$plotElements$legend_CNV$legend[["-2"]]$baseline<- baseline
    res$plotElements$legend_CNV$legend[["-2"]]$size <- 15
    res$plotElements$legend_CNV$legend[["-2"]]$color<-DEL_color
    res$plotElements$legend_CNV$legend[["-2"]]$left<-legendLeft+lenOfTitle
    res$plotElements$legend_CNV$legend[["-2"]]$text<-"DEL"
    
	legendOffset<-legendOffset+getStrLen("DEL",0.5)*300+15
    res$plotElements$legend_CNV$legend[["-1"]]$baseline<-baseline
    res$plotElements$legend_CNV$legend[["-1"]]$size<-15
    res$plotElements$legend_CNV$legend[["-1"]]$color<-LOSS_color
    res$plotElements$legend_CNV$legend[["-1"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend_CNV$legend[["-1"]]$text<-"LOSS"
    
	legendOffset<-legendOffset+getStrLen("LOSS",0.5)*300+15
    res$plotElements$legend_CNV$legend[["0"]]$baseline<-baseline
    res$plotElements$legend_CNV$legend[["0"]]$size<-15
    res$plotElements$legend_CNV$legend[["0"]]$color<-NOR_color
    res$plotElements$legend_CNV$legend[["0"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend_CNV$legend[["0"]]$text<-"NOR"
    
    
    legendOffset<-legendOffset+getStrLen("NOR",0.5)*300+15
    res$plotElements$legend_CNV$legend[["1"]]$baseline<-baseline
    res$plotElements$legend_CNV$legend[["1"]]$size<-15
    res$plotElements$legend_CNV$legend[["1"]]$color<-GAIN_color
    res$plotElements$legend_CNV$legend[["1"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend_CNV$legend[["1"]]$text<-"GAIN"
    
	legendOffset<-legendOffset+getStrLen("GAIN",0.5)*300+15
    res$plotElements$legend_CNV$legend[["2"]]$baseline<-baseline
    res$plotElements$legend_CNV$legend[["2"]]$size<-15
    res$plotElements$legend_CNV$legend[["2"]]$color<-AMP_color
    res$plotElements$legend_CNV$legend[["2"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend_CNV$legend[["2"]]$text<-"AMP"
    
    baseline<-baseline+buckHeight+buckInterv
  }
  if("mutation" %in% Omics){
    res$plotElements$legend_mutation$plotType <-"drawLegend"
    res$plotElements$legend_mutation$title$text<-"Mutation type:"
    res$plotElements$legend_mutation$title$baseline<-baseline
    res$plotElements$legend_mutation$title$left<-legendLeft
    lenOfTitle<-150
    
	legendOffset<-0
    res$plotElements$legend_mutation$legend[["Mutation"]]$baseline<-baseline
    res$plotElements$legend_mutation$legend[["Mutation"]]$size<-15
    res$plotElements$legend_mutation$legend[["Mutation"]]$color<- Mutation_color
    res$plotElements$legend_mutation$legend[["Mutation"]]$left<-legendLeft+lenOfTitle
    res$plotElements$legend_mutation$legend[["Mutation"]]$text<-"Mutation"
    
	legendOffset<-legendOffset+getStrLen("Mutation",0.5)*300+15
    res$plotElements$legend_mutation$legend[["None"]]$baseline<-baseline
    res$plotElements$legend_mutation$legend[["None"]]$size<-15
    res$plotElements$legend_mutation$legend[["None"]]$color<- No_mutation_color
    res$plotElements$legend_mutation$legend[["None"]]$left<-legendLeft+lenOfTitle+legendOffset
    res$plotElements$legend_mutation$legend[["None"]]$text<-"No alteration"
    
    
    baseline<-baseline+buckHeight+buckInterv
  }
  baseline<-baseline+buckHeight*0.5
  
  #-----------绘制正文的参数
  #定义正文参数（组学）======================================================================================
  
  if(reference[2] == "CNV"){
    lab<-c("-2","-1","0","1","2")
    color_tmp<-c(DEL_color,LOSS_color,NOR_color,GAIN_color,AMP_color)
    for(i in 1:length(lab)){
      res$plotElements[[paste0("reference_",i)]]$plotType<-"numericalBarplot"
      res$plotElements[[paste0("reference_",i)]]$pillarColor<-color_tmp[i]
      res$plotElements[[paste0("reference_",i)]]$baseline<-baseline
      CNV_data_tmp<-rep(0,length(multiomics[,reference_gene]))
      CNV_data_tmp[which(multiomics[,reference_gene]==lab[i])]<-0.8
      res$plotElements[[paste0("reference_",i)]]$data <- CNV_data_tmp
      res$plotElements[[paste0("reference_",i)]]$buckHeight<-buckHeight
      res$plotElements[[paste0("reference_",i)]]$pillarClassName<-"CellAndDrugResistTags"
      res$plotElements[[paste0("reference_",i)]]$resize<-"false"
      if(i==1){
        res$plotElements[[paste0("reference_",i)]]$leftText <- paste0(reference[1]," ",reference[2])
      }
    }
  }else if(reference[2] =="mutation"){
    lab<-c("Mutation","None")
    color_tmp<-c(Mutation_color,No_mutation_color)
    for(i in 1:length(lab)){
      res$plotElements[[paste0("reference_",i)]]$plotType<-"numericalBarplot"
      res$plotElements[[paste0("reference_",i)]]$pillarColor<-color_tmp[i]
      res$plotElements[[paste0("reference_",i)]]$baseline<-baseline
      mutation_data_tmp<-rep(0,length(multiomics[,reference_gene]))
      mutation_data_tmp[which(multiomics[,reference_gene]==lab[i])]<-0.8
      res$plotElements[[paste0("reference_",i)]]$data <- mutation_data_tmp
      res$plotElements[[paste0("reference_",i)]]$buckHeight<-buckHeight
      res$plotElements[[paste0("reference_",i)]]$pillarClassName<-"CellAndDrugResistTags"
      res$plotElements[[paste0("reference_",i)]]$resize<-"false"
      if(i==1){
        res$plotElements[[paste0("reference_",i)]]$leftText <- paste0(reference[1]," ",reference[2])
      }
    }
  }else{
    res$plotElements[["reference"]]$plotType<-"numericalBarplot"
    res$plotElements[["reference"]]$pillarColor<- reference_pillar_color
    res$plotElements[["reference"]]$baseline<-baseline
    if (reference[2] %in% c("methylation_site","microRNA","Chromatin","Metabolomics")) {
      res$plotElements[["reference"]]$leftText <- reference[1]
    }else if(reference[2] =="methylation_promoter"){
      res$plotElements[["reference"]]$leftText <- paste0(reference[1]," promoter") 
    }else if(reference[2] =="methylation_body"){
      res$plotElements[["reference"]]$leftText <- paste0(reference[1]," body") 
    }else if(reference[2] =="CNV_continuous"){
      res$plotElements[["reference"]]$leftText <- paste0(reference[1]," CNV") 
    }else if(reference[2] =="Metastatic"){
      res$plotElements[["reference"]]$leftText <- paste0(reference[1]," metastatic potential") 
    }else if(reference[2] =="RNAi_combine"){
      res$plotElements[["reference"]]$leftText <- paste0(reference[1]," RNAi dependency") 
    }else if(reference[2] =="CRISPR"){
      res$plotElements[["reference"]]$leftText <- paste0(reference[1]," CRISPR-Cas9 dependency") 
    }else{
      res$plotElements[["reference"]]$leftText <- paste0(reference[1]," ",reference[2]) 
    }
    res$plotElements[["reference"]]$data <-multiomics[,reference_gene]
    res$plotElements[["reference"]]$buckHeight<-buckHeight
    res$plotElements[["reference"]]$pillarClassName<-"CellAndDrugResistTags"
    res$plotElements[["reference"]]$resize<-"false"
  }
  if (reference[2] =="methylation_site") {
    methylation_site_annotation <- paste0("./data/GDSC/methylation_site_annotation/",substr(reference[1],1,6),"_ann.RData")
    load(methylation_site_annotation)
    final_annotation <- list[[reference[1]]]
    names(final_annotation)<-list[["item"]]
    final_annotation[2] <- gsub(" ","",final_annotation[2] )
    library(rjson)
    res$plotElements[["reference"]]$leftTextAdditionInfos <- toJSON(final_annotation)
    res$plotElements[["reference"]]$leftTextClassName <- "MethSiteTags"
  }
  baseline<-baseline+buckInterv+buckHeight
  
  for (j in row.names(pro_new)) {
    if (strsplit(j,"&")[[1]][2] =="CNV") {
      lab<-c("-2","-1","0","1","2","No_alteration")
      color_tmp<-c(DEL_color,LOSS_color,NOR_color,GAIN_color,AMP_color,"white")
      for(i in 1:length(lab)){
        
        res$plotElements[[paste0(j,"_",i)]]$plotType<-"numericalBarplot"
        res$plotElements[[paste0(j,"_",i)]]$pillarColor <- color_tmp[i]
        res$plotElements[[paste0(j,"_",i)]]$baseline<-baseline
        CNV_data_tmp <- rep(0,length(multiomics[,j]))
        CNV_data_tmp[which(multiomics[,j]==lab[i])]<-0.8
        res$plotElements[[paste0(j,"_",i)]]$data <- CNV_data_tmp
        res$plotElements[[paste0(j,"_",i)]]$buckHeight<-buckHeight
        res$plotElements[[paste0(j,"_",i)]]$pillarClassName<-"CellAndDrugResistTags"
        res$plotElements[[paste0(j,"_",i)]]$resize<-"false"
        if(i==1){
          res$plotElements[[paste0(j,"_",i)]]$leftText <- gsub("&"," ",j)
          res$plotElements[[paste0(j,"_",i)]]$rightText <- pro_new$wenben[row.names(pro_new)==j]
          res$plotElements[[paste0(j,"_",i)]]$rightTextColor <- pro_new$col[row.names(pro_new)==j]
        }}}else if(strsplit(j,"&")[[1]][2] =="mutation"){
          lab<-c("Mutation","None","No_alteration")
          color_tmp<-c(Mutation_color,No_mutation_color,"white")
          for(i in 1:length(lab)){
            res$plotElements[[paste0(j,"_",i)]]$plotType<-"numericalBarplot"
            res$plotElements[[paste0(j,"_",i)]]$pillarColor<-color_tmp[i]
            res$plotElements[[paste0(j,"_",i)]]$baseline<-baseline
            mutationdata_tmp<-rep(0,length(multiomics[,j]))
            mutationdata_tmp[which(multiomics[,j]==lab[i])]<-0.8
            res$plotElements[[paste0(j,"_",i)]]$data <- mutationdata_tmp
            res$plotElements[[paste0(j,"_",i)]]$buckHeight<-buckHeight
            res$plotElements[[paste0(j,"_",i)]]$pillarClassName<-"CellAndDrugResistTags"
            res$plotElements[[paste0(j,"_",i)]]$resize<-"false"
            if(i==1){
              res$plotElements[[paste0(j,"_",i)]]$leftText <- gsub("&"," ",j)
              res$plotElements[[paste0(j,"_",i)]]$rightText <- pro_new$wenben[row.names(pro_new)==j]
              res$plotElements[[paste0(j,"_",i)]]$rightTextColor <- pro_new$col[row.names(pro_new)==j]
            }}}else{
              res$plotElements[[j]]$plotType<-"numericalBarplot"
              res$plotElements[[j]]$data <- multiomics[,j]
              if(strsplit(j,"&")[[1]][2] == "methylation_site"){
                res$plotElements[[j]]$leftText <- gsub("&methylation_site","",j)
                
              }else if(strsplit(j,"&")[[1]][2] %in% c("methylation_promoter","methylation_body")){
                res$plotElements[[j]]$leftText <- gsub("&methylation_"," ",j)
              }else if(strsplit(j,"&")[[1]][2] == "microRNA"){
                res$plotElements[[j]]$leftText <- gsub("&microRNA","",j)
              }else if(strsplit(j,"&")[[1]][2] == "Chromatin"){
                res$plotElements[[j]]$leftText <- gsub("&Chromatin","",j)
              }else if(strsplit(j,"&")[[1]][2] == "Metabolomics"){
                res$plotElements[[j]]$leftText <- gsub("&Metabolomics","",j)
              }else if(strsplit(j,"&")[[1]][2] == "CNV_continuous"){
                res$plotElements[[j]]$leftText <- paste0(gsub("&CNV_continuous","",j)," CNV")
              }else if(strsplit(j,"&")[[1]][2] == "Metastatic"){
                library(stringr)
                res$plotElements[[j]]$leftText <- str_to_title(paste0(gsub("&"," ",j)," potential"))
              }else if(strsplit(j,"&")[[1]][2] == "RNAi_combine"){
                res$plotElements[[j]]$leftText <- paste0(gsub("&RNAi_combine"," ",j)," RNAi dependency")
              }else if(strsplit(j,"&")[[1]][2] == "CRISPR"){
                res$plotElements[[j]]$leftText <- paste0(gsub("&CRISPR"," ",j)," CRISPR-Cas9 dependency")
              }else{
                res$plotElements[[j]]$leftText <- gsub("&"," ",j)
              }
              
              res$plotElements[[j]]$rightText <-  pro_new$wenben[row.names(pro_new)==j]
              res$plotElements[[j]]$rightTextColor <- pro_new$col[row.names(pro_new)==j]
              res$plotElements[[j]]$pillarColor<- Omics_pillar_color
              res$plotElements[[j]]$baseline<- baseline
              res$plotElements[[j]]$buckHeight<- buckHeight
              if(strsplit(j,"&")[[1]][2] %in% c("methylation_promoter","methylation_body","methylation_site")){
                res$plotElements[[j]]$resize<-"false"
              }
              if(strsplit(j,"&")[[1]][2] =="methylation_site"){
                res$plotElements[[j]]$leftTextClassName <- "MethSiteTags"
                res$plotElements[[j]]$pillarClassName<-"CellAndDrugResistTags"
                methylation_site_annotation_file_name <- paste0("./data/GDSC/methylation_site_annotation/",substr(j,1,6),"_ann.RData")
                load(methylation_site_annotation_file_name)
                CpG <- gsub("&methylation_site","",j)
                methylation_site_annotation <- list[[CpG]]
                names(methylation_site_annotation)<-list[["item"]]
                methylation_site_annotation[2] <- gsub(" ","",methylation_site_annotation[2])
                res$plotElements[[j]]$leftTextAdditionInfos <- toJSON(methylation_site_annotation)
              }else{
                res$plotElements[[j]]$pillarClassName<-"CellAndDrugResistTags"
              }}
    baseline<-baseline+buckInterv+buckHeight
  }
  #收尾获得文档高度，并返回======================================================================================
  res$commomParams$DocHeight<- baseline+topMargin
  res$cell_line <- row.names(multiomics)
  return(res)
}


############
Drug_response_model_calculation <- function(Dataset="GDSC",Omics = "mRNA",
                                            genesets=  c("ACOT7", "ADM", "ALDOA", "CDKN3", "ENO1", "LDHA", "MIF",
                                                         "MRPS17", "NDRG1", "P4HA1", "PGAM1", "SLC2A1", "TPI1", 
                                                         "TUBB6" , "VEGFA"),Algorithm="Logistic",
                                            cell_lines=c("A253", "BB49-HNC","Detroit562"),
                                            drug_input	= "AZD5363"){
  options(stringsAsFactors = F)
  
  if(Omics =="methylation_site"){
    gene_information_path <- paste0("data/",Dataset,"/",Dataset,"_methylation_site_infomation.RData")
    load(gene_information_path)
    gene_information <- get(ls()[grep(paste0(Dataset,"_methylation_site_infomation"),ls())])
    chaji <- setdiff(genesets,gene_information) 
    genesets <- intersect(genesets,gene_information)
  }else{
    gene_information_path <- paste0("data/",Dataset,"/",Dataset,"_gene_information.RData")
    load(gene_information_path)
    gene_information <- get(ls()[grep(paste0(Dataset,"_gene_information"),ls())])
    chaji <- setdiff(genesets,gene_information[[Omics]]$Gene) 
    genesets <- intersect(genesets,gene_information[[Omics]]$Gene)
  }
  name_end <- ifelse(Omics =="methylation_site",6,2)
  Omics_file_name <- paste0("./data/",Dataset,"/",Omics,"/",substr(genesets[1],1,name_end),".RData")
  load(Omics_file_name)
  Final_Omics <- as.data.frame(list[c("cell_lines",genesets[1])])
  row.names(Final_Omics) <- Final_Omics$cell_lines
  names(Final_Omics)[names(Final_Omics)=="cell_lines"] <- 
    paste0("cell_lines_",genesets[1])
  if (length(genesets) >1 ) {
    for (i in 2:length(genesets)) {
      Omics_file_name <- paste0("./data/",Dataset,"/",Omics,"/",substr(genesets[i],1,name_end),".RData")
      load(Omics_file_name)
      Final_Omics_tmp <- as.data.frame(list[c("cell_lines",genesets[i])])
      row.names(Final_Omics_tmp) <- Final_Omics_tmp$cell_lines
      names(Final_Omics_tmp)[names(Final_Omics_tmp)=="cell_lines"] <- 
        paste0("cell_lines_",genesets[i])
      Final_Omics <- cbind(Final_Omics,Final_Omics_tmp)
    }}else{
      Final_Omics <- Final_Omics
    }
  Final_Omics <- Final_Omics[cell_lines,-grep("cell_lines",names(Final_Omics))]
  Final_Omics_scale <- as.data.frame(scale(Final_Omics))
  if(sum(is.na(Final_Omics_scale))>0){
    library("impute")
    mat=as.matrix(Final_Omics_scale)
    mat=impute.knn(mat)
    Final_Omics_scale= as.data.frame(mat$data)
  }
  load("data/GDSC/GDSC_drug_response_Binary.RData")
  GDSC_drug_response_Binary$cell_lines <- row.names(GDSC_drug_response_Binary)
  GDSC_drug_response_Binary <- GDSC_drug_response_Binary[cell_lines,c("cell_lines",drug_input)]
  tmp <- as.data.frame.array(table(GDSC_drug_response_Binary[,drug_input] ))
  #names(tmp)[1] <- "Freq"
  if(dim(tmp)[1]==1){
    res <- list()
    res$error <- paste0("There is olny one response of ", drug_input, " in selected cell lines; Please try another combination of pharmacogenomic parameters")
    return(res)
  }
  
  GDSC_drug_response_Binary[,drug_input] <- as.numeric(ifelse(GDSC_drug_response_Binary[,drug_input] =="sensitivity",1,0))
  Final_Omics_scale <- cbind(GDSC_drug_response_Binary,Final_Omics_scale)
  Final_Omics_scale <- Final_Omics_scale[,-which(names(Final_Omics_scale)=="cell_lines")]
  
  Final_Omics_scale <- Final_Omics_scale[!is.na(Final_Omics_scale[,drug_input]),]
  names(Final_Omics_scale) <- gsub("\\.","-",names(Final_Omics_scale))
  genesets1 <- genesets
  index_crossbar  <- grep("-",genesets1)
  if(length(index_crossbar)>0){
    for (i in index_crossbar) {
      genesets1[i] <- paste0("`",genesets1[i]  ,"`")
      
    }
  }
  if (length(grep("-",drug_input)) ==1) {
    fmla <- as.formula(paste0("`",drug_input ,"`","~", paste0(genesets1, collapse = "+")))
  }else{
    fmla <- as.formula(paste0(drug_input ,"~", paste0(genesets1, collapse = "+")))
  }
  library(caret)
  inValidation <- createDataPartition(y=Final_Omics_scale[,drug_input],p=0.3,list=F)
  
  train <- Final_Omics_scale[-inValidation,]
  vallidation <- Final_Omics_scale[inValidation,]
  
  res <- list()
    if (length(table(train[,drug_input]))==1) {
      res$error <- paste0("There is only one drug response of selected cell lines in train cohort; Please try another combination of pharmacogenomic parameters")
	  return(res)
    }
    if (length(table(vallidation[,drug_input]))==1) {
      res$error <- paste0("There is only one drug response of selected cell lines in validation cohort; Please try another combination of pharmacogenomic parameters")
	  return(res)
    }
  fc<-as.numeric()
  mod_pre<-as.numeric()
  if(Algorithm=="randomForest"){
    library(randomForest)
    Final_Omics_scale[,drug_input] <- as.factor(Final_Omics_scale[,drug_input])
    train <- Final_Omics_scale[-inValidation,]
    vallidation <- Final_Omics_scale[inValidation,]
    model<-randomForest(fmla,data=train,proximity=T,importance=T)
    model_pre<-predict(model,newdata = vallidation,type="prob")
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,model_pre[,1])
    }else if(Algorithm=="Logistic"){
    model<-glm(fmla,family=binomial(link=logit),data=train)
    model_pre<-predict(model,type='response',newdata=vallidation)
    fc<-append(fc,vallidation[,drug_input])
    mod_pre<-append(mod_pre,as.numeric(model_pre))
   }else if(Algorithm=="SVM"){
    library(e1071)
    model<-svm(fmla,data=train,probability=T)
    model_pre<-predict(model,newdata = vallidation,decision.values = TRUE, probability = TRUE)
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,as.numeric(model_pre))
    }else if(Algorithm=="KNN"){
    library(kknn)
    model<-kknn(fmla,train,vallidation,distance = 1, kernel = "triangular",
                  k=round(sqrt(dim(train)[1])))
    model_pre <- fitted(model)  
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,as.numeric(model_pre))
    }else if(Algorithm=="DecisionTree"){
    library(rpart)
    model <- rpart(fmla, method = "class", data = train,control=rpart.control(minsplit = 1))
    model_pre <- predict(model, newdata = vallidation,type = 'prob')
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,as.numeric(model_pre[,2]))
    }else if(Algorithm=="GBM"){
	library(gbm)
    library(caret)
    fitControl <- trainControl( method = "repeatedcv", number = 4, repeats = 4)
    model <- caret::train(fmla, data = train, method = "gbm", trControl = fitControl,verbose = FALSE)
    model_pre = predict(model, newdata = vallidation)
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,as.numeric(model_pre))
    }else if(Algorithm=="Adaboost"){
    library(adabag)
    Final_Omics_scale[,drug_input] <- as.factor(Final_Omics_scale[,drug_input])
    train <- Final_Omics_scale[-inValidation,]
    vallidation <- Final_Omics_scale[inValidation,]
    model <- boosting(fmla,data = train,boos=TRUE, mfinal=100 )
    model_pre <- predict(model,newdata = vallidation)$class
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,as.numeric(model_pre))
    }else if(Algorithm=="Ridge"){
    library(glmnet)
    lambdas <- seq(0,5, length.out = 200)
    index <- which(names(train)==drug_input)
    X =as.matrix(train[,-index])
    Y=train[,index]
    ridge_model <- cv.glmnet(X,Y,alpha = 0,lambda = lambdas,nfolds =3)
    ridge_min <- ridge_model$lambda.min
    ridge_best <- glmnet(X,Y,alpha = 0,lambda = ridge_min)
    model_pre <- predict(ridge_best,as.matrix(vallidation[,-index]))
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,as.numeric(model_pre))
    }else if(Algorithm=="Lasso"){
    library(glmnet)
    lambdas <- seq(0,5, length.out = 200)
    index <- which(names(train)==drug_input)
    X =as.matrix(train[,-index])
    Y=train[,index]
    lasso_model <- cv.glmnet(X,Y,alpha = 1,lambda = lambdas,nfolds =3)
    lasso_min <- lasso_model$lambda.min
    lasso_best <- glmnet(X,Y,alpha = 1,lambda = lasso_min)
    model_pre <- predict(lasso_best,as.matrix(vallidation[,-index]))
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,as.numeric(model_pre))
    }else if(Algorithm=="NNet"){
    library(nnet)
    model <- nnet(fmla,data = train,size = 2,rang = 0.1,decay = 5e-4,maxit = 200)
    model_pre = predict(model, newdata = vallidation)
    fc<-append(fc,as.numeric(vallidation[,drug_input]))
    mod_pre<-append(mod_pre,as.numeric(model_pre))
    }
  library(pROC)
  df<-cbind(fc,as.numeric(mod_pre))
  modelroc <-roc(df[,1],df[,2],ci=T)
  modelRocData<-list()
  modelRocData$specificities<-modelroc$specificities
  modelRocData$sensitivities<-modelroc$sensitivities
  modelRocData$ci<-as.numeric(modelroc$ci)
  modelRocData$model<-Algorithm
  return(modelRoc_translate(modelRocData))
}
###########################
DepMap_MultiOmics_calculation <- function (Dataset="GDSC", gene= "TP53",cell_lines=c("A253", "BB49-HNC","Detroit562") ,
                                    mRNA_pillar_color = "#34495e",methylation_site_pillar_color="#928a97",
                                    TMB_pillar_color="#34495e",ploidy_pillar_color="#928a97",
                                    Mutation_color="#000091",No_mutation_color ="grey",
                                    MSI_color="#000091",MSS_color ="grey",DEL_color = "#000091",
                                    LOSS_color = "#4343b5",NOR_color = "#8e8e8e",GAIN_color = "#ce5555",
                                    AMP_color = "#ba0000",protein_pillar_color = "#928a97",
                                    cut_off=0.2,cor_method="spearman",mutation_method="wilcoxon",statistic_method="FDR"){
  options(stringsAsFactors = F)
  library(dplyr)
  #-----mRNA
  mRNA_file_name <- paste0("./data/",Dataset,"/mRNA/",substr(gene,1,2),".RData")
  load(mRNA_file_name)
  mRNA_Omics <- as.data.frame(list[c("cell_lines",gene)])
  names(mRNA_Omics) <- gsub("\\.","-",names(mRNA_Omics))
  row.names(mRNA_Omics) <- mRNA_Omics$cell_lines
  names(mRNA_Omics)[names(mRNA_Omics)==gene] <- paste0(gene,"_mRNA")
  #-----CNV
  CNV_file_name <- paste0("./data/",Dataset,"/CNV_continuous/",substr(gene,1,2),".RData")
  load(CNV_file_name)
  CNV_Omics <- as.data.frame(list[c("cell_lines",gene)])
  names(CNV_Omics) <- gsub("\\.","-",names(CNV_Omics))
  row.names(CNV_Omics) <- CNV_Omics$cell_lines
  names(CNV_Omics)[names(CNV_Omics)==gene] <- paste0(gene,"_CNV")
  #-------mutation
  mutation_file_name <- paste0("./data/",Dataset,"/mutation/",substr(gene,1,2),".RData")
  load(mutation_file_name)
  mutation_Omics <- as.data.frame(list[c("cell_lines",gene)])
  names(mutation_Omics) <- gsub("\\.","-",names(mutation_Omics))
  row.names(mutation_Omics) <- mutation_Omics$cell_lines
  names(mutation_Omics)[names(mutation_Omics)==gene] <- paste0(gene,"_mutation")
  
  
    #methylation_promoter
    gene_information_path <- paste0("data/",Dataset,"/",Dataset,"_gene_information.RData")
    load(gene_information_path)
    gene_information <- get(ls()[grep(paste0(Dataset,"_gene_information"),ls())])
    if(gene %in% gene_information$methylation_promoter$Gene){
      methylation_promoter_file_name <- paste0("./data/",Dataset,"/methylation_promoter/",substr(gene,1,2),".RData")
      load(methylation_promoter_file_name)
      methylation_promoter_Omics <- as.data.frame(list[c("cell_lines",gene)])
      names(methylation_promoter_Omics) <- gsub("\\.","-",names(methylation_promoter_Omics))
      row.names(methylation_promoter_Omics) <- methylation_promoter_Omics$cell_lines
      names(methylation_promoter_Omics)[names(methylation_promoter_Omics)==gene] <- paste0(gene,"_methylation_promoter")
    }
    #CRISPR
    if(gene %in% gene_information$CRISPR$Gene){
      CRISPR_file_name <- paste0("./data/",Dataset,"/CRISPR/",substr(gene,1,2),".RData")
      load(CRISPR_file_name)
      CRISPR_Omics <- as.data.frame(list[c("cell_lines",gene)])
      names(CRISPR_Omics) <- gsub("\\.","-",names(CRISPR_Omics))
      row.names(CRISPR_Omics) <- CRISPR_Omics$cell_lines
      names(CRISPR_Omics)[names(CRISPR_Omics)==gene] <- paste0(gene,"CRISPR")
    }
    #RNAi
    if(gene %in% gene_information$RNAi_combine$Gene){
      RNAi_file_name <- paste0("./data/",Dataset,"/CRISPR/",substr(gene,1,2),".RData")
      load(RNAi_file_name)
      RNAi_Omics <- as.data.frame(list[c("cell_lines",gene)])
      names(RNAi_Omics) <- gsub("\\.","-",names(RNAi_Omics))
      row.names(RNAi_Omics) <- RNAi_Omics$cell_lines
      names(RNAi_Omics)[names(RNAi_Omics)==gene] <- paste0(gene,"RNAi")
    }
  #----添加蛋白组学
    if(gene %in% gene_information$protein$Gene){
      protein_file_name <- paste0("./data/",Dataset,"/protein/",substr(gene,1,2),".RData")
      load(protein_file_name)
      protein_Omics <- as.data.frame(list[c("cell_lines",gene)])
      names(protein_Omics) <- gsub("\\.","-",names(protein_Omics))
      row.names(protein_Omics) <- protein_Omics$cell_lines
      names(protein_Omics)[names(protein_Omics)==gene] <- paste0(gene,"_protein")
    }

    load("data/DepMap/DepMap_demo_data.RData")
    multiomics <- left_join(left_join(left_join(DepMap_demo_data,mutation_Omics),CNV_Omics),mRNA_Omics)
    if(length(grep("methylation_promoter_Omics",ls()))==1){
      multiomics <- left_join(multiomics,methylation_promoter_Omics)
    }
    if(length(grep("CRISPR_Omics",ls()))==1){
      multiomics <- left_join(multiomics,CRISPR_Omics)
    }
    if(length(grep("RNAi_Omics",ls()))==1){
      multiomics <- left_join(multiomics,RNAi_Omics)
    }
    if(length(grep("protein_Omics",ls()))==1){
      multiomics <- left_join(multiomics,protein_Omics)
    }
  row.names(multiomics) <- multiomics$cell_lines
  multiomics <- multiomics[cell_lines,]
  
  
  multiomics <- multiomics[,-grep("cell_line",names(multiomics))]
  multiomics <- multiomics[,-grep("demo",names(multiomics))]
  #------------计算相关系数
  pro <- matrix(NA,dim(multiomics)[2]-1,2)
  row.names(pro) <- names(multiomics)[-grep("mRNA",names(multiomics))]
  colnames(pro) <- c("cor","P")
  pro <- as.data.frame(pro)
  corrlation_variable <- row.names(pro)[-grep("mutation",row.names(pro))]
  for (i in corrlation_variable) {
    tmp <- cor.test(as.numeric(multiomics[,i]),
                    as.numeric(multiomics[,grep("mRNA",names(multiomics))]),
                    method = cor_method)
    pro[i,1] <-tmp$estimate
    pro[i,2]<- tmp$p.value
  }
  Mutation_cell <- row.names(multiomics)[multiomics[,grep("mutation",names(multiomics))]=="Mutation"]
  None_Mutation_cell <- row.names(multiomics)[multiomics[,grep("mutation",names(multiomics))]=="None"]
  if( length(Mutation_cell) == 0 ){
    res$warning <- "There is no mutation of this gene in your selected cell lines"
    pro[grep("mutation",row.names(pro)),2] <- NA
  }else{
    test_method <- ifelse(mutation_method =="wilcoxon","wilcox.test","t.test")
    tmp <- get(test_method)(multiomics[Mutation_cell,grep("mRNA",names(multiomics))],
                            multiomics[None_Mutation_cell,grep("mRNA",names(multiomics))])
    pro[grep("mutation",row.names(pro)),2] <- tmp$p.value
  }
  
  pro$cor1 <- sprintf("%0.3f", pro$cor)
  pro_new <- pro[order(pro$P),]
  pro_new$FDR <- p.adjust(pro_new$P)
  pro_new$P_new  <- sprintf("%0.3f",  pro_new$P)
  pro_new$FDR_new  <- sprintf("%0.3f",  pro_new$FDR)
  pro_new$col <- ifelse(abs(pro_new$cor) < cut_off,"#928a97","black")
  if (statistic_method =="P") {
    pro_new$biaoji <- ifelse(pro_new$P < 0.0001,"****",ifelse(pro_new$P < 0.001,"***",
                                                              ifelse(pro_new$P < 0.01,"**",ifelse(pro_new$P < 0.05,"*",""))))
    pro_new[grep("mutation",row.names(pro_new)),"col"] <-  ifelse(pro_new[grep("mutation",row.names(pro_new)),"P"] < 0.05,
                                                                  "black","#928a97")
   
  }else{
    pro_new$biaoji <- ifelse(pro_new$FDR < 0.0001,"****",ifelse(pro_new$FDR < 0.001,"***",
                                                                ifelse(pro_new$FDR < 0.01,"**",ifelse(pro_new$FDR < 0.05,"*",""))))
    pro_new[grep("mutation",row.names(pro_new)),"col"] <-  ifelse(pro_new[grep("mutation",row.names(pro_new)),"FDR"] < 0.05,
                                                                  "black","#928a97")
   
  }
  
  
  
  mutaion_location <- grep("mutation",row.names(pro_new))
  
  
  pro_new$qianzui <- NA
  pro_new$wenben <- NA
  
      if(is.na(pro_new[grep("mutation",row.names(pro_new)),"P"])){
        pro_new[grep("mutation",row.names(pro_new)),"wenben"] <- ""
      }else{
        pro_new[-mutaion_location,"qianzui"] <- ifelse(cor_method=="spearman","Rho = ","r = ")
        if (statistic_method =="P") {
          pro_new[mutaion_location,"wenben"] <- 
            paste0("P = ",pro_new[mutaion_location,"P_new"],
                   pro_new[mutaion_location,"biaoji"])
        }else{pro_new[mutaion_location,"wenben"] <- 
          paste0("FDR = ",pro_new[mutaion_location,"FDR_new"],
                 pro_new[mutaion_location,"biaoji"])}
        pro_new[-mutaion_location,"wenben"] <- 
          paste0(pro_new[-mutaion_location,"qianzui"],pro_new[-mutaion_location,"cor1"],
                 pro_new[-mutaion_location,"biaoji"])
        
        pro_new$wenben <- gsub("= 0.000","< 0.001",pro_new$wenben)
      }

  #----------mRNA 赋值
  multiomics <- multiomics[order(multiomics[,paste0(gene,"_mRNA")],na.last = F),]
  multiomics[!is.na(multiomics[,grep("mRNA",names(multiomics))]),grep("mRNA",names(multiomics))] <- 
    seq(0.3,0.7,length.out =length(multiomics[!is.na(multiomics[,grep("mRNA",names(multiomics))]),
                                              grep("mRNA",names(multiomics))] ))
  multiomics[is.na(multiomics[,grep("mRNA",names(multiomics))]),grep("mRNA",names(multiomics))] <- 0
  multiomics[,-c(grep("mutation",names(multiomics)),grep("MSI",names(multiomics)))] <- 
                               round(multiomics[,-c(grep("mutation",names(multiomics)),
                                                    grep("MSI",names(multiomics)))],3)
  pro_new <- pro_new[names(multiomics)[-grep("mRNA",names(multiomics))],]
  #---存储计算结果
  res<-list()
  buckHeight <- 25
  buckInterv <- buckHeight/10
  if(length(cell_lines)>600){
	barWidth = 2.25*600/length(cell_lines)
	barInterval = 0.75*600/length(cell_lines)
  }else{
	barWidth = 2.25
	barInterval = 0.75
  }
  rightMargin=180
  bottomMargin=10
  topMargin=25
  baseline=0+bottomMargin
  dpi <- 300
  leftMargin<-max(sapply(names(multiomics),function(x)getStrLen(x,0.5)))*dpi+5
  #当画幅宽度小于570，则不再继续减小画幅宽度
  minWidth<-round(getStrLen("* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001",0.5)*dpi)
  
  minWidth<-700
	currWidth<-length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin
  if(currWidth<minWidth){
    leftMargin<-leftMargin+(minWidth-currWidth)/2
    rightMargin<-rightMargin+(minWidth-currWidth)/2
  }
  #定义公共画图参数======================================================================================
  res$commomParams$barWidth<-barWidth
  res$commomParams$barInterval<-barInterval
  res$commomParams$leftMargin<-leftMargin
  res$commomParams$rightMargin<-rightMargin
  res$commomParams$width<-(barWidth+barInterval)*length(cell_lines)
  #定义legend参数======================================================================================
  if(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin==minWidth){
    legendLeft=-leftMargin+30
  }else if(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin<870){
    legendLeft=(length(cell_lines)*(barWidth+barInterval)+rightMargin+leftMargin-minWidth)/2-leftMargin+30
  }else{
    legendLeft=0
  }
  
  
  
  ######
  
  res$plotElements$legend_text$plotType <-"drawLegend"
  
  res$plotElements$legend_text$left<-legendLeft
  res$plotElements$legend_text$comment$baseline<- baseline
  res$plotElements$legend_text$comment$text<-"* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001"
  
  baseline<-baseline+buckHeight+buckInterv
  
  
  res$plotElements$legend_mutation$plotType <-"drawLegend"
  res$plotElements$legend_mutation$title$text<-"Mutation type:"
  res$plotElements$legend_mutation$title$baseline<-baseline
  res$plotElements$legend_mutation$title$left<-legendLeft
  lenOfTitle<-150
  
  legendOffset<-0
  res$plotElements$legend_mutation$legend[["Mutation"]]$baseline<-baseline
  res$plotElements$legend_mutation$legend[["Mutation"]]$size<-15
  res$plotElements$legend_mutation$legend[["Mutation"]]$color<- Mutation_color
  res$plotElements$legend_mutation$legend[["Mutation"]]$left<-legendLeft+lenOfTitle
  res$plotElements$legend_mutation$legend[["Mutation"]]$text<-"Mutation"
  
  legendOffset<-legendOffset+getStrLen("Mutation",0.5)*300+15
  res$plotElements$legend_mutation$legend[["None"]]$baseline<-baseline
  res$plotElements$legend_mutation$legend[["None"]]$size<-15
  res$plotElements$legend_mutation$legend[["None"]]$color<- No_mutation_color
  res$plotElements$legend_mutation$legend[["None"]]$left<-legendLeft+lenOfTitle+legendOffset
  res$plotElements$legend_mutation$legend[["None"]]$text<-"No alteration"
  
  
  baseline<-baseline+buckHeight*1.5+buckInterv
  
  #-----------绘制正文的参数
  
  #mRNA
  res$plotElements$mRNA$plotType<-"numericalBarplot"
  res$plotElements$mRNA$pillarColor<- mRNA_pillar_color
  res$plotElements$mRNA$baseline<-baseline
  res$plotElements$mRNA$leftText <- paste0(gene," mRNA")
  res$plotElements$mRNA$data <- multiomics[,grep("mRNA",names(multiomics))]
  res$plotElements$mRNA$buckHeight<-buckHeight
  res$plotElements$mRNA$rightText <- ""
  res$plotElements$mRNA$pillarClassName<-"CellAndDrugResistTags"
  res$plotElements$mRNA$resize<-"false"
  
  baseline<-baseline+buckInterv+buckHeight
  
      res$plotElements$CNV$plotType<-"numericalBarplot"
      res$plotElements$CNV$pillarColor<- protein_pillar_color
      res$plotElements$CNV$baseline<-baseline
      res$plotElements$CNV$leftText <- paste0(gene," CNV")
      res$plotElements$CNV$data <- multiomics[,grep("CNV",names(multiomics))]
      res$plotElements$CNV$buckHeight<- buckHeight
      res$plotElements$CNV$rightText <- pro_new$wenben[grep("CNV",row.names(pro_new))]
      res$plotElements$CNV$rightTextColor <-pro_new$col[grep("CNV",row.names(pro_new))]
      res$plotElements$CNV$pillarClassName<-"CellAndDrugResistTags"
      #res$plotElements$protein$resize<-Ture
      baseline<-baseline+buckInterv+buckHeight
      
 
  
  #Mutation
  lab<-c("Mutation","None")
  color_tmp<-c(Mutation_color,No_mutation_color)
  for(i in 1:length(lab)){
    res$plotElements[[paste0("mutation_",i)]]$plotType<-"numericalBarplot"
    res$plotElements[[paste0("mutation_",i)]]$pillarColor<-color_tmp[i]
    res$plotElements[[paste0("mutation_",i)]]$baseline<-baseline
    mutationdata_tmp<-rep(0,length(multiomics[,grep("mutation",names(multiomics))]))
    mutationdata_tmp[which(multiomics[,grep("mutation",names(multiomics))]==lab[i])]<-0.8
    res$plotElements[[paste0("mutation_",i)]]$data <- mutationdata_tmp
    res$plotElements[[paste0("mutation_",i)]]$buckHeight<-buckHeight
    res$plotElements[[paste0("mutation_",i)]]$pillarClassName<-"CellAndDrugResistTags"
    res$plotElements[[paste0("mutation_",i)]]$resize<-"false"
    if(i==1){
      res$plotElements[[paste0("mutation_",i)]]$leftText <- paste0(gene," Mutation")
      res$plotElements[[paste0("mutation_",i)]]$rightText <- pro_new$wenben[grep("mutation",row.names(pro_new))]
      res$plotElements[[paste0("mutation_",i)]]$rightTextColor <- pro_new$col[grep("mutation",row.names(pro_new))]
    }
  }
  baseline<-baseline+buckInterv+buckHeight
  
  #--protien
  if(length(grep("protein_Omics",ls()))==1){
    res$plotElements$protein$plotType<-"numericalBarplot"
    res$plotElements$protein$pillarColor<- protein_pillar_color
    res$plotElements$protein$baseline<-baseline
    res$plotElements$protein$leftText <- paste0(gene," protein")
    res$plotElements$protein$data <- multiomics[,grep("protein",names(multiomics))]
    res$plotElements$protein$buckHeight<- buckHeight
    res$plotElements$protein$rightText <- pro_new$wenben[grep("protein",row.names(pro_new))]
    res$plotElements$protein$rightTextColor <-pro_new$col[grep("protein",row.names(pro_new))]
    res$plotElements$protein$pillarClassName<-"CellAndDrugResistTags"
    #res$plotElements$protein$resize<-Ture
    baseline<-baseline+buckInterv+buckHeight
  }
  
      #methylation_promoter
  if(length(grep("methylation_promoter_Omics",ls()))==1){
      res$plotElements$methylation_promoter$plotType<-"numericalBarplot"
      res$plotElements$methylation_promoter$pillarColor<- protein_pillar_color
      res$plotElements$methylation_promoter$baseline<-baseline
      res$plotElements$methylation_promoter$leftText <- paste0(gene," methylation promoter")
      res$plotElements$methylation_promoter$data <- multiomics[,grep("methylation_promoter",names(multiomics))]
      res$plotElements$methylation_promoter$buckHeight<- buckHeight
      res$plotElements$methylation_promoter$rightText <- pro_new$wenben[grep("methylation_promoter",row.names(pro_new))]
      res$plotElements$methylation_promoter$rightTextColor <-pro_new$col[grep("methylation_promoter",row.names(pro_new))]
      res$plotElements$methylation_promoter$pillarClassName<-"CellAndDrugResistTags"
      #res$plotElements$protein$resize<-Ture
      baseline<-baseline+buckInterv+buckHeight
  }    
      ##########CRISPR
  if(length(grep("CRISPR_Omics",ls()))==1){
      res$plotElements$CRISPR$plotType<-"numericalBarplot"
      res$plotElements$CRISPR$pillarColor<- protein_pillar_color
      res$plotElements$CRISPR$baseline<-baseline
      res$plotElements$CRISPR$leftText <- paste0(gene," CRISPR")
      res$plotElements$CRISPR$data <- multiomics[,grep("CRISPR",names(multiomics))]
      res$plotElements$CRISPR$buckHeight<- buckHeight
      res$plotElements$CRISPR$rightText <- pro_new$wenben[grep("CRISPR",row.names(pro_new))]
      res$plotElements$CRISPR$rightTextColor <-pro_new$col[grep("CRISPR",row.names(pro_new))]
      res$plotElements$CRISPR$pillarClassName<-"CellAndDrugResistTags"
      #res$plotElements$protein$resize<-Ture
      baseline<-baseline+buckInterv+buckHeight
  }
      ##########RNAi
  if(length(grep("RNAi_Omics",ls()))==1){
      res$plotElements$RNAi$plotType<-"numericalBarplot"
      res$plotElements$RNAi$pillarColor<- protein_pillar_color
      res$plotElements$RNAi$baseline<-baseline
      res$plotElements$RNAi$leftText <- paste0(gene," RNAi")
      res$plotElements$RNAi$data <- multiomics[,grep("RNAi",names(multiomics))]
      res$plotElements$RNAi$buckHeight<- buckHeight
      res$plotElements$RNAi$rightText <- pro_new$wenben[grep("RNAi",row.names(pro_new))]
      res$plotElements$RNAi$rightTextColor <-pro_new$col[grep("RNAi",row.names(pro_new))]
      res$plotElements$RNAi$pillarClassName<-"CellAndDrugResistTags"
      #res$plotElements$protein$resize<-Ture
      baseline<-baseline+buckInterv+buckHeight
  }
  
  #收尾获得文档高度，并返回======================================================================================
  res$commomParams$DocHeight<- baseline+topMargin
  res$cell_line <- row.names(multiomics)
  return(res)
}