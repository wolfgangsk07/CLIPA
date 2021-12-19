

function addSVGEvent(Model_ID){
	
	var Model=document.getElementById(Model_ID);
	var plotArea=Model.getElementsByClassName("plotArea")[0]
	var Tags=["drugTags","CellAndDrugResistTags","MethSiteTags","summBarplotInfo","warnings"]
	for(var t=0;t<Tags.length;t++){
		var drugTags=plotArea.getElementsByClassName(Tags[t]);
		for(var i=0;i<drugTags.length;i++){

			drugTags[i].addEventListener("mouseenter",showInfo,false)
			drugTags[i].addEventListener("mouseleave",closeInfo,false)
		}
	}
	
	
}


function addOmics(Model_ID){
	
	var Model=document.getElementById(Model_ID)
	var omicsSelector=Model.getElementsByClassName("omicsSelector")[0]
	var genename=Model.getElementsByClassName("genename")[0]
	if(genename.value==""){
		genename.style.background="red"
		return
	}else{
		genename.style.background=""
	}
	var OmicResultItems=Model.getElementsByClassName("OmicResultItem")
	for(var o=0;o<OmicResultItems.length;o++){
		if(omicsSelector.value+","+genename.value==OmicResultItems[o].getAttribute("Omic")+","+OmicResultItems[o].getAttribute("gene")){
			var addButton=Model.getElementsByClassName("addButton")[0]
			addButton.setAttribute("errorInfo",'<p style="font-size:12px"><span style="color:red">'+omicsSelector.value+":"+genename.value+"</span> has already been added!</p>")
			showInfo_mani(addButton)
			return;
		}
	}
	var OmicResultBox=Model.getElementsByClassName("OmicResultBox")[0]
	var item=document.createElement("div")
		item.setAttribute("class","OmicResultItem")
		item.setAttribute("Omic",omicsSelector.value)
		item.setAttribute("gene",genename.value)
		item.setAttribute("rank",OmicResultItems.length+1)
		item.setAttribute("onclick","moveToFirst(this)")
		
		var p=document.createElement("p")
			var selectedIndex=omicsSelector.selectedIndex
			var options=omicsSelector.getElementsByTagName("option")
			p.innerHTML="("+options[selectedIndex].innerHTML+'): <span style="font-weight:bold">'+genename.value+"</span>"
		var clear_div=document.createElement("div")
			var clear=SVG(clear_div)
				clear.viewbox(0,0,20,20)
				var group = clear.group().attr({"onclick":"removeOmics(this)",class:"justPointer"})
				var circle=clear.circle(18).attr({
					cx:10,
					cy:10,
					fill:"transparent"
				})
				var line1=clear.line(5, 5, 15, 15).attr({
					"stroke-width": 1.5,
					'stroke-linecap':"round",
					'stroke-linejoin':"round",
					'stroke':"white"
				})
				var line2=clear.line(15, 5, 5, 15).attr({
					"stroke-width": 1.5,
					'stroke-linecap':"round",
					'stroke-linejoin':"round",
					'stroke':"white"
				})
				group.add(circle)
				group.add(line1)
				group.add(line2)
		item.appendChild(p)
		item.appendChild(clear_div)
	OmicResultBox.appendChild(item)
	reArrangeOmics(OmicResultBox)
	
	//checkStartButton(Model_ID)
}


function applyFilter(Model_ID){
	
	var dataBaseName=Model_ID.split(',')[0]
	var ModelName=Model_ID.split(',')[1]
	var Model=document.getElementById(Model_ID);
	
	var filterSwitch=Model.getElementsByClassName("filterSwitch")[0];
	var CellLinesFilterContent=Model.getElementsByClassName("CellLinesFilterContent")[0];
	var OmicsSelector=Model.getElementsByClassName("omicsSelector")[0];
	var omicsSelectorAndMultigenes=Model.getElementsByClassName("omicsSelectorAndMultigenes")[0]
	

	dealingCellLines[Model_ID]=JSON.parse(initParams_str).cellLines[dataBaseName];
	var filteredCellLine=0;
	if(!filterSwitch.checked){
		CellLinesFilterContent.setAttribute("style","height:0px");
	}
	
	var MSI_cellLines=CellLinesFilterContent.getElementsByClassName("MSI_cellLines")[0]
		var MSI_cellLines_Status=MSI_cellLines.getElementsByClassName("subFilterCheckbox")[0].checked;
		if(MSI_cellLines.getElementsByTagName("input")[1].checked){
		MSI_status="MSS"
		}else{
		MSI_status="MSI"
		}
	var ploidy_cellLines=CellLinesFilterContent.getElementsByClassName("ploidy_cellLines")[0]
		var ploidy_cellLines_Status=ploidy_cellLines.getElementsByClassName("subFilterCheckbox")[0].checked;
		var ploidy_cellLines_type=ploidy_cellLines.getElementsByTagName("select")[0].value
		var ploidy_cellLines_value=ploidy_cellLines.getElementsByClassName("subFilterPercInput")[0].value
	
	var mutation_cellLines=CellLinesFilterContent.getElementsByClassName("mutation_cellLines")[0]
		var smutation_cellLines_Status=mutation_cellLines.getElementsByClassName("subFilterCheckbox")[0].checked;
		var mutation_cellLines_type=mutation_cellLines.getElementsByTagName("select")[0].value
		var mutation_cellLines_value=mutation_cellLines.getElementsByClassName("subFilterPercInput")[0].value
		
	if(["Cancer Types Summary","Omics-Omics (cis-regulation)","Cancer Subtypes","Omics-Omics (trans-regulation)","Fusions"].indexOf(ModelName)>-1){
		var shouldHasDrug=false
	}else{
		var shouldHasDrug=true
	}
		
	
	var cancerTypes_Names=Object.keys(dealingCellLines[Model_ID]);
	for(var i=0;i<cancerTypes_Names.length;i++){
		var cellLines_Names=Object.keys(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines);
		for(var j=0;j<cellLines_Names.length;j++){
			var cellLine_enabled=true; //初始化结果
			
			//判断条件
			//console.log(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]][OmicsSelector.value.split("_")[0]]);
			if(OmicsSelector!=undefined&&["Omics-Omics (trans-regulation)"].indexOf(ModelName)==-1){
				if(!dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]][OmicsSelector.value.split("_")[0]]){
					cellLine_enabled=false
					filteredCellLine++
				}
			}else if(ModelName=="Drug-Pathway"){
				if(!dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]]["mRNA"]){
					cellLine_enabled=false
					filteredCellLine++
				}
			}else if(ModelName=="Fusions"){
				if(!dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]]["Fusions"]){
					cellLine_enabled=false
					filteredCellLine++
				}
			}
			
			if(omicsSelectorAndMultigenes!=undefined){
				if(!dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]][omicsSelectorAndMultigenes.value.split("_")[0]]){
					cellLine_enabled=false
					filteredCellLine++
				}
			}
			
			
			if(shouldHasDrug&!dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]].hasDrug){
				
				cellLine_enabled=false
				filteredCellLine++
			}
			
			
			
			
			
			
			if(filterSwitch.checked){
				if(MSI_cellLines_Status){
					MSI_cellLines.getElementsByTagName("input")[1].disabled=false;
					MSI_cellLines.getElementsByTagName("input")[2].disabled=false;
					if(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]].MSI_Status!=MSI_status){
					cellLine_enabled=false
					filteredCellLine++
					}
				}else{
					MSI_cellLines.getElementsByTagName("input")[1].disabled=true;
					MSI_cellLines.getElementsByTagName("input")[2].disabled=true;
				}
				
				
				if(ploidy_cellLines_Status){
					ploidy_cellLines.getElementsByTagName("input")[1].disabled=false;
					ploidy_cellLines.getElementsByTagName("select")[0].disabled=false;
					if(isPercentage(ploidy_cellLines_value)){
					if(ploidy_cellLines_type=="Above"){
						if(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]].ploidyPerc<ploidy_cellLines_value.replace("%","")/100){
						cellLine_enabled=false
						filteredCellLine++
						}
					}else if(ploidy_cellLines_type=="Below"){
						if(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]].ploidyPerc>ploidy_cellLines_value.replace("%","")/100){
						cellLine_enabled=false
						filteredCellLine++
						}
					}
					ploidy_cellLines.getElementsByClassName("subFilterPercInput")[0].style.background="";
					}else{
					ploidy_cellLines.getElementsByClassName("subFilterPercInput")[0].style.background="red";
					ploidy_cellLines.getElementsByClassName("subFilterCheckbox")[0].checked=false;

					}
				}else{
					ploidy_cellLines.getElementsByTagName("input")[1].disabled=true;
					ploidy_cellLines.getElementsByTagName("select")[0].disabled=true;
				}
				
				
				if(smutation_cellLines_Status){
					mutation_cellLines.getElementsByTagName("input")[1].disabled=false;
					mutation_cellLines.getElementsByTagName("select")[0].disabled=false;
					if(isPercentage(mutation_cellLines_value)){
					if(mutation_cellLines_type=="Above"){
						if(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]].mutational_burdenPerc<mutation_cellLines_value.replace("%","")/100){
						cellLine_enabled=false
						filteredCellLine++
						}
					}else if(mutation_cellLines_type=="Below"){
						if(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]].mutational_burdenPerc>mutation_cellLines_value.replace("%","")/100){
						cellLine_enabled=false
						filteredCellLine++
						}
					}
					mutation_cellLines.getElementsByClassName("subFilterPercInput")[0].style.background="";
					}else{
					mutation_cellLines.getElementsByClassName("subFilterPercInput")[0].style.background="red";
					mutation_cellLines.getElementsByClassName("subFilterCheckbox")[0].checked=false;

					}
				}else{
					mutation_cellLines.getElementsByTagName("input")[1].disabled=true;
					mutation_cellLines.getElementsByTagName("select")[0].disabled=true;
				}
				CellLinesFilterContent.removeAttribute("style");
			}
			
			
			
			
			
			
			
			dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLines_Names[j]].enabled=cellLine_enabled
				
		}
	}
	
	
}


function asArray(x){
	
	if(typeof x=="string"|typeof x=="number"){
		return([x])
	}else{
		return(x)
	}
	
}
function clearNodelist(nodeList){
	
	var postChild = [];
    for (var i = 0; i < nodeList.length; i++) {
        if(!(nodeList[i].nodeType == '3' && nodeList[i].nodeName == '#text' && !/\S/.test(nodeList[i].nodeValue))){ //文本节点并且是空的文本节点时，将空文本节点删除
            postChild.push(nodeList[i]);
        }
    }
	
    return(postChild);

}



function checkGenes(obj){
	
	//for gene list only
	if(vennList==undefined)return
	var Model=document.getElementById(currentModel);
	var omicsSelectorAndMultigenes=Model.getElementsByClassName("omicsSelectorAndMultigenes")[0];
	if(omicsSelectorAndMultigenes==undefined){
		var vennTag="all"
	}else{
		var vennTag=currentModel.split(",")[0]+"_"+omicsSelectorAndMultigenes.value
	}
	getGeneList(vennTag)
	if(obj.value=="")return
	var genes=obj.value.toUpperCase()

	//genes="CEBPB,gata2,gata3\na1bg,,,,sdfsdg"
	//obj.value=genes;
	genes=genes.replace(/ /g,"")
	genes=genes.replace(/, /g,",")
	genes=genes.replace(/\n/g,",")
	genes=genes.replace(/\t/g,",")
	while(genes.indexOf(",,")>-1){
		genes=genes.replace(/,,/g,",")
	}

	genes=genes.split(",")
	
	var ourGenes=[]
	var invalid=[]
	for(var g=0;g<genes.length;g++){
		var gene=genes[g]
		if(geneList[vennTag].indexOf(gene)>-1){
			ourGenes.push(gene)
			continue
		}else{
			invalid.push(gene)
		}
		//console.log(gene)
	}
	var ourGenes_col="green"
	var invalid_col="green"
	if(ourGenes.length==0){
		ourGenes_col="red"
	}
	if(invalid.length!=0){
		invalid_col="red"
	}
	var text='<p style="color:'+ourGenes_col+';font-weight:bold">'+ourGenes.length +' genes passed the inspection!</p><p style="color:'+invalid_col+';font-weight:bold"><br>'+invalid.length+' gene(s) are invalid:<br>'+invalid.join(", ")+'</p>'
	var color="red"
	if(ourGenes.length==1){
		text='<p style="color:red;font-weight:bold">Less than 2 genes passed the inspection!</p>'
		ourGenes=[]
	}


	obj.setAttribute("statusText",text)
	var index=0;
	obj.value=ourGenes.join(",")
	showInfo_mani(obj);
	checkStartButton(currentModel)
	

	

}


function checkOneDrug(obj){
	

	var Model=document.getElementById(currentModel)
	var dataBaseName=currentModel.split(',')[0]
	var ModelName=currentModel.split(',')[1]
	//folderList_type="foldList_drugs";
	var drug_selector_brothers=obj.parentNode.parentNode.children;
	//console.log(drug_selector_brothers[0].getElementsByClassName("drug_selector")[0].checked);
	//var selectorSubALLs=document.getElementsByClassName("selectorSubALL");
	//userParams[currentModel][folderList_type][itemName]=object.checked;

	//console.log(JSON.stringify(userParams[currentModel]));
	
	//var keys= Object.keys(initParams[folderList_type]);
	brotherStatus=0;
	for(var i=0;i< drug_selector_brothers.length;i++){
		if(drug_selector_brothers[i].getElementsByClassName("drug_selector")[0].checked){
			brotherStatus++;
			//console.log(brotherStatus);
		}
	}
	
	var parentNode=drug_selector_brothers[0].parentNode.parentNode;
	
	if(brotherStatus==0){
		parentNode.getElementsByClassName("selectorSubALL")[0].checked=false;
		parentNode.getElementsByClassName("selectorSubALL")[0].indeterminate =false;
	}else if(brotherStatus==drug_selector_brothers.length){
		parentNode.getElementsByClassName("selectorSubALL")[0].checked=true;
		parentNode.getElementsByClassName("selectorSubALL")[0].indeterminate =false;
	}else{
		parentNode.getElementsByClassName("selectorSubALL")[0].indeterminate =true;
	}
	var item_tag=parentNode.getElementsByClassName("countTag")[0]
	if(!initParams.ModelInfo[dataBaseName][ModelName].includes("foldList_drugs_single")){
		item_tag.innerHTML="("+brotherStatus+"/"+drug_selector_brothers.length+")"
	}
	if(brotherStatus!=0){
		item_tag.setAttribute("style","color:white;")
	}else{
		item_tag.setAttribute("style","color:grey;")
	}
	
	var count_drugs=countFoldList_drugs(currentModel);
	if(count_drugs[0].length==0){

		Model.getElementsByClassName("drug_selectorALL")[0].checked=false;
		Model.getElementsByClassName("drug_selectorALL")[0].indeterminate =false;
	}else if(count_drugs[0].length==count_drugs[1]){

		Model.getElementsByClassName("drug_selectorALL")[0].checked=true;
		Model.getElementsByClassName("drug_selectorALL")[0].indeterminate =false;
	}else{

		Model.getElementsByClassName("drug_selectorALL")[0].indeterminate =true;
	}

	checkStartButton(currentModel)
	
	
}

function countFoldList_drugs(Model_ID){
	
	var Model=document.getElementById(Model_ID)
	var dataBaseName=Model_ID.split(',')[0]
	var ModelName=Model_ID.split(',')[1]
	var folderList_type="foldList_drugs";
	var foldList_drugs=Model.getElementsByClassName(folderList_type)[0];
	var drugs=foldList_drugs.getElementsByClassName("drug_selector");

	
	
	var selectedDrugs=[];
	for(var i=0;i< drugs.length;i++){
		if(drugs[i].checked){
			selectedDrugs.push(drugs[i].getAttribute("drugname"));
		}
	}
	
	var checkedCount=selectedDrugs.length;
	if(initParams.ModelInfo[dataBaseName][ModelName].includes("foldList_drugs_single")){
		if(selectedDrugs==0)selectedDrugs=["None"]
		Model.getElementsByClassName("totalCount_Drug")[0].innerHTML="("+selectedDrugs.join("")+")";
	}else{
		Model.getElementsByClassName("totalCount_Drug")[0].innerHTML="("+checkedCount+"/"+drugs.length+")";
	}
	
	return([selectedDrugs,drugs.length])
}

function countFoldList_cellLines(Model_ID){
	var selectedCellLines=[];
	var Model=document.getElementById(Model_ID)
	//var dealingCellLines[Model_ID]=applyCellLineFilter(Model_ID)
	var cellLines=Model.getElementsByClassName("cellLine_selector");
	for(var i=0;i<cellLines.length;i++){
		
		if(cellLines[i].checked){
			selectedCellLines.push(cellLines[i].getAttribute("cellname"));
		}
	}
		
		
		
	var checkedCount=selectedCellLines.length;
	var totalCount=0;
	var cancerTypes_Names= Object.keys(dealingCellLines[Model_ID]);
	
	for(var i=0;i< cancerTypes_Names.length;i++){	
		var cellLineNames=Object.keys(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines)
		for(var j=0;j<cellLineNames.length;j++){
			if(dealingCellLines[Model_ID][cancerTypes_Names[i]].cellLines[cellLineNames[j]].enabled==false){
				
			}else{
				totalCount++;
			}
		}
	}
	
	
	Model.getElementsByClassName("totalCount_CellLines")[0].innerHTML="("+checkedCount+"/"+totalCount+")";
	return([selectedCellLines,totalCount])
	
}



function closeInfo(event){
	isOutElement=true;
	setTimeout(function(){
		if(isOutBox&isOutElement){
			popupInfoBox.style.opacity=0
			popupInfoBox.style.pointerEvents ="none";

		}
	},500)
	
}

function checkOneCancer(obj){
	var Model=document.getElementById(currentModel)
	var dataBaseName=currentModel.split(',')[0]
	var ModelName=currentModel.split(',')[1]
	//folderList_type="foldList_drugs";
	var cellline_selector_brothers=obj.parentNode.parentNode.children;
	//console.log(cellline_selector_brothers[0].getElementsByClassName("drug_selector")[0].checked);
	//var selectorSubALLs=document.getElementsByClassName("selectorSubALL");
	//userParams[currentModel][folderList_type][itemName]=object.checked;

	//console.log(JSON.stringify(userParams[currentModel]));
	
	//var keys= Object.keys(initParams[folderList_type]);
	brotherStatus=0;
	for(var i=0;i< cellline_selector_brothers.length;i++){
		if(cellline_selector_brothers[i].getElementsByClassName("cellLine_selector")[0].checked){
			brotherStatus++;
			//console.log(brotherStatus);
		}
	}
	
	var parentNode=cellline_selector_brothers[0].parentNode.parentNode;
	
	if(brotherStatus==0){
		parentNode.getElementsByClassName("selectorSubALL")[0].checked=false;
		parentNode.getElementsByClassName("selectorSubALL")[0].indeterminate =false;
	}else if(brotherStatus==cellline_selector_brothers.length){
		parentNode.getElementsByClassName("selectorSubALL")[0].checked=true;
		parentNode.getElementsByClassName("selectorSubALL")[0].indeterminate =false;
	}else{
		parentNode.getElementsByClassName("selectorSubALL")[0].indeterminate =true;
	}
	var item_tag=parentNode.getElementsByClassName("countTag")[0]
	if(!initParams.ModelInfo[dataBaseName][ModelName].includes("foldList_drugs_single")){
		item_tag.innerHTML="("+brotherStatus+"/"+cellline_selector_brothers.length+")"
	}
	if(brotherStatus!=0){
		item_tag.setAttribute("style","color:white;")
	}else{
		item_tag.setAttribute("style","color:grey;")
	}
	
	var count_drugs=countFoldList_cellLines(currentModel);
	if(count_drugs[0].length==0){

		Model.getElementsByClassName("cellLine_selectorALL")[0].checked=false;
		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =false;
	}else if(count_drugs[0].length==count_drugs[1]){

		Model.getElementsByClassName("cellLine_selectorALL")[0].checked=true;
		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =false;
	}else{

		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =true;
	}

	checkStartButton(currentModel)
	
	/*
	//countFoldList_cellLines(currentModel);
	var Model=document.getElementById(currentModel)
	var count_cellLines=countFoldList_cellLines(currentModel);
	if(count_cellLines[0].length==0){
		Model.getElementsByClassName("cellLine_selectorALL")[0].checked=false;
		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =false;
	}else if(count_cellLines[0].length==count_cellLines[1]){

		Model.getElementsByClassName("cellLine_selectorALL")[0].checked=true;
		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =false;
	}else{

		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =true;
	}
	checkStartButton(currentModel);
	*/
	
}


function changeFoldListCheckbox(obj,Model_ID){
	var Model=document.getElementById(Model_ID);
	var foldList_folderInputs=Model.getElementsByClassName("foldList_folderInput")
	var balloon=Model.getElementsByClassName("balloon")[0]
	if(obj.checked){
		var count=0
		var myVar = setInterval(function(){
			if(balloon.clientHeight==0){
				for(var f=0;f<foldList_folderInputs.length;f++){
					if(foldList_folderInputs[f]!=obj&foldList_folderInputs[f].checked){
						foldList_folderInputs[f].checked=false;
						
						clearInterval(myVar)
						break
					}
				}
			}
			count++
			if(count>1000/10)clearInterval(myVar)
		},10)
		

	}
}

function changeInput(id,promotInfo){
	var inputBox=document.getElementById(id);
	inputBox.value=promotInfo;
	var promotBox=document.getElementById("promotBox");
	promotBox.style.display="none";
	var Model=document.getElementById(currentModel);
	var inputAndSearch_Search=Model.getElementsByClassName("inputAndSearch_Search")[0];
	if(inputAndSearch_Search!=undefined){
		inputAndSearch_Search.style.display="none";
	}
	checkStartButton(currentModel);
}
function closePromot(obj){
	isOnInput=false
	var Model=document.getElementById(currentModel);
	var omicsSelector=Model.getElementsByClassName("omicsSelector")[0];

	var promotBox=document.getElementById("promotBox");
	if(isOutPromotBox){
		promotBox.style.display="none";
		if(obj.getAttribute("class").split(" ").indexOf("genename")>-1&(omicsSelector==undefined||omicsSelector.value!="methylation_site")){
			var gene=obj.value;
			if(omicsSelector==undefined){
				var vennTag="all"
			}else{
				var vennTag=currentModel.split(",")[0]+"_"+omicsSelector.value
			}
			if(geneList[vennTag]!=undefined){
				if(geneList[vennTag].indexOf(gene)>-1){
					obj.style.background="white"
					
				}else{
					obj.style.background="red"
				}
			}else{
				obj.style.background="white"
			}
		}
		checkStartButton(currentModel);

	}
	

}



function closeSearchBox(obj){
	inputAndSearch_input_isOn=false
	setTimeout(function(){
		if(!(inputAndSearch_input_isOn|inputAndSearch_Search_isOn)){
			var inputAndSearch_Search=obj.parentNode.parentNode.getElementsByClassName("inputAndSearch_Search")[0];
			inputAndSearch_Search.style.display="none";
			checkMSigDB(currentModel)
			checkStartButton(currentModel);
		}
	},200)
}



function checkMSigDB(Model_ID){
	var Model=document.getElementById(Model_ID);
	var obj=Model.getElementsByClassName("inputAndSearch_Input_input")[0]
	var correct=false;
	var genes=obj.value.toUpperCase()
	var parents=Object.keys(initParams.foldList_geneSets);
	for(var i=0;i<parents.length;i++){
		var geneSet_Names=Object.keys(initParams.foldList_geneSets[parents[i]])
		for(var j=0;j<geneSet_Names.length;j++){
			if(parents[i].toUpperCase()+"_"+geneSet_Names[j]==obj.value){
				correct=true;
				break;
			}
		}
	}
	if(correct|genes==""){
		obj.style.background="";
	}else{
		obj.style.background="red";
	}
	
}
function checkStartButton(Model_ID){
	
	var Model=document.getElementById(Model_ID);
	
	var startButtonStatus=true;
	var dataBaseName=Model_ID.split(',')[0]
	var ModelName=Model_ID.split(',')[1]
	if(dataBaseName=="NATIVE")return;
	var Children_names=initParams.ModelInfo[dataBaseName][ModelName];
	disableReasons=[];
	if(!Children_names.includes("plotButton")){
		return;
	}
	if(Children_names.includes("omicsSelector")){
		if(Model.getElementsByClassName("genename")[0].style["background-color"]=="red"|Model.getElementsByClassName("genename")[0].value==""){
			startButtonStatus=false;
			disableReasons.push("No gene is assigned!")
		}
	}
	if(Children_names.includes("drugMethod")){
		if(Model.getElementsByClassName("pcutoff")[0].style["background-color"]=="red"){
			startButtonStatus=false;
			disableReasons.push("Cut-off of P value is NOT assigned!")
		}
	}
	
	if(Children_names.includes("onlyGene")){
		if(Model.getElementsByClassName("genename")[0].style["background-color"]=="red"|Model.getElementsByClassName("genename")[0].value==""){
			startButtonStatus=false;
			disableReasons.push("No gene is assigned!")
		}
	}
	if(Children_names.includes("geneSetSelector")){
		if(geneSets_Textarea){
			if(Model.getElementsByClassName("geneSetsTextarea")[0].value==""){
				startButtonStatus=false;
				disableReasons.push("No gene set are assigned!")
			}
		}else if(geneSets_Preset){
			if(Model.getElementsByClassName("inputAndSearch_Input_input")[0].value==""|Model.getElementsByClassName("inputAndSearch_Input_input")[0].style["background-color"]=="red"){
				startButtonStatus=false;
				disableReasons.push("No gene set are assigned!")
			}
		}
	}
	if(Children_names.includes("omicsSelectorAndMultigenes")){
		var geneSetsTextarea=Model.getElementsByClassName("geneSetsTextarea")[0].value;
		if(geneSetsTextarea==""){
			startButtonStatus=false;
			disableReasons.push("No gene set are assigned!")
		}
	}
	if(Children_names.includes("foldList_drugs")){
		if(countFoldList_drugs(Model_ID)[0]==0){
			startButtonStatus=false;
			disableReasons.push("No drug(s) are selected!")
		}
	}
	if(Children_names.includes("foldList_drugs_single")){
		if(countFoldList_drugs(Model_ID)[0][0]=="None"){
			startButtonStatus=false;
			disableReasons.push("No drug are selected!")
		}
	}
	if(Children_names.includes("foldList_cellLines")){
		if(countFoldList_cellLines(Model_ID)[0].length<3){
			startButtonStatus=false;
			disableReasons.push("Less than 3 cell line(s) are selected!")
		}
	}
	
	if(Children_names.includes("multiOmicsSelector")){
		var OmicResultItem=Model.getElementsByClassName("OmicResultItem")
		if(OmicResultItem.length<2){
			startButtonStatus=false;
			disableReasons.push("Less than 2 target(s) are selected!")
		}
	}
	
	
	
			
	
	

		
	
	if(!startButtonStatus){
		Model.getElementsByClassName("startplot")[0].style.background="grey";
		Model.getElementsByClassName("startplot")[0].style.cursor="default";
		Model.getElementsByClassName("startplot")[0].removeAttribute("onclick")
		Model.getElementsByClassName("startplot")[0].addEventListener("mouseenter",showInfo,false)
		Model.getElementsByClassName("startplot")[0].addEventListener("mouseleave",closeInfo,false)
	}else{
		Model.getElementsByClassName("startplot")[0].style.background="";
		Model.getElementsByClassName("startplot")[0].style.cursor="";
		Model.getElementsByClassName("startplot")[0].setAttribute("onclick","startplot()")
		Model.getElementsByClassName("startplot")[0].removeEventListener("mouseenter",showInfo,false)
		Model.getElementsByClassName("startplot")[0].removeEventListener("mouseleave",closeInfo,false)
	}

}

function clearInput(obj){
	var input=obj.parentNode.getElementsByTagName("input")
	input.value=""
}

function clusterOmicChanged(Model_ID){
	var Model=document.getElementById(Model_ID);
	var omicsSelectorAndMultigenes=Model.getElementsByClassName("omicsSelectorAndMultigenes")[0];
	var geneSetsTextarea=Model.getElementsByClassName("geneSetsTextarea")[0];
	var genenameTag=Model.getElementsByClassName("genenameTag")[0];
	if(omicsSelectorAndMultigenes.value=="methylation_site"){
		geneSetsTextarea.placeholder='Paste your methylation site list here. They should be seperated by "," or "line break" like below:\ncg00009944,cg00009943';
		genenameTag.innerHTML="Please input methylation sites:"
		geneSetsTextarea.setAttribute("onblur","")
		geneSetsTextarea.value=""

	}else{
		genenameTag.innerHTML="Please input gene symbols:"
		geneSetsTextarea.setAttribute("onblur","checkGenes(this)")
		geneSetsTextarea.placeholder='Paste your gene list here. Genes should be seperated by "," or "line break" like below:\nTP53, GATA3\nGATA2,A1CF,A1BG';
		geneSetsTextarea.value=""
	}
}
function disableIndividual(Model_ID){
	var Model=document.getElementById(Model_ID);
	var individualCheckBox=Model.getElementsByClassName("individualCheckBox")[0]
	var checkbox=individualCheckBox.getElementsByTagName("input")[0]
	var select=individualCheckBox.getElementsByTagName("select")[0]
	if(checkbox.checked){
		select.removeAttribute("disabled")
	}else{
		select.setAttribute("disabled",true)
	}
		
}
function expandModelTips(obj){
	obj.setAttribute("expand",true);

	for(var i=1;i<obj.children.length;i++){
		obj.children[i].style.top=10+i*45+"px";
		obj.children[i].style.boxShadow="3px 3px 8px black"
	}
	setTimeout(function(){
		if(obj.getAttribute("expand")=="true"){
			for(var i=1;i<obj.children.length;i++){

				obj.children[i].style.width=obj.children[i].children[0].getBoundingClientRect().width+20+"px"
			}
		}
	},200)
	obj.style.height=obj.children.length*45+"px"
}

function expandHelp(obj){
	var helpInfoTitle=obj.getElementsByClassName("helpInfoTitle")[0]
	var helpTitle_width=helpInfoTitle.getBoundingClientRect().width
	var helpTitle_height=helpInfoTitle.getBoundingClientRect().height
	var helpInfoContent=obj.getElementsByClassName("helpInfoContent")[0]
	var helpContent_width=helpInfoContent.getBoundingClientRect().width
	var helpContent_height=helpInfoContent.getBoundingClientRect().height
	obj.style.width=Math.max(helpTitle_width,helpContent_width)+10+"px";
	obj.style.height=helpTitle_height+helpContent_height+20+"px"
	obj.style.boxShadow="3px 3px 5px grey"
}
function geneSets_PresetChecked(Model_ID){
	var Model=document.getElementById(Model_ID);
	var geneSetsTextarea=Model.getElementsByClassName("geneSetsTextarea")[0];
	geneSetsTextarea.style.display="none"
	var inputAndSearchBox=Model.getElementsByClassName("inputAndSearchBox")[0];
	inputAndSearchBox.style.display=""
	geneSets_Preset=true;
	geneSets_Textarea=false;
	checkMSigDB(currentModel)
	checkStartButton(currentModel);
}

function geneSets_TextareaChecked(Model_ID){
	var Model=document.getElementById(Model_ID);
	var geneSetsTextarea=Model.getElementsByClassName("geneSetsTextarea")[0];
	geneSetsTextarea.style.display=""
	var inputAndSearchBox=Model.getElementsByClassName("inputAndSearchBox")[0];
	inputAndSearchBox.style.display="none"
	geneSets_Textarea=true;
	geneSets_Preset=false;
	
}

function get_chevronsup(){
	var svg=document.createElementNS('http://www.w3.org/2000/svg',"svg");

		svg.setAttribute("width","10");
		svg.setAttribute("height","10");
		svg.setAttribute("viewBox","0 0 10 10");
		svg.setAttribute("class","chevrons-up");
		svg.innerHTML='<polyline points="1 7 5 3 9 7" color="white" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"></polyline>'
	return(svg);
}

function get_searchButton(){
	var svg=document.createElementNS('http://www.w3.org/2000/svg',"svg");


		svg.setAttribute("viewBox","0 0 15 15");
		svg.setAttribute("class","searchButton");
		
		svg.innerHTML='<circle cx="7.5" cy="7.5" r="4" color="grey" fill="none" stroke="currentColor" stroke-width="1" stroke-linecap="round" stroke-linejoin="round">'
		svg.innerHTML+='<line x1='+(7.5+4/Math.sqrt(2)) +' y1='+(7.5+4/Math.sqrt(2))+' x2='+(7.5+7/Math.sqrt(2)) +' y2='+(7.5+7/Math.sqrt(2))+ ' stroke="grey" stroke-width="1" stroke-linecap="round" stroke-linejoin="round">'
	return(svg);
}
function get_clearButton(){
	var svg=document.createElementNS('http://www.w3.org/2000/svg',"svg");


		svg.setAttribute("viewBox","0 0 15 15");
		svg.setAttribute("class","clearButton");
		svg.setAttribute("onclick","clearInput(this)");
		svg.innerHTML+='<line x1="4" y1="4" x2="11" y2="11" stroke="#d25252" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round">'
		svg.innerHTML+='<line x1="4" y1="11" x2="11" y2="4" stroke="#d25252" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round">'
	return(svg);
}


function getGeneList(vennTag){
	console.log(vennTag)
	var realVennTag=PrefixInteger((vennList.index.indexOf(vennTag)+1),3)
	var Time1 = new Date()
	if(geneList[vennTag]==undefined||geneList[vennTag].length==0){
		if(vennTag=="all"){
			geneList[vennTag]=[]
			var vennList_Name=Object.keys(vennList)
			for(var v=0;v<vennList_Name.length;v++){
				geneList[vennTag].push.apply(geneList[vennTag],asArray(vennList[vennList_Name[v]].data))
			}
		}else if(vennTag.indexOf("_methylation_site")>-1){
			var vennList_Name=Object.keys(cgvennList)
			geneList[vennTag]=[]
			for(var v=0;v<vennList_Name.length;v++){
				if(cgvennList[vennList_Name[v]].group.indexOf(vennTag)>-1){
					var orgCgs=asArray(cgvennList[vennList_Name[v]].data)
					for(var o=0;o<orgCgs.length;o++){
						geneList[vennTag].push("cg"+PrefixInteger(orgCgs[o],8))
					}
					//geneList[vennTag].push.apply(geneList[vennTag],orgCgs)
				}
			}
			
		}else{
			var vennList_Name=Object.keys(vennList)
			geneList[vennTag]=[]
			for(var v=0;v<vennList_Name.length;v++){
				if(vennList_Name[v]=="index")continue
				if(vennList[vennList_Name[v]].group.indexOf(realVennTag)>-1){
					var orgGenes=asArray(vennList[vennList_Name[v]].data)
					for(var o=0;o<orgGenes.length;o++){
						geneList[vennTag].push(orgGenes[o])
					}
					
				}
			}
		}
		console.log("合并cgvenn:"+(new Date()-Time1))
		var Time2 = new Date()
		geneList[vennTag]=geneList[vennTag].sort()
		console.log("venn排序："+(new Date()-Time2))
	}else{
		
	}
}


function get_chevronsdown(){
	var svg=document.createElementNS('http://www.w3.org/2000/svg',"svg");

		svg.setAttribute("width","10");
		svg.setAttribute("height","10");
		svg.setAttribute("viewBox","0 0 10 10");
		svg.setAttribute("class","chevrons-down");
		svg.innerHTML='<polyline points="1 3 5 7 9 3" color="white" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"></polyline>'
	return(svg);
}


function get_plus(){
	var svg=document.createElementNS('http://www.w3.org/2000/svg',"svg");

		svg.setAttribute("width","10");
		svg.setAttribute("height","10");
		svg.setAttribute("viewBox","0 0 10 10");
		svg.setAttribute("class","plus");
		svg.innerHTML='<rect width=8 height=8 x=1 y=1 stroke="black" stroke-width=0.5 fill="none"></rect>'
		svg.innerHTML+='<line x1=2 y1=5 x2=8 y2=5  fill="none" stroke="black" stroke-width="0.5"  ></line>'
		svg.innerHTML+='<line x1=5 y1=2 x2=5 y2=8  fill="none" stroke="black" stroke-width="0.5"  ></line>'
	return(svg);
}


function get_minus(){
	var svg=document.createElementNS('http://www.w3.org/2000/svg',"svg");

		svg.setAttribute("width","10");
		svg.setAttribute("height","10");
		svg.setAttribute("viewBox","0 0 10 10");
		svg.setAttribute("class","minus");
		svg.innerHTML='<rect width=8 height=8 x=1 y=1 stroke="black" stroke-width=0.5 fill="none"></rect>'
		svg.innerHTML+='<line x1=2 y1=5 x2=8 y2=5  fill="none" stroke="black" stroke-width="0.5"  ></line>'
	return(svg);
}


function getlength(obj){
	return(Object.keys(obj).length)
	
}


function guid() {
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
    var r = Math.random()*16|0, v = c == 'x' ? r : (r&0x3|0x8);
    return v.toString(16);
  });
}

function getOffsetTop(el){
 return el.offsetParent
  ? el.offsetTop + getOffsetTop(el.offsetParent)
  : el.offsetTop
}

function getOffsetLeft(el){
 return el.offsetParent
  ? el.offsetLeft + getOffsetLeft(el.offsetParent)
  : el.offsetLeft
}

function getScrollTop(obj){
	var scrollTop=0;
	if(obj&&obj.scrollTop){
		scrollTop=obj.scrollTop;
	}
	return scrollTop;
}

function initHTML(){
	funcDebug=[]

	var xmlhttp=new XMLHttpRequest();
		//document.getElementById("plotarea").innerHTML +="chrome1";

	xmlhttp.onreadystatechange=function(){
		//document.getElementById("plotarea").innerHTML =genename+"got:"+xmlhttp.responseText;
		if (xmlhttp.readyState==4 && xmlhttp.status==200){
			initParams_str =xmlhttp.responseText.toString();
			initParams = JSON.parse(initParams_str);
			
			initModels();
			initModelChildren();
			
			popupInfoBox=document.getElementById("popupInfoBox");
			popupInfoBoxBg=document.getElementById("popupInfoBoxBg");
			popupInfoBoxContent=document.getElementById("popupInfoBoxContent");
			popupInfoBox.style.opacity=0;
			popupInfoBoxContent.addEventListener("mouseleave",isOnInfoBox,false);
			popupInfoBoxContent.addEventListener("mouseenter",isOnInfoBox,false);
			promotBox=document.getElementById("promotBox");
			promotBox.addEventListener("mouseleave",isOnPromotBox,false);
			promotBox.addEventListener("mouseenter",isOnPromotBox,false);
			isOutPromotBox=true
			isOnInput=false
			var contactus=document.getElementsByClassName("Contact")[0]
			shrinkHelp(contactus)
			document.body.removeAttribute("style");
			userParams=new Object();
			isWarnings=new Object();
			Params=new Object();
			var dataBaseNames=Object.keys(initParams.cellLines)
			model_ids=new Object();
			for(var d=0;d<dataBaseNames.length;d++){
				var cancerTypes=Object.keys(initParams.cellLines[dataBaseNames[d]])
				for(var c=0;c<cancerTypes.length;c++){
					var cellLinePool=Object.keys(initParams.cellLines[dataBaseNames[d]][cancerTypes[c]].cellLines)
					for(var l=0;l<cellLinePool.length;l++){
						model_ids[cellLinePool[l]]=initParams.cellLines[dataBaseNames[d]][cancerTypes[c]].cellLines[cellLinePool[l]].model_id
					}
				}
			}
			updateHomeDiagram()
			updateDownload()
			helpDone=false
			updateVisits()
			recordvisit("recordVisit.php")
		}else{
			
		}
	}
	//document.getElementById("plotarea").innerHTML= +"start";
	xmlhttp.open("GET","php/getInitParams.php",true);
	xmlhttp.send();
	

	var vennList_xmlhttp=new XMLHttpRequest();
		vennList_xmlhttp.onreadystatechange=function(){
		//document.getElementById("plotarea").innerHTML =genename+"got:"+xmlhttp.responseText;
		if (vennList_xmlhttp.readyState==4 && vennList_xmlhttp.status==200){
			var vennList_str =vennList_xmlhttp.responseText.toString();
			vennList = JSON.parse(vennList_str);
			geneList=new Object();
		}else{	
		}
	}
	//document.getElementById("plotarea").innerHTML= +"start";
	vennList_xmlhttp.open("GET","php/getvennList.php",true);
	vennList_xmlhttp.send();
	
	var cgvennList_xmlhttp=new XMLHttpRequest();
		cgvennList_xmlhttp.onreadystatechange=function(){
		//document.getElementById("plotarea").innerHTML =genename+"got:"+xmlhttp.responseText;
		if (cgvennList_xmlhttp.readyState==4 && cgvennList_xmlhttp.status==200){
			var cgvennList_str =cgvennList_xmlhttp.responseText.toString();
			cgvennList = JSON.parse(cgvennList_str);
		}else{	
		}
	}
	//document.getElementById("plotarea").innerHTML= +"start";
	cgvennList_xmlhttp.open("GET","php/getcgvennList.php",true);
	cgvennList_xmlhttp.send();
	 
	var modifyData_xmlhttp=new XMLHttpRequest();
		modifyData_xmlhttp.onreadystatechange=function(){
		//document.getElementById("plotarea").innerHTML =genename+"got:"+xmlhttp.responseText;
		if (modifyData_xmlhttp.readyState==4 && modifyData_xmlhttp.status==200){
			var ModifyData_str =modifyData_xmlhttp.responseText.toString();
			//document.getElementById("modifyData").appendChild(modifyData);
			var modifyData=JSON.parse(ModifyData_str);
			var output="";
			for(var i=0;i<modifyData.length;i++){
				output+=modifyData[i]+"<br/>";
			}
			document.getElementById("modifyData").innerHTML=output;

		}else{	
		}
	}
	//document.getElementById("plotarea").innerHTML= +"start";
	modifyData_xmlhttp.open("GET","php/getModifyData.php",true);
	modifyData_xmlhttp.send();
	gain=new Object();
	originalWidth=new Object();
	originalHeight=new Object();
}


function initModels(){
	
	var dataBaseNames= Object.keys(initParams.ModelInfo);
	var modelsTips=document.getElementById("modelsTips");
	var main=document.getElementById("main");
	for(var d=0;d<dataBaseNames.length;d++){
		var dataBase_div=document.createElement("div")
			dataBase_div.setAttribute("class","dataBaseTip")
			dataBase_div.setAttribute("onmouseenter","expandModelTips(this)")
			dataBase_div.setAttribute("onmouseleave","shrinkModelTips(this)")
			var label=document.createElement("label");
				var p=document.createElement("p");
					p.innerHTML=dataBaseNames[d];
				label.appendChild(p);
			dataBase_div.appendChild(label)
		var ModelNames= Object.keys(initParams.ModelInfo[dataBaseNames[d]]);
		for(var m=0;m<ModelNames.length;m++){
			
			var label=document.createElement("label");
				var p=document.createElement("p");
					p.innerHTML=ModelNames[m];

				label.appendChild(p);
				label.setAttribute("name",dataBaseNames[d]+","+ModelNames[m]);
				label.setAttribute("onclick","switchModel(this)");
			dataBase_div.appendChild(label)
			
			//插入Model
			var div=document.createElement("div");
			div.setAttribute("class","Model");
			div.setAttribute("id",dataBaseNames[d]+","+ModelNames[m]);
			div.style.opacity=0;
			div.style.pointerEvents="none";
			div.style.zIndex="1";

			
			var childNodes=clearNodelist(main.childNodes);
			
			main.insertBefore(div,childNodes[childNodes.length-2])
		}
		var childNodes=clearNodelist(modelsTips.childNodes);
			
		modelsTips.insertBefore(dataBase_div,childNodes[childNodes.length-2])
		
	}
	
	
	
	currentModel="NATIVE,HOME";
	//initModelsTipBg()
	
}

function initModelChildren(){
	dealingCellLines=new Object();
	var dataBaseNames= Object.keys(initParams.ModelInfo);
	for(var d=0;d<dataBaseNames.length;d++){
		var ModelNames= Object.keys(initParams.ModelInfo[dataBaseNames[d]]);
		
		for(var i=0;i<ModelNames.length;i++){
			var isFirstFoldList=true;
			var Model=document.getElementById(dataBaseNames[d]+","+ModelNames[i]);
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("sidebar")){
				var div=document.createElement("div");
				div.setAttribute("class","sidebar");
				Model.appendChild(div)
			}
			var foldPointlessID=pointlessID()
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("omicsSelector")){
		
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						
						foldList_folderInput.setAttribute("style","display:none;");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var p1=document.createElement("p");
								p1.innerHTML="Select Omics:";
							foldList_label.appendChild(p1);
						foldList_bg_title.appendChild(foldList_label)
						
					var foldList=document.createElement("div");
						foldList.setAttribute("class","foldList");
						var omicsSelectorDiv=document.createElement("div");
						var p=document.createElement("p");
							p.innerHTML="Please select a omics:";
						omicsSelectorDiv.appendChild(p);
						var omicsSelector=document.createElement("select");
							if(dataBaseNames[d]=="GDSC"){
								var option1=document.createElement("option");
									option1.value="mRNA";
									option1.innerHTML="mRNA";
								var option2=document.createElement("option");
									option2.value="mutation";
									option2.innerHTML="Somatic Mutation";
								var option3=document.createElement("option");
									option3.value="methylation_promoter";
									option3.innerHTML="DNA Methylation (Promoter)";
								var option4=document.createElement("option");
									option4.value="methylation_body";
									option4.innerHTML="DNA Methylation (Gene body)";
								var option5=document.createElement("option");
									option5.value="methylation_site";
									option5.innerHTML="DNA Methylation (site)";
								var option6=document.createElement("option");
									option6.value="CNV";
									option6.innerHTML="Copy Number Variation";
								omicsSelector.appendChild(option1);
								omicsSelector.appendChild(option2);
								omicsSelector.appendChild(option3);
								omicsSelector.appendChild(option4);
								omicsSelector.appendChild(option5);
								omicsSelector.appendChild(option6);

							
							}else if(dataBaseNames[d]=="DepMap"){
								var option1=document.createElement("option");
									option1.value="mRNA";
									option1.innerHTML="mRNA";
								var option2=document.createElement("option");
									option2.value="mutation";
									option2.innerHTML="Somatic Mutation";
								var option3=document.createElement("option");
									option3.value="methylation_promoter";
									option3.innerHTML="DNA Methylation (Promoter)";
								var option4=document.createElement("option");
									option4.value="Chromatin";
									option4.innerHTML="Global Chromatin Profiling";
								var option5=document.createElement("option");
									option5.value="Metabolomics";
									option5.innerHTML="Metabolomics";
								var option6=document.createElement("option");
									option6.value="CNV_continuous";
									option6.innerHTML="Copy Number Variation";
								var option7=document.createElement("option");
									option7.value="protein";
									option7.innerHTML="Protein";
								var option8=document.createElement("option");
									option8.value="microRNA";
									option8.innerHTML="microRNA";
								var option9=document.createElement("option");
									option9.value="Metastatic";
									option9.innerHTML="Metastatic Potential";
								var option10=document.createElement("option");
									option10.value="CRISPR";
									option10.innerHTML="CRISPR";
								var option11=document.createElement("option");
									option11.value="RNAi_combine";
									option11.innerHTML="RNAi";
								omicsSelector.appendChild(option1);
								omicsSelector.appendChild(option2);
								omicsSelector.appendChild(option3);
								omicsSelector.appendChild(option4);
								omicsSelector.appendChild(option5);
								omicsSelector.appendChild(option6);
								omicsSelector.appendChild(option7);
								omicsSelector.appendChild(option8);
								omicsSelector.appendChild(option9);
								omicsSelector.appendChild(option10);
								omicsSelector.appendChild(option11);
								
							}else if(dataBaseNames[d]=="NCI60"){
								var option1=document.createElement("option");
									option1.value="mRNA";
									option1.innerHTML="mRNA";
								var option2=document.createElement("option");
									option2.value="mutation";
									option2.innerHTML="Somatic Mutation";
								var option3=document.createElement("option");
									option3.value="methylation_promoter";
									option3.innerHTML="DNA Methylation (Promoter)";
								var option4=document.createElement("option");
									option4.value="methylation_body";
									option4.innerHTML="DNA Methylation (Gene body)";
								var option5=document.createElement("option");
									option5.value="methylation_site";
									option5.innerHTML="DNA Methylation (site)";
								var option6=document.createElement("option");
									option6.value="CNV";
									option6.innerHTML="Copy Number Variation";
								var option7=document.createElement("option");
									option7.value="protein";
									option7.innerHTML="Protein";
								var option8=document.createElement("option");
									option8.value="microRNA";
									option8.innerHTML="microRNA";
								omicsSelector.appendChild(option1);
								omicsSelector.appendChild(option2);
								omicsSelector.appendChild(option3);
								omicsSelector.appendChild(option4);
								omicsSelector.appendChild(option5);
								omicsSelector.appendChild(option6);
								omicsSelector.appendChild(option7);
								omicsSelector.appendChild(option8);
							}
							
							omicsSelector.setAttribute("class","omicsSelector");
							omicsSelector.setAttribute("onchange","initList_foldList_cellLines(currentModel);omicsChanged(currentModel);checkStartButton(currentModel)")
						omicsSelectorDiv.appendChild(omicsSelector);
						foldList.appendChild(omicsSelectorDiv);
							var geneInputDiv=document.createElement("div")
								geneInputDiv.setAttribute("class","inputgene");
								var p=document.createElement("p")
									p.innerHTML="Please input a gene symbol:"
									p.setAttribute("class","genenameTag")
								var inputComp=document.createElement("div")
									inputComp.setAttribute("class","inputComp")
									var input=document.createElement("input");
										input.setAttribute("class","genename searchInput");
										
										input.setAttribute("value","TP53");
										input.setAttribute("placeholder","Gene symbol here");
										input.setAttribute("onblur","closePromot(this)");
										input.setAttribute("onfocus","showPromot(this)");
										input.setAttribute("oninput","updatePromot(this)");
										input.setAttribute("id",pointlessID());
									var searchButton=get_searchButton()
									inputComp.appendChild(input)
									inputComp.appendChild(searchButton)
									var clearButton=get_clearButton()
										clearButton.setAttribute("onclick",'changeInput("'+pointlessID(0)+'","")')
									inputComp.appendChild(clearButton)
							geneInputDiv.appendChild(p);
							geneInputDiv.appendChild(inputComp);
						foldList.appendChild(geneInputDiv);
							var individualCheckBox_div=document.createElement("div");
								individualCheckBox_div.setAttribute("class","individualCheckBox")
								individualCheckBox_div.style.display="none"
								var div=document.createElement("div")
									var individualCheckBox=document.createElement("input");
										individualCheckBox.setAttribute("type","checkbox")
										individualCheckBox.setAttribute("onchange","disableIndividual(currentModel)")
									var p=document.createElement("p");
										p.innerHTML="Show Individual Cell Lines"
									div.appendChild(individualCheckBox)
									div.appendChild(p)
								var select=document.createElement("select")
									select.setAttribute("disabled",true)
									var cancerType=Object.keys(initParams.cellLines[dataBaseNames[d]])
									for(var c=0;c<cancerType.length;c++){
										var option=document.createElement("option")
											option.value=cancerType[c];
											option.innerHTML=cancerType[c];
											select.appendChild(option)
									}
								individualCheckBox_div.appendChild(div)
								
								individualCheckBox_div.appendChild(select)
						foldList.appendChild(individualCheckBox_div);
					foldList_bg.appendChild(foldList_bg_title)
					foldList_bg.appendChild(foldList)
				sidebar.appendChild(foldList_bg)
					
			}
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("omicsSelectorAndMultigenes")){
		
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						
						foldList_folderInput.setAttribute("style","display:none;");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var p1=document.createElement("p");
								p1.innerHTML="Select Omics:";
							foldList_label.appendChild(p1);
						foldList_bg_title.appendChild(foldList_label)
						
					var foldList=document.createElement("div");
						foldList.setAttribute("class","foldList");
						var omicsSelectorDiv=document.createElement("div");
						var p=document.createElement("p");
							p.innerHTML="Please select a omics:";
						omicsSelectorDiv.appendChild(p);
						var omicsSelectorAndMultigenes=document.createElement("select");

							if(dataBaseNames[d]=="GDSC"){
								var option1=document.createElement("option");
									option1.value="mRNA";
									option1.innerHTML="mRNA";

								var option3=document.createElement("option");
									option3.value="methylation_promoter";
									option3.innerHTML="DNA Methylation (Promoter)";
								var option4=document.createElement("option");
									option4.value="methylation_body";
									option4.innerHTML="DNA Methylation (Gene body)";
								var option5=document.createElement("option");
									option5.value="methylation_site";
									option5.innerHTML="DNA Methylation (site)";

								omicsSelectorAndMultigenes.appendChild(option1);

								omicsSelectorAndMultigenes.appendChild(option3);
								omicsSelectorAndMultigenes.appendChild(option4);
								omicsSelectorAndMultigenes.appendChild(option5);


							}else if(dataBaseNames[d]=="DepMap"){
								var option1=document.createElement("option");
									option1.value="mRNA";
									option1.innerHTML="mRNA";
								var option3=document.createElement("option");
									option3.value="methylation_promoter";
									option3.innerHTML="DNA Methylation (Promoter)";
								var option4=document.createElement("option");
									option4.value="Chromatin";
									option4.innerHTML="Global Chromatin Profiling";
								var option5=document.createElement("option");
									option5.value="Metabolomics";
									option5.innerHTML="Metabolomics";
								var option6=document.createElement("option");
									option6.value="CNV_continuous";
									option6.innerHTML="Copy Number Variation";
								var option7=document.createElement("option");
									option7.value="protein";
									option7.innerHTML="Protein";
								var option8=document.createElement("option");
									option8.value="microRNA";
									option8.innerHTML="microRNA";
								var option9=document.createElement("option");
									option9.value="Metastatic";
									option9.innerHTML="Metastatic Potential";
								var option10=document.createElement("option");
									option10.value="CRISPR";
									option10.innerHTML="CRISPR";
								var option11=document.createElement("option");
									option11.value="RNAi_combine";
									option11.innerHTML="RNAi";
								omicsSelectorAndMultigenes.appendChild(option1);
								omicsSelectorAndMultigenes.appendChild(option3);
								omicsSelectorAndMultigenes.appendChild(option4);
								omicsSelectorAndMultigenes.appendChild(option5);
								omicsSelectorAndMultigenes.appendChild(option6);
								omicsSelectorAndMultigenes.appendChild(option7);
								omicsSelectorAndMultigenes.appendChild(option8);
								omicsSelectorAndMultigenes.appendChild(option9);
								omicsSelectorAndMultigenes.appendChild(option10);
								omicsSelectorAndMultigenes.appendChild(option11);
							}else if(dataBaseNames[d]=="NCI60"){
								var option1=document.createElement("option");
									option1.value="mRNA";
									option1.innerHTML="mRNA";
								var option2=document.createElement("option");
									option2.value="mutation";
									option2.innerHTML="Somatic Mutation";
								var option3=document.createElement("option");
									option3.value="methylation_promoter";
									option3.innerHTML="DNA Methylation (Promoter)";
								var option4=document.createElement("option");
									option4.value="methylation_body";
									option4.innerHTML="DNA Methylation (Gene body)";
								var option5=document.createElement("option");
									option5.value="methylation_site";
									option5.innerHTML="DNA Methylation (site)";
								var option6=document.createElement("option");
									option6.value="CNV";
									option6.innerHTML="Copy Number Variation";
								var option7=document.createElement("option");
									option7.value="protein";
									option7.innerHTML="Protein";
								var option8=document.createElement("option");
									option8.value="microRNA";
									option8.innerHTML="microRNA";
								omicsSelectorAndMultigenes.appendChild(option1);
								omicsSelectorAndMultigenes.appendChild(option2);
								omicsSelectorAndMultigenes.appendChild(option3);
								omicsSelectorAndMultigenes.appendChild(option4);
								omicsSelectorAndMultigenes.appendChild(option5);
								omicsSelectorAndMultigenes.appendChild(option6);
								omicsSelectorAndMultigenes.appendChild(option7);
								omicsSelectorAndMultigenes.appendChild(option8);
							}
							
							
							omicsSelectorAndMultigenes.setAttribute("class","omicsSelectorAndMultigenes");
							omicsSelectorAndMultigenes.setAttribute("onchange","initList_foldList_cellLines(currentModel);clusterOmicChanged(currentModel);checkStartButton(currentModel)")
						omicsSelectorDiv.appendChild(omicsSelectorAndMultigenes);
						foldList.appendChild(omicsSelectorDiv);
							var geneInputDiv=document.createElement("div")
								geneInputDiv.setAttribute("class","inputgene");
								var p=document.createElement("p")
									p.innerHTML="Please input gene symbols:"
									p.setAttribute("class","genenameTag")
								var inputComp=document.createElement("div")
									inputComp.setAttribute("class","inputComp")
									var input=document.createElement("textarea");
										input.setAttribute("class","geneSetsTextarea");
										
										input.value='A1BG,A1CF,A2M,A2ML1,A4GALT,A4GNT,AAAS,AACS';
										input.setAttribute("placeholder",'Paste your gene list here. Genes should be seperated by "," or "line break" like below:\nTP53, GATA3\nGATA2,A1CF,A1BG')
										input.setAttribute("onblur","checkGenes(this)");
										input.setAttribute("id",pointlessID());
									inputComp.appendChild(input)
							geneInputDiv.appendChild(p);
							geneInputDiv.appendChild(inputComp);
						foldList.appendChild(geneInputDiv);
					foldList_bg.appendChild(foldList_bg_title)
					foldList_bg.appendChild(foldList)
				sidebar.appendChild(foldList_bg)
			}
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("multiOmicsSelector")){
			
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						
						foldList_folderInput.setAttribute("style","display:none;");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var p1=document.createElement("p");
								p1.innerHTML="Select Targets:";
							foldList_label.appendChild(p1);
						foldList_bg_title.appendChild(foldList_label)
						
					var foldList=document.createElement("div");
						foldList.setAttribute("class","foldList");

							var OmicInputBox=document.createElement("div")
								var p=document.createElement("p");
									p.innerHTML="Please select a omics:";
								var omicsSelector=document.createElement("select");
								if(dataBaseNames[d]=="GDSC"){
									var option1=document.createElement("option");
										option1.value="mRNA";
										option1.innerHTML="mRNA";
									var option2=document.createElement("option");
										option2.value="mutation";
										option2.innerHTML="Somatic Mutation";
									var option3=document.createElement("option");
										option3.value="CNV";
										option3.innerHTML="Copy Number Variation";
									var option4=document.createElement("option");
										option4.value="methylation_promoter";
										option4.innerHTML="DNA Methylation (Promoter)";
									var option5=document.createElement("option");
										option5.value="methylation_body";
										option5.innerHTML="DNA Methylation (Gene body)";
									var option6=document.createElement("option");
										option6.value="methylation_site";
										option6.innerHTML="DNA Methylation (site)";
									
									omicsSelector.appendChild(option1);
									omicsSelector.appendChild(option2);
									omicsSelector.appendChild(option3);
									omicsSelector.appendChild(option4);
									omicsSelector.appendChild(option5);
									omicsSelector.appendChild(option6);

								
								}else if(dataBaseNames[d]=="DepMap"){
									var option1=document.createElement("option");
										option1.value="mRNA";
										option1.innerHTML="mRNA";
									var option2=document.createElement("option");
										option2.value="mutation";
										option2.innerHTML="Somatic Mutation";
									var option3=document.createElement("option");
										option3.value="CNV_continuous";
										option3.innerHTML="Copy Number Variation";
									var option4=document.createElement("option");
										option4.value="methylation_promoter";
										option4.innerHTML="DNA Methylation (Promoter)";
									var option5=document.createElement("option");
										option5.value="Chromatin";
										option5.innerHTML="Global Chromatin Profiling";
									var option6=document.createElement("option");
										option6.value="Metabolomics";
										option6.innerHTML="Metabolomics";
									
									var option7=document.createElement("option");
										option7.value="protein";
										option7.innerHTML="Protein";
									var option8=document.createElement("option");
										option8.value="microRNA";
										option8.innerHTML="microRNA";
									var option9=document.createElement("option");
										option9.value="Metastatic";
										option9.innerHTML="Metastatic Potential";
									var option10=document.createElement("option");
										option10.value="CRISPR";
										option10.innerHTML="CRISPR";
									var option11=document.createElement("option");
										option11.value="RNAi_combine";
										option11.innerHTML="RNAi";
									omicsSelector.appendChild(option1);
									omicsSelector.appendChild(option2);
									omicsSelector.appendChild(option3);
									omicsSelector.appendChild(option4);
									omicsSelector.appendChild(option5);
									omicsSelector.appendChild(option6);
									omicsSelector.appendChild(option7);
									omicsSelector.appendChild(option8);
									omicsSelector.appendChild(option9);
									omicsSelector.appendChild(option10);
									omicsSelector.appendChild(option11);
									
								}else if(dataBaseNames[d]=="NCI60"){
									var option1=document.createElement("option");
										option1.value="mRNA";
										option1.innerHTML="mRNA";
									var option2=document.createElement("option");
										option2.value="mutation";
										option2.innerHTML="Somatic Mutation";
									var option3=document.createElement("option");
										option3.value="CNV";
										option3.innerHTML="Copy Number Variation";
									var option4=document.createElement("option");
										option4.value="methylation_promoter";
										option4.innerHTML="DNA Methylation (Promoter)";
									var option5=document.createElement("option");
										option5.value="methylation_body";
										option5.innerHTML="DNA Methylation (Gene body)";
									var option6=document.createElement("option");
										option6.value="methylation_site";
										option6.innerHTML="DNA Methylation (site)";
									
									var option7=document.createElement("option");
										option7.value="protein";
										option7.innerHTML="Protein";
									var option8=document.createElement("option");
										option8.value="microRNA";
										option8.innerHTML="microRNA";
									omicsSelector.appendChild(option1);
									omicsSelector.appendChild(option2);
									omicsSelector.appendChild(option3);
									omicsSelector.appendChild(option4);
									omicsSelector.appendChild(option5);
									omicsSelector.appendChild(option6);
									omicsSelector.appendChild(option7);
									omicsSelector.appendChild(option8);
								}
								omicsSelector.setAttribute("class","omicsSelector");
								omicsSelector.setAttribute("onchange","omicsChanged(currentModel)")
							OmicInputBox.appendChild(p)
							OmicInputBox.appendChild(omicsSelector)
							
								var inputComp=document.createElement("div")
									inputComp.setAttribute("class","inputComp")
										var p=document.createElement("p")
										p.innerHTML="Please input a gene symbol:"
										p.setAttribute("class","genenameTag")
									var input=document.createElement("input");
										input.setAttribute("class","genename searchInput");
										
										input.setAttribute("value","");
										input.setAttribute("placeholder","Gene symbol here");
										input.setAttribute("onblur","closePromot(this)");
										input.setAttribute("onfocus","showPromot(this)");
										input.setAttribute("oninput","updatePromot(this)");
										input.setAttribute("id",pointlessID());
									var searchButton=get_searchButton()
									inputComp.appendChild(input)
									inputComp.appendChild(searchButton)
									var clearButton=get_clearButton()
										clearButton.setAttribute("onclick",'changeInput("'+pointlessID(0)+'","")')
									inputComp.appendChild(clearButton)
							OmicInputBox.appendChild(p)
							OmicInputBox.appendChild(inputComp)
							
								var addButton_div=document.createElement("div")
									
									addButton_div.style.width="40px"
									addButton_div.style.height="40px"
									addButton_div.style.margin="auto"
									var addButton=SVG(addButton_div)
										var group = addButton.group().attr({"onclick":"addOmics(currentModel);checkStartButton(currentModel)",class:"justPointer addButton"})
										addButton.viewbox(0,0,40,40)
										var circle=addButton.circle(30).fill('white').attr({
											cx:20,
											cy:20
										})
										var line1=addButton.line(10, 20, 30, 20).attr({
											"stroke-width": 2,
											"stroke-linecap":"round",
											'stroke-linejoin':"round",
											'stroke':"grey"
										})
										var line2=addButton.line(20, 10, 20, 30).attr({
											"stroke-width": 2,
											'stroke-linecap':"round",
											'stroke-linejoin':"round",
											'stroke':"grey"
										})
										group.add(circle)
										group.add(line1)
										group.add(line2)
									
								
							OmicInputBox.appendChild(addButton_div)
						foldList.appendChild(OmicInputBox)
						var OmicResultBox=document.createElement("div")
							OmicResultBox.setAttribute("class","OmicResultBox")
								var p=document.createElement("p")
								p.innerHTML="<p style='text-align: center'>Selected compare components:"+"<br><span style='font-size:12px;font-style:italic;'>(The red component is the reference)</span><p>"
							OmicResultBox.appendChild(p)
							
						foldList.appendChild(OmicResultBox)
					foldList_bg.appendChild(foldList_bg_title)
					foldList_bg.appendChild(foldList)
				sidebar.appendChild(foldList_bg)
				//设置3个默认值
			
				var omics=[0,1,2]
				var genes=["TP53","A1BG","A2ML1"]
				for(var index=0;index<omics.length;index++){
					//console.log(index)
					omicsSelector[omics[index]].selected=true;
					input.value=genes[index];
					addOmics(dataBaseNames[d]+","+ModelNames[i])
				}
			
			}
			
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("onlyGene")){
		
				var sidebar=Model.getElementsByClassName("sidebar")[0];
					var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						
						foldList_folderInput.setAttribute("style","display:none;");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var p1=document.createElement("p");
								p1.innerHTML="Input Gene:";
							foldList_label.appendChild(p1);
						foldList_bg_title.appendChild(foldList_label)
						
					var foldList=document.createElement("div");
						foldList.setAttribute("class","foldList");
							var geneInputDiv=document.createElement("div")
								geneInputDiv.setAttribute("class","inputgene");
								var p=document.createElement("p")
									p.innerHTML="Please input a gene symbol:"
									p.setAttribute("class","genenameTag")
								var inputComp=document.createElement("div")
									inputComp.setAttribute("class","inputComp")
									var input=document.createElement("input");
										input.setAttribute("class","genename searchInput");
										
										input.setAttribute("value","TP53");
										input.setAttribute("placeholder","Gene symbol here");
										input.setAttribute("onblur","closePromot(this)");
										input.setAttribute("onfocus","showPromot(this)");
										input.setAttribute("oninput","updatePromot(this)");
										input.setAttribute("id",pointlessID());
									var searchButton=get_searchButton()
									inputComp.appendChild(input)
									inputComp.appendChild(searchButton)
									var clearButton=get_clearButton()
										clearButton.setAttribute("onclick",'changeInput("'+pointlessID(0)+'","")')
									inputComp.appendChild(clearButton)
								geneInputDiv.appendChild(p);
								geneInputDiv.appendChild(inputComp);

						foldList.appendChild(geneInputDiv);
					foldList_bg.appendChild(foldList_bg_title)
					foldList_bg.appendChild(foldList)
				sidebar.appendChild(foldList_bg)
				
			}
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("geneSetSelector")){
		
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					foldList_bg.style.flexShrink=0;
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						
						
						foldList_folderInput.setAttribute("style","display:none;");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
							
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var p1=document.createElement("p");
								p1.innerHTML="Input molecular signature:";
							foldList_label.appendChild(p1);
						foldList_bg_title.appendChild(foldList_label);
						foldList_bg_title.appendChild(get_chevronsup());
						foldList_bg_title.appendChild(get_chevronsdown());
					foldList_bg.appendChild(foldList_bg_title);

					var foldList_geneSets=document.createElement("div");
						foldList_geneSets.setAttribute("class","foldList foldList_geneSets");
							var op1=document.createElement("div");
								op1.setAttribute("class","geneSetsInputType")
								var inputName=pointlessID();
								var input1=document.createElement("input");
									input1.setAttribute("type","radio")
									input1.setAttribute("name",inputName);
									input1.setAttribute("checked",true);
									input1.setAttribute("onchange","geneSets_TextareaChecked(currentModel)");
									input1.setAttribute("id",pointlessID());
								var label1=document.createElement("label")
									label1.setAttribute("for",pointlessID(0))
									label1.innerHTML="<p>Based on user-defined molecular signature:</p>"

								op1.appendChild(input1)
								op1.appendChild(label1)
							var op2=document.createElement("div");
								op2.setAttribute("class","geneSetsInputType")
								var input2=document.createElement("input");
									input2.setAttribute("type","radio")
									input2.setAttribute("name",inputName);
									input2.setAttribute("onchange","geneSets_PresetChecked(currentModel)");
									input2.setAttribute("id",pointlessID());
								var label2=document.createElement("label")
									label2.setAttribute("for",pointlessID(0))
									label2.innerHTML="<p>From MSigDB:</p>"
		
								op2.appendChild(input2)
								op2.appendChild(label2)
							var opALL=document.createElement("div");
								opALL.appendChild(op1)
								opALL.appendChild(op2)
								opALL.setAttribute("class","geneSetsInputTypeALL")
							var textArea_geneSets=document.createElement("textarea");
								textArea_geneSets.setAttribute("class","geneSetsTextarea")
								textArea_geneSets.setAttribute("placeholder",'Paste your gene list here. Genes should be seperated by "," or "line break" like below:\nTP53, GATA3\nGATA2,A1CF,A1BG')
								textArea_geneSets.setAttribute("onblur","checkGenes(this)")
								textArea_geneSets.value='A1BG,A1CF,A2M,A2ML1,A4GALT,A4GNT,AAAS,AACS';
							var inputAndSearchBox=document.createElement("div")
								inputAndSearchBox.setAttribute("class","inputAndSearchBox");
								inputAndSearchBox.style.display="none"
								var inputAndSearch_Input=document.createElement("div")
									inputAndSearch_Input.setAttribute("class","inputAndSearch_Input inputComp");
									var inputAndSearch_Input_input=document.createElement("input");
										inputAndSearch_Input_input.setAttribute("class","inputAndSearch_Input_input searchInput");
										inputAndSearch_Input_input.setAttribute("onfocus","showSearchBox(this)");
										inputAndSearch_Input_input.setAttribute("oninput","updateSearchBox(this)");
										inputAndSearch_Input_input.setAttribute("onblur","closeSearchBox(this)");
										
										inputAndSearch_Input_input.setAttribute("id",pointlessID());
										inputAndSearch_Input_input.setAttribute("placeholder","Search for Gene Set")
									var searchButton=get_searchButton()
									inputAndSearch_Input.appendChild(inputAndSearch_Input_input)
									inputAndSearch_Input.appendChild(searchButton)
								var inputAndSearch_Search=document.createElement("div");
									inputAndSearch_Search.setAttribute("class","inputAndSearch_Search");
									inputAndSearch_Search.setAttribute("onfocus","inputAndSearch_Search_onfocus(this)");
									inputAndSearch_Search.setAttribute("onblur","inputAndSearch_Search_onblur(this)");
									inputAndSearch_Search.setAttribute("tabindex","0");
									inputAndSearch_Search.setAttribute("hidefocus","true");
									inputAndSearch_Search.style.display="none";
								inputAndSearchBox.appendChild(inputAndSearch_Input)
								inputAndSearchBox.appendChild(inputAndSearch_Search)
						foldList_geneSets.appendChild(opALL);
						foldList_geneSets.appendChild(textArea_geneSets);
						foldList_geneSets.appendChild(inputAndSearchBox);
						updateSearchBox(inputAndSearch_Input_input)
					foldList_bg.appendChild(foldList_geneSets);
				sidebar.appendChild(foldList_bg);
				geneSets_Textarea=true;
			}
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("foldList_drugs")|initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("foldList_drugs_single")){
				if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("foldList_drugs_single")){
					var isSingle=true;
					var title="Select Drug:";
				}else{
					var isSingle=false;
					var title="Select Drug(s):";
				}
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());

						foldList_folderInput.setAttribute("style","display:none;");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var foldlist_inputAll_label=document.createElement("label")
								foldlist_inputAll_label.setAttribute("for",pointlessID());
								foldlist_inputAll_label.setAttribute("class","selectALL")
								var foldlist_inputAll=document.createElement("input")
									foldlist_inputAll.type="checkbox";
									if(!isSingle){
										foldlist_inputAll.setAttribute("onchange","selectAllDrugs(this)")
									}else{
										foldlist_inputAll.style.display="none"
									}
									foldlist_inputAll.setAttribute("id",pointlessID(0))
									foldlist_inputAll.setAttribute("class","drug_selectorALL")
									foldlist_inputAll.style.display="none"
								var p=document.createElement("p")
									
								
									p.innerHTML="Select<br>All"
								foldlist_inputAll_label.appendChild(p)
								foldlist_inputAll_label.appendChild(foldlist_inputAll);
							var p1=document.createElement("p");
								p1.innerHTML=title;
							var p2=document.createElement("p");
								p2.setAttribute("class","totalCount_Drug");
							foldList_label.appendChild(foldlist_inputAll_label);
							foldList_label.appendChild(p1);
							foldList_label.appendChild(p2);
						
						foldList_bg_title.appendChild(foldList_label);
						foldList_bg_title.appendChild(get_chevronsup());
						foldList_bg_title.appendChild(get_chevronsdown());
					foldList_bg.appendChild(foldList_bg_title);
					var foldList_drugs=document.createElement("div");
						foldList_drugs.setAttribute("class","foldList foldList_drugs");
					foldList_bg.appendChild(foldList_drugs);
				sidebar.appendChild(foldList_bg);
				
				initList_foldList_drugs(dataBaseNames[d]+","+ModelNames[i],isSingle=isSingle);
				countFoldList_drugs(dataBaseNames[d]+","+ModelNames[i]);
			}
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("foldList_cellLines")|initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("foldList_cellLines_single")){
				if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("foldList_cellLines_single")){
					var isSingle=true;
					var title="Select Cell Line:";
				}else{
					var isSingle=false;
					var title="Select Cell Line(s):";
				}
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						
						foldList_folderInput.setAttribute("style","display:none;");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var foldlist_inputAll_label=document.createElement("label")
								foldlist_inputAll_label.setAttribute("for",pointlessID());
								foldlist_inputAll_label.setAttribute("class","selectALL")
								var foldlist_inputAll=document.createElement("input")
									foldlist_inputAll.type="checkbox";
									if(!isSingle){
										foldlist_inputAll.setAttribute("onchange","selectAllCellLines(this)")
									}else{
										foldlist_inputAll.style.display="none"
									}
									foldlist_inputAll.setAttribute("id",pointlessID(0))
									foldlist_inputAll.setAttribute("class","cellLine_selectorALL")
									foldlist_inputAll.style.display="none"
								var p=document.createElement("p")
									
								
									p.innerHTML="Select<br>All"
								foldlist_inputAll_label.appendChild(p)
								foldlist_inputAll_label.appendChild(foldlist_inputAll);
							var p1=document.createElement("p");
								p1.innerHTML=title
							var p2=document.createElement("p");
								p2.setAttribute("class","totalCount_CellLines");
							foldList_label.appendChild(foldlist_inputAll_label);
							foldList_label.appendChild(p1);
							foldList_label.appendChild(p2);
						foldList_bg_title.appendChild(foldList_label);
						foldList_bg_title.appendChild(get_chevronsup());
						foldList_bg_title.appendChild(get_chevronsdown());
					foldList_bg.appendChild(foldList_bg_title);
					
					var foldList=document.createElement("div");
						foldList.setAttribute("class","foldList");
						
						
						//CellLinesFilter
						var CellLinesFilter=document.createElement("div");
							CellLinesFilter.setAttribute("class","CellLinesFilter");
							var CellLinesFilterTitle=document.createElement("div");
								CellLinesFilterTitle.setAttribute("class","CellLinesFilterTitle");
								var CellLinesFilterTitle_label=document.createElement("label")
									CellLinesFilterTitle_label.setAttribute("for",pointlessID());
									var CellLinesFilterSwitch=document.createElement("input")
										CellLinesFilterSwitch.setAttribute("type","checkbox")
										CellLinesFilterSwitch.setAttribute("id",pointlessID(0))
										CellLinesFilterSwitch.setAttribute("class","filterSwitch")
										if(dataBaseNames[d]=="DepMap"){
											CellLinesFilterSwitch.setAttribute("disabled",true)
										}
										CellLinesFilterSwitch.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
									var CellLinesFilterTitle_p=document.createElement("p");
										CellLinesFilterTitle_p.innerHTML="Apply Filter";
									
									CellLinesFilterTitle_label.appendChild(CellLinesFilterSwitch);
									CellLinesFilterTitle_label.appendChild(CellLinesFilterTitle_p);
								CellLinesFilterTitle.appendChild(CellLinesFilterTitle_label);
							var CellLinesFilterContent=document.createElement("div");
								CellLinesFilterContent.setAttribute("class","CellLinesFilterContent");
								CellLinesFilterContent.setAttribute("style","height:0px;");
								
								var MSI_cellLines=document.createElement("div")
									MSI_cellLines.setAttribute("class","subFilterContent MSI_cellLines");
									var MSI_checkbox=document.createElement("input");
										MSI_checkbox.setAttribute("type","checkbox")
										MSI_checkbox.setAttribute("class","subFilterCheckbox")
										MSI_checkbox.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
									MSI_cellLines.appendChild(MSI_checkbox);
									var ptitle=document.createElement("p");
										ptitle.innerHTML="By MSI Status:";
									MSI_cellLines.appendChild(ptitle);
									var MSSbox=document.createElement("div");
										var checkbox1=document.createElement("input");
											checkbox1.setAttribute("type","radio");
											checkbox1.setAttribute("checked","checked");
											checkbox1.setAttribute("name",pointlessID());
											checkbox1.setAttribute("value","MSS");
											checkbox1.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
										var p4=document.createElement("p");
											p4.innerHTML="MSS";
									MSSbox.appendChild(checkbox1);
									MSSbox.appendChild(p4);
									var MSIbox=document.createElement("div");
										var checkbox2=document.createElement("input");
											checkbox2.setAttribute("type","radio");
											checkbox2.setAttribute("name",pointlessID(0));
											checkbox2.setAttribute("value","MSI");
											checkbox2.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
										var p5=document.createElement("p");
											p5.innerHTML="MSI";
									MSIbox.appendChild(checkbox2);
									MSIbox.appendChild(p5);
									MSI_cellLines.appendChild(MSSbox);
									MSI_cellLines.appendChild(MSIbox);
								CellLinesFilterContent.appendChild(MSI_cellLines);
								
								var ploidy_cellLines=document.createElement("div")
									ploidy_cellLines.setAttribute("class","subFilterContent ploidy_cellLines");
									var ploidy_checkbox=document.createElement("input");
										ploidy_checkbox.setAttribute("type","checkbox")
										ploidy_checkbox.setAttribute("class","subFilterCheckbox")
										ploidy_checkbox.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
									ploidy_cellLines.appendChild(ploidy_checkbox);
									var ptitle=document.createElement("p");
										ptitle.innerHTML="By ploidy:";
									var AbvOrBelw_sp=document.createElement('select');
										AbvOrBelw_sp.setAttribute("class","AbvOrBelw_sp")
										AbvOrBelw_sp.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
										var option0=document.createElement("option");
											option0.value="Above";
											option0.innerHTML="Above";
										var option1=document.createElement("option");
											option1.value="Below";
											option1.innerHTML="Below";
										AbvOrBelw_sp.appendChild(option0);
										AbvOrBelw_sp.appendChild(option1);
									var ploidy_input=document.createElement("input");
										ploidy_input.setAttribute("placeholder","e.g. 5%");
										ploidy_input.setAttribute("class","subFilterPercInput");
										ploidy_input.setAttribute("value","5%");
										ploidy_input.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
									ploidy_cellLines.appendChild(ptitle);
									ploidy_cellLines.appendChild(AbvOrBelw_sp);
									ploidy_cellLines.appendChild(ploidy_input);
								CellLinesFilterContent.appendChild(ploidy_cellLines);
								
								var mutation_cellLines=document.createElement("div")
									mutation_cellLines.setAttribute("class","subFilterContent mutation_cellLines");
									var MT_checkbox=document.createElement("input");
										MT_checkbox.setAttribute("type","checkbox")
										MT_checkbox.setAttribute("class","subFilterCheckbox")
										MT_checkbox.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
									mutation_cellLines.appendChild(MT_checkbox);
									var ptitle=document.createElement("p");
										ptitle.innerHTML="By mutation burden:";
									var AbvOrBelw_mt=document.createElement('select');
										AbvOrBelw_mt.setAttribute("class","AbvOrBelw_mt")
										AbvOrBelw_mt.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
										var option0=document.createElement("option");
											option0.value="Above";
											option0.innerHTML="Above";
										var option1=document.createElement("option");
											option1.value="Below";
											option1.innerHTML="Below";
										AbvOrBelw_mt.appendChild(option0);
										AbvOrBelw_mt.appendChild(option1);
									var mutation_input=document.createElement("input");
										mutation_input.setAttribute("placeholder","e.g. 5%");
										mutation_input.setAttribute("value","5%");
										mutation_input.setAttribute("class","subFilterPercInput");
										mutation_input.setAttribute("onchange","initList_foldList_cellLines(currentModel);checkStartButton(currentModel)")
									mutation_cellLines.appendChild(ptitle);
									mutation_cellLines.appendChild(AbvOrBelw_mt);
									mutation_cellLines.appendChild(mutation_input);
								CellLinesFilterContent.appendChild(mutation_cellLines);
								
							CellLinesFilter.appendChild(CellLinesFilterTitle);
							CellLinesFilter.appendChild(CellLinesFilterContent);
						foldList.appendChild(CellLinesFilter);
						

						var foldList_cellLines=document.createElement("div");
							foldList_cellLines.setAttribute("class","foldList_cellLines");
		

							
							
				
						foldList.appendChild(foldList_cellLines);
					foldList_bg.appendChild(foldList);
				sidebar.appendChild(foldList_bg);
				initList_foldList_cellLines(dataBaseNames[d]+","+ModelNames[i],isSingle=isSingle);
				countFoldList_cellLines(dataBaseNames[d]+","+ModelNames[i]);
				
				
			}
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("clusterOptions")){
		
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						foldList_folderInput.setAttribute("style","display:none;");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var p1=document.createElement("p");
								p1.innerHTML="Options:";
							foldList_label.appendChild(p1);
						foldList_bg_title.appendChild(foldList_label)
					var foldList=document.createElement("div");
						foldList.setAttribute("class","foldList");
						foldList.style.maxHeight="110px"
						foldList.style.flexShrink=0;
							var groupNum_div=document.createElement("div");
								groupNum_div.setAttribute("class","selecterInFoldList")
								var groupNum_p=document.createElement("p");
									groupNum_p.innerHTML="Number of Groups:"
								var groupNum_select=document.createElement("select");
									groupNum_select.setAttribute("class","groupNum_select")
									var option0=document.createElement("option");
										option0.value="2";
										option0.innerHTML="2";
									var option1=document.createElement("option");
										option1.value="3";
										option1.innerHTML="3";
									var option2=document.createElement("option");
										option2.value="4";
										option2.innerHTML="4";
									groupNum_select.appendChild(option0);
									groupNum_select.appendChild(option1);
									groupNum_select.appendChild(option2);
									groupNum_select.selectedIndex=1;
								groupNum_div.appendChild(groupNum_p);
								groupNum_div.appendChild(groupNum_select);
						foldList.appendChild(groupNum_div);
							var clusterAlgorithm_div=document.createElement("div");
								clusterAlgorithm_div.setAttribute("class","selecterInFoldList")
								var clusterAlgorithm_p=document.createElement("p");
									clusterAlgorithm_p.innerHTML="Cluster Algorithm:"
								var clusterAlgorithm_select=document.createElement("select");
									clusterAlgorithm_select.setAttribute("class","clusterAlgorithm_select")
									var option0=document.createElement("option");
										option0.value="km";
										option0.innerHTML="k-means";
									var option1=document.createElement("option");
										option1.value="hc";
										option1.innerHTML="hclust";
									var option2=document.createElement("option");
										option2.value="pam";
										option2.innerHTML="paritioning around medoids";
									clusterAlgorithm_select.appendChild(option0);
									clusterAlgorithm_select.appendChild(option1);
									clusterAlgorithm_select.appendChild(option2);
								clusterAlgorithm_div.appendChild(clusterAlgorithm_p);
								clusterAlgorithm_div.appendChild(clusterAlgorithm_select);
						foldList.appendChild(clusterAlgorithm_div);
					foldList_bg.appendChild(foldList_bg_title)
					foldList_bg.appendChild(foldList)
				sidebar.appendChild(foldList_bg)
			}
			
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("MachineLearning")){
		
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						foldList_folderInput.setAttribute("style","display:none;");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var p1=document.createElement("p");
								p1.innerHTML="Options:";
							foldList_label.appendChild(p1);
						foldList_bg_title.appendChild(foldList_label)
					var foldList=document.createElement("div");
						foldList.setAttribute("class","foldList");
						foldList.style.maxHeight="60px"
						foldList.style.flexShrink=0;
							var Algorithm_div=document.createElement("div");
								Algorithm_div.setAttribute("class","selecterInFoldList")
								var Algorithm_p=document.createElement("p");
									Algorithm_p.innerHTML="Deep Learning Algorithm Options:"
								var Algorithm_select=document.createElement("select");
									Algorithm_select.setAttribute("class","Algorithm_select")
									var option0=document.createElement("option");
										option0.value="Logistic";
										option0.innerHTML="Logistic Regression";
									var option1=document.createElement("option");
										option1.value="Lasso";
										option1.innerHTML="Lasso Regression";
									var option2=document.createElement("option");
										option2.value="Ridge";
										option2.innerHTML="Ridge Regression";
									var option3=document.createElement("option");
										option3.value="KNN";
										option3.innerHTML="k-Nearest Neighbor";
									var option4=document.createElement("option");
										option4.value="DecisionTree";
										option4.innerHTML="Decision Tree";
									var option5=document.createElement("option");
										option5.value="SVM";
										option5.innerHTML="Support vector machines";
									var option6=document.createElement("option");
										option6.value="GBM";
										option6.innerHTML="Gradient Boosting Machine";
									var option7=document.createElement("option");
										option7.value="Adaboost";
										option7.innerHTML="Adaboost";
									var option8=document.createElement("option");
										option8.value="NNet";
										option8.innerHTML="Neural networks";
									var option9=document.createElement("option");
										option9.value="randomForest";
										option9.innerHTML="Random Forest";
									Algorithm_select.appendChild(option0);
									Algorithm_select.appendChild(option1);
									Algorithm_select.appendChild(option2);
									Algorithm_select.appendChild(option3);
									Algorithm_select.appendChild(option4);
									Algorithm_select.appendChild(option5);
									Algorithm_select.appendChild(option6);
									Algorithm_select.appendChild(option7);
									Algorithm_select.appendChild(option8);
									Algorithm_select.appendChild(option9);
									Algorithm_select.selectedIndex=0;
								Algorithm_div.appendChild(Algorithm_p);
								Algorithm_div.appendChild(Algorithm_select);
						foldList.appendChild(Algorithm_div);
					foldList_bg.appendChild(foldList_bg_title)
					foldList_bg.appendChild(foldList)
				sidebar.appendChild(foldList_bg)
			}
			
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("drugMethod")){
		
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var foldList_bg=document.createElement("div");
					foldList_bg.setAttribute("class","foldList_bg");
					var foldList_folderInput=document.createElement("input");
						foldList_folderInput.setAttribute("type","checkbox");
						
						foldList_folderInput.setAttribute("class","foldList_folderInput");
						foldList_folderInput.setAttribute("id",pointlessID());
						foldList_folderInput.setAttribute("style","display:none;");
						foldList_folderInput.setAttribute("onclick","changeFoldListCheckbox(this,currentModel)");
						if(isFirstFoldList){
							isFirstFoldList=false;
							foldList_folderInput.checked=true;
						}
					foldList_bg.appendChild(foldList_folderInput);
					var foldList_bg_title=document.createElement("div");
						foldList_bg_title.setAttribute("class","foldList_bg_title");
						var foldList_label=document.createElement("label");
							foldList_label.setAttribute("class","foldList_label");
							foldList_label.setAttribute("for",pointlessID(0));
							var p1=document.createElement("p");
								p1.innerHTML="More options:";
							foldList_label.appendChild(p1);
						foldList_bg_title.appendChild(foldList_label)
					var foldList=document.createElement("div");
						foldList.setAttribute("class","foldList");
						foldList.style.maxHeight="150px"
						var methodBox=document.createElement("div");
							methodBox.setAttribute("class","methodBox");
							var cor_method=document.createElement("div")
								cor_method.setAttribute("class","cor_method")
								var cor_method_title=document.createElement("p");
									cor_method_title.innerHTML="Statistical  Method:";
								var cor_method_op1=document.createElement("input");
									cor_method_op1.setAttribute("type","radio");
									cor_method_op1.setAttribute("checked","true");
									cor_method_op1.setAttribute("value","spearman");
									cor_method_op1.setAttribute("name",pointlessID());
								var cor_method_p1=document.createElement("p");
									cor_method_p1.innerHTML="Spearman"
								var cor_method_op2=document.createElement("input");
									cor_method_op2.setAttribute("type","radio");
									cor_method_op2.setAttribute("value","pearson");
									cor_method_op2.setAttribute("name",pointlessID(0));
								var cor_method_p2=document.createElement("p");
									cor_method_p2.innerHTML="Pearson"
								cor_method.appendChild(cor_method_title)
								cor_method.appendChild(cor_method_op1)
								cor_method.appendChild(cor_method_p1)
								cor_method.appendChild(cor_method_op2)
								cor_method.appendChild(cor_method_p2)
							var mutation_method=document.createElement("div")
								mutation_method.setAttribute("class","mutation_method")
								
								var mutation_method_title=document.createElement("p");
									mutation_method_title.innerHTML="Statistical  Method:";
								var mutation_method_op1=document.createElement("input");
									mutation_method_op1.setAttribute("type","radio");
									mutation_method_op1.setAttribute("checked","true");
									mutation_method_op1.setAttribute("value","wilcoxon");
									mutation_method_op1.setAttribute("name",pointlessID());
								var mutation_method_p1=document.createElement("p");
									mutation_method_p1.innerHTML="Wilcoxon"
								var mutation_method_op2=document.createElement("input");
									mutation_method_op2.setAttribute("type","radio");
									mutation_method_op2.setAttribute("value","t.test");
									mutation_method_op2.setAttribute("name",pointlessID(0));
								var mutation_method_p2=document.createElement("p");
									mutation_method_p2.innerHTML="t.test"
								mutation_method.appendChild(mutation_method_title)
								mutation_method.appendChild(mutation_method_op1)
								mutation_method.appendChild(mutation_method_p1)
								mutation_method.appendChild(mutation_method_op2)
								mutation_method.appendChild(mutation_method_p2)
							if(["Omics-Omics (cis-regulation)","Omics-Omics (trans-regulation)"].indexOf(ModelNames[i])==-1){
								mutation_method.style.display="none"
							}else{
								cor_method_title.innerHTML="Statistical  Method 1:";
								mutation_method_title.innerHTML="Statistical  Method 2:";
							}
							var pValue_div=document.createElement("div")
								pValue_div.setAttribute("class","MHT")
								var p=document.createElement("p")
									p.innerHTML="Multiple Testing:"
								var op1=document.createElement("input")
									op1.setAttribute("class","FDR")
									op1.setAttribute("type","radio");
									op1.setAttribute("checked","true");
									op1.setAttribute("value","FDR");
									op1.setAttribute("name",pointlessID());
								var p1=document.createElement("p")
									p1.innerHTML="FDR"
								var op2=document.createElement("input")
									op2.setAttribute("type","radio");
									op2.setAttribute("value","P");
									op2.setAttribute("name",pointlessID(0));
								var p2=document.createElement("p")
									p2.innerHTML="P"
								pValue_div.appendChild(p)
								pValue_div.appendChild(op1)
								pValue_div.appendChild(p1)
								pValue_div.appendChild(op2)
								pValue_div.appendChild(p2)
							var pCufoff_div=document.createElement("div")
								var p=document.createElement("p")
									p.innerHTML="Cut-off: r/Rho ="
								var inputbox=document.createElement("input")
									inputbox.placeholder="e.g. 0.2"
									inputbox.value="0.2"
									inputbox.style.width="80px"
									inputbox.setAttribute("oninput","isPvalue(this)")
									inputbox.setAttribute("class","pcutoff")
								pCufoff_div.appendChild(p)
								pCufoff_div.appendChild(inputbox)
							methodBox.appendChild(cor_method)
							methodBox.appendChild(mutation_method)
							methodBox.appendChild(pValue_div)
							methodBox.appendChild(pCufoff_div)
						foldList.appendChild(methodBox)
							//cor_method="spearman",mutation_method="wilcoxon"
					foldList_bg.appendChild(foldList_bg_title)
					foldList_bg.appendChild(foldList)
				sidebar.appendChild(foldList_bg)
			}
			
			if(true){//设置气球
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var balloon=document.createElement("div");
					balloon.setAttribute("class","balloon")
					balloon.style.height="50vh"
					balloon.style.flexShrink=99999
					
				sidebar.appendChild(balloon);
			}
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("plotButton")){
				var sidebar=Model.getElementsByClassName("sidebar")[0];
				var startplot=document.createElement("div");
					startplot.setAttribute("class","startplot");
					startplot.setAttribute("style","background:grey;cursor:default;");
					startplot.innerHTML='<h3>START</h1>';
				sidebar.appendChild(startplot);
			}
			
			
			if(initParams.ModelInfo[dataBaseNames[d]][ModelNames[i]].includes("plotArea")){
				
				var plotArea=document.createElement("div");
					plotArea.setAttribute("class","plotArea");
					var toolbar=document.createElement("div")
						toolbar.setAttribute("class","toolbar")
						var helpInfoBox=document.createElement("div")
							helpInfoBox.setAttribute("class","helpInfoBox")
							helpInfoBox.setAttribute("onmouseenter","expandHelp(this)")
							helpInfoBox.setAttribute("onmouseleave","shrinkHelp(this)")
							var helpTitle_p=document.createElement("p")
								
								
								helpTitle_p.innerHTML=ModelNames[i]
							var helpInfoTitle=document.createElement("div")
								helpInfoTitle.setAttribute("class","helpInfoTitle")
								helpInfoTitle.appendChild(helpTitle_p)
								var helpTitle_svg=SVG(helpInfoTitle)
									helpTitle_svg.attr({"viewBox":"0 0 30 30"})
									helpTitle_svg.circle(22).attr({
										cx:15,
										cy:15,
										stroke:"none",
										fill:"white"
									})
									helpTitle_svg.plain("?").attr({
										x:15,
										y:15,
										"text-anchor":"middle",
										"dominant-baseline":  "central",
										"font-size":16,
										"fill":"#aaa"
									})
							var helpInfoContent=document.createElement("div")
								helpInfoContent.setAttribute("class","helpInfoContent")
							var p=document.createElement("p")
								p.innerHTML=initParams.ModelHelp[ModelNames[i]]
								helpInfoContent.appendChild(p)
								
							helpInfoBox.appendChild(helpInfoTitle)
							helpInfoBox.appendChild(helpInfoContent)
						var tweakbar=document.createElement("div");
							tweakbar.setAttribute("class","tweakbar");
							var zoombox=document.createElement("div");
								zoombox.setAttribute("class","zoombox");
								var shrink=document.createElement("div");
									shrink.setAttribute("onclick","shrinkSVG()")
									shrink.innerHTML="-"
								var zoomindex=document.createElement("div")
									zoomindex.innerHTML="100%";
								var enlarge=document.createElement("div");
									enlarge.setAttribute("onclick","enlargeSVG()")
									enlarge.innerHTML="+"
								var init=document.createElement("div");
									init.setAttribute("onclick","initSVG()")
									init.innerHTML="◎"
							zoombox.appendChild(shrink);
							zoombox.appendChild(zoomindex);
							zoombox.appendChild(enlarge);
							zoombox.appendChild(init);
							var typeSelector=document.createElement("select");
								typeSelector.setAttribute("onchange","downloadFig(this)");
								var option0=document.createElement("option");
									option0.value="title";
									option0.innerHTML="Download Figure:";
								var option1=document.createElement("option");
									option1.value=".eps";
									option1.innerHTML="*.eps";
								var option2=document.createElement("option");
									option2.value=".pdf";
									option2.innerHTML="*.pdf";
								var option3=document.createElement("option");
									option3.value=".tiff";
									option3.innerHTML="*.tiff";
								var option4=document.createElement("option");
									option4.value=".svg";
									option4.innerHTML="*.svg";
								typeSelector.appendChild(option0);
								typeSelector.appendChild(option1);
								typeSelector.appendChild(option2);
								typeSelector.appendChild(option3);
								typeSelector.appendChild(option4);
								typeSelector.selectedIndex=0;	
							tweakbar.appendChild(zoombox)
							tweakbar.appendChild(typeSelector)
						toolbar.appendChild(helpInfoBox)
						toolbar.appendChild(tweakbar)
					var svgBackground=document.createElement("div")
						svgBackground.setAttribute("class","svgBackground");
					var svgContainer=document.createElement("div")
						svgContainer.setAttribute("class","svgContainer");
					plotArea.appendChild(toolbar);
					plotArea.appendChild(svgBackground);
					plotArea.appendChild(svgContainer);
				Model.appendChild(plotArea);

				plotted=false;
				shrinkHelp(helpInfoBox)
			}
			
			
			checkStartButton(dataBaseNames[d]+","+ModelNames[i]);
		}
	}
}
function initList_foldList_drugs(Model_ID,isSingle=false){
	var Model=document.getElementById(Model_ID);
	var dataBaseName=Model_ID.split(",")[0]
	module="foldList_drugs";
	if(initParams[module][dataBaseName]==undefined)return
	if(isSingle){
		var checkSubAllAttr='style="display:none"'
		var checkSubAll_labelAttr='style="border-right:none;width:10px;"'
		var checkboxType="radio"
		var countTagNum=""
		var countTagCol='style="color:white"'
	}else{
		
		var checkSubAllAttr='onchange=selectAllSubDrugs(this)'
		var checkSubAll_labelAttr=''
		var checkboxType="checkbox"
		var countTagNum="0/"
		var countTagCol=''
	}
	drugListHtml='';
	//console.log(initParams[module][dataBaseName].length)
	var pathways= Object.keys(initParams[module][dataBaseName]);
	var radioGroup=pointlessID()
	for(var i=0;i< pathways.length;i++){
		var drugs= Object.keys(initParams[module][dataBaseName][pathways[i]].drugs);
		drugListHtml+=
			'<div><input type="checkbox" id="'+
			pointlessID()+
			'" class="drug_folder_input" checked><div class="foldlist_folder_comp"><label '+checkSubAll_labelAttr+' for="'+
			pointlessID()+
			'"><input type="checkbox" class="selectorSubALL" '+checkSubAllAttr+' id="'+
			pointlessID(0)+
			'" name="'+initParams[module][dataBaseName][pathways[i]].typename+
			'"></label><label for="'+
			pointlessID(-1)+
			'"><p>'+
			initParams[module][dataBaseName][pathways[i]].typename+
			'</p><p '+countTagCol+' class="countTag"> ('+countTagNum+drugs.length+')</p></label>'+
			'<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10" viewBox="0 0 10 10"  class="chevrons-up">'+
			'<polyline points="1 7 5 3 9 7" color="white" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"></polyline>'+
			'</svg>'+
			'<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10" viewBox="0 0 10 10"  class="chevrons-down">'+
			'<polyline points="1 3 5 7 9 3" color="white" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"></polyline>'+
			'</svg>'+		
			'</div>'+
			'<div class="foldlist_items" style=height:'+drugs.length*31+'px;>'
		
		for(var j=0;j< drugs.length;j++){
			//userParams[currentModel][module][initParams[module][dataBaseName][pathways[i]].children[j]]=false;
			drugListHtml+=
				'<div><input type="'+checkboxType+'" onchange=checkOneDrug(this) id="'+ pointlessID()+
				'" name='+radioGroup+' drugname="'+initParams[module][dataBaseName][pathways[i]].drugs[drugs[j]].drugName+
				'" class="drug_selector"><label for="'+
				pointlessID(0)+
				'">'+
				'<p class="drugTags_L">'+
				initParams[module][dataBaseName][pathways[i]].drugs[drugs[j]].drugName+
				'</p></label></div>'
		}
		drugListHtml+='</div></div>'
	}
	//drugListHtml="";
	var foldList=Model.getElementsByClassName(module)[0];
	foldList.innerHTML=drugListHtml;
	//checkStartButton();
	var drugTags_L=foldList.getElementsByClassName("drugTags_L")
	for(var i=0;i<drugTags_L.length;i++){
		drugTags_L[i].addEventListener("mouseenter",showInfo,false)
		drugTags_L[i].addEventListener("mouseleave",closeInfo,false)
	}
	
}




function initList_foldList_cellLines(Model_ID,isSingle=false){
	
	var Model=document.getElementById(Model_ID);
	module="foldList_cellLines"
	var foldList=Model.getElementsByClassName(module)[0];
	if(foldList==undefined)return;
	if(isSingle){
		var checkSubAllAttr='style="display:none"'
		var checkSubAll_labelAttr='style="border-right:none;width:10px;"'
		var checkboxType="radio"
		var countTagNum=""
		var countTagCol='style="color:white"'
	}else{
		
		var checkSubAllAttr='onchange=selectAllSubCellLines(this)'
		var checkSubAll_labelAttr=''
		var checkboxType="checkbox"
		var countTagNum="0/"
		var countTagCol=''
	}
	applyFilter(Model_ID)
	var CancerListHtml=""

	var cancerTypes= Object.keys(dealingCellLines[Model_ID]);
	var radioGroup=pointlessID()
	for(var i=0;i< cancerTypes.length;i++){
		var cellLineCount=0;
		var cellLineNames=Object.keys(dealingCellLines[Model_ID][cancerTypes[i]].cellLines)
		
		
		
		var cellListHtml=""
		var cellCount=0
		for(var j=0;j< cellLineNames.length;j++){
			if(dealingCellLines[Model_ID][cancerTypes[i]].cellLines[cellLineNames[j]].enabled){
				cellCount++
				cellListHtml+=
					'<div><input type="'+checkboxType+'" onchange=checkOneCancer(this) id="'+ pointlessID()+
					'" name='+radioGroup+' cellname="'+cellLineNames[j]+
					'" class="cellLine_selector"><label for="'+
					pointlessID(0)+
					'">'+
					'<p class="CellAndDrugResistTags" cellLineName='+cellLineNames[j]+'>'+
					cellLineNames[j]+
					'</p></label></div>'
			}
		}
		if(cellCount==0)continue
		CancerListHtml+=
			'<div><input type="checkbox" id="'+
			pointlessID()+
			'" class="drug_folder_input" checked><div class="foldlist_folder_comp"><label '+checkSubAll_labelAttr+' for="'+
			pointlessID()+
			'"><input type="checkbox" class="selectorSubALL" '+checkSubAllAttr+' id="'+
			pointlessID(0)+
			'" name="'+cancerTypes[i]+
			'"></label><label for="'+
			pointlessID(-1)+
			'"><p>'+
			cancerTypes[i]+
			'</p><p '+countTagCol+' class="countTag"> ('+countTagNum+cellCount+')</p></label>'+
			'<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10" viewBox="0 0 10 10"  class="chevrons-up">'+
			'<polyline points="1 7 5 3 9 7" color="white" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"></polyline>'+
			'</svg>'+
			'<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10" viewBox="0 0 10 10"  class="chevrons-down">'+
			'<polyline points="1 3 5 7 9 3" color="white" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"></polyline>'+
			'</svg>'+		
			'</div>'+
			'<div class="foldlist_items" style=height:'+cellCount.length*31+'px;>'
			
		CancerListHtml+=cellListHtml+'</div></div>'

	}
	
	
	foldList.innerHTML=CancerListHtml;
	Model.getElementsByClassName("cellLine_selectorALL")[0].checked=false;
	Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate=false;
	var CellAndDrugResistTags=foldList.getElementsByClassName("CellAndDrugResistTags")
	for(var i=0;i<CellAndDrugResistTags.length;i++){
		CellAndDrugResistTags[i].addEventListener("mouseenter",showInfo,false)
		CellAndDrugResistTags[i].addEventListener("mouseleave",closeInfo,false)
	}

}


function inputAndSearch_Search_onfocus(obj){
	inputAndSearch_Search_isOn=true
}
function inputAndSearch_Search_onblur(obj){
	inputAndSearch_Search_isOn=false
	setTimeout(function(){
		if(!(inputAndSearch_input_isOn|inputAndSearch_Search_isOn)){
			var inputAndSearch_Search=obj.parentNode.getElementsByClassName("inputAndSearch_Search")[0];
			inputAndSearch_Search.style.display="none";
			checkStartButton(currentModel);
		}
	},200)
}


function isPercentage(str){
   var patt1 = new RegExp(/^(100|[1-9]?\d(\.\d\d\d\d\d\d\d?)?)%$/);
 
   return(patt1.test(str));

}

function isOnInfoBox(event){

	if(event.type=="mouseleave"){
		isOutBox=true;
		if(isOutElement){
			popupInfoBox.style.opacity=0
			popupInfoBox.style.pointerEvents ="none";
	
		}
	}else if(event.type=="mouseenter"){
		isOutBox=false;
	}
}

function isOnPromotBox(event){

	if(event.type=="mouseleave"){
		isOutPromotBox=true;
		if(!isOnInput){
			var promotBox=document.getElementById("promotBox");
			promotBox.style.display="none";
		}
	}else if(event.type=="mouseenter"){
		isOutPromotBox=false;
	}
}

function isPvalue(obj){
	console.log(obj.value)
	obj.value=obj.value.replace(/^\D*(\d*(?:\.\d{0,2})?).*$/g, '$1')
	if(parseFloat(obj.value)>1){
		obj.value=1
	}
	if(isNaN(parseFloat(obj.value))){
		obj.style.background="red"
	}else{
		obj.style.background=""
	}
	checkStartButton(currentModel)
}

function moveToFirst(obj){
	
	var orgRank=parseInt(obj.getAttribute("rank"))
	var OmicResultBox=obj.parentNode
	var OmicResultItems=OmicResultBox.getElementsByClassName("OmicResultItem")
	for(var o=0;o<OmicResultItems.length;o++){
		var currRank=parseInt(OmicResultItems[o].getAttribute("rank"))
		if(orgRank>currRank){
			OmicResultItems[o].setAttribute("rank",currRank+1)
		}
		
	}
	obj.setAttribute("rank",1)
	
	reArrangeOmics(OmicResultBox)
}

function omicsChanged(Model_ID){
	var Model=document.getElementById(Model_ID);
	var omicsSelector=Model.getElementsByClassName("omicsSelector")[0];
	var genename=Model.getElementsByClassName("genename")[0];
	var genenameTag=Model.getElementsByClassName("genenameTag")[0];
	if(omicsSelector.value=="methylation_site"){
		genename.value="cg00009944"
		genename.placeholder="Methylation site name here";
		genenameTag.innerHTML="Please input a methylation site:"
	}else if(omicsSelector.value=="microRNA"){
		genename.value="hsa-let-7b"
		genename.placeholder="microRNA name here";
		genenameTag.innerHTML="Please input a microRNA symbol:"
	}else if(omicsSelector.value=="Chromatin"){
		genename.value="H3K4me0"
		genename.placeholder="Histone modification here";
		genenameTag.innerHTML="Please input a histone modification:"
	}else if(omicsSelector.value=="Metabolomics"){
		genename.value="2-aminoadipate"
		genename.placeholder="Metabolite here";
		genenameTag.innerHTML="Please input a metabolite:"
	}else if(omicsSelector.value=="Metastatic"){
		genename.value="all"
		genename.placeholder="Metastatic site here";
		genenameTag.innerHTML="Please input a metastatic site:"
	}else{
		genenameTag.innerHTML="Please input a gene symbol:"
		genename.value="TP53"
		genename.placeholder="Gene symbol here";
	}
	if(omicsSelector.value=="mutation"){
		var pcutoff=Model.getElementsByClassName("pcutoff")[0];
		if(pcutoff!=undefined){
			pcutoff.parentNode.style.display="none"
		}
		
	}else{
		var pcutoff=Model.getElementsByClassName("pcutoff")[0];
		if(pcutoff!=undefined){
			pcutoff.parentNode.style.display=""
		}
	}
	genename.style.background="white"
	/*
	if(Model_ID=="DepMap,Cancer Types Summary"){
		var individualCheckBox=Model.getElementsByClassName("individualCheckBox")[0]
		if(omicsSelector.value=="Metastatic"){
			individualCheckBox.style.display=""
		}else{
			individualCheckBox.style.display="none"
		}
	}
	*/
	var cor_method=Model.getElementsByClassName("cor_method")[0];
	if(cor_method!=undefined){
		var mutation_method=Model.getElementsByClassName("mutation_method")[0];
		if(["CNV","mutation"].indexOf(omicsSelector.value)>-1){
			cor_method.style.display="none"
			mutation_method.style.display=""
		}else{
			cor_method.style.display=""
			mutation_method.style.display="none"
		}
	}
	//closePromot(genename)
}



function pointlessID(type){
	if(!window.poID){
		poID=0;
	};
	if(type!=undefined){
		return("pointlessID"+(poID+type));
	}else{		
		poID++;
		return("pointlessID"+poID);
	}
}

function PrefixInteger(num, n) {
    return (Array(n).join(0) + num).slice(-n);
}

function removeOmics(obj){
	event.cancelBubble=true;
	var orgRank=parseInt(obj.parentNode.parentNode.parentNode.getAttribute("rank"))
	
	var OmicResultBox=obj.parentNode.parentNode.parentNode.parentNode
	obj.parentNode.parentNode.parentNode.remove()
	var OmicResultItems=OmicResultBox.getElementsByClassName("OmicResultItem")
	for(var o=0;o<OmicResultItems.length;o++){
		var currRank=parseInt(OmicResultItems[o].getAttribute("rank"))
		if(orgRank<currRank){
			OmicResultItems[o].setAttribute("rank",currRank-1)
		}
		
	}
	reArrangeOmics(OmicResultBox)
	checkStartButton(currentModel)
}

function reArrangeOmics(OmicResultBox){

	var OmicResultItems=OmicResultBox.getElementsByClassName("OmicResultItem")
	for(var o=0;o<OmicResultItems.length;o++){
		var currRank=parseInt(OmicResultItems[o].getAttribute("rank"))
		if(currRank==1){
			OmicResultItems[o].style.backgroundImage="linear-gradient(#bc4e4e,#6c2a2a)";
		}else{
			OmicResultItems[o].style.backgroundImage="";
		}
		OmicResultItems[o].style.top=50+42*(currRank-1)+"px"
	}
	OmicResultBox.style.height=(OmicResultItems.length+1)*42+30+"px"
}

function showPromot(obj){
	isOnInput=true
	var x=obj.getBoundingClientRect().x;
	var y=obj.getBoundingClientRect().y;
	var targetWidth=obj.getBoundingClientRect().width;
	var targetHeight=obj.getBoundingClientRect().height;
	var promotBox=document.getElementById("promotBox");
	promotBox.innerHTML=""
	promotBox.style.minWidth=targetWidth+"px";
	promotBox.style.left=x+"px";
	promotBox.style.top=y+targetHeight+"px";
	promotBox.style.display="flex";
	updatePromot(obj)
}

function showSearchBox(obj){
	inputAndSearch_Search_isOn=false
	var inputAndSearch_Search=obj.parentNode.parentNode.getElementsByClassName("inputAndSearch_Search")[0];
	inputAndSearch_Search.style.display="";
	inputAndSearch_input_isOn=true
}


function switchModel(obj){
	var fromModel=currentModel;
	var toModel=obj.getAttribute("name");
	if(toModel=="NATIVE,HELP")initHelp()
	if(toModel==fromModel)return;
	var allModelNames=[];
	var allModels=document.getElementsByClassName("Model");
	for(var m =0;m<allModels.length;m++){
		allModelNames.push(allModels[m].id);
	}
	
	
	if(allModelNames.indexOf(fromModel)<allModelNames.indexOf(toModel)){
		startLeft=document.body.clientWidth;
	}else{
		startLeft=-document.body.clientWidth;
	}
	var fromModel_obj=document.getElementById(fromModel);
	//fromModel_obj.style.opacity=0;
	fromModel_obj.style.display="none";
	//fromModel_obj.style.zIndex="1";
	//toModel=currentModel
	var toModel_obj=document.getElementById(toModel);
	toModel_obj.setAttribute("style","transition:left 0s ease;left:"+startLeft+"px;");
	toModel_obj.style.display="";
	setTimeout(function(){
		toModel_obj.removeAttribute("style");
	//toModel_obj.style.opacity=1
	
	toModel_obj.style.left=0
	//toModel_obj.style.zIndex="2";
	},1)
	currentModel=toModel;
	var dataBaseName=currentModel.split(',')[0]
	var ModelName=currentModel.split(',')[1]
	if(["NATIVE,HOME","NATIVE,DOWNLOAD","NATIVE,HELP"].indexOf(currentModel)==-1){
		if(initParams.ModelInfo[dataBaseName][ModelName].indexOf("plotArea")){
			plotSvgBackground();
		}
	}
	var dataBaseTip=document.getElementsByClassName("dataBaseTip")
	for(var d=0;d<dataBaseTip.length;d++){
		dataBaseTip[d].style.borderTop=""
		var labels=dataBaseTip[d].getElementsByTagName("label")
		for(var l=0;l<labels.length;l++){
			if(["NATIVE,HOME","NATIVE,DOWNLOAD","NATIVE,HELP"].indexOf(toModel)==-1){
				labels[l].style.background=""
			}
		}
	}
	obj.parentNode.style.borderTop="5px solid #ff8888"
	if(["NATIVE,HOME","NATIVE,DOWNLOAD","NATIVE,HELP"].indexOf(toModel)==-1){
		obj.style.background="#ccc"
	}
	checkStartButton(toModel);

}


function showInfo(event){
	isOutElement=false;
	isOutBox=true;
	var dataBaseName=currentModel.split(",")[0]
	
	var popupTarget=event.target
	var className=popupTarget.getAttribute("class");
	if(className=="drugTags"|className=="drugTags_L"){
		if(className=="drugTags_L"){
			var indicator="left";
		}else{
			var indicator="right";
		}
		var drugName=popupTarget.innerHTML;
		var pathwayNames=Object.keys(initParams.foldList_drugs[dataBaseName])
		for(var i=0;i<pathwayNames.length;i++){
			var drugNames=Object.keys(initParams.foldList_drugs[dataBaseName][pathwayNames[i]].drugs)
			for(var j=0;j<drugNames.length;j++){
				if(drugNames[j]==drugName){
					var class_name=initParams.foldList_drugs[dataBaseName][pathwayNames[i]].drugs[drugName].class
					pubchem=initParams.foldList_drugs[dataBaseName][pathwayNames[i]].drugs[drugName].pubchem
					if(pubchem!="-"){
						var pubchems=pubchem.split(", ")
						
						for(var p=0;p<pubchems.length;p++){
							pubchems[p]='<a target="view_window" href=https://pubchem.ncbi.nlm.nih.gov/compound/'+pubchems[p]+'>'+pubchems[p]+'</a>'
						}
						var pubchem_linkComp=pubchems.join(",")

					}else{
						var pubchem_linkComp=pubchem
					}
					var synonyms=initParams.foldList_drugs[dataBaseName][pathwayNames[i]].drugs[drugName].synonyms
					var targets=initParams.foldList_drugs[dataBaseName][pathwayNames[i]].drugs[drugName].targets
					break;
				}
			}
		}
		popupInfoBoxContent_tmp.innerHTML='<table>'+
								'<tr><td><p>Drug name: </p></td><td><p>'+drugName+'</p></td></tr>'+
								'<tr><td><p>Class name: </p></td><td><p>'+class_name+'</p></td></tr>'+
								'<tr><td><p>pubchem: </p></td><td><p>'+pubchem_linkComp+'</p></td></tr>'+
								'<tr><td><p>Synonyms: </p></td><td><p>'+synonyms+'</p></td></tr>'+
								'<tr><td><p>Targets: </p></td><td><p>'+targets+'</p></td></tr>'+
								'</table>'
		
		
		var x=event.target.getBoundingClientRect().x;
		var y=event.target.getBoundingClientRect().y;
		var targetWidth=event.target.getBoundingClientRect().width;
	}else if(className=="genename"){

		var indicator="left";
	
		var statusText=popupTarget.getAttribute("statusText");
		var statusColor=popupTarget.getAttribute("statusColor");
		popupInfoBoxContent_tmp.innerHTML='<p color="'+statusColor+'">'+statusText+'</p>'
		
		
		var x=event.target.getBoundingClientRect().x;
		var y=event.target.getBoundingClientRect().y;
		var targetWidth=event.target.getBoundingClientRect().width;
		
	}else if(className=="cellLinesTag"){
		var cancerTypeName=popupTarget.parentNode.children[0].getAttribute("name");
		var cancerTypes=Object.keys(dealingCellLines[currentModel])
		//var str='<table><tr><th><p>Cell Name</p></th><th><p>Tissue</p></th><th><p>cancer_type_detail</p></th><th><p>gender</p></th><th><p>age_at_sampling</p></th><th><p>Growth Properties</p></th><th><p>smoking_status</p></th><th><p>MSI_Status</p></th><th><p>cancer_type_ncit_id</p></th><th><p>ploidyPerc</p></th><th><p>mutational_burden</p></th></tr>'
		for(var i=0;i<cancerTypes.length;i++){
			if(cancerTypes[i]==cancerTypeName){
				var cellLineNames=Object.keys(dealingCellLines[currentModel][cancerTypes[i]].cellLines)
				var cellLineNames_enabled=[];
				for(var j=0;j<cellLineNames.length;j++){
						if(dealingCellLines[currentModel][cancerTypes[i]].cellLines[cellLineNames[j]].enabled==true){
							var model_id=dealingCellLines[currentModel][cancerTypes[i]].cellLines[cellLineNames[j]].model_id
							var DepMap_ID=dealingCellLines[currentModel][cancerTypes[i]].cellLines[cellLineNames[j]].DepMap_ID
							if(DepMap_ID!=undefined){
								cellLineNames_enabled.push('<a target="view_window" href=https://depmap.org/portal/cell_line/'+DepMap_ID+'>'+cellLineNames[j]+'</a>')
							}else{
								cellLineNames_enabled.push('<a target="view_window" href=https://cellmodelpassports.sanger.ac.uk/passports/'+model_id+'>'+cellLineNames[j]+'</a>')
							}
						}
						//str+='<tr><th><p>'+cellLineNames[j]+'</p></th><th><p>'+Tissue+'</p></th><th><p>'+cancer_type_detail+'</p></th><th><p>'+gender+'</p></th><th><p>'+age_at_sampling+'</p></th><th><p>'+Growth_Properties+'</p></th><th><p>'+smoking_status+'</p></th><th><p>'+MSI_Status+'</p></th><th><p>'+cancer_type_ncit_id+'</p></th><th><p>'+ploidyPerc+'</p></th><th><p>'+mutational_burden+'</p></th></tr>'
				}
				break;
			}
		}
		//str+='</table>'
		

		
		
		popupInfoBoxContent_tmp.innerHTML='<table>'+
								'<tr><td><p>Cancer Type: </p></td><td><p>'+cancerTypeName+'</p></td></tr>'+
								'<tr><td><p>Included Cell Lines: </p></td><td><p>'+cellLineNames_enabled.join(", ")+'</p></td></tr>'+
								'</table>'
		var indicator="left";
		var x=getOffsetLeft(event.target);
		var scrolltop=getScrollTop(event.target.parentNode.parentNode.parentNode.parentNode)
		var y=getOffsetTop(event.target)-scrolltop;
		var targetWidth=event.target.clientWidth
	}else if(className=="CellAndDrugResistTags"){
		var cellLineName=popupTarget.getAttribute("cellLineName");
		var dataBaseName=currentModel.split(",")[0]
		var currCellLines=initParams.cellLines[dataBaseName]
		var cancerTypes=Object.keys(currCellLines)
		//var str='<table><tr><th><p>Cell Name</p></th><th><p>Tissue</p></th><th><p>cancer_type_detail</p></th><th><p>gender</p></th><th><p>age_at_sampling</p></th><th><p>Growth Properties</p></th><th><p>smoking_status</p></th><th><p>MSI_Status</p></th><th><p>cancer_type_ncit_id</p></th><th><p>ploidyPerc</p></th><th><p>mutational_burden</p></th></tr>'
		for(var i=0;i<cancerTypes.length;i++){
			
			var cellLineNames=Object.keys(currCellLines[cancerTypes[i]].cellLines)
			for(var j=0;j<cellLineNames.length;j++){
				if(currCellLines[cancerTypes[i]].cellLines[cellLineNames[j]].Cell_Name==cellLineName){
					var cellInfos=currCellLines[cancerTypes[i]].cellLines[cellLineNames[j]]

					break;
				}
					
			}
		}
		//console.log("cell:"+i+","+j+","+cellLineNames[j]+","+cellLineName)
		//str+='</table>'
		if(popupTarget.getAttribute("AdditionInfo")!=undefined){
			var drugResistance=JSON.parse(popupTarget.getAttribute("AdditionInfo"));
			var oriValue=popupTarget.getAttribute("oriValue");
			var drugResistanceHTML='<tr><td><p>Response to <span style="font-style:italic;">'+drugResistance[0]+'</span>:</p></td><td><p>'+drugResistance[1]+ '</p></td></tr>'+
				'<tr><td><p>IC50: </p></td><td><p>'+Math.exp(oriValue).toFixed(2)+' uM</p></td></tr>'
		}else{
			var drugResistanceHTML=""
		}
		var notForShow=["Cell_Name","MSI_Status","ploidyPerc","mutational_burdenPerc","CNV","methylation","mRNA","mutation","model_id","microRNA","protein","CRISPR","Metabolomics","Metastatic","RNAi","Fusions","Chromatin","hasDrug","enabled","DepMap_ID"]
		html_tmp='<table><tr><td><p>Cell line: </p></td><td><p>'+cellLineName+'</p></td></tr>'
		if(cellInfos.model_id!=""){
			html_tmp+='<tr><td><p>Cell Model Link: </p></td><td><p><a target="view_window" href=https://cellmodelpassports.sanger.ac.uk/passports/'+cellInfos.model_id+'>'+cellInfos.model_id+'</a></p></td></tr>'
		}
		if(cellInfos.DepMap_ID!=undefined&cellInfos.DepMap_ID!=""&cellInfos.DepMap_ID!="NA"){
			html_tmp+='<tr><td><p>DepMap Link: </p></td><td><p><a target="view_window" href=https://depmap.org/portal/cell_line/'+cellInfos.DepMap_ID+'>'+cellInfos.DepMap_ID+'</a></p></td></tr>'
		}
		
		
		var infoNames=Object.keys(cellInfos);
		for(var i=0;i<infoNames.length;i++){
			if(notForShow.indexOf(infoNames[i])>-1)continue
			if(cellInfos[infoNames[i]]==""|cellInfos[infoNames[i]]=="Unknown"|cellInfos[infoNames[i]]=="NA")continue
			
			html_tmp+='<tr><td><p>'+infoNames[i]+': </p></td><td><p>'+cellInfos[infoNames[i]]+'</p></td></tr>'
			
		}
		html_tmp+=drugResistanceHTML+'</table>'
		popupInfoBoxContent_tmp.innerHTML=html_tmp
		/*
		popupInfoBoxContent_tmp.innerHTML='<table>'+
								'<tr><td><p>Cell line: </p></td><td><p><a target="view_window" href=https://cellmodelpassports.sanger.ac.uk/passports/'+model_id+'>'+cellLineName+'</a></p></td></tr>'+
								'<tr><td><p>Tissue: </p></td><td><p>'+Tissue+'</p></td></tr>'+
								'<tr><td><p>Cancer type detail: </p></td><td><p>'+cancer_type_detail+'</p></td></tr>'+
								'<tr><td><p>Gender: </p></td><td><p>'+gender+'</p></td></tr>'+
								'<tr><td><p>Age at sampling: </p></td><td><p>'+age_at_sampling+'</p></td></tr>'+
								'<tr><td><p>Growth Properties: </p></td><td><p>'+Growth_Properties+'</p></td></tr>'+
								'<tr><td><p>Smoking status: </p></td><td><p>'+smoking_status+'</p></td></tr>'+
								'<tr><td><p>MSI Status: </p></td><td><p>'+MSI_Status+'</p></td></tr>'+
								'<tr><td><p>NCIT id: </p></td><td><p><a target="view_window" href=https://ncit.nci.nih.gov/ncitbrowser/ConceptReport.jsp?dictionary=NCI_Thesaurus&code='+cancer_type_ncit_id+'>'+cancer_type_ncit_id+'</a></p></td></tr>'+
								'<tr><td><p>Ploidy: </p></td><td><p>'+ploidy+'</p></td></tr>'+
								'<tr><td><p>Mutational burden: </p></td><td><p>'+mutational_burden+'</p></td></tr>'+
								drugResistanceHTML+
								'</table>'
		*/
		
		var indicator="left";
		var x=event.target.getBoundingClientRect().x;
		var y=event.target.getBoundingClientRect().y;
		var targetWidth=event.target.getBoundingClientRect().width;
	}else if(className=="MethSiteTags"){
		var MethSiteName=popupTarget.innerHTML;
		var MethSiteAdditionInfo=JSON.parse(popupTarget.getAttribute("AdditionInfo"));
		var MethSiteAdditionInfo_item=Object.keys(MethSiteAdditionInfo)
		var str='<table><tr><td><p>Methylation Site: </p></td><td><p>'+MethSiteName+'</p></td></tr>'
		for(var i =0;i<MethSiteAdditionInfo_item.length;i++){
			str+='<tr><td><p>'+MethSiteAdditionInfo_item[i]+': </p></td><td><p>'+MethSiteAdditionInfo[MethSiteAdditionInfo_item[i]]+'</p></td></tr>'
			
		}
		str+='</table>'
		popupInfoBoxContent_tmp.innerHTML=str
		var indicator="right";
		var x=event.target.getBoundingClientRect().x;
		var y=event.target.getBoundingClientRect().y;
		var targetWidth=event.target.getBoundingClientRect().width;
	}else if(className=="startplot"){
		
		var str='<p style="font-weight:bold;">Things to check: <p><table>'
		for(var i =0;i<disableReasons.length;i++){
			str+='<tr style="color:red;"><td><p>'+(i+1)+':</p></td><td><p>'+disableReasons[i]+'</p></td></tr>'
			
		}
		str+='</table>'
		popupInfoBoxContent_tmp.innerHTML=str
		var indicator="left";
		var x=event.target.getBoundingClientRect().x;
		var y=event.target.getBoundingClientRect().y;
		var targetWidth=event.target.getBoundingClientRect().width;
		
	}else if(className=="summBarplotInfo"){
		var str=popupTarget.getAttribute("AdditionInfo");
		popupInfoBoxContent_tmp.innerHTML=str
		var indicator="right";
		var x=event.target.getBoundingClientRect().x;
		var y=event.target.getBoundingClientRect().y;
		var targetWidth=event.target.getBoundingClientRect().width;
	}else if(className=="warnings"){
		var str=popupTarget.getAttribute("AdditionInfo");
		var warnings_json=JSON.parse(str)
		var str_content=['<p style="font-size:12px;font-weight:bold"> Warning(s):</p>']
		for(var w=0;w<warnings_json.length;w++){
			str_content.push('<p style="color:red;font-size:12px;text-indent:2em">'+(w+1)+". "+warnings_json[w]+'</p>')
		}
		popupInfoBoxContent_tmp.innerHTML=str_content.join("")
		var indicator="left";
		var x=event.target.getBoundingClientRect().x;
		var y=event.target.getBoundingClientRect().y;
		var targetWidth=event.target.getBoundingClientRect().width;
	}else{
		return;
	}
		
	
	//popupInfoBoxContent.style.transform="scale(0.5)"
	popupInfoBoxContent.innerHTML=popupInfoBoxContent_tmp.innerHTML
	
	var Margin=10;
	popupInfoBoxContent_tmp.style.width="";
	var boxWidth=popupInfoBoxContent_tmp.clientWidth;
	popupInfoBoxContent_tmp.style.width=boxWidth+"px";
	var boxHeight=popupInfoBoxContent_tmp.clientHeight+Margin*3;
	
	popupInfoBoxContent.style.margin=Margin+"px"
	popupInfoBox.style.width=boxWidth+2*Margin+"px";
	popupInfoBox.style.height=boxHeight+"px";
	popupInfoBoxBg.innerHTML="";
	popupInfoBoxBg.style.width=boxWidth+2*Margin+"px";
	popupInfoBoxBg.style.height=boxHeight+"px";
	popupInfoBoxBg.style.marginBottom=-boxHeight+"px";
	popupInfoBoxContent.style.margin=Margin+"px";
	popupInfoBoxContent.style.width=boxWidth+"px";
	popupInfoBoxContent.style.height=boxHeight-Margin*3+"px";
	var drawing=SVG(popupInfoBoxBg);
	drawing.attr({"viewBox":"0 0 "+(boxWidth+2*Margin)+" "+boxHeight,"style":"background:transparent"})
	if(indicator=="right"){
		drawing.polygon(
			Margin/5+","+Margin/2+" "+
			(boxWidth+2*Margin-Margin/5)+","+Margin/2+" "+
			(boxWidth+2*Margin-Margin/5)+","+(boxHeight-Margin/2-Margin)+" "+
			(boxWidth+2*Margin-Margin/5-Margin/2)+","+(boxHeight-Margin/2-Margin)+" "+
			(boxWidth+2*Margin-Margin/5)+","+(boxHeight-Margin/5)+" "+
			(boxWidth+2*Margin-Margin/5-2*Margin)+","+(boxHeight-Margin/2-Margin)+" "+
			Margin/5+","+(boxHeight-Margin/2-Margin)
		).attr({
			fill:"white",
			stroke:"black",
			"stroke-width":0.5
		})
			

		
		popupInfoBox.style.left=x-boxWidth-2*Margin+"px";
		popupInfoBox.style.top=y-boxHeight+"px";
	}else if(indicator=="left"){
		drawing.polygon(
			Margin/5+","+Margin/2+" "+
			(boxWidth+2*Margin-Margin/5)+","+Margin/2+" "+
			(boxWidth+2*Margin-Margin/5)+","+(boxHeight-Margin/2-Margin)+" "+
			(Margin/5+2*Margin)+","+(boxHeight-Margin/2-Margin)+" "+
			(Margin/5)+","+(boxHeight-Margin/5)+" "+
			(Margin/5+Margin/2)+","+(boxHeight-Margin/2-Margin)+" "+
			Margin/5+","+(boxHeight-Margin/2-Margin)
		).attr({
			fill:"white",
			stroke:"black",
			"stroke-width":0.5
		})


		
		popupInfoBox.style.left=x+targetWidth+"px";
		popupInfoBox.style.top=y-boxHeight+"px";
	}
	popupInfoBox.style.opacity=1;
	popupInfoBox.style.pointerEvents ="all";
}

function showInfo_mani(popupTarget){
	isOutElement=true;
	isOutBox=true;
	
	

	var className=popupTarget.getAttribute("class");

	if(className.indexOf("geneSetsTextarea")>-1){

		var indicator="left";
	
		var statusText=popupTarget.getAttribute("statusText");
		var statusColor=popupTarget.getAttribute("statusColor");
		popupInfoBoxContent_tmp.innerHTML=statusText
		
		
		var x=popupTarget.getBoundingClientRect().x;
		var y=popupTarget.getBoundingClientRect().y;
		var targetWidth=popupTarget.getBoundingClientRect().width;

	}
	if(className.indexOf("addButton")>-1){

		var indicator="left";
	
		var statusText=popupTarget.getAttribute("errorInfo");

		popupInfoBoxContent_tmp.innerHTML=statusText
		
		
		var x=popupTarget.getBoundingClientRect().x;
		var y=popupTarget.getBoundingClientRect().y;
		var targetWidth=popupTarget.getBoundingClientRect().width;

	}
	if(className.indexOf("OmicResultItem")>-1){

		var indicator="left";
	
		var statusText=popupTarget.getAttribute("errorInfo");

		popupInfoBoxContent_tmp.innerHTML=statusText
		
		
		var x=popupTarget.getBoundingClientRect().x;
		var y=popupTarget.getBoundingClientRect().y;
		var targetWidth=popupTarget.getBoundingClientRect().width;

	}
	popupInfoBoxContent.innerHTML=popupInfoBoxContent_tmp.innerHTML
	
	var Margin=10;
	popupInfoBoxContent_tmp.style.width="";
	var boxWidth=popupInfoBoxContent_tmp.clientWidth;
	popupInfoBoxContent_tmp.style.width=boxWidth+"px";
	var boxHeight=popupInfoBoxContent_tmp.clientHeight+Margin*3;
	
	popupInfoBoxContent.style.margin=Margin+"px"
	popupInfoBox.style.width=boxWidth+2*Margin+"px";
	popupInfoBox.style.height=boxHeight+"px";
	popupInfoBoxBg.innerHTML="";
	popupInfoBoxBg.style.width=boxWidth+2*Margin+"px";
	popupInfoBoxBg.style.height=boxHeight+"px";
	popupInfoBoxBg.style.marginBottom=-boxHeight+"px";
	popupInfoBoxContent.style.margin=Margin+"px";
	popupInfoBoxContent.style.width=boxWidth+"px";
	popupInfoBoxContent.style.height=boxHeight-Margin*3+"px";
	var drawing=SVG(popupInfoBoxBg);
	drawing.attr({"viewBox":"0 0 "+(boxWidth+2*Margin)+" "+boxHeight,"style":"background:transparent"})
	if(indicator=="right"){
		drawing.polygon(
			Margin/5+","+Margin/5+" "+
			(boxWidth+2*Margin-Margin/5)+","+Margin/2+" "+
			(boxWidth+2*Margin-Margin/5)+","+(boxHeight-Margin/2-Margin)+" "+
			(boxWidth+2*Margin-Margin/5-Margin/2)+","+(boxHeight-Margin/2-Margin)+" "+
			(boxWidth+2*Margin-Margin/5)+","+(boxHeight-Margin/5)+" "+
			(boxWidth+2*Margin-Margin/5-2*Margin)+","+(boxHeight-Margin/2-Margin)+" "+
			Margin/5+","+(boxHeight-Margin/2-Margin)
		).attr({
			fill:"white",
			stroke:"black",
			"stroke-width":0.5
		})
			

		
		popupInfoBox.style.left=x-boxWidth-2*Margin+"px";
		popupInfoBox.style.top=y-boxHeight+"px";
	}else if(indicator=="left"){
		drawing.polygon(
			Margin/5+","+Margin/5+" "+
			(boxWidth+2*Margin-Margin/5)+","+Margin/2+" "+
			(boxWidth+2*Margin-Margin/5)+","+(boxHeight-Margin/2-Margin)+" "+
			(Margin/5+2*Margin)+","+(boxHeight-Margin/2-Margin)+" "+
			(Margin/5)+","+(boxHeight-Margin/5)+" "+
			(Margin/5+Margin/2)+","+(boxHeight-Margin/2-Margin)+" "+
			Margin/5+","+(boxHeight-Margin/2-Margin)
		).attr({
			fill:"white",
			stroke:"black",
			"stroke-width":0.5
		})


		
		popupInfoBox.style.left=x+targetWidth+"px";
		popupInfoBox.style.top=y-boxHeight+"px";
	}
	popupInfoBox.style.opacity=1;
	popupInfoBox.style.pointerEvents ="all";
	setTimeout(function(){
		if(isOutElement&isOutBox){
			popupInfoBox.style.opacity=0;
			popupInfoBox.style.pointerEvents ="none";
		}
	},3000)
}

function startplot(){
	var Model_ID=currentModel
	var dataBaseName=Model_ID.split(',')[0]
	var ModelName=Model_ID.split(',')[1]
	var Model=document.getElementById(Model_ID);
	
	var plottingStatus="plotting";
	userParams[Model_ID]=new Object();
	
	var plotButton=Model.getElementsByClassName("startplot")[0];
	plotButton.style.animationName="myfirst";
	plotButton.setAttribute("onclick","");
	plotButton.style.cursor="default";
	//var DebugTips=document.getElementById("DebugTips");
	//DebugTips.style.background="";
	
	var Children_names=initParams.ModelInfo[dataBaseName][ModelName];

	
	userParams[Model_ID].Model=Model_ID;
	userParams[Model_ID].use="jsShow";
	var requestGUID=guid();
	userParams[Model_ID].requestGUID=requestGUID;
	var plotGUID=guid();
	userParams[Model_ID].plotGUID=plotGUID;
	userParams[Model_ID].path="../plotCache/data_jsShow_"+plotGUID+".json";
	if(Children_names.includes("foldList_drugs")|Children_names.includes("foldList_drugs_single")){
		userParams[Model_ID].foldList_drugs=new Object();
		userParams[Model_ID].foldList_drugs=countFoldList_drugs(Model_ID)[0];
	}
	if(Children_names.includes("foldList_cellLines")|Children_names.includes("foldList_cellLines_single")){
		userParams[Model_ID].foldList_cellLines=new Object();
		userParams[Model_ID].foldList_cellLines=countFoldList_cellLines(Model_ID)[0]
		
	}
	
	if(Children_names.includes("geneSetSelector")){
		userParams[Model_ID].foldList_geneSets=new Object();
		if(Children_names.includes("geneSetSelector")){
			if(geneSets_Textarea){
				userParams[Model_ID].foldList_geneSets.type="customGeneList"
				userParams[Model_ID].foldList_geneSets.genes=Model.getElementsByClassName("geneSetsTextarea")[0].value.split(",")
			}else if(geneSets_Preset){
				userParams[Model_ID].foldList_geneSets.type="MSigDB"
				userParams[Model_ID].foldList_geneSets.MSigDB_Name=Model.getElementsByClassName("inputAndSearch_Input_input")[0].value
			}
		}
		
	}
	
	
	if(Children_names.includes("omicsSelectorAndMultigenes")){
		userParams[Model_ID].omics=Model.getElementsByClassName("omicsSelectorAndMultigenes")[0].value;
		userParams[Model_ID].genes=Model.getElementsByClassName("geneSetsTextarea")[0].value.split(",");
	}
 
	if(Children_names.includes("clusterOptions")){
		userParams[Model_ID].groupNum=Model.getElementsByClassName("groupNum_select")[0].value;
		userParams[Model_ID].clusterAlgorithm=Model.getElementsByClassName("clusterAlgorithm_select")[0].value;
	}
	
	if(Children_names.includes("MachineLearning")){
		userParams[Model_ID].Algorithm=Model.getElementsByClassName("Algorithm_select")[0].value;
	}
	if(Children_names.includes("drugMethod")){
		var cor_method=Model.getElementsByClassName("cor_method")[0]
			var input=cor_method.getElementsByTagName("input")
			for(var i=0;i<input.length;i++){
				if(input[i].checked){
					userParams[Model_ID].cor_method=input[i].value
				}
			}
		var mutation_method=Model.getElementsByClassName("mutation_method")[0]
			var input=mutation_method.getElementsByTagName("input")
			for(var i=0;i<input.length;i++){
				if(input[i].checked){
					userParams[Model_ID].mutation_method=input[i].value
				}
			}
		var MHT=Model.getElementsByClassName("MHT")[0]
			var input=MHT.getElementsByTagName("input")
			for(var i=0;i<input.length;i++){
				if(input[i].checked){
					userParams[Model_ID].MHT=input[i].value
				}
			}
		var pcutoff=Model.getElementsByClassName("pcutoff")[0]
			userParams[Model_ID].pcutoff=pcutoff.value

	}
	if(Children_names.includes("omicsSelector")){
		userParams[Model_ID].omicsSelector=Model.getElementsByClassName("omicsSelector")[0].value;
		userParams[Model_ID].singleGeneInput=Model.getElementsByClassName("genename")[0].value;
		if(userParams[Model_ID].omicsSelector=="CNV"){
			userParams[Model_ID].colors=new Object();
			userParams[Model_ID].colors.DEL_color = "#000091";
			userParams[Model_ID].colors.LOSS_color = "#4343b5"
			userParams[Model_ID].colors.NOR_color = "#8e8e8e"
			userParams[Model_ID].colors.GAIN_color = "#ce5555"
			userParams[Model_ID].colors.AMP_color = "#ba0000"
		}
	}
	
	if(Children_names.includes("onlyGene")){
		userParams[Model_ID].singleGeneInput=Model.getElementsByClassName("genename")[0].value;
	}
	
	if(Children_names.includes("multiOmicsSelector")){
		var OmicResultItem=Model.getElementsByClassName("OmicResultItem")
		userParams[Model_ID].multiOmicsSelector=new Object();
		userParams[Model_ID].multiOmicsSelector.Omics=[];
		userParams[Model_ID].multiOmicsSelector.Genes=[];
		userParams[Model_ID].multiOmicsSelector.Rank=[];
		for(var o=0;o<OmicResultItem.length;o++){
			userParams[Model_ID].multiOmicsSelector.Omics.push(OmicResultItem[o].getAttribute("Omic"));
			userParams[Model_ID].multiOmicsSelector.Genes.push(OmicResultItem[o].getAttribute("gene"));
			userParams[Model_ID].multiOmicsSelector.Rank.push(OmicResultItem[o].getAttribute("rank"));
		}
	}
	
	jsDebug=[];
	xmlhttp=new XMLHttpRequest();
	//var Debug=document.getElementById("NATIVE,HELP");
	xmlhttp.onreadystatechange=function(){
		//document.getElementById("plotarea").innerHTML =genename+"got:"+xmlhttp.responseText;
		if (xmlhttp.readyState==4 && xmlhttp.status==200){
			
			
			Params_str=xmlhttp.responseText;
			
			//console.log(Params_str)
			//document.getElementById("svgContainer").innerHTML=null;
			jsDebug["jsDrawingStart"]=new Date().getTime()
			
			try {
				Params[Model_ID]=JSON.parse(Params_str);
			} catch (error) {
				//Debug.innerHTML=Params_str+'<br><a href="/CLIPA/plotCache/userParams_'+requestGUID+'.json" target="_blank" style="height:20px;">下载请求参数</a>'
				
				
				plotButton.style.animationName="";
				plotButton.setAttribute("onclick","startplot()")
				plotButton.style.cursor="";
				return;
			}
			if(Params[Model_ID].requestGUID!=requestGUID){
				console.log("requestGUID not match")
				plotButton.style.animationName="";
				plotButton.setAttribute("onclick","startplot()")
				plotButton.style.cursor="";
				return
			}
			
			
			if(Params[Model_ID].error){
				//Debug.innerHTML=new Date().toLocaleTimeString()+":<br>"+Params[Model_ID].error+'<br><a href="/CLIPA/plotCache/userParams_'+requestGUID+'.json" target="_blank" style="height:20px;">下载请求参数</a>'
				//DebugTips.style.background="red";
				plotButton.style.animationName="";
				plotButton.setAttribute("onclick","startplot()")
				plotButton.style.cursor="";
				return
			}
			if(Params[Model_ID].res.error){
				plotButton.style.animationName="";
				plotButton.setAttribute("onclick","startplot()")
				plotButton.style.cursor="";
				alert(Params[Model_ID].res.error);
				return;
			}
			plotButton.style.animationName="";
			plotButton.setAttribute("onclick","startplot()")
			plotButton.style.cursor="";
			createSVG(Params,Model_ID)
			
			var plottingStatus="done";
			jsDebug["TotalDelay"]=(new Date().getTime()-jsDebug["TotalStart"])/1000;
			jsDebug["jsDrawingDelay"]=(new Date().getTime()-jsDebug["jsDrawingStart"])/1000;
			jsDebug["transferDelay"]=(jsDebug["TotalDelay"]-jsDebug["jsDrawingDelay"]-jsDebug["DrawingDelay"]);
			//Debug.innerHTML=jsDebug+'<br><a href="/CLIPA/plotCache/userParams_'+requestGUID+'.json" target="_blank" style="height:20px;">下载请求参数</a>'
			document.body.style.background="white"
	
			//DebugTips.style.background="green";
		}else{
			
		}
			
	}
	xmlhttp.onerror=function(e){
		plotButton.style.animationName="";
		plotButton.setAttribute("onclick","startplot()")
		plotButton.style.cursor="";
		//Debug.innerHTML=new Date().toLocaleTimeString()+":<br>"+e
		//DebugTips.style.background="red";
	}
	//document.getElementById("plotarea").innerHTML= +"start";
	userParams_uri="userParams="+btoa(JSON.stringify(userParams[Model_ID]));
	//console.log(userParams[Model_ID]_uri)
	xmlhttp.open("POST","php/sendToR.php",true);
	xmlhttp.setRequestHeader("Content-Type","application/x-www-form-urlencoded");
	jsDebug["TotalStart"]=new Date().getTime();

	xmlhttp.send(userParams_uri);
	recordvisit("recordPlot.php")
		


}


function selectAllSubDrugs(obj){
	targetStatus=obj.checked;
	obj.indeterminate =false;
	var drugs=obj.parentNode.parentNode.parentNode.getElementsByClassName("drug_selector");
	for(var i=0;i<drugs.length;i++){
		drugs[i].checked=targetStatus;
	}
	var item_tag=obj.parentNode.parentNode.getElementsByClassName("countTag")[0]
	item_tag.innerHTML="("+drugs.length*targetStatus+"/"+drugs.length+")"
	if(targetStatus){
		item_tag.setAttribute("style","color:white;")
	}else{
		item_tag.setAttribute("style","color:grey;")
	}
	//countFoldList_drugs(currentModel);
	var Model=document.getElementById(currentModel)
	var count_drugs=countFoldList_drugs(currentModel);
	if(count_drugs[0].length==0){
		console.log("set0")
		Model.getElementsByClassName("drug_selectorALL")[0].checked=false;
		Model.getElementsByClassName("drug_selectorALL")[0].indeterminate =false;
	}else if(count_drugs[0].length==count_drugs[1]){

		Model.getElementsByClassName("drug_selectorALL")[0].checked=true;
		Model.getElementsByClassName("drug_selectorALL")[0].indeterminate =false;
	}else{

		Model.getElementsByClassName("drug_selectorALL")[0].indeterminate =true;
	}
	checkStartButton(currentModel);
}

function selectAllSubCellLines(obj){
	targetStatus=obj.checked;
	obj.indeterminate =false;
	var cellLines=obj.parentNode.parentNode.parentNode.getElementsByClassName("cellLine_selector");
	for(var i=0;i<cellLines.length;i++){
		cellLines[i].checked=targetStatus;
	}
	var item_tag=obj.parentNode.parentNode.getElementsByClassName("countTag")[0]
	item_tag.innerHTML="("+cellLines.length*targetStatus+"/"+cellLines.length+")"
	if(targetStatus){
		item_tag.setAttribute("style","color:white;")
	}else{
		item_tag.setAttribute("style","color:grey;")
	}
	//countFoldList_drugs(currentModel);
	var Model=document.getElementById(currentModel)
	var count_cellLines=countFoldList_cellLines(currentModel);
	if(count_cellLines[0].length==0){
		Model.getElementsByClassName("cellLine_selectorALL")[0].checked=false;
		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =false;
	}else if(count_cellLines[0].length==count_cellLines[1]){

		Model.getElementsByClassName("cellLine_selectorALL")[0].checked=true;
		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =false;
	}else{

		Model.getElementsByClassName("cellLine_selectorALL")[0].indeterminate =true;
	}
	checkStartButton(currentModel);
}

function selectAllDrugs(obj){
	var targetStatus=obj.checked;

	obj.indeterminate =false;
	var drugs=obj.parentNode.parentNode.parentNode.parentNode.getElementsByClassName("drug_selector");
	for(var i=0;i<drugs.length;i++){
		drugs[i].checked=targetStatus;
	}
	var selectorSubALL=obj.parentNode.parentNode.parentNode.parentNode.getElementsByClassName("selectorSubALL")
	for(var i=0;i<selectorSubALL.length;i++){
		selectorSubALL[i].checked=targetStatus;
		selectorSubALL[i].indeterminate=false;
	}
	//countFoldList_drugs(currentModel);
	
	checkStartButton(currentModel);
}

function selectAllCellLines(obj){
	var targetStatus=obj.checked;
	
	obj.indeterminate =false;
	var Model=document.getElementById(currentModel)
	var cellLine_selector=Model.getElementsByClassName("cellLine_selector");
	for(var i=0;i<cellLine_selector.length;i++){
		cellLine_selector[i].checked=targetStatus;
	}
	var selectorSubALL=obj.parentNode.parentNode.parentNode.parentNode.getElementsByClassName("selectorSubALL")
	for(var i=0;i<selectorSubALL.length;i++){
		selectorSubALL[i].checked=targetStatus;
		selectorSubALL[i].indeterminate=false;
	}
	//countFoldList_drugs(currentModel);
	
	checkStartButton(currentModel);
}

function shrinkModelTips(obj){
	obj.setAttribute("expand",false);

	for(var i=1;i<obj.children.length;i++){
		obj.children[i].style.width="120px"
		
	} 
	setTimeout(function(){
		if(obj.getAttribute("expand")=="false"){
			for(var i=1;i<obj.children.length;i++){

				obj.children[i].style.top="0px";
				obj.children[i].style.boxShadow=""
			}
		}
	},200)
	obj.style.height=40+"px"
}
function shrinkHelp(obj){
	helpInfoTitle=obj.getElementsByClassName("helpInfoTitle")[0]
	var helpTitle_width=helpInfoTitle.getBoundingClientRect().width
	var helpTitle_height=helpInfoTitle.getBoundingClientRect().height
	var helpInfoContent=obj.getElementsByClassName("helpInfoContent")[0]
	var helpContent_width=helpInfoContent.getBoundingClientRect().width
	var helpContent_height=helpInfoContent.getBoundingClientRect().height
	obj.style.width=helpTitle_width+10+"px";
	obj.style.height="100%"
	obj.style.boxShadow=""
	
}
function updatePromot(obj){
	var cellHeight=28
	var Model=document.getElementById(currentModel);
	var dataBaseName=currentModel.split(",")[0]
	var ModelName=currentModel.split(",")[1]
	var promotBox=document.getElementById("promotBox");
	if(vennList==undefined){
		promotBox.innerHTML='<div><p style="font-style:italic;color:grey;font-size:12px;margin:2px;height:20px;">"Gene list has not been prepared. Please retry!</p></div>'
		promotBox.style.height=cellHeight+'px'
	}else{
		if(obj.getAttribute("class").indexOf("genename")>-1){
			var omicsSelector=Model.getElementsByClassName("omicsSelector")[0];
			if(omicsSelector==undefined){
				
				if(ModelName=="Omics-Omics (cis-regulation)"){
					var vennTag=dataBaseName+"_cis"
				}else{
					var vennTag="all"
				}
			}else{
				var vennTag=currentModel.split(",")[0]+"_"+omicsSelector.value
			}
			getGeneList(vennTag)
			if(geneList[vennTag].length!=0){
				var currInput=obj.value;
				//obj.value=currInput;
				var divHeight=0;
				var promotInfo=[];
				for(var i=0;i<geneList[vennTag].length;i++){
					if(geneList[vennTag][i].toUpperCase().indexOf(currInput.toUpperCase())==0){
						
						promotInfo.push(geneList[vennTag][i])
						
					}
				}
				
				promotHtml_tmp="";
				if(promotInfo.length>0){
					for(var p=0;p<promotInfo.length;p++){
						promotHtml_tmp+='<div class="promotBoxItem" onclick=changeInput("'+obj.id+'","'+promotInfo[p]+'")><p>'+promotInfo[p].substr(0,promotInfo[p].toUpperCase().indexOf(currInput.toUpperCase()))+'<span style="font-weight:900">'+promotInfo[p].substr(promotInfo[p].toUpperCase().indexOf(currInput.toUpperCase()),promotInfo[p].toUpperCase().indexOf(currInput.toUpperCase())+currInput.length)+"</span>"+promotInfo[p].substr(promotInfo[p].toUpperCase().indexOf(currInput.toUpperCase())+currInput.length,promotInfo[p].length)+'</p></div>'
						divHeight+=cellHeight;
						if(p>8&p!=promotInfo.length-1){
							divHeight+=cellHeight;
							promotHtml_tmp+='<div style="font-style:italic;color:grey;font-size:12px;margin:2px;height:20px;"><p>'+(promotInfo.length-p-1)+' more...</p></div>'
							break;
						}
					}
					obj.style.background="white"
				}else{
					promotHtml_tmp+='<div><p style="font-style:italic;color:grey;margin:2px;height:20px;">"'+currInput+'" is invalid! Please check!</p></div>'
					divHeight+=cellHeight;
					obj.style.background="red"
				}

				
				promotBox.innerHTML=promotHtml_tmp;
				promotBox.style.height=divHeight+"px"
			
			}else{
				promotBox.style.display='none'
				obj.style.background="";
		
			}
		}
	}
	
}


function updateSearchBox(obj){
	inputAndSearch_Search_isOn=false
	var currInput=obj.value.toUpperCase();
	obj.value=currInput;
	


	var geneSets_Parent_Name= Object.keys(initParams["foldList_geneSets"]);
	newList=new Object();
	for(var i=0;i< geneSets_Parent_Name.length;i++){
		newList[geneSets_Parent_Name[i]]=[];
		var geneSets_Name= Object.keys(initParams["foldList_geneSets"][geneSets_Parent_Name[i]]);
		for(var j=0;j<geneSets_Name.length;j++){
			if(geneSets_Name[j].indexOf(currInput)>-1){
				newList[geneSets_Parent_Name[i]].push(geneSets_Name[j])
			}
		}
	}
	
	var inputAndSearch_Search=obj.parentNode.parentNode.getElementsByClassName("inputAndSearch_Search")[0];
	inputAndSearch_Search.innerHTML=""
	var newList_Parent_Names= Object.keys(newList);
	for(var i=0;i< newList_Parent_Names.length;i++){
		var geneSet_Names= newList[newList_Parent_Names[i]];
		if(geneSet_Names.length==0)continue;
		var div=document.createElement("div")
			var folderInput=document.createElement("input")
				folderInput.style.display="none";
				folderInput.setAttribute("id",pointlessID())
				folderInput.setAttribute("type","checkbox")
				folderInput.setAttribute("class","newList_fold_input")
			var label=document.createElement("label")
				label.setAttribute("for",pointlessID(0))
				label.appendChild(get_plus())
				label.appendChild(get_minus())
				label.innerHTML+="<p>"+newList_Parent_Names[i]+"</p>"
				label.setAttribute("class","Parent_title_comp")
			var newList_foldlist=document.createElement("div")
				newList_foldlist.setAttribute("class","newList_foldlist")
			
			for(var j=0;j< geneSet_Names.length;j++){
				var geneSet_Name=document.createElement("label")
					geneSet_Name.setAttribute("class","newList_foldlist_items")
					geneSet_Name.setAttribute("onclick",'changeInput("'+obj.getAttribute("id")+'","'+newList_Parent_Names[i].toUpperCase()+"_"+geneSet_Names[j]+'"),checkMSigDB(currentModel)')
					var boldedStr='<p>'+geneSet_Names[j].substr(0,geneSet_Names[j].indexOf(currInput))+'<span style="font-weight:900">'+
					currInput+'</span>'+geneSet_Names[j].substr(geneSet_Names[j].indexOf(currInput)+currInput.length,geneSet_Names[j].length)+'<span style="font-weight:900">('+initParams["foldList_geneSets"][newList_Parent_Names[i]][geneSet_Names[j]].num+' genes)</span></p>'
					geneSet_Name.innerHTML=boldedStr
				newList_foldlist.appendChild(geneSet_Name)
			}
			//newList_foldlist.style.height=geneSet_Names.length*20+"px"
			div.appendChild(folderInput)
			
			div.appendChild(label)
			div.appendChild(newList_foldlist)
		inputAndSearch_Search.appendChild(div)
	}
		

}


window.onresize = function(){
	//initModelsTipBg();
	plotSvgBackground();
	
}

function getFUNname (str){
	
   
    var lastIndex = str.search(/\(/);
    var funName = str.substring(9,lastIndex);
    return(funName);
}
function recordvisit(php){
	var modifyData_xmlhttp=new XMLHttpRequest();
	modifyData_xmlhttp.open("GET","php/"+php,true);
	modifyData_xmlhttp.send();
}