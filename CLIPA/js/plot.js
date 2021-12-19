


plotSvgBackground=function(){
	
	var Model=document.getElementById(currentModel);
	var dataBaseName=currentModel.split(',')[0]

	if(initParams.ModelInfo[dataBaseName]==undefined)return;
	var svgBackground=Model.getElementsByClassName("svgBackground")[0];
	vDPI=svgBackground.parentNode.clientWidth*0.9/8.25;
	
	var A4Width=8.25*vDPI;
	var A4Height=11.75*vDPI;
	A4margin=5;
	var svgBackgroundWidth=A4Width+A4margin*2;
	var svgBackgroundHeight=A4Height+A4margin*2;
	svgBackground.setAttribute("style","width:"+svgBackgroundWidth+"px;height:"+svgBackgroundHeight+"px;margin-bottom:"+(-svgBackgroundHeight-50)+"px;")
	
	svgBackground.innerHTML="";
	if (SVG.supported) {
	  var drawing=SVG(svgBackground);
	  drawing.viewbox(0,0,svgBackgroundWidth,svgBackgroundHeight)

	} else {
	  alert('svg.js not supported');
	  return;
	}
	var foldcorner=0.03
	
	if(isWarnings[currentModel]){
		var disp=""
	}else{
		var disp="none"
	}
	var triA=A4Width*0.03
	var tipsfontsize=0.02*Math.min((svgBackgroundWidth-A4margin*2),svgBackgroundHeight-A4margin*2);
	var warnings=drawing.group().attr({class:"warnings",style:'display:'+disp})
	var warnings_triA=drawing.polygon((A4margin+A4Width*0.02)+","+(A4margin+A4Height*0.03)+" "+(A4margin+A4Width*0.02+triA)+","+(A4margin+A4Height*0.03)+" "+(A4margin+A4Width*0.02+triA/2)+","+(A4margin+A4Height*0.03-triA*Math.cos(Math.PI/6))).attr({
		fill:"red",
		stroke:"none",
		
		"stroke-linejoin":"round"
		
	})
	var warnings_text=drawing.plain("!").attr({
		"font-size":tipsfontsize,
		x:(A4margin+A4Width*0.02+triA/2),
		y:(A4margin+A4Height*0.03)-triA*0.35,
		stroke:"none",
		fill:"white",
		"font-weight":"bold",
		"text-anchor":"middle" ,
		"dominant-baseline":"central",
		"pointer-events":"none"
	})
	warnings.add(warnings_triA)
	warnings.add(warnings_text)
	plotA4=false
	
	drawing.polygon(A4margin+","+A4margin+" "+A4margin+","+(A4margin+A4Height)+" "+(A4margin+A4Width)+","+(A4margin+A4Height)+" "+(A4margin+A4Width)+","+(A4margin+A4Width*foldcorner)+" "+(A4margin+A4Width*(1-foldcorner))+","+A4margin).attr({
		fill:"white",
		stroke:"black",
		"stroke-width":"1px"
		
	})
	drawing.polyline((A4margin+A4Width*(1-foldcorner))+","+A4margin+" "+(A4margin+A4Width*(1-foldcorner))+","+(A4margin+A4Width*foldcorner)+" "+(A4margin+A4Width)+","+(A4margin+A4Width*foldcorner)).attr({
		stroke:"black",
		fill:"#ddd",
		"stroke-width":"1px"
		
	})
	if(plotA4){
		drawing.plain("Size of A4: 8.25 inch * 11.75 inch").attr({
			"font-size":tipsfontsize,
			x:A4margin+A4Width*0.5,
			y:A4margin*1.1+tipsfontsize,
			stroke:"none",
			fill:"#d1d1d1",
			"text-anchor":"middle" ,
			"pointer-events":"none"
			

		})
	}
	if(plotted){
		var svgContainer=Model.getElementsByClassName("svgContainer")[0];
		if(originalWidth[currentModel]!=undefined){
			svgContainer.setAttribute("style","width:"+originalWidth[currentModel]*gain[currentModel]/100/printDPI*vDPI+"px;height:"+originalHeight[currentModel]*gain[currentModel]/100/printDPI*vDPI+"px;margin-top:"+(50+A4margin+margin_top_inch*vDPI)+"px;");
			drawing.rect(originalWidth[currentModel]/printDPI*vDPI,originalHeight[currentModel]/printDPI*vDPI).attr({
				stroke:"#aaa",
				fill:"none",
				y:margin_top_inch*vDPI+A4margin,
				x:(A4Width-originalWidth[currentModel]/printDPI*vDPI)/2+A4margin,
				"stroke-dasharray":'4 2'
			})
			drawing.plain("Figure Width: "+(originalWidth[currentModel]/printDPI).toFixed(2)+" inch").attr({
				y:margin_top_inch*vDPI+A4margin-A4Width*0.01,
				x:A4Width/2+A4margin,
				"text-anchor":"middle",
				"dominant-baseline":  "central",
				"font-size":A4Width*0.012,
				fill:"#aaa",
				stroke:"none"
			})
			if(originalHeight[currentModel]/printDPI>0.84){
				drawing.plain("Figure Height: "+(originalHeight[currentModel]/printDPI).toFixed(2)+" inch").attr({
					y:margin_top_inch*vDPI+A4margin+originalHeight[currentModel]/printDPI*vDPI/2,
					x:(A4Width-originalWidth[currentModel]/printDPI*vDPI)/2+A4margin-A4Width*0.01,
					"text-anchor":"middle",
					"dominant-baseline":  "central",
					"font-size":A4Width*0.012,
					fill:"#aaa",
					stroke:"none",
					transform:'rotate(-90,'+((A4Width-originalWidth[currentModel]/printDPI*vDPI)/2+A4margin-A4Width*0.01)+','+(margin_top_inch*vDPI+A4margin+originalHeight[currentModel]/printDPI*vDPI/2)+')'
				})
			}
		}
	}
	
}


createSVG=function(Params,Model_ID){
	//Params=JSON.parse(Params_str);
	
	printDPI=300;
	var res=Params[Model_ID].res;
	plotted=true;
	var Model=document.getElementById(Model_ID);
	var svgContainer=Model.getElementsByClassName("svgContainer")[0];
	svgContainer.innerHTML=null;
	var sizeinfo_top_div=document.createElement("div")
	var sizeinfo_left_div=document.createElement("div")
	svgContainer.appendChild(sizeinfo_top_div)
	svgContainer.appendChild(sizeinfo_left_div)
	//document.getElementById("blankSVG").setAttribute("display","none")
	if (SVG.supported) {
	  var drawing=SVG(svgContainer)
	 
	} else {
	  alert('svg.js not supported');
	  return;
	}
	//var lwd=res.commomParams.barWidth;
	var lineInterval=res.commomParams.barInterval;
	var leftMargin=res.commomParams.leftMargin;
	var rightMargin=res.commomParams.rightMargin;
	var width=res.commomParams.width;
	var DocHeight=res.commomParams.DocHeight;
	if(res.table!=undefined){
		DocHeight+=res.table.height
		if(res.table.plotElements!=undefined){
			var tableElementsNames=Object.keys(res.table.plotElements)
			for(var t=0;t<tableElementsNames.length;t++){
				Params[Model_ID].res.plotElements[tableElementsNames[t]]=res.table.plotElements[tableElementsNames[t]]
			}
		}
		
	}
	var warnings=Model.getElementsByClassName("warnings")[0]
	if(res.warnings!=undefined){
		isWarnings[Params[Model_ID].Model]=true;
		if(typeof res.warnings=="string"){
			var warings_array=[res.warnings]
		}else{
			var warings_array=res.warnings
		}
		warnings.setAttribute("AdditionInfo",JSON.stringify(warings_array))
		warnings.style.display=""
		
		
	}else{
		isWarnings[Model_ID]=false;
		
		warnings.style.display="none"
		warnings.removeAttribute("AdditionInfo")
	}
	jsDebug["DrawingDelay"]=Params[Model_ID].debug.Delay
	var sizeinfo_top=30
	var sizeinfo_left=30

	var svgHeight=DocHeight;
	var svgWidth=width+leftMargin+rightMargin;
	//console.log(svgWidth);
	//drawing.size(svgWidth,svgHeight);
	//drawing.attr({"shape-rendering":"crispEdges"})
	drawing.attr({"text-rendering":"optimizeLegibility"})
	
	drawing.attr({"viewBox":"0 0 "+svgWidth+" "+svgHeight,"style":"background:white"})
	margin_top_inch=11.75/2-svgHeight/printDPI/2;
	if(margin_top_inch<0.1)margin_top_inch=0.1;
	if(margin_top_inch>0.2)margin_top_inch=0.2;
	margin_top_inch=0.2
	svgContainer.setAttribute("style","width:"+svgWidth/printDPI*vDPI+"px;height:"+svgHeight/printDPI*vDPI+"px;margin-top:"+(50+A4margin+margin_top_inch*vDPI)+"px;");
	var plotElements=Params[Model_ID].res.plotElements;
	var plotElementNames=Object.keys(plotElements)
	for(var i=0;i<plotElementNames.length;i++){
		if(plotElements[plotElementNames[i]].plotType=="numericalBarplot"){
			numericalBarplot(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="categoryBarplot"){
			categoryBarplot(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="drawLegend"){
			drawLegend(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="Boxplot"){
			Boxplot(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="drawAxis"){
			drawAxis(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="SummBarplot"){
			SummBarplot(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="drawRect"){
			drawRect(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="drawText"){
			drawText(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="drawPolyline"){
			drawPolyline(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="drawRawSVG"){
			drawRawSVG(drawing,plotElementNames[i],Model_ID)
		}else if(plotElements[plotElementNames[i]].plotType=="drawColBar"){
			drawColBar(drawing,plotElementNames[i],Model_ID)
		}
		
	}
	
	if(Params[Model_ID].res.table!=undefined){
		if(Params[Model_ID].res.table.url!=undefined){
			drawing.rect(Params[Model_ID].res.table.downloadBarWidth,Params[Model_ID].res.table.downloadBarHeight).attr({
				x:Params[Model_ID].res.table.downloadBarLeft,
				y:DocHeight-Params[Model_ID].res.table.downloadBarBaseline-Params[Model_ID].res.table.height-Params[Model_ID].res.table.downloadBarHeight,
				fill:"grey",
				rx:5,
				stroke:""
			})
			drawing.plain("Download cluster result(*.csv)").attr({
				x:Params[Model_ID].res.table.downloadBarLeft+Params[Model_ID].res.table.downloadBarWidth/2,
				y:DocHeight-Params[Model_ID].res.table.downloadBarBaseline-Params[Model_ID].res.table.height-Params[Model_ID].res.table.downloadBarHeight/2,
				"font-size":22+"px",
				"dominant-baseline":  "central",
				"text-anchor":"middle",
				fill:"white",
				onclick:'download_url("cluster_result.csv",'+'"/CLIPA/plotCache/clust_result_'+userParams[currentModel].plotGUID+'.csv")',
				class:"justPointer"
			})
		}
		
	}
	//addColorPicker(Params.Model)
	
	Model.getElementsByClassName("tweakbar")[0].style.display="flex";
	gain[Model_ID]=100;
	//enlargeSVG(100)
	originalWidth[Model_ID]=svgWidth;
	originalHeight[Model_ID]=svgHeight;
	addSVGEvent(Model_ID);
	if(svgWidth/printDPI<6){
		console.log(svgWidth/printDPI)
		setTimeout(function(){
			enlargeSVG(Math.floor(10*6/(svgWidth/printDPI))*10)
		},1000);
	}
	plotSvgBackground()
	/*
	if(svgWidth/printDPI>8){
		setTimeout(function(){
			shrinkSVG(8/(svgWidth/printDPI)*100)
		},1000);
	}
	*/
}

numericalBarplot=function(drawing,plotElementName,Model_ID){

	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var lwd=Params[Model_ID].res.commomParams.barWidth;
	var lineInterval=Params[Model_ID].res.commomParams.barInterval;
	var leftMargin=Params[Model_ID].res.commomParams.leftMargin;
	var width=Params[Model_ID].res.commomParams.width;
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	var baseline=DocHeight-plotElement.baseline
	var minOfData=Math.min.apply(null,NAomit(plotElement.data));
	var maxOfData=Math.max.apply(null,NAomit(plotElement.data));
	for(i=0;i<plotElement.data.length;i++){
		if(plotElement.data[i]=="NA"|isNaN(plotElement.data[i])){
			continue
		}
		if(plotElement.resize!="false"){
			if(maxOfData==minOfData){
				var normed=0
			}else{
				var normed=(plotElement.data[i]-minOfData)/(maxOfData-minOfData)
			}
		}else{
			var normed=plotElement.data[i]
		}
		
		if(plotElement.pillarColor instanceof Object){
			var pillarColor=plotElement.pillarColor[i];
		}else{
			var pillarColor=plotElement.pillarColor;
		}
		if(plotElement.pillarAdditionInfos!=undefined){
			var AdditionInfo=plotElement.pillarAdditionInfos[i]
		}
		/*
		drawing.line(i*(lwd+lineInterval)+leftMargin,baseline,i*(lwd+lineInterval)+leftMargin,baseline-normed*plotElement.buckHeight).attr({
			"stroke":pillarColor,
			"stroke-width":lwd+"px",
			"class":plotElement.pillarClassName,
			"name":Params[Model_ID].res.cell_line[i],
			"AdditionInfo":AdditionInfo,
			"oriValue":plotElement.data[i]
		})
		*/
		drawing.rect(lwd,normed*plotElement.buckHeight).attr({
			x:i*(lwd+lineInterval)+leftMargin,
			y:baseline-normed*plotElement.buckHeight,
			"fill":pillarColor,
			"stroke":"transparent",
			"stroke-width":lineInterval+"px",//stroke是在图像边线居中
			"class":plotElement.pillarClassName,
			"cellLineName":Params[Model_ID].res.cell_line[i],
			"AdditionInfo":AdditionInfo,
			"oriValue":plotElement.data[i]
		})
	}
	if(plotElement.leftText!=undefined){
		drawing.plain(plotElement.leftText).attr({
			x:leftMargin-5,
			y:baseline,
			"font-size":"22px",
			"text-anchor":"end" ,
			"font-weight":"bold",
			"class":plotElement.leftTextClassName,
			"AdditionInfo":plotElement.leftTextAdditionInfos

		})
	}
	if(plotElement.rightText!=undefined){
		drawing.plain(plotElement.rightText).attr({
			x:leftMargin+width+5,
			y:baseline,
			"font-size":"22px",
			"text-anchor":"start" ,
			"fill":plotElement.rightTextColor,
			"font-weight":"bold"
			
		})
	}
}

categoryBarplot=function(drawing,plotElementName,Model_ID){

	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var lwd=Params[Model_ID].res.commomParams.barWidth;
	var lineInterval=Params[Model_ID].res.commomParams.barInterval;
	var leftMargin=Params[Model_ID].res.commomParams.leftMargin;
	var width=Params[Model_ID].res.commomParams.width;
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	var baseline=DocHeight-plotElement.baseline
	for(i=0;i<plotElement.data.length;i++){
		if(plotElement.data[i]=="NA"){
			console.log("skipped")
			continue
		}
		drawing.line(i*(lwd+lineInterval)+leftMargin,baseline,i*(lwd+lineInterval)+leftMargin,baseline-0.8*plotElement.buckHeight).attr({
			stroke:plotElement.color[plotElement.data[i]],
			"class":plotElement.pillarClassName,
			"stroke-width":lwd+"px",
			"class":"cellLineTag",
			"name":Params[Model_ID].res.cell_line[i]
		})
	}
	if(plotElement.leftText!=undefined){
		drawing.plain(plotElement.leftText).attr({
			x:leftMargin-5,
			y:baseline,
			"font-size":"22px",
			"text-anchor":"end" ,
			"font-weight":"bold",
			"class":plotElement.className

		})
	}
	if(plotElement.rightText!=undefined){
		drawing.plain(plotElement.rightText).attr({
			x:leftMargin+width+5,
			y:baseline,
			"font-size":"22px",
			"text-anchor":"start" ,
			"fill":plotElement.rightTextColor,
			"font-weight":"bold"
			
		})
	}
}

drawRawSVG=function(drawing,plotElementName,Model_ID){
	//console.log(plotElementName)
	//debug_drawing=drawing
	var Model=document.getElementById(Model_ID)
	var svgContainer=Model.getElementsByClassName("svgContainer")[0];
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	svgContainer.innerHTML=plotElement.data
	//drawing.node.appendChild(plotElement.data)


}
drawLegend=function(drawing,plotElementName,Model_ID){
	
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var lwd=Params[Model_ID].res.commomParams.barWidth;
	var lineInterval=Params[Model_ID].res.commomParams.barInterval;
	var leftMargin=Params[Model_ID].res.commomParams.leftMargin;
	var width=Params[Model_ID].res.commomParams.width;
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	if(plotElement.comment!=undefined){
		drawing.plain(plotElement.comment.text).attr({
			x:plotElement.left+leftMargin,
			y:DocHeight-plotElement.comment.baseline,
			"font-size":"22px",
			"text-anchor":"start" ,
			"font-weight":"bold"

		})
	}
	if(plotElement.title!=undefined){
		drawing.plain(plotElement.title.text).attr({
			x:plotElement.title.left+leftMargin,
			y:DocHeight-plotElement.title.baseline,
			"font-size":"22px",
			"text-anchor":"start" ,
			"font-weight":"bold"

		})
	}
	if(plotElement.legend!=undefined){
		var lengendNames=Object.keys(plotElement.legend);

		for(var l=0;l<lengendNames.length;l++){
			drawing.rect(plotElement.legend[lengendNames[l]].size,plotElement.legend[lengendNames[l]].size).attr({
				x:plotElement.legend[lengendNames[l]].left+leftMargin,
				y:DocHeight-plotElement.legend[lengendNames[l]].baseline-plotElement.legend[lengendNames[l]].size,
				fill:plotElement.legend[lengendNames[l]].color,
				stroke:"none"
			
			})
			
			drawing.plain(plotElement.legend[lengendNames[l]].text).attr({
				x:plotElement.legend[lengendNames[l]].left+plotElement.legend[lengendNames[l]].size+plotElement.legend[lengendNames[l]].size/4+leftMargin,
				y:DocHeight-plotElement.legend[lengendNames[l]].baseline,
				"font-size":"22px",
				"text-anchor":"start" ,
				"font-weight":"bold"

			})
		}
	}


}
drawRect=function(drawing,plotElementName,Model_ID){
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var leftMargin=Params[Model_ID].res.commomParams.leftMargin;
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	
	var className=plotElement.className
	if(className==undefined){
		className=""
	}
	if(plotElement.border==undefined){
		plotElement.border=""
	}
	if(plotElement.color==undefined){
		plotElement.color="transparent"
	}
	drawing.rect(plotElement.width,plotElement.height).attr({
		x:plotElement.left+leftMargin,
		y:DocHeight-plotElement.bottom-plotElement.height,
		fill:plotElement.color,
		stroke:plotElement.border,
		class:className
	})
}
drawPolyline=function(drawing,plotElementName,Model_ID){
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var leftMargin=Params[Model_ID].res.commomParams.leftMargin;
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	if(plotElement.lty!=undefined){
		if(plotElement.lty=="dashed"){
			var lty='10'
		}else{
			var lty=''
		}
	}else{
		var lty=''
	}
	var pos=[]
	for(var a=0;a<plotElement.xs.length;a++){
		pos.push(plotElement.xs[a]+leftMargin)
		pos.push(DocHeight-plotElement.ys[a])
	}
	console.log(pos)
	if(plotElement.className!=undefined){
		var className=plotElement.className
	}else{
		var className="";
	}
	drawing.polyline(pos).attr({
		stroke:plotElement.color,
		'stroke-width':plotElement.lwd*3,
		fill:"none",
		'stroke-dasharray':lty,
		class:className
		
	})
}
Boxplot=function(drawing,plotElementName,Model_ID){
	
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	var baseline=plotElement.baseline
	drawing.rect(plotElement.boxWidth,plotElement.Qs[3]-plotElement.Qs[1]).attr({
		x:plotElement.x-plotElement.boxWidth/2+Params[Model_ID].res.commomParams.leftMargin,
		y:DocHeight-plotElement.Qs[3]-baseline,
		fill:plotElement.color,
		"stroke-width":"1.5px",
		stroke:"black"
	})
	drawing.line(plotElement.x-plotElement.boxWidth/2+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[0]-baseline,plotElement.x+plotElement.boxWidth/2+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[0]-baseline).attr({
		stroke:"black",
		'stroke-width':"1.5px"
	})
	drawing.line(plotElement.x-plotElement.boxWidth/2+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[2]-baseline,plotElement.x+plotElement.boxWidth/2+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[2]-baseline).attr({
		stroke:"black",
		'stroke-width':"1.5px"
	})
	drawing.line(plotElement.x-plotElement.boxWidth/2+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[4]-baseline,plotElement.x+plotElement.boxWidth/2+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[4]-baseline).attr({
		stroke:"black",
		'stroke-width':"1.5px"
	})
	drawing.line(plotElement.x+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[0]-baseline,plotElement.x+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[1]-baseline).attr({
		stroke:"black",
		'stroke-width':"1.5px"
	})
	drawing.line(plotElement.x+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[3]-baseline,plotElement.x+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.Qs[4]-baseline).attr({
		stroke:"black",
		'stroke-width':"1.5px"
	})
	if(plotElement.plotData!=undefined){
		var g=drawing.group().attr({
			class:"showBoxPlot"
		});
		var rectBG=drawing.rect(plotElement.plotData.width,Params[Model_ID].res.plotElements.axis.y.startEnd[1]-Params[Model_ID].res.plotElements.axis.y.startEnd[0]).attr({
			x:plotElement.x-plotElement.plotData.width/2+Params[Model_ID].res.commomParams.leftMargin,
			y:DocHeight-Params[Model_ID].res.plotElements.axis.y.startEnd[1],
			fill:"#ffffff88",
			stroke:"none"
		})
		g.add(rectBG)
		for(var p=0;p<plotElement.plotData.ys.length;p++){
			var c=drawing.circle(plotElement.plotData.r).attr({
				cx:plotElement.plotData.xs[p]+Params[Model_ID].res.commomParams.leftMargin,
				cy:DocHeight-plotElement.plotData.ys[p]-baseline,
				fill:plotElement.color,
				//"stroke-width":"0.5px",
				//stroke:"black",
				cellLineName:plotElement.plotData.cellNames[p],
				class:"CellAndDrugResistTags"
			})
			g.add(c)
		}
		
	}
}
drawText=function(drawing,plotElementName,Model_ID){
	//text(plotElement$text,x=plotElement$x,y=plotElement$y,cex=plotElement$cex,col=plotElement$color,adj=plotElement$adj,srt=plotElement$srt)
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var leftMargin=Params[Model_ID].res.commomParams.leftMargin;
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	if(plotElement.adj[0]==0){
		var xadj="start"
	}else if(plotElement.adj[0]==1){
		var xadj="end"
	}else{
		var xadj="middle"
	}
	
	if(plotElement.adj[1]==0){
		var yadj="text-before-edge"
	}else if(plotElement.adj[1]==1){
		var yadj="text-after-edge"
	}else{
		var yadj="central"
	}
	
	drawing.plain(plotElement.text).attr({
		x:plotElement.x+Params[Model_ID].res.commomParams.leftMargin,
		y:DocHeight-plotElement.y,
		"font-size":plotElement.cex*50+"px",
		"dominant-baseline":  yadj,
		"text-anchor":xadj,
		fill:plotElement.color,
		transform:"rotate("+(-plotElement.srt)+" "+(plotElement.x+Params[Model_ID].res.commomParams.leftMargin)+","+(DocHeight-plotElement.y)+")"
	})
}
SummBarplot=function(drawing,plotElementName,Model_ID){
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	bardebug=plotElement
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	var baseline=plotElement.baseline
	var color=[]
	var data=[]
	var subText=[]
	var cell_lines_text=[]
	var mainText=plotElement.popupText.main
	var typeNames=Object.keys(plotElement.data);
	for(var t=0;t<typeNames.length;t++){
		color.push(plotElement.color[typeNames[t]])
		data.push(plotElement.data[typeNames[t]])
		if(plotElement.popupText.Sub_Text!=undefined){
			subText.push(plotElement.popupText.Sub_Text[typeNames[t]])
		}
		cell_lines_text.push(plotElement.popupText.cell_lines_text[typeNames[t]])
	}
	
	
	var baseline_tmp=0;
	for(var i=0;i<data.length;i++){
		baseline_tmp+=data[i]
		if(plotElement.popupText.Sub_Text!=undefined){
			var popupText='<p style="margin-bottom:2px;border-bottom:1px solid black;font-size:18px;text-align: center;">'+mainText+"</p><table>"
			for(var t=0;t<typeNames.length;t++){
				
				if(typeNames.length-t-1==i){
					var sty='style="background:#c6dff9;"';
				}else{
					var sty="";
				}
				if(typeof cell_lines_text[typeNames.length-t-1]=="string"){
					cell_lines_text_tmp='<a target="view_window" href=https://cellmodelpassports.sanger.ac.uk/passports/'+model_ids[cell_lines_text[typeNames.length-t-1]]+'>'+cell_lines_text[typeNames.length-t-1]+'</a>'
				}else{
					if(cell_lines_text[typeNames.length-t-1].length==0){

						continue;
					}
					//var cell_lines_text_tmp=cell_lines_text[typeNames.length-t-1].join(", ")
					var cell_lines_text_tmp=[]
					for(var j=0;j<cell_lines_text[typeNames.length-t-1].length;j++){
						cell_lines_text_tmp.push('<a target="view_window" href=https://cellmodelpassports.sanger.ac.uk/passports/'+model_ids[cell_lines_text[typeNames.length-t-1][j]]+'>'+cell_lines_text[typeNames.length-t-1][j]+'</a>')
					}
					var cell_lines_text_tmp=cell_lines_text_tmp.join(", ")
				}
				popupText+='<tr '+sty+'><td><p>'+subText[typeNames.length-t-1]+':</p></td><td><p>'+cell_lines_text_tmp+"</p></td></tr>"
			}
			popupText+="</table>"
		}else{
			var popupText='<p style="font-size:12px">'+mainText+":<br>"
			var cell_lines_text_tmp=[]
			if(typeof cell_lines_text[i]=="string"){
				cell_lines_text_tmp='<a target="view_window" href=https://cellmodelpassports.sanger.ac.uk/passports/'+model_ids[cell_lines_text[i]]+'>'+cell_lines_text[i]+'</a>'
			}else{
				if(cell_lines_text[i].length==0){
						continue;
				}
				//var cell_lines_text_tmp=cell_lines_text[typeNames.length-t-1].join(", ")
				var cell_lines_text_tmp=[]
				for(var j=0;j<cell_lines_text[i].length;j++){
					cell_lines_text_tmp.push('<a target="view_window" href=https://cellmodelpassports.sanger.ac.uk/passports/'+model_ids[cell_lines_text[i][j]]+'>'+cell_lines_text[i][j]+'</a>')
				}
				var cell_lines_text_tmp=cell_lines_text_tmp.join(", ")
			}
			popupText+=cell_lines_text_tmp+"</p>"
		}
		
		drawing.rect(plotElement.barWidth,data[i]).attr({
			x:plotElement.x-plotElement.barWidth/2+Params[Model_ID].res.commomParams.leftMargin,
			y:DocHeight-baseline_tmp-baseline,
			fill:color[i],
			stroke:"none",
			class:"summBarplotInfo",
			AdditionInfo:popupText
		})
	}
}

drawAxis=function(drawing,plotElementName,Model_ID){
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	if(plotElement.box!=undefined){
		if(plotElement.box){
			var isBox=true
		}else{
			var isBox=false
		}
	}else{
		var isBox=false
	}
	if(isBox){
		drawing.rect(plotElement.x.startEnd[1]-plotElement.x.startEnd[0],plotElement.y.startEnd[1]-plotElement.y.startEnd[0]).attr({
			x:Params[Model_ID].res.commomParams.leftMargin,
			y:DocHeight-plotElement.y.startEnd[1],
			fill:"none",
			stroke:"black",
			'stroke-width':"1.5px"
		})
		
	}
	if(plotElement.x!=undefined){
		if(!isBox){
			drawing.line(plotElement.x.startEnd[0]+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.x.baseline,plotElement.x.startEnd[1]+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.x.baseline).attr({
				stroke:"black",
				'stroke-width':"1.5px"
			})
		}
		if(plotElement.x.srt!=0){
			var textAnchor="end"
			var transformAng=-plotElement.x.srt
			
		}else{
			var textAnchor="middle"
			var transformAng=0
		}
		for(var a=0;a<plotElement.x.ats.length;a++){
			drawing.line(plotElement.x.ats[a]+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.x.baseline,plotElement.x.ats[a]+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.x.baseline+10).attr({
				stroke:"black",
				'stroke-width':"1.5px"
			})
			drawing.plain(plotElement.x.labs[a]).attr({
				x:plotElement.x.ats[a]+Params[Model_ID].res.commomParams.leftMargin,
				y:DocHeight-plotElement.x.baseline+20,
				"font-size":"24px",
				"dominant-baseline":  "middle",
				"text-anchor":textAnchor,
				//"font-family":"Helvetica",
				transform:"rotate("+transformAng+" "+(plotElement.x.ats[a]+Params[Model_ID].res.commomParams.leftMargin)+","+(DocHeight-plotElement.x.baseline+10)+")"
			})
			
		}
		if(plotElement.x.xlab!=undefined){
			drawing.plain(plotElement.x.xlab).attr({
				x:(plotElement.x.startEnd[1]+plotElement.x.startEnd[0])/2+Params[Model_ID].res.commomParams.leftMargin,
				y:DocHeight-plotElement.x.baseline+75,
				"dominant-baseline":  "middle",
				"font-size":"36px",
				"text-anchor":"middle",
				"font-weight":"bold"
			})
			
		}
	}
	if(plotElement.y!=undefined){
		if(!isBox){
			drawing.line(plotElement.y.x+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.y.startEnd[0],plotElement.y.x+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.y.startEnd[1]).attr({
				stroke:"black",
				'stroke-width':"1.5px"
			})
		}
		for(var a=0;a<plotElement.y.ats.length;a++){
			drawing.line(plotElement.y.x+Params[Model_ID].res.commomParams.leftMargin,DocHeight-plotElement.y.ats[a],plotElement.y.x+Params[Model_ID].res.commomParams.leftMargin-10,DocHeight-plotElement.y.ats[a]).attr({
				stroke:"black",
				'stroke-width':"1.5px"
			})
			drawing.plain(plotElement.y.labs[a]).attr({
				x:plotElement.y.x+Params[Model_ID].res.commomParams.leftMargin-15,
				y:DocHeight-plotElement.y.ats[a],
				"dominant-baseline":  "middle",
				"font-size":"24px",
				"text-anchor":"end"
			})
		}
		if(plotElement.y.ylab!=undefined){
			if(plotElement.y.ylabRescale!=undefined){
				var ylabRescale=plotElement.y.ylabRescale
			}else{
				var ylabRescale=1
			}
			drawing.plain(plotElement.y.ylab).attr({
				x:plotElement.y.x+Params[Model_ID].res.commomParams.leftMargin-95,
				y:DocHeight-(plotElement.y.startEnd[1]+plotElement.y.startEnd[0])/2,
				"dominant-baseline":  "middle",
				"font-size":36*ylabRescale+"px",
				"text-anchor":"middle",
				"font-weight":"bold",
				transform:"rotate(-90 "+(plotElement.y.x+Params[Model_ID].res.commomParams.leftMargin-95)+","+(DocHeight-(plotElement.y.startEnd[1]+plotElement.y.startEnd[0])/2)+")"
			})
			
		}
		
	}
	
}
drawColBar=function(drawing,plotElementName,Model_ID){
	var plotElement=Params[Model_ID].res.plotElements[plotElementName]
	var leftMargin=Params[Model_ID].res.commomParams.leftMargin;
	var DocHeight=Params[Model_ID].res.commomParams.DocHeight;
	var gradient = drawing.gradient('linear', function (stop) {
		stop.at(0, plotElement.color_range[0]);
		stop.at(0.5, plotElement.color_range[1]);
		stop.at(1, plotElement.color_range[2]);
	}).from(0, 1).to(0, 0);
	
	if(plotElement.border==undefined){
		plotElement.border=""
	}
	if(plotElement.color==undefined){
		plotElement.color="transparent"
	}
	drawing.rect(plotElement.width,plotElement.height).attr({
		x:plotElement.left+leftMargin,
		y:DocHeight-plotElement.bottom-plotElement.height,
		fill:gradient,
		stroke:plotElement.border
	})
	drawing.plain("High").attr({
		x:plotElement.left+leftMargin+plotElement.width/2,
		y:DocHeight-plotElement.bottom-plotElement.height-10,
		"font-size":plotElement.cex*50+"px",
		"dominant-baseline":  "text-after-edge",
		"text-anchor":"middle",
		fill:"black"
	})
	drawing.plain("Low").attr({
		x:plotElement.left+leftMargin+plotElement.width/2,
		y:DocHeight-plotElement.bottom+10,
		"font-size":plotElement.cex*50+"px",
		"dominant-baseline":  "text-before-edge",
		"text-anchor":"middle",
		fill:"black"
	})

}
NAomit=function(arr){
	var newarr=[];
	for(var i =0;i<arr.length;i++){
		if(arr[i]!="NA"&!isNaN(arr[i])){
			newarr.push(arr[i]);
		}
	}
	return(newarr)
}
addColorPicker=function(Model_Name){
	if(Params.res.plotElements.legend.legend!=undefined){
		var legends=Params.res.plotElements.legend.legend;
		var Model=document.getElementById(Model_Name);
		
		if(Model.getElementsByClassName("colorPickers").length==0){
			var tweakbar=Model.getElementsByClassName("tweakbar")[0];
			var colorPickers=document.createElement("div");
			tweakbar.appendChild(colorPickers);
		}else{
			var colorPickers=Model.getElementsByClassName("colorPickers")[0]
			colorPickers.innerHTML="";
		}
			
			colorPickers.setAttribute("class","colorPickers")
			colorPickers.innerHTML="Pick your own color:"
			var legendsNames=Object.keys(legends);
			for(var i=0;i<legendsNames.length;i++){
				var color=document.createElement("input")
				color.type="color";
				color.value=legends[legendsNames[i]].color
				color.setAttribute("valueOfColor",legends[legendsNames[i]].text)
				color.setAttribute("onchange","changeColor(this)")
				var p=document.createElement("p")
				p.innerHTML=legends[legendsNames[i]].text
				colorPickers.appendChild(color);
				colorPickers.appendChild(p);
			}
			
		
	}
}

changeColor=function(obj){
	userParams[currentModel].colors[obj.getAttribute("valueOfColor")+"_color"]=obj.value
	createSVG(Params_str);
}
shrinkSVG=function(tar=""){
	var gain_former=gain[currentModel];
	if(tar!=""){
		gain[currentModel]=Math.round(tar)
	}else{
		gain[currentModel]=(Math.floor(gain[currentModel]/1.1/10))*10;
	}
	
	
	if(gain[currentModel]<10){
		gain[currentModel]=10;
	}
	var Model=document.getElementById(currentModel);
	var svgContainer=Model.getElementsByClassName("svgContainer")[0];
	svgContainer.setAttribute("style","width:"+originalWidth[currentModel]*gain[currentModel]/100/printDPI*vDPI+"px;height:"+originalHeight*gain[currentModel]/100/printDPI*vDPI+"px;margin-top:"+(50+A4margin+margin_top_inch*vDPI)+"px;transition:"+(gain_former-gain[currentModel])*0.005+"s width ease,"+(gain_former-gain[currentModel])*0.005+"s height ease;");
	Model.getElementsByClassName("zoombox")[0].children[1].innerHTML=gain[currentModel]+"%";
}

enlargeSVG=function(tar=""){
	var gain_former=gain[currentModel];
	if(tar!=""){
		gain[currentModel]=Math.round(tar)
	}else{
		gain[currentModel]=(Math.ceil(1.1*gain[currentModel]/10))*10;
	}
	
	if(gain[currentModel]>=1000){
		gain[currentModel]=1000;
	}
	var Model=document.getElementById(currentModel);
	var svgContainer=Model.getElementsByClassName("svgContainer")[0];
	svgContainer.setAttribute("style","width:"+originalWidth[currentModel]*gain[currentModel]/100/printDPI*vDPI+"px;height:"+originalHeight*gain[currentModel]/100/printDPI*vDPI+"px;margin-top:"+(50+A4margin+margin_top_inch*vDPI)+"px;transition:"+(gain[currentModel]-gain_former)*0.005+"s width ease,"+(gain[currentModel]-gain_former)*0.005+"s height ease;");
	Model.getElementsByClassName("zoombox")[0].children[1].innerHTML=gain[currentModel]+"%";
}

initSVG=function(){
	var gain_former=gain[currentModel];
	gain[currentModel]=100;
	var Model=document.getElementById(currentModel);
	var svgContainer=Model.getElementsByClassName("svgContainer")[0];
	svgContainer.setAttribute("style","width:"+originalWidth[currentModel]*gain[currentModel]/100/printDPI*vDPI+"px;height:"+originalHeight*gain[currentModel]/100/printDPI*vDPI+"px;margin-top:"+(50+A4margin+margin_top_inch*vDPI)+"px;transition:"+Math.abs((gain[currentModel]-gain_former))*0.005+"s width ease,"+Math.abs((gain[currentModel]-gain_former))*0.005+"s height ease;");
	Model.getElementsByClassName("zoombox")[0].children[1].innerHTML=gain[currentModel]+"%";
}

downloadFig=function(obj){
	var Model_ID=currentModel
	if(obj.value==".svg"){
		console.log("直接blob，待做。。。。");
	}
	var filetype=obj.value;
	
	obj.selectedIndex=0;
	var Model=document.getElementById(Model_ID);
	//var DebugTips=document.getElementById("NATIVE,HELP");
	var dataBaseName=Model_ID.split(',')[0]
	var ModelName=Model_ID.split(',')[1]
	var Children_names=initParams.ModelInfo[dataBaseName][ModelName];
	//var Debug=document.getElementById("NATIVE,HELP");
	userParams[Model_ID].Model=Model_ID;
	userParams[Model_ID].use="download";
	var requestGUID=guid();
	userParams[Model_ID].requestGUID=requestGUID;
	userParams[Model_ID].path="../plotCache/data_download_"+userParams[Model_ID].plotGUID+".json"
	userParams[Model_ID].filetype=filetype;
	
	
	jsDebug=[];
	xmlhttp=new XMLHttpRequest();

	xmlhttp.onreadystatechange=function(){
		//document.getElementById("plotarea").innerHTML =genename+"got:"+xmlhttp.responseText;
		if (xmlhttp.readyState==4 && xmlhttp.status==200){
			
			var Params_str=xmlhttp.responseText;
			
			//console.log(Params_str)
			//document.getElementById("svgContainer").innerHTML=null;
			jsDebug["jsDrawingStart"]=new Date().getTime()
			try {
				Params_dl=JSON.parse(Params_str);
			} catch (error) {
				alert(Params_str);
				//DebugTips.style.background="red";
				return;
			}
			if(Params_dl.requestGUID!=requestGUID)return
			
			if(Params_dl.error){
				alert(new Date().toLocaleTimeString()+":<br>"+Params_dl.error)
				//DebugTips.style.background="red";
				return
			}
			//DebugTips.style.background="green";
			
			
			download_url(userParams[Model_ID].plotGUID+filetype,Params_dl.Figpath);
		}else{
			
		}
			
	}
	//document.getElementById("plotarea").innerHTML= +"start";
	userParams_uri="userParams="+btoa(JSON.stringify(userParams[Model_ID]));
	//console.log(userParams[Model_ID]_uri)
	xmlhttp.open("POST","php/sendToR.php",true);
	xmlhttp.setRequestHeader("Content-Type","application/x-www-form-urlencoded");
	jsDebug["TotalStart"]=new Date().getTime();
	xmlhttp.send(userParams_uri);
	

}

function download_url(name, url) {
	console.log(url)
	function fake_click(obj) {
		var ev = document.createEvent("MouseEvents");
		ev.initMouseEvent(
			"click", true, false, window, 0, 0, 0, 0, 0, false, false, false, false, 0, null
		);
		obj.dispatchEvent(ev);
	}
  

    var save_link = document.createElementNS("http://www.w3.org/1999/xhtml", "a")
    save_link.href = url;
    save_link.download = name;
    fake_click(save_link);
	//调用方法
	//download("save.txt","内容");
}


function download_blob(name, data) {
	function fake_click(obj) {
		var ev = document.createEvent("MouseEvents");
		ev.initMouseEvent(
			"click", true, false, window, 0, 0, 0, 0, 0, false, false, false, false, 0, null
		);
		obj.dispatchEvent(ev);
	}
    var urlObject = window.URL || window.webkitURL || window;

    var downloadData = new Blob([data]);

    var save_link = document.createElementNS("http://www.w3.org/1999/xhtml", "a")
    save_link.href = urlObject.createObjectURL(downloadData);
    save_link.download = name;
    fake_click(save_link);
	//调用方法
	//download("save.txt","内容");
}

getlength=function(obj){
	return(Object.keys(obj).length)
	
}

