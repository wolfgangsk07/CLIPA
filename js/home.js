function updateHomeDiagram(){
	var diagram=document.getElementsByClassName("diagram")[0]
	var plotWidth=300
	var databaseNames=["GDSC","DepMap","NCI60"]
	var colorPool=['#E41A1C', '#A73C52', '#6B5F88', '#3780B3', '#3F918C', '#47A266',
		'#53A651', '#6D8470', '#87638F', '#A5548D', '#C96555', '#ED761C',
		'#FF9508', '#FFC11A', '#FFEE2C', '#EBDA30', '#CC9F2C', '#AD6428',
		'#BB614F', '#D77083', '#F37FB8', '#DA88B3', '#B990A6', '#999999']
	for(var d=0;d<databaseNames.length;d++){
		var database_div=document.createElement("div")
			database_div.setAttribute("class","database_div")
			var p=document.createElement("p")
				p.innerHTML=databaseNames[d]
			database_div.appendChild(p)
		//开始画PIE
		var pieR=plotWidth/2
		var subdiv=document.createElement("div")
			subdiv.setAttribute("class","subdiagram")
			/*
			var infoDiv=document.createElement("div")
				infoDiv.setAttribute("class","subdiagram_info")
					var ul=document.createElement("ul")
						for(var l=0;l<7;l++){
							var li=document.createElement("li")
								li.innerHTML="xxxxxxxxxxx"
								ul.appendChild(li)
						}
					infoDiv.appendChild(ul)
					subdiv.appendChild(infoDiv)
					*/
			var pieDiv=document.createElement("div")
				pieDiv.setAttribute("class","subdiagram_gram")
			
			
			subdiv.appendChild(pieDiv)
		database_div.appendChild(subdiv)
		var pie=SVG(pieDiv)
		var Vmargin=0.1*pieR
		var Hmargin=0.1*pieR
		pie.attr({"text-rendering":"optimizeLegibility"})
		//pie.width(pieR*2+Hmargin+"px")
		//pie.height(pieR*2+Vmargin+"px")
		pie.width("100%")
		pie.height("100%")
		pie.viewbox(0,0,pieR*2+Hmargin,pieR*2+Vmargin)
		var cancerTypes=Object.keys(initParams.cellLines[databaseNames[d]])
		var pieData=[]
		var totalCells=0
		for(var c=0;c<cancerTypes.length;c++){
			var cellNames=Object.keys(initParams.cellLines[databaseNames[d]][cancerTypes[c]].cellLines)
			pieData.push(cellNames.length)
			totalCells+=cellNames.length
		}
		var perimeter=Math.PI*2*pieR
		var plottedLength=0
		var pieBlankRatio=0.7
		for(var p=0;p<pieData.length;p++){
			
			pie.path(drawPie(pieR+Hmargin/2,pieR+Vmargin/2,pieR,pieR*pieBlankRatio,pieData[p]/totalCells*2*Math.PI,plottedLength/totalCells*2*Math.PI)).attr({
				stroke:"white",
				"stroke-width":"0.5px",
				fill:colorPool[p%colorPool.length],
				class:"PieInfo"
			})
			var g=pie.group().attr({
				class:"PieInfo_txt"
			})
			var txt1=pie.plain(cancerTypes[p]+": ").attr({
				x:pieR+Hmargin/2,
				y:pieR+Vmargin/2-7,
				"text-anchor":"middle",
				"dominant-baseline":  "central",
				"font-size":13,
				fill:colorPool[p%colorPool.length],
				stroke:"none"
			})
			var txt2=pie.plain(pieData[p]+" Cell Lines ("+(pieData[p]/totalCells*100).toFixed(2)+"%)").attr({
				x:pieR+Hmargin/2,
				y:pieR+Vmargin/2+7,
				"text-anchor":"middle",
				"dominant-baseline":  "central",
				"font-size":14,
				
				fill:colorPool[p%colorPool.length],
				stroke:"none"
			})
			g.add(txt1)
			g.add(txt2)
			plottedLength+=pieData[p]
		}
		pie.plain(toThousands(totalCells)).attr({
			y:pieR+Vmargin/2-12,
			x:pieR+Hmargin/2,
			"text-anchor":"middle",
			"dominant-baseline":  "central",
			"font-size":24,
			"font-weight":"bold",
			fill:"#888",
			stroke:"none",
			class:"PieInfo_maintxt"
		})
		pie.plain("Cell Lines").attr({
			y:pieR+Vmargin/2+12,
			x:pieR+Hmargin/2,
			"text-anchor":"middle",
			"dominant-baseline":  "central",
			"font-size":24,
			"font-weight":"bold",
			fill:"#888",
			stroke:"none",
			class:"PieInfo_maintxt"
		})
		//开始画细胞图
		var histWidth=plotWidth
		var histHeight=plotWidth
		var subdiv=document.createElement("div")
			subdiv.setAttribute("class","subdiagram")
			/*
			var infoDiv=document.createElement("div")
				infoDiv.setAttribute("class","subdiagram_info")
					var ul=document.createElement("ul")
						for(var l=0;l<7;l++){
							var li=document.createElement("li")
								li.innerHTML="xxxxxxxxxxx"
								ul.appendChild(li)
						}
					infoDiv.appendChild(ul)
		
			subdiv.appendChild(infoDiv)
			*/
				var histADiv=document.createElement("div")
				histADiv.setAttribute("class","subdiagram_gram")
			subdiv.appendChild(histADiv)
		database_div.appendChild(subdiv)
		var counts=initParams.homeInfo.cellInfo[databaseNames[d]].counts
		var Omics=initParams.homeInfo.cellInfo[databaseNames[d]].Omics
		var Ylab="Count of Cell Lines"
		drawHomeMultiBar(Omics,counts,histADiv,histWidth,histHeight,colorPool,Ylab)
		
		//开始画gene柱状图
		var histWidth=plotWidth
		var histHeight=plotWidth
		var subdiv=document.createElement("div")
			subdiv.setAttribute("class","subdiagram")
			/*
			var infoDiv=document.createElement("div")
				infoDiv.setAttribute("class","subdiagram_info")
					var ul=document.createElement("ul")
						for(var l=0;l<7;l++){
							var li=document.createElement("li")
								li.innerHTML="xxxxxxxxxxx"
								ul.appendChild(li)
						}
					infoDiv.appendChild(ul)
			
			subdiv.appendChild(infoDiv)
			*/
			var histADiv=document.createElement("div")
				histADiv.setAttribute("class","subdiagram_gram")
			subdiv.appendChild(histADiv)
		database_div.appendChild(subdiv)
		var counts=initParams.homeInfo.genesInfo[databaseNames[d]].counts
		var Omics=initParams.homeInfo.genesInfo[databaseNames[d]].Omics
		var Ylab="Count of Cases"
		drawHomeBar(Omics,counts,histADiv,histWidth,histHeight,colorPool,Ylab)
		
		
		diagram.appendChild(database_div)
	}
	
}

function drawPie(x,y,r,r0,por,start){
	//x,y代表圆心位置，r代表扇形圆的半径，r0代表缺损扇形的半径，por代表扇形弧度，start代表扇形开始的位置（从6点钟方向逆时针）
	var d="M "+(x+Math.sin(start)*r0)+" "+(y+Math.cos(start)*r0)
	d+=" L "+(x+Math.sin(start)*r)+" "+(y+Math.cos(start)*r)
	if(por>Math.PI){
		large=1
	}else{
		large=0
	}
	d+=" A "+r+" "+r+" "+0+" "+large+" "+0+" "+(x+Math.sin(start+por)*r)+" "+(y+Math.cos(start+por)*r)
	d+=" L "+(x+Math.sin(start+por)*r0)+" "+(y+Math.cos(start+por)*r0)
	d+=" A "+r0+" "+r0+" "+0+" "+large+" "+1+" "+(x+Math.sin(start)*r0)+" "+(y+Math.cos(start)*r0)
	d+=" Z"
	return(d)
}


function median(arr){
  const mid = Math.floor(arr.length / 2),
    nums = [...arr].sort((a, b) => a - b);
  return arr.length % 2 !== 0 ? nums[mid] : (nums[mid - 1] + nums[mid]) / 2;
};

function toThousands(num = 0){
   return num.toString().replace(/\d+/, function(n) {
      return n.replace(/(\d)(?=(?:\d{3})+$)/g, '$1,');
   });
};

function drawHomeBar(Omics,counts,histADiv,histWidth,histHeight,colorPool,Ylab){
	
	var marginBottom=100
	var marginLeft=70
	var barVsBlank=2
	var hist=SVG(histADiv)
		hist.attr({"text-rendering":"optimizeLegibility"})
		//hist.width(histWidth+"px")
		//hist.height(histHeight+"px")
		hist.width("100%")
		hist.height("100%")
		hist.viewbox(0,0,histWidth,histHeight)
		/*
		if(Omics.indexOf("DNA methylation (Site)")>-1&counts[Omics.indexOf("DNA methylation (Site)")]!=0){
			var cglength=counts[Omics.indexOf("DNA methylation (Site)")]
			var counts_remain=counts
			counts_remain[Omics.indexOf("DNA methylation (Site)")]=0
			counts[Omics.indexOf("DNA methylation (Site)")]=Math.max.apply(null,counts_remain)*1.2
			
		}
		*/
		var countTags=[]
		for(var c=0;c<counts.length;c++){
			if(counts[c]>100000){
				countTags.push(counts[c])
				var counts_remain=counts
				counts_remain[c]=0
				counts[c]=Math.max.apply(null,counts_remain)*1.2
			}else{
				countTags.push(counts[c])
			}
		}
		var scale=(histHeight-marginBottom)/Math.max.apply(null,counts)*0.8
		var barWidth=(histWidth-marginLeft)/Omics.length*barVsBlank/(1+barVsBlank)
		for(var O=0;O<Omics.length;O++){
			var g=hist.group().attr({
				class:"showHistTag"
			})
			if(countTags[O]!=counts[O]){
				var gap=5
				
				var x1=marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2
				var x2=marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2+barWidth
				var po1=hist.polyline([x1,histHeight-marginBottom-scale*counts[O]*0.6,x2,histHeight-marginBottom-scale*counts[O]*0.6-barWidth,x2,histHeight-marginBottom-scale*counts[O],x1,histHeight-marginBottom-scale*counts[O]]).attr({
					fill:colorPool[O%colorPool.length],
					stroke:"none",
					"stroke-width":"1px",
					class:"Hist"
				})
				var po2=hist.polyline([x1,histHeight-marginBottom-scale*counts[O]*0.6+gap,x1,histHeight-marginBottom,x2,histHeight-marginBottom,x2,histHeight-marginBottom-scale*counts[O]*0.6+gap-barWidth]).attr({
					fill:colorPool[O%colorPool.length],
					stroke:"none",
					"stroke-width":"1px",
					class:"Hist"
				})
				g.add(po1)
				g.add(po2)
			}else{
				
				var r1=hist.rect(barWidth,scale*counts[O]).attr({
					x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2,
					y:histHeight-marginBottom-scale*counts[O],
					fill:colorPool[O%colorPool.length],
					stroke:"none",
					"stroke-width":"1px",
					class:"Hist"
					
				})
				g.add(r1)
			}
			
			var bg=hist.rect(barWidth,histHeight-marginBottom).attr({
				x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2,
				y:0,
				fill:"transparent",
				stroke:"transparent"
				
			})
			var tag=hist.plain(toThousands(countTags[O])).attr({
				x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2+barWidth/2,
				y:histHeight-marginBottom-scale*counts[O],
				fill:colorPool[O%colorPool.length],
				stroke:"none",
				"font-size":barWidth/4,
				"text-anchor":"middle",
				"font-weight":"bold",
				"dominant-baseline":  "text-after-edge"
				
			})
			g.add(bg)
			g.add(tag)
			hist.plain(Ylab).attr({
				x:marginLeft-16,
				y:(histHeight-marginBottom)/2,
				fill:"black",
				stroke:"none",
				"font-size":16,
				"text-anchor":"middle",
				"dominant-baseline":  "central",
				"transform":"rotate(-90 "+(marginLeft-16)+","+(histHeight-marginBottom)/2+")"
			})
			
			hist.plain(Omics[O]).attr({
				x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2+barWidth/2,
				y:histHeight-marginBottom+4,
				fill:"black",
				stroke:"none",
				"font-size":10,
				"text-anchor":"end",
				"dominant-baseline":  "central",
				"transform":"rotate(-45 "+(marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2+barWidth/2)+","+(histHeight-marginBottom+4)+")"
			})
			
		}
		hist.polyline([marginLeft,0,marginLeft,histHeight-marginBottom,histWidth,histHeight-marginBottom]).attr({
			stroke:"black",
			'stroke-width':"1px",
			fill:"none"
			
		})
}


function drawHomeMultiBar(Omics,counts,histADiv,histWidth,histHeight,colorPool,Ylab){
	
	var marginBottom=100
	var marginLeft=70
	var barVsBlank=2
	var hist=SVG(histADiv)
		hist.attr({"text-rendering":"optimizeLegibility"})
		//hist.width(histWidth+"px")
		//hist.height(histHeight+"px")
		hist.width("100%")
		hist.height("100%")
		hist.viewbox(0,0,histWidth,histHeight)

		var countTags=[]
		var countsSum=[]
		for(var c=0;c<Omics.length;c++){
				countTags.push([counts.a[c],counts.b[c]])
				countsSum[c]=counts.a[c]+counts.b[c]

		}
		var scale=(histHeight-marginBottom)/Math.max.apply(null,countsSum)*0.7
		var barWidth=(histWidth-marginLeft)/Omics.length*barVsBlank/(1+barVsBlank)
		for(var O=0;O<Omics.length;O++){

			var g=hist.group()
			var r1=hist.rect(barWidth,scale*counts.a[O]).attr({
				x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2,
				y:histHeight-marginBottom-scale*counts.a[O],
				fill:colorPool[O%colorPool.length],
				stroke:"none",
				"stroke-width":"1px",
				class:"Hist"
				
			})
			
			var tagT1=hist.plain("Number of cell lines").attr({
				x:(histWidth+marginLeft)/2,
				y:(histHeight-marginBottom)*0.15-8,
				fill:colorPool[O%colorPool.length],
				stroke:"none",
				"font-size":16,
				"text-anchor":"middle",
				"dominant-baseline":  "text-after-edge",
				opacity:0
			})
			var tagT2=hist.plain("with drug screen:").attr({
				x:(histWidth+marginLeft)/2,
				y:(histHeight-marginBottom)*0.15+8,
				fill:colorPool[O%colorPool.length],
				stroke:"none",
				"font-size":16,
				"text-anchor":"middle",
				"dominant-baseline":  "text-after-edge",
				opacity:0
			})
			var tagN=hist.plain(counts.a[O]-counts.b[O]).attr({
				x:(histWidth+marginLeft)/2,
				y:(histHeight-marginBottom)*0.15+24,
				fill:colorPool[O%colorPool.length],
				stroke:"none",
				"font-size":16,
				"text-anchor":"middle",
				"dominant-baseline":  "text-after-edge",
				"font-weight":"bold",
				opacity:0
			})
			g.add(r1)
			g.add(tagT1)
			g.add(tagT2)
			g.add(tagN)
			var g=hist.group()
			var r2=hist.rect(barWidth,scale*counts.b[O]).attr({
				x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2,
				y:histHeight-marginBottom-scale*counts.a[O]-scale*counts.b[O],
				fill:enhanceCol(colorPool[O%colorPool.length],0.4),
				stroke:"none",
				"stroke-width":"1px",
				class:"Hist"
				
			})
			
			var tagT1=hist.plain("Number of cell lines").attr({
				x:(histWidth+marginLeft)/2,
				y:(histHeight-marginBottom)*0.15-8,
				fill:enhanceCol(colorPool[O%colorPool.length],0.4),
				stroke:"none",
				"font-size":16,
				"text-anchor":"middle",
				"dominant-baseline":  "text-after-edge",
				opacity:0
			})
			var tagT2=hist.plain("without drug screen:").attr({
				x:(histWidth+marginLeft)/2,
				y:(histHeight-marginBottom)*0.15+8,
				fill:enhanceCol(colorPool[O%colorPool.length],0.4),
				stroke:"none",
				"font-size":16,
				"text-anchor":"middle",
				"dominant-baseline":  "text-after-edge",
				opacity:0
			})
			var tagN=hist.plain(counts.b[O]).attr({
				x:(histWidth+marginLeft)/2,
				y:(histHeight-marginBottom)*0.15+24,
				fill:enhanceCol(colorPool[O%colorPool.length],0.4),
				stroke:"none",
				"font-size":16,
				"text-anchor":"middle",
				"dominant-baseline":  "text-after-edge",
				"font-weight":"bold",
				opacity:0
			})
			g.add(r2)
			g.add(tagT1)
			g.add(tagT2)
			g.add(tagN)
			/*
			var bg=hist.rect(barWidth,histHeight-marginBottom).attr({
				x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2,
				y:0,
				fill:"transparent",
				stroke:"transparent"
				
			})
			var tag=hist.plain(toThousands(countTags[O][0])).attr({
				x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2+barWidth/2,
				y:histHeight-marginBottom-scale*countsSum[O],
				fill:"grey",
				stroke:"none",
				"font-size":barWidth/4,
				"text-anchor":"middle",
				"dominant-baseline":  "text-after-edge"
			})

			g.add(bg)
			g.add(tag)
			*/
			hist.plain(Ylab).attr({
				x:marginLeft-16,
				y:(histHeight-marginBottom)/2,
				fill:"black",
				stroke:"none",
				"font-size":16,
				"text-anchor":"middle",
				"dominant-baseline":  "central",
				"transform":"rotate(-90 "+(marginLeft-16)+","+(histHeight-marginBottom)/2+")"
			})
			
			hist.plain(Omics[O]).attr({
				x:marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2+barWidth/2,
				y:histHeight-marginBottom+4,
				fill:"black",
				stroke:"none",
				"font-size":10,
				"text-anchor":"end",
				"dominant-baseline":  "central",
				"transform":"rotate(-45 "+(marginLeft+O*(barWidth*(1+barVsBlank)/barVsBlank)+(barWidth/(1+barVsBlank))/2+barWidth/2)+","+(histHeight-marginBottom+4)+")"
			})
			
		}
		hist.polyline([marginLeft,0,marginLeft,histHeight-marginBottom,histWidth,histHeight-marginBottom]).attr({
			stroke:"black",
			'stroke-width':"1px",
			fill:"none"
			
		})
}

function enhanceCol(str,por){
	if(por>1)return(str)
	var r=str.substr(1,2)
	var g=str.substr(3,2)
	var b=str.substr(5,2)
	var adjusted=[]
	for(var i=0;i<3;i++){
		var c=str.substr(1+i*2,2)
		var ten=parseInt(c,16)
		var re=255-ten
		adjusted.push(parseInt(ten+re*por).toString(16))
	}
	return("#"+adjusted.join(""))
	
}


function initHelp(){
	if(helpDone)return
	var help_main=document.getElementById("NATIVE,HELP")
	help_main.innerHTML=""
	var xmlhttp=new XMLHttpRequest();
	xmlhttp.onreadystatechange=function(){
		//document.getElementById("plotarea").innerHTML =genename+"got:"+xmlhttp.responseText;
		if (xmlhttp.readyState==4 && xmlhttp.status==200){
			help_str =xmlhttp.responseText.toString();
			help_obj = JSON.parse(help_str);
			var helpCatalog_div=document.createElement("div")
				helpCatalog_div.setAttribute("class","helpCatalog_div")
			help_main.append(helpCatalog_div)
			var index=document.createElement("p")
				index.innerHTML="INDEX"
				index.setAttribute("style","font-size:25px;font-weight:bold;color:var(--morelightBG);margin-left:auto;margin-right:auto;margin-top:30px;")
			helpCatalog_div.appendChild(index)
			parseHelpCata(help_obj,helpCatalog_div,[])
			
			var helpContent=document.createElement("div")
				helpContent.setAttribute("class","helpContent")
			help_main.append(helpContent)
			var title=document.createElement("p")
				title.innerHTML="Help Document for CLIPA"
				title.setAttribute("style","font-size:25px;font-weight:bold;margin-left:auto;margin-right:auto;margin-top:30px;")
			helpContent.appendChild(title)
			parseHelpContent(help_obj,helpContent,[])
			helpDone=true
		}else{
			
		}
	}
	//document.getElementById("plotarea").innerHTML= +"start";
	xmlhttp.open("GET","php/getHelp.php",true);
	xmlhttp.send();
}

function parseHelpCata(help_tmp,div,heading){
	var keys=Object.keys(help_tmp).sort()
	var folderCount=0
	for(var k=0;k<keys.length;k++){
		if(keys[k]=="00_title"){
			var p=document.createElement("p")
				if(heading.length==1){
					var heading_flat=heading[0]+"."
				}else{
					var heading_flat=heading.join(".")
				}
				p.innerHTML='<span style="font-weight:bold">'+heading_flat+"</span>  "+help_tmp[keys[k]]
				p.setAttribute("class","helpCatalogTitle")
			var a=document.createElement("a")
				a.setAttribute("href","#"+"helpID_"+heading_flat)
			a.appendChild(p)
			div.appendChild(a)
		}else if(keys[k].indexOf("_plain")>-1){
		}else if(keys[k].indexOf("_italic")>-1){
		}else if(keys[k].indexOf("_bold")>-1){
		}else if(keys[k].indexOf("_table")>-1){
		}else if(keys[k].indexOf("_png")>-1){
			
		}else{
			folderCount++
			var subdiv=document.createElement("div")
			subdiv.setAttribute("class","helpTitleDiv")
			var newheading=JSON.parse(JSON.stringify(heading))
			newheading.push(folderCount)
			parseHelpCata(help_tmp[keys[k]],subdiv,newheading)
			div.appendChild(subdiv)
		}
	}
	
}


function parseHelpContent(help_tmp,div,heading){
	var toItalics=["GDSC","DepMap","NCI60","GSVA","HGNC","MSigDB","IlluminaHumanMethylation450kanno.ilmn12.hg19","pubchem","Cell Model Passports"]
	var keys=Object.keys(help_tmp).sort()

	var folderCount=0
	for(var k=0;k<keys.length;k++){
		if(keys[k]=="00_title"){
			var p=document.createElement("p")
				if(heading.length==1){
					var heading_flat=heading[0]+"."
				}else{
					var heading_flat=heading.join(".")
				}
				p.innerHTML=heading_flat+"  "+help_tmp[keys[k]]
				p.setAttribute("class","helpContentTitle")
				p.style.fontSize=(12+12/heading.length)+"px"
				p.setAttribute("id","helpID_"+heading_flat)
			div.appendChild(p)
		}else if(keys[k].indexOf("_plain")>-1){
			
			var plains=help_tmp[keys[k]]
			
			for(var pl=0;pl<plains.length;pl++){
				for(var toIt=0;toIt<toItalics.length;toIt++){

					plains[pl]=plains[pl].replace(new RegExp(toItalics[toIt],"gm"),'<span style="font-style:italic">'+toItalics[toIt]+'</span>')

				}
				var p=document.createElement("p")
					p.innerHTML=plains[pl]
					p.setAttribute("class","helpContentPlain")
				div.appendChild(p)
			}
		}else if(keys[k].indexOf("_italic")>-1){
			var italic=help_tmp[keys[k]]
			for(var pl=0;pl<italic.length;pl++){
				var p=document.createElement("p")
					p.innerHTML=italic[pl]
					p.setAttribute("class","helpContentItalic")
				div.appendChild(p)
			}
		}else if(keys[k].indexOf("_bold")>-1){
			
			var bold=help_tmp[keys[k]]
			
			for(var pl=0;pl<bold.length;pl++){
				var p=document.createElement("p")
					p.innerHTML=bold[pl]
					p.setAttribute("class","helpContentBold")
				div.appendChild(p)
			}
		}else if(keys[k].indexOf("_table")>-1){
			var tables=help_tmp[keys[k]]
			var tabTitle=tables[0].split("\t")
			var txttmp="<table><tr>"
			for(var t=0;t<tabTitle.length;t++){
				txttmp+='<th>'+tabTitle[t]+'</th>'
			}
			txttmp+='</tr>'
			
			for(var tabRow=1;tabRow<tables.length;tabRow++){
				var onerow=tables[tabRow].split("\t")
				txttmp+='<tr>'
				for(var acell=0;acell<onerow.length;acell++){
					txttmp+='<td>'+onerow[acell]+'</td>'
				}
				txttmp+='</tr>'
			}
			txttmp+='</table>'
			console.log(txttmp)
			var table=document.createElement("table")
				table.innerHTML=txttmp
				table.setAttribute("class","helpContentTable")
			div.appendChild(table)
		}else if(keys[k].indexOf("_png")>-1){
			var img=document.createElement("img")
				img.setAttribute("src","/CLIPA/help/"+help_tmp[keys[k]])
				var imgwidth=help_tmp[keys[k]].split("_")[1]
				var scale=2.5
				img.style.width=imgwidth/scale+"px"
			div.appendChild(img)
		}else{
			folderCount++
			var subdiv=document.createElement("div")
			subdiv.setAttribute("class","helpContentDiv")
			var newheading=JSON.parse(JSON.stringify(heading))
			newheading.push(folderCount)
			parseHelpContent(help_tmp[keys[k]],subdiv,newheading)
			div.appendChild(subdiv)
		}
	}
}

function updateDownload(){
	var downloadblock_div=document.getElementsByClassName("downloadblock_div")[0]
	downloadblock_div.innerHTML=""
	var databaseNames=Object.keys(initParams.dl)
	for(var d=0;d< databaseNames.length;d++){
		var downloadblock=document.createElement("div")
			downloadblock.setAttribute("class","downloadblock")
		var p=document.createElement("p");
		p.innerHTML="Downloads for "+databaseNames[d]
		downloadblock.appendChild(p)
		var table=document.createElement("table")
		var files=Object.keys(initParams.dl[databaseNames[d]])
		var tableTmp=""
		for(var f=0;f<files.length;f++){
			tableTmp+='<tr><td><a href=/CLIPA/download/'+databaseNames[d]+'/'+files[f]+'>'+files[f]+'</a></td><td>'+initParams.dl[databaseNames[d]][files[f]].date+'</td><td>'+initParams.dl[databaseNames[d]][files[f]].size+'</td></tr>'
		}
		table.innerHTML=tableTmp
		downloadblock.appendChild(table)
		downloadblock_div.appendChild(downloadblock)
	}
	
}

function updateVisits(){
	
	xmlhttp=new XMLHttpRequest();
	//var Debug=document.getElementById("NATIVE,HELP");
	xmlhttp.onreadystatechange=function(){
		//document.getElementById("plotarea").innerHTML =genename+"got:"+xmlhttp.responseText;
		if (xmlhttp.readyState==4 && xmlhttp.status==200){
			var visitRecords=JSON.parse(xmlhttp.responseText);
			var visitRecords_div=document.getElementsByClassName("visitRecords")[0]
				var p1=document.createElement("p")
					p1.innerHTML='Total Visits: <span style="font-weight:bold;color:red;">'+visitRecords.visitips.times+'</span> times <span style="color:grey;font-style:italic;">since '+visitRecords.visitips.from+"</span>"
				var p2=document.createElement("p")
					p2.innerHTML='Total Plots: <span style="font-weight:bold;color:red;">'+visitRecords.plotips.times+'</span> times <span style="color:grey;font-style:italic;">since '+visitRecords.plotips.from+"</span>"
			visitRecords_div.appendChild(p1)
			visitRecords_div.appendChild(p2)
		}
			
	}

	xmlhttp.open("POST","php/getvisitRecords.php",true);
	xmlhttp.setRequestHeader("Content-Type","application/x-www-form-urlencoded");
	xmlhttp.send();
}