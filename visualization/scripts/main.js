require.config({
	paths: {
		'Kinetic': "./kinetic-v5.1.0.min",
		"jquery" : "./jquery-2.1.3.min",
		"ui":"./jquery-ui.min"
	}
});

require(['Kinetic', 'jquery', 'ui'], function(Kinetic, jquery){
	
    $( "#tabs" ).tabs();
	if (window.File && window.FileReader && window.FileList && window.Blob) {
		// Great success! All the File APIs are supported.
	} else {
		alert('The File APIs are not fully supported in this browser.');
	}
	var matrix = {};
	var matrixKeys = [];
	var stage;
	var layer;
	var canvasWidth = 1024;
	var canvasHeight = 768;
	var margin = [80, 10, 0, 100];
	var currentKey = "";
	var Nucs = ["A","C","G","T"];
	var NucColors = ["green","red","blue","orange"];
	var CurrentViewData = {};
	var yScaleFactor = 1;
	var mode = 0; //0: display full model, 1:
	var inspectSeqArray = [];
	//var canvas = document.getElementById('canavs');
	stage = new Kinetic.Stage({
	    container: 'canvas',
		width: canvasWidth,
		height: canvasHeight
	});
	$("#yFactor").text(yScaleFactor);
	$("#yFactorSlider").slider({
		max: 3.1,
		min:0.1,
		step:0.1,
		value: yScaleFactor,
		slide:function( event, ui ) {
			yScaleFactor = ui.value;
			$("#yFactor").text(yScaleFactor);
		},
		change: function( event, ui ) {
			yScaleFactor = ui.value;
			$("#yFactor").text(yScaleFactor);
			if(currentKey){
				drawMatrices(currentKey);
			}
		}
	});

	function handleFileSelect(evt) {
		evt.stopPropagation();
		evt.preventDefault();
		//matrix = {};
		//$("#lists").empty();

		var files = evt.dataTransfer.files; // FileList object.
		// files is a FileList of File objects. List some properties.
		var output = [];
		for (var i = 0, f; f = files[i]; i++) {
			//output.push('<li><strong>', escape(f.name), '</strong> (', f.type || 'n/a', ') - ', f.size, ' bytes, last modified: ', f.lastModifiedDate ? f.lastModifiedDate.toLocaleDateString() : 'n/a','</li>');
			if(f.type == ""){
			}else if (!f.type.match('text.*')){
	    		continue;
	    	}
			var reader = new FileReader();
			reader.onload = (function(file, index, total){
				var key;
				if(matrix[file.name] != undefined){
					var c = 1;
					key = file.name + "("+c+")";
					while(matrix[key] !== undefined){
						key = file.name + "("+c+")";
						c ++;
					}
				}else{
					key = file.name;
				}
				matrix[key] = {AF: {}, SF:{}, s: 0, intercept: 0, motifLength : 0};
				$("#lists").append("<span class='matrixName' title='"+key+"'>"+key+"</span><br>");
				return function(e){
					var text = e.target.result;
					var lines = text.split("\n");
					var sfCount = 0;
					for(var j=0, l; l = lines[j];j++){
						var elements = l.split("\t");
						if(elements[0] == "(Intercept)"){
							matrix[key]["intercept"] = elements[1];
						}else if(codedSF = elements[0].match(/^Z(\d+)/)){
							//Single-base feature
							var pos = parseInt((codedSF[1] - 1) / 3) +1;
							var nuc = Nucs[(codedSF[1]-1)%3];
							var SFterm = nuc + pos;
							matrix[key]["SF"][SFterm] = parseFloat(elements[1]);
							sfCount ++;
							if(nuc == "G"){
								matrix[key]["SF"]["T"+pos] = -1 * (matrix[key]["SF"]["A"+pos] + matrix[key]["SF"]["C"+pos] + matrix[key]["SF"]["G"+pos]);
								sfCount ++;
							}
						}else if(elements[0].match(/([ATCG])(\d+)/g)){
							//Association features
							if(elements[1] != 0){
								matrix[key]["AF"][elements[0]] = parseFloat(elements[1]);
							}
						}else if(elements[0] == "s"){
							matrix[key]["s"] = parseFloat(elements[1]);
						}
					}
					matrix[key]["motifLength"] = sfCount /4 ;
					//drawMatrices(file.name);
					/*					
					if(index === total -1){
						//After reading all files, draw the first matrix
						drawMatrices(0);
						$($("#lists").find(".matrixName").get(0)).css("font-weight","bold");
						$("#controls").html("<span class='hideInter' style='margin:0 5px;text-decoration:underline;cursor:pointer;'>Hide Interactions</span>");
					}
					*/
				};
			}(f,i, files.length));
			reader.readAsText(f);
			$(evt.target).removeClass("dragEnter");
		}
	}
	function handleDragOver(evt) {
		evt.stopPropagation();
		evt.preventDefault();
		$(evt.target).addClass("dragEnter");
		evt.dataTransfer.dropEffect = 'copy'; // Explicitly show this is a copy.
	}
	function handleDragLeave(evt){
		$(evt.target).removeClass("dragEnter");
	}
	// Setup the dnd listeners.
	var dropZone = document.getElementById('drop_zone');
	dropZone.addEventListener('dragover', handleDragOver, false);
	dropZone.addEventListener('drop', handleFileSelect, false);
	dropZone.addEventListener('dragleave', handleDragLeave, false);
	function drawMatrices(key){
		var i, j, k, l;
		var monoTerms = [];
		var interTerms = {};
		var maxValue = null;
		currentKey = key;
		var motifLength = matrix[key]["motifLength"];
		var SFs = Object.keys(matrix[key]["SF"]).sort(sortTerms);
		var AFs = Object.keys(matrix[key]["AF"]).sort(sortTerms);
		var CoefToSpecificity = {};
		if(layer){
			layer.destroyChildren();
			layer.destroy();
		}

		layer  = new Kinetic.Layer();
		CurrentViewData = {
			layer: layer,
			SFTerms: [], //Use position as the index (start from 0)
			AFlinks: {}  
		};
		var background = new Kinetic.Rect({
			x: 0, y: 0, width: canvasWidth, height:canvasHeight, fill:'white'
		});
		layer.add(background);
		title = new Kinetic.Text({
			text: "The ELRM for "+key,
			fontSize: 18,
			padding: 0,
			fill:"black",
			align:"center"
		});
		title.position({
			x: canvasWidth/2 - title.getWidth()/2,
			y: 2
		});
		layer.add(title);
		var hr = new Kinetic.Line({
			points: [0,canvasHeight*2/3, canvasWidth,canvasHeight*2/3],
			stroke: 'silver',
			strokeWidth: 1
		});
		layer.add(hr);
		hr = new Kinetic.Line({
			points: [0,canvasHeight*2/3 + 20, canvasWidth,canvasHeight*2/3 + 20],
			stroke: 'silver',
			strokeWidth: 1
		});
		layer.add(hr);

		hr = new Kinetic.Line({
			points: [0,canvasHeight-80, canvasWidth,canvasHeight-80],
			stroke: 'silver',
			strokeWidth: 1
		});
		layer.add(hr);
		
		var tempText = new Kinetic.Text({
			text: "Association feature",
			fontSize: 16,
			padding: 0,
			fill:"black",
			align:"center",
			rotation: -90
		});
		tempText.offsetX(tempText.getWidth()/2);
		tempText.offsetY(tempText.getHeight()/2);
		tempText.position({
			x: 50,
			y: ((canvasHeight - margin[0])/3)  + margin[0]
		});
		layer.add(tempText);
		var tempText = new Kinetic.Text({
			text: "Position",
			fontSize: 16,
			padding: 0,
			fill:"black",
			align:"center"
		});
		tempText.position({
			x: 5,
			y: ((canvasHeight - margin[1] - margin[2])*2/3) + (tempText.getHeight()/2)
		});
		layer.add(tempText);
		var tempText = new Kinetic.Text({
			text: "Independent\nsingle-base feature",
			fontSize: 16,
			padding: 0,
			fill:"black",
			align:"center",
			rotation: -90
		});
		tempText.offsetX(tempText.getWidth()/2);
		tempText.offsetY(tempText.getHeight()/2);
		tempText.position({
			x: 50,
			y: ((canvasHeight - margin[0])*2/3) + tempText.getWidth()/2  + margin[0] 
		});
		layer.add(tempText);

		var tempText = new Kinetic.Text({
			text: "Irrelavent\nfeature",
			fontSize: 16,
			padding: 0,
			fill:"black",
			align:"center",
			rotation:-90
		});
		tempText.offsetX(tempText.getWidth()/2);
		tempText.offsetY(tempText.getHeight()/2);
		tempText.position({
			x: 50,
			y: (canvasHeight - margin[1] - margin[2] - 80  + 25 + tempText.getWidth()/2)
		});
		layer.add(tempText);

		var eachNucWidth = (canvasWidth- margin[3] - margin[1])/motifLength;

		for(i=0; i < motifLength;i++){
			var temp = new Kinetic.Line({
				points: [ margin[3]+eachNucWidth*i, margin[0],  margin[3]+eachNucWidth*i, canvasHeight - margin[2]],
				stroke: 'silver',
				strokeWidth:0.5
			});
			layer.add(temp);
			if(i >= 0){
				var temp = new Kinetic.Text({
					text: i+1,
					fontSize:16,
					fill:"black",
					padding:5,
					fontStyle:"bold"
				});
				temp.position({
					x: margin[3] + eachNucWidth*(i+0.5) - temp.getWidth()/2,
					y: ((canvasHeight - margin[1] - margin[2])*2/3) + 5
				});
				layer.add(temp);
			}
		}
		
		for(var i = 0; i < AFs.length; i ++){
			var SFinAF = AFs[i].split(":");
			var AFcoef = matrix[key]["AF"][AFs[i]];
			for(var j = 0; j < SFinAF.length; j++){
				temp = SFinAF[j].match(/([ATCG])(\d+)/);
				var nuc = temp[1];
				var motifPos = temp[2];
				var thePos = {
					x: margin[3] + eachNucWidth *(motifPos-1),
					y: margin[0] + canvasHeight*1/3 + 20
				};
				if(CurrentViewData.SFTerms[motifPos - 1] == undefined){
					CurrentViewData.SFTerms[motifPos - 1] = {};
				}
				if(CurrentViewData.SFTerms[motifPos - 1][nuc] == undefined){
					var SFLabel = createSFLabel(thePos, SFinAF[j], matrix[key]["SF"][SFinAF[j]]);
					CurrentViewData.SFTerms[motifPos - 1][nuc] = SFLabel;
					//Check if there are other SF terms of the same position
					var nucsAtSamePos = Object.keys(CurrentViewData.SFTerms[motifPos - 1]);
					for(var t = 0; t < nucsAtSamePos.length; t++){
						var theLabel = CurrentViewData.SFTerms[motifPos - 1][nucsAtSamePos[t]];
						theLabel.position({
							x: margin[3] + eachNucWidth *(motifPos-0.5) - theLabel.getWidth()/2 + t*15,
							y: canvasHeight*1/3 - theLabel.getHeight()/2 + t * 10
						});
					}
					layer.add(SFLabel);
				}
			}

			//Draw association links
			for(var j = 0; j < SFinAF.length; j++){
				temp1 = SFinAF[j].match(/([ATCG])(\d+)/);
				if(j+1 == SFinAF.length){
					if(SFinAF.length > 2){
						temp2 = SFinAF[0].match(/([ATCG])(\d+)/);
						var AFkey = SFinAF[0]+"_"+SFinAF[j];
					}else{
						continue;
					}
				}else{
					temp2 = SFinAF[j+1].match(/([ATCG])(\d+)/);
					var AFkey = SFinAF[j]+"_"+SFinAF[j+1];
				}
				
				var nuc1 = temp1[1];
				var motifPos1 = temp1[2];
				var nuc2 = temp2[1];
				var motifPos2 = temp2[2];
				
				var node1 = CurrentViewData.SFTerms[motifPos1 - 1][nuc1];
				var node2 = CurrentViewData.SFTerms[motifPos2 - 1][nuc2];
				var existed = 0;
				if(CurrentViewData.AFlinks[AFkey] == undefined){
					CurrentViewData.AFlinks[AFkey] = [];
					existed = 0;
				}else{
					existed = CurrentViewData.AFlinks[AFkey].length;
				}
				var sceneFunction = (function(node1, node2, AFcoef){
					var node1Pos = node1.position();
					var node2Pos = node2.position();
					var controlPoint, startX, startY, endX, endY;
					startX = node1Pos.x + node1.getWidth()/2;
					endX = node2Pos.x + node2.getWidth()/2;
					if(AFcoef > 0){
						startY = node1Pos.y;
						endY = node2Pos.y;
						controlPoint = {
							x: (node1Pos.x + (node1.getWidth()/2) + node2Pos.x + (node2.getWidth()/2) )/2,
							y: (canvasHeight/3) * (1 - yScaleFactor*((Math.abs(node2Pos.x - node1Pos.x)/eachNucWidth)/(motifLength-1))) - existed*20
						};

					}else{
						startY = node1Pos.y + node1.getHeight();
						endY = node2Pos.y + node2.getHeight();
						controlPoint = {
							x: (node1Pos.x + (node1.getWidth()/2) + node2Pos.x + (node2.getWidth()/2) )/2,
							y: (canvasHeight/3) * (1 + yScaleFactor*((Math.abs(node2Pos.x - node1Pos.x)/eachNucWidth)/(motifLength-1))) + existed*20
						};
					}
					return function(context){
						context.beginPath();
						context.moveTo(startX, startY);
						context.quadraticCurveTo(controlPoint.x, controlPoint.y, endX, endY);
						//context.closePath();
						//KineticJS specific context method
						context.fillStrokeShape(this);
					};
				})(node1, node2, AFcoef);
				var associationLink = new Kinetic.Shape({
					sceneFunc: sceneFunction,
					stroke: coef2ColorLink(AFcoef),
					strokeWidth: coef2LinkWidth(AFcoef)
				});
				associationLink.setZIndex(2);
				CurrentViewData.AFlinks[AFkey].push(associationLink);
				layer.add(associationLink);
			}
		}
		//Draw single-base features
		for(motifPos = 1; motifPos <= motifLength; motifPos++){
			var posNuc = {};
			var negNuc = {};
			var nouseSFHash = [];
			for(i = 0; i < 4; i++){
				nuc = Nucs[i];
				tempTerm = nuc + motifPos;
				tempCoef = matrix[key]["SF"][tempTerm];
				if(tempCoef > 0){
					posNuc[nuc] = tempCoef;
				}else if(tempCoef < 0){
					negNuc[nuc] = tempCoef;
				}else{
					nouseSFHash.push(nuc);
				}
			}
			var posSFs = Object.keys(posNuc).sort(function(a, b){
				return posNuc[b] - posNuc[a];
			});
			var negSFs = Object.keys(negNuc).sort(function(a, b){
				return negNuc[a] - negNuc[b];
			});
			if(CurrentViewData.SFTerms[motifPos - 1] == undefined){
				CurrentViewData.SFTerms[motifPos - 1] = {};
			}
			var tempCount = 0;
			for(i = 0; i < posSFs.length; i ++){
				nuc = posSFs[i];
				SFTerm = nuc+ motifPos;
				var SFcoef = matrix[key]["SF"][SFTerm];
				var thePos = {
					x: margin[3] + eachNucWidth *(motifPos-0.5),
					y: (canvasHeight*2/3) + 25
				};

				if(CurrentViewData.SFTerms[motifPos - 1][nuc] == undefined){
					var SFLabel = createSFLabel(thePos, SFTerm, SFcoef);
					SFLabel.position({
						x: thePos.x - SFLabel.getWidth()/2,
						y: thePos.y + tempCount*SFLabel.getHeight()
					});
					CurrentViewData.SFTerms[motifPos - 1][nuc] = SFLabel;
					//SFLabel.getText().setSize({width: eachNucWidth});
					layer.add(SFLabel);
					tempCount ++;
				}
			}
			var tempCount = 0;
			for(i = 0; i < negSFs.length; i ++){
				nuc = negSFs[i];
				SFTerm = nuc+ motifPos;
				var SFcoef = matrix[key]["SF"][SFTerm];
				var thePos = {
					x: margin[3] + eachNucWidth *(motifPos-0.5),
					y: canvasHeight - 80 - margin[2] - 10
				};
				if(CurrentViewData.SFTerms[motifPos - 1][nuc] == undefined){
					var SFLabel = createSFLabel(thePos, SFTerm, SFcoef);
					SFLabel.position({
						x: thePos.x - SFLabel.getWidth()/2,
						y: thePos.y - SFLabel.getHeight() - tempCount*SFLabel.getHeight()
					});
					CurrentViewData.SFTerms[motifPos - 1][nuc] = SFLabel;
					//SFLabel.getText().setSize({width: eachNucWidth});
					layer.add(SFLabel);
					tempCount ++;
				}
			}
			var tempCount = 0;

			for(i = 0; i < nouseSFHash.length; i++){
				if(CurrentViewData.SFTerms[motifPos - 1][nouseSFHash[i]] == undefined){
					var thePos = {
						x: margin[3] + eachNucWidth *(motifPos-0.5),
						y: canvasHeight - 80 - margin[2] + 5
					};
					var SFLabel = createSFLabel(thePos, nouseSFHash[i]+motifPos, 0);
					SFLabel.position({
						x: thePos.x - SFLabel.getWidth()/2,
						y: thePos.y + tempCount * SFLabel.getHeight()
					});
					CurrentViewData.SFTerms[motifPos - 1][nouseSFHash[i]] = SFLabel;
					layer.add(SFLabel);
					tempCount ++;
				}
			}
		}

		var SFtoSpecificity = {};
		for(i = 0; i < SFs.length; i++){
			var sfterm = SFs[i];
			totalCoef = 0;
			totalCoef += matrix[currentKey]["SF"][sfterm];
			for(j = 0; j < AFs.length; j++){
				
				//if(afterm.indexOf(sfterm) >= 0){
				if(checkSFinAF(sfterm, AFs[j], ":")){
					totalCoef += matrix[currentKey]["AF"][AFs[j]];
				}
			}
			SFtoSpecificity[sfterm] = totalCoef;
		}
		var PosSFs2Spec = Object.keys(SFtoSpecificity).sort(function(a, b){
			return parseFloat(SFtoSpecificity[b]) - parseFloat(SFtoSpecificity[a]);
		});
		
		var PosToSpec = [];
		var NegToSpec = [];
		for(i = 0; i < PosSFs2Spec.length; i++){
			if(SFtoSpecificity[PosSFs2Spec[i]] > 0){
				PosToSpec.push(PosSFs2Spec[i] + "("+parseInt(1000*SFtoSpecificity[PosSFs2Spec[i]])/1000+")");
			}
		}
		var NegSFs2Spec = PosSFs2Spec.reverse();
		for(i = 0; i < NegSFs2Spec.length; i++){
			if(SFtoSpecificity[NegSFs2Spec[i]] < 0){
				NegToSpec.push(NegSFs2Spec[i] + "("+parseInt(1000*SFtoSpecificity[NegSFs2Spec[i]])/1000+")");
			}
		}
		/*
		var text = "Positive contributions: "+ PosToSpec.join(",") + "\nNegative contributions: "+NegToSpec.join(",");
		var InvestigationText = new Kinetic.Text({
			x: 10, y: 25,
			text: text,
			fill: "black",
			fontSize: 16
		});
		if(CurrentViewData.investigationText){
			CurrentViewData.investigationText.destroy();
		}		
		CurrentViewData.layer.add(InvestigationText);
		CurrentViewData.investigationText = InvestigationText;
		*/

		console.log(PosToSpec, NegToSpec);
		//Attach event handler
		for(i = 1; i <= CurrentViewData.SFTerms.length; i++){
			for(j = 0; j < 4; j++){
				var sf = Nucs[j] + i;
				CurrentViewData.SFTerms[i-1][Nucs[j]].on("click", (function(pos, nuc){
					return function(event){
						//highlight the all Associated single-base features;
						event.evt.stopImmediatePropagation();
						event.evt.preventDefault();
						highlightAssociateLinks(pos, nuc);
					}
				}(i, Nucs[j])));
				//CurrentViewData.SFTerms[i-1][Nucs[j]].on("mouseleave", resetHighlight);
			}
		}
		stage.add(layer);
		
	}
	
	function createSFLabel(Point, label, SFcoef){
		var SFLabel = new Kinetic.Label({
			x: Point.x,
			y: Point.y
		});
		var SFText = new Kinetic.Text({
			text: label,
			fontSize: 18,
			padding: 0,
			fill:"black",
			align:"center",
			width:40
		});
		var SFTag = new Kinetic.Tag({
			stroke: "#333",
			fill: coef2Color(SFcoef)
		});
		SFLabel.add(SFTag).add(SFText);
		return SFLabel;
	}
	function unbindAllEvents(){
		for(sf in CurrentViewData.SFTerms){
			CurrentViewData.SFTerms[sfterm].off("mouseover");
			CurrentViewData.SFTerms[sfterm].off("mouseout");
		}
		return true;
	}
	function createToolTip(Point, Head,Content){
		var layer = CurrentViewData.layer;
		var tooltip = new Kinetic.Group({
			x: Point.x,
			y: Point.y,
		});

		var tooltipTag = new Kinetic.Tag({
			width: 250,
			fill: "#fff",
			opacity: 0.6,
			stroke: 'silver',
			strokeWidth: 1,
			pointerDirection: "up",
			pointerWidth: 15,
			pointerHeight: 10,
			shadowColor: 'black',
	        shadowBlur: 10,
	        shadowOffset: {x: 5, y:5},
	        shadowOpacity: 0.2,
		});
		var ContentHead = new Kinetic.Text({
			text: Head,
			fill: "black",
			fontStyle: "bold",
			fontSize: 18,
			padding:10
		});
		var ContentText = new Kinetic.Text({
			text: Content,
			fill: "black",
			padding: 10,
			fontSize: 16
		});
		tooltip.add(tooltipTag).add(ContentHead).add(ContentText);
		
		var newX = Point.x - tooltipTag.getWidth()/2;
		var newY = Point.y + 5;
		
		tooltipTag.height(ContentText.getHeight() + ContentHead.getHeight() + 20);

		if(newY + tooltipTag.getHeight() > canvasHeight){
			newY = Point.y  - 30- tooltipTag.getHeight(); 
			tooltipTag.setPointerDirection("down");
		}
		if(newX + tooltipTag.getWidth() > canvasWidth){
			newX = canvasWidth - tooltipTag.getWidth() - 10;
			tooltipTag.setPointerDirection("none");
		}
		ContentText.position({
			x: ContentText.position().x,
			y: ContentText.position().y + 20
		});
		tooltip.position({
			x: newX,
			y: newY
		});
		tooltip.on("click", resetHighlight);
		return tooltip;
	}
	function resetHighlight(){
		drawMatrices(currentKey);
	}
	function checkSFinAF(SF, AF, sep){
		var isIn = false;
		var i=0, temp = AF.split(sep);
		for(i = 0; i < temp.length; i++){
			if(temp[i] == SF){
				isIn = true;
				break;
			}
		}
		return isIn;
	}
	function highlightAssociateLinks(pos, nuc){
		var sfTerm = nuc+pos;
		var sfTerms = [];
		var associatedSFs = {};

		for(i = 1; i <= CurrentViewData.SFTerms.length; i++){
			for(j = 0; j < 4; j ++){
				if(pos == i && nuc == Nucs[j]){
					CurrentViewData.SFTerms[i-1][Nucs[j]].opacity(1);
					CurrentViewData.SFTerms[i-1][Nucs[j]].setZIndex(100);
				}else{
					CurrentViewData.SFTerms[i-1][Nucs[j]].opacity(0.1);
					CurrentViewData.SFTerms[i-1][Nucs[j]].setZIndex(-1);
				}
			}
		}
		for(aflink in CurrentViewData.AFlinks){
			if(checkSFinAF(sfTerm, aflink, "_")){
			//if(aflink.indexOf(sfTerm)>= 0){
				for(i = 0; i < CurrentViewData.AFlinks[aflink].length; i++){
					CurrentViewData.AFlinks[aflink][i].opacity(1);
					CurrentViewData.AFlinks[aflink][i].setZIndex(99);
				}
				//CurrentViewData.valueLabels[iT].hide();
				sfTerms = aflink.split(/\_/);
				for(i = 0; i < sfTerms.length;i ++){
					associatedSFs[sfTerms[i]] = 1;
				}
			}else{
				for(i = 0; i < CurrentViewData.AFlinks[aflink].length; i++){
					CurrentViewData.AFlinks[aflink][i].opacity(0.1);
					CurrentViewData.AFlinks[aflink][i].setZIndex(-1);
				}
				//CurrentViewData.valueLabels[iT].show();
			}
		}
		for(sf in associatedSFs){
			temp = sf.match(/([ATCG])(\d+)/);
			var nucleotide = temp[1];
			var motifPos = temp[2];
			CurrentViewData.SFTerms[motifPos-1][nucleotide].opacity(1);
			CurrentViewData.SFTerms[motifPos-1][nucleotide].setZIndex(99);
			//Draw the value label
			//CurrentViewData.valueLabels[eachTerm[i]].show();
		}
		//CurrentViewData.interTerms[iT].draw();
		var theSF =  CurrentViewData.SFTerms[pos-1][nuc];
		var tooltipPos = {
			x: theSF.position().x + theSF.getWidth()/2,
			y: theSF.position().y + theSF.getHeight() + 5
		};

		if(CurrentViewData.tooltip !== undefined){
			CurrentViewData.tooltip.destroyChildren();
			CurrentViewData.tooltip.destroy();
		}
		var SFterm = nuc+pos;
		var AFterms = {};
		console.log("XDD",SFterm,matrix[currentKey]["AF"]);

		for(af in matrix[currentKey]["AF"]){
			if(checkSFinAF(SFterm, af, ":")){
			//if(af.indexOf(SFterm) >= 0){
				AFterms[af] = parseInt(1000*matrix[currentKey]["AF"][af])/1000;
			}
		}
		var tooltipHead = SFterm + ": " + parseInt(1000*matrix[currentKey]["SF"][SFterm])/1000;
		var AFKeys = [];
		AFKeys = Object.keys(AFterms).sort(function(a, b){
			return AFterms[b] - AFterms[a]
		});
		var tooltipText = "";
		if(matrix[currentKey]["SF"][SFterm] == 0 && AFKeys.length == 0){
			tooltipText = "- Not considered in the model\n";
		}else{

			if(AFKeys.length == 0){
				tooltipText += "- No invovled associations\n";
			}else{
				tooltipText += "- Invovled associations:\n";	
			}
		}
		for(i = 0; i < AFKeys.length; i++){
			tooltipText += "-- "+AFKeys[i]+" => "+AFterms[AFKeys[i]]+"\n";
		}
		
		var tooltip = createToolTip(tooltipPos, tooltipHead ,tooltipText);
		


		CurrentViewData.layer.add(tooltip);
		CurrentViewData.tooltip = tooltip;
		stage.draw();
	}
	
	function highlightForSeq(seq){
		mode = 1;
		var seqArray = seq.split("");
		inspectSeqArray = seqArray;
		var terms = {};
		var posNuc = [];
		var usedSFNumber = 0;
		var totalSFNumber = 0;
		var usedAFNumber = 0;
		var i = 0, j=0, k=0;
		if(CurrentViewData.tooltip !== undefined){
			CurrentViewData.tooltip.destroyChildren();
			CurrentViewData.tooltip.destroy();
			CurrentViewData.tooltip = undefined;
		}
		for(i = 0; i < seqArray.length; i++){
			var term = seqArray[i] + (i+1);
			terms[term] = matrix[currentKey]["SF"][term];
			if(matrix[currentKey]["SF"][term] !== 0){
				usedSFNumber ++;
			}
		}
		for(i = 1; i <= CurrentViewData.SFTerms.length; i++){
			for(j = 0; j < 4; j ++){
				CurrentViewData.SFTerms[i-1][Nucs[j]].opacity(0.1);
				CurrentViewData.SFTerms[i-1][Nucs[j]].setZIndex(-1);
				if(matrix[currentKey]["SF"][Nucs[j]+(i)] !== 0){
					totalSFNumber ++;
					console.log(Nucs[j]+(i));
				}
			}
		}
		for(i = 1; i <= CurrentViewData.SFTerms.length; i++){
			for(k =1; k <= seqArray.length; k++){
				for(j = 0; j < 4; j ++){
					if(k == i && seqArray[k-1] == Nucs[j]){
						CurrentViewData.SFTerms[i-1][Nucs[j]].opacity(1);
						CurrentViewData.SFTerms[i-1][Nucs[j]].setZIndex(100);
					}
				}
			}
		}
		for(aflink in CurrentViewData.AFlinks){
			for(i = 0; i < CurrentViewData.AFlinks[aflink].length; i++){
				CurrentViewData.AFlinks[aflink][i].opacity(0.1);
				CurrentViewData.AFlinks[aflink][i].setZIndex(-1);
			}
		}
		
		for(aflink in CurrentViewData.AFlinks){
			var isHighLight = 1;
			var termsInAf = aflink.split("_");
			for(i = 0; i < termsInAf.length; i++){
				if(terms[termsInAf[i]] == undefined){
					isHighLight = 0;
					break;
				}
			}
			if(isHighLight == 1){
				for(i = 0; i < CurrentViewData.AFlinks[aflink].length; i++){
					CurrentViewData.AFlinks[aflink][i].opacity(1);
					CurrentViewData.AFlinks[aflink][i].setZIndex(99);
				}
			}
		}
		//Finding used Association features
		var AFs = Object.keys(matrix[currentKey]["AF"]).sort(sortTerms);
		
		for(i = 0; i < AFs.length; i++){
			var temp = AFs[i].split(":");
			var containAF = 0;
			for(j = 0; j < temp.length; j++){
				for(k = 1; k <= inspectSeqArray.length; k++){
					if(temp[j] == inspectSeqArray[k-1]+k){
						containAF ++;
					}
				}
			}
			if(containAF == temp.length){
				terms[AFs[i]] = matrix[currentKey]["AF"][AFs[i]];
				usedAFNumber ++;
			}
		}
		var score = parseFloat(matrix[currentKey]["intercept"]);

		for(term in terms){
			score += terms[term];
		}
		var probability = 1/(1+ Math.exp(-1*score));
		probability = parseInt(10000 * probability)/10000;

		var InvestigationText = new Kinetic.Text({
			x: 10, y: 25,
			text: "Test sequence: "+seq + "\nProbability: "+probability + "\nNo. of used single-base features: "+usedSFNumber+"/"+totalSFNumber+"\nNo. of used association features:"+usedAFNumber+"/"+AFs.length,
			fill: "black",
			fontSize: 16
		});
		if(CurrentViewData.investigationText){
			CurrentViewData.investigationText.destroy();
		}

		CurrentViewData.layer.add(InvestigationText);
		CurrentViewData.investigationText = InvestigationText;
		stage.draw();
	}

	function sortTerms(a, b){
		var A = a.split(":");
		var B = b.split(":");
		var ref = [];
		if(A.length < B.length){
			ref = A;
		}else{
			ref = B;
		}
		for(var i=0; i< ref.length; i++){
			var arrayA = A[i].match(/([ATCG])(\d+)/);
			var arrayB = B[i].match(/([ATCG])(\d+)/);
			if(parseInt(arrayA[2]) > parseInt(arrayB[2])){
				return 1;
			}else if(parseInt(arrayA[2]) < parseInt(arrayB[2])){
				return -1;
			}else{
				if(A[i] > B[i]){
					return 1;
				}else{
					return -1;
				}
			}
		}
	}
	function hideInterTermGruop(){
		CurrentViewData.interTermGroup.hide();
		//Show the values?
		stage.draw();
	}
	function coef2LinkWidth(coef){
		var maxWidth = 5;
		var minWidth = 1;
		var returnWidth = 0;
		if(coef > 0 ){
			if(coef <= 3){
				returnWidth = ((maxWidth - minWidth)/2)*coef + minWidth;
			}else{
				returnWidth = maxWidth;
			}
		}else{
			if(coef >= -3){
				returnWidth = ((maxWidth - minWidth)/2)*Math.abs(coef) + minWidth;
			}else{
				returnWidth = maxWidth;
			}
		}
		return returnWidth;
	}
	function coef2Color(coef){
		var returnColor = "";
		var R, G, B;
		if(coef > 0){
			if(coef <= 3){
				R = parseInt(255 - (coef*34));
				G = parseInt(255 - (coef*85));
				B = parseInt(255 - (coef*85));
				returnColor = "rgb("+R+","+G+","+B+")";
			}else{
				returnColor = "rgb(255, 0, 0)";
			}
		}else if(coef < 0){
			if(coef >= -3){
				R = parseInt(85 *( coef + 3));
				G = parseInt(34*(coef + 3) + 153);
				B = parseInt(85 *( coef + 3));
				returnColor = "rgb("+R+","+G+","+B+")";
			}else{
				returnColor = "rgb(0, 255, 0)";
			}
		}else{
			returnColor = "rgb(200,200,200)";
		}
		return returnColor;
	}
	function coef2ColorLink(coef){
		var returnColor = "";
		var R, G, B;
		if(coef > 0){
			if(coef <= 3){
				R = parseInt(coef*153/2);
				returnColor = "rgb("+R+",0,0)";
			}else{
				returnColor = "rgb(255, 0, 0)";
			}
		}else if(coef < 0){
			if(coef >= -3){
				G = -1*parseInt(coef*153/2);
				returnColor = "rgb(0,"+G+",0)";
			}else{
				returnColor = "rgb(0, 255, 0)";
			}
		}else{
			returnColor = "rgb(200,200,200)";
		}
		return returnColor;
	}
	
	$("#lists").delegate(".matrixName","click", function(event){		
		var nameSpan = $(event.target);
		$("#leftContainer").css("left", -180);
		$("#leftControlBar").children("div").text(">");
		$('#Info').hide();
		$('#canvas').show();
		drawMatrices(nameSpan.prop("title"));
		$("#lists").find(".matrixName").css("font-weight","normal");
		nameSpan.css("font-weight","bold");
		
	});
	$('#leftControlBar').click(function(eventj){
		var pos = $("#leftContainer").position();
		if(pos.left == 0){
			$("#leftContainer").css("left", -180);
			$("#leftControlBar").children("div").text(">");
		}else{
			$("#leftContainer").css("left", 0);
			$("#leftControlBar").children("div").text("<");
		}
	});
	$('#leftControlBar').hover(function(){
		$("#leftControlBar").css('background-color', '#ddd');
	}, function(){
		$("#leftControlBar").css('background-color', '#aaa');
	});
	$('#toggleControls').hover(function(){
		$("#toggleControls").css('background-color', '#eee');
	}, function(){
		$("#toggleControls").css('background-color', '#aaa');
	});
	$('#toggleControls').click(function(){
		$("#controlElement").toggle();
	});
	$("#testSeq").keyup(function(event){
		event.preventDefault();
		if(event.which == 13){
			var currentSeq = $( this ).val();
			var currentLength = currentSeq.length;
			var motifLength = matrix[currentKey]["motifLength"];
			$(this).parent().find(".warning").remove();
			
			var newSeq = currentSeq.replace(/[^ACGTacgt]+/, "");
			var newLength = newSeq.length;
			newSeq = newSeq.toUpperCase();
			if(currentLength !== motifLength){
				$(this).after("<span class='warning'>&nbsp*The sequence length should be "+motifLength+" bps which is the same as the motif length.</span>");
				return false;
			}
			if(newLength !== currentLength){
				$(this).after("<span class='warning'>&nbsp*Only accept characters, ACGT.</span>");
				return false;
			}else{
				//Valid sequence
				highlightForSeq(newSeq);
			}
		}
	});
});
	
