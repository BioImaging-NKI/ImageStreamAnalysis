#@ File (label = "Input image stack", style = "file") inputImageStack
#@ File (label = "Output folder", style = "directory") output
#@ Integer (label = "Brightfield channel", value = 1, min=1) ch_BF
#@ Integer (label = "Tumor cell marker", value = 3, min=1) ch_tumor
#@ Integer (label = "T-cell marker", value = 4, min=1) ch_Tcell
#@ Integer (label = "Marker Of Interest (MOI)", value = 2, min=1) ch_MOI
#@ String (label = "Tumor / T-cell marker localization", choices={"nuclear","membrane"}, style="radioButtonHorizontal", value="nuclear") markerLoc
//#@ String (label = "Threshold method for determining cell classes", choices = {"Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, style="list", value="Huang") autoThresholdMethod
#@ String (label = "MOI expressed in cell type", choices={"T cell","Tumor cell"}, style="radioButtonHorizontal", value="T cell") CellOfInterest
#@ Double (label = "cell diameter", value = 30, min=1) cellDiameter
#@ Integer (label = "Shrink cells with (pixels)", value = 2, min=0) erosionSteps
#@ Integer (label = "Membrane thickness (pixels)", value = 3, min=1) membraneThickness
#@ Integer (label = "Membrane thickness for straightening (pixels)", value = 7, min=1) membraneThicknessForStraightening
#@ Boolean (label = "Limit number of analyzed cells", value = false) bool_limitNrCells
#@ Integer (label = "...to", value = 100, min=1) limitNrCells
#@ Boolean (label = "Load labelmap instead of performing Cellpose", value = false) bool_loadLabelmap
#@ File (label = "Labelmap file", style = "file") labelmapFile
#@ File (label = "Cellpose environment path", style="Directory", value="D:\\Software\\Python\\Cellpose3\\venv") env_path
#@ String (label = "Cellpose environment type", choices={"conda","venv"}, style="listBox") env_type
#@ String (label = "Cellpose model", default="cyto3") CellposeModel

List.setCommands;
if (List.get("Cellpose ...")!="") oldCellposeWrapper = false;
else if (List.get("Cellpose Advanced")!="") oldCellposeWrapper = true;
else exit("ERROR: Cellpose wrapper not found! Update Fiji and activate the PT-BIOP Update site.");

output = output + File.separator + File.getNameWithoutExtension(inputImageStack);
if(!File.exists(output)) File.makeDirectory(output);

/* ChangeLog
 *  v1.0.2: Store all script parameters
 *  v1.0.3: Store graphs as png, don't store the montage
 *  v1.0.4: Get ch_array with channel names
 *  v1.1: Added analysis of intensity along the membrane, for synapse detection
 *  v1.2: Updated saving results and some other small things (not affecting the analysis)
 *  v1.3: Adaptations for the updated BIOP Cellpose wrapper
 */


// ---- CONSTANTS ---- //
VERSION = "1.5";
CELL_COUNT_TABLE = "All cell combinations";
CELL_INTENSITY_TABLE = "Cell classification";
DOUBLETS_COUNT_TABLE = "Single Tumor Single T cell";
CellposeFlowThreshold = 0.5;
OPACITY = 100;		//Of the membrane overlays
minCellSize = 100;	//In pixels
maxCellSize = 10000;
RESULT_DECIMALS = 3;
MONTAGE_ASPECT_RATIO = 2.8;
MINIMUM_CLASS_FRACTION = 1; //in percentage
MEASURE_STRAIGHTENED_MEMBRANE_METHOD = "MIP";	//MIP or Mean
N_CLUSTERS = 2;  //Nr of celltypes

// -----STORE SETTINGS ----- //
List.clear();
List.set("Version",VERSION);
List.set("inputImageStack", inputImageStack);
List.set("output",output);
List.set("ch_BF",ch_BF);
List.set("ch_tumor",ch_tumor);
List.set("ch_Tcell",ch_Tcell);
//List.set("ch_TC",ch_TCR);
List.set("cellDiameter",cellDiameter);
List.set("erosionSteps",erosionSteps);
List.set("ch_MOI",ch_MOI);
List.set("CellOfInterest",CellOfInterest);
List.set("membraneThickness",membraneThickness);
List.set("minCellSize",minCellSize);
List.set("maxCellSize",maxCellSize);
List.set("markerLoc",markerLoc);
List.set("bool_limitNrCells",bool_limitNrCells);
List.set("limitNrCells",limitNrCells);

MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
List.set("Date",DayNames[dayOfWeek]+" "+dayOfMonth +"-"+MonthNames[month]+"-"+year);
List.set("Time",hour +":"+minute+":"+second);
settings = List.getList;
File.saveString(settings, output+File.separator+"settings.txt");

//interfaceDilationSteps = round(interfaceThickness/2)

//Close all images and tables
run("Close All");
if(isOpen(CELL_INTENSITY_TABLE)) close(CELL_INTENSITY_TABLE);
if(isOpen(CELL_COUNT_TABLE)) close(CELL_COUNT_TABLE);
if(isOpen(DOUBLETS_COUNT_TABLE)) close(DOUBLETS_COUNT_TABLE);
roiManager("reset");
run("Clear Results");

saveSettings();

run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();
run("Set Measurements...", "area mean min median stack redirect=None decimal="+RESULT_DECIMALS);
run("Conversions...", " ");

open(inputImageStack);
original = getTitle();
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);
ch_array = newArray(channels);
for (i = 1; i <= channels; i++) {
	Stack.setPosition(i, 1, 1);
	ch_array[i-1] = getMetadata("Label");
}
Overlay.remove;
metadata = split(getImageInfo(), "\n");
medians = split(metadata[1], ',');
if(medians.length == 1) medians = newArray(channels);

if(bool_limitNrCells == true && frames > limitNrCells) {
	run("Duplicate...", "duplicate title=[First "+limitNrCells+" cells] frames=1-"+limitNrCells);
	original = getTitle();
	getDimensions(width, height, channels, slices, frames);
}

//Set display scaling & colors
setBatchMode("hide");
selectWindow(original);

wait(10);	//Somehow this is necessary (Why ImageJ, why?!)
run("32-bit");
run("Grays");
Stack.setChannel(ch_Tcell);
run("Biop-Azure");
run("Enhance Contrast", "saturated=0.35");
getMinAndMax(min, max);
setMinAndMax(medians[ch_Tcell-1], max);
run("Set Label...", "label=T cell marker");
Stack.setChannel(ch_tumor);
run("Green");
run("Enhance Contrast", "saturated=0.35");
getMinAndMax(min, max);
setMinAndMax(medians[ch_tumor-1], max);
run("Set Label...", "label=Tumor marker");
Stack.setChannel(ch_MOI);
if(ch_MOI != ch_Tcell && ch_MOI != ch_tumor) {
	run("Red");
	run("Enhance Contrast", "saturated=0.35");
	getMinAndMax(min, max);
	setMinAndMax(medians[ch_MOI-1], max);
	run("Set Label...", "label=Marker of Interest");
}
Stack.setChannel(ch_BF);
run("Grays");
setMinAndMax(medians[ch_BF-1]-100, medians[ch_BF-1]+400);
run("Set Label...", "label=Brightfield");

Stack.setDisplayMode("composite");
Stack.setActiveChannels("1111000");
setBatchMode("show");

//Create montage
montageWidth = floor(sqrt(MONTAGE_ASPECT_RATIO)*floor(sqrt(frames))+1);
montageHeight = floor(frames/montageWidth + round(0.4999 + (frames/montageWidth)%1));
run("Make Montage...", "columns="+montageWidth+" rows="+montageHeight+" scale=1");
rename("Montage of "+original);
montage = getTitle();
montageID = getImageID();
Stack.setActiveChannels("1111000");
setBatchMode("exit and display");

selectWindow(original);
//Preprocess and segment cells
if(markerLoc == "nuclear") {
	selectWindow(original);
	Stack.setChannel(ch_BF);
	run("Duplicate...", "title=BF_variance duplicate channels="+ch_BF);
	run("Variance...", "radius=2 stack");
	
	selectWindow(original);
	Stack.setChannel(ch_tumor);
	tumor = "tumor";
	run("Duplicate...", "title="+tumor+ " duplicate channels="+ch_tumor);
	normalizeImage("tumor");
	Ext.CLIJ2_push(tumor);
	
	selectWindow(original);
	Stack.setChannel(ch_Tcell);
	Tcell = "Tcell";
	run("Duplicate...", "title="+Tcell+ " duplicate channels="+ch_Tcell);
	normalizeImage("Tcell");
	Ext.CLIJ2_push(Tcell);
	
	imageCalculator("Add create stack", "tumor","Tcell");
	rename("nuclei");
	//close("Tcell");
	
	run("Merge Channels...", "c1=BF_variance c2=nuclei create");
	Stack.setChannel(2);
	run("Green");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(1);
	run("Magenta");
	run("Enhance Contrast", "saturated=0.35");
	getDimensions(width, height, channels, slices, frames);
	setBatchMode("show");
	
	nrColumns = frames;
	run("Make Montage...", "columns="+nrColumns+" rows=1 scale=1");
	setBatchMode("show");
	
	if(!bool_loadLabelmap) {
		if(oldCellposeWrapper) {
			run("Cellpose Advanced", "diameter="+cellDiameter+" cellproba_threshold=0.0 flow_threshold="+CellposeFlowThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=2 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
			print("[INFO] Using old CellPose wrapper. Update Fiji / PT-BIOP Update site to use the new version.")
		}
		else run("Cellpose ...", "env_path="+env_path+" env_type="+env_type+" model="+CellposeModel+" model_path=path\\to\\own_cellpose_model diameter="+cellDiameter+" ch1=1 ch2=2 additional_flags=[--use_gpu, --flow_threshold, "+CellposeFlowThreshold+", --cellprob_threshold, 0.0]");
	else open(labelmapFile);
	}
}
else if(markerLoc == "membrane") {	//Sum the membrane marker and the Brightfield image
	selectWindow(original);
	Stack.setChannel(ch_BF);
	run("Duplicate...", "title=BF_variance duplicate channels="+ch_BF);
	run("Variance...", "radius=2 stack");
	normalizeImage("BF_variance");
	
	selectWindow(original);
	Stack.setChannel(ch_tumor);
	tumor = "tumor";
	run("Duplicate...", "title="+tumor+ " duplicate channels="+ch_tumor);
	normalizeImage("tumor");
	Ext.CLIJ2_push(tumor);
	
	selectWindow(original);
	Stack.setChannel(ch_Tcell);
	Tcell = "Tcell";
	run("Duplicate...", "title="+Tcell+ " duplicate channels="+ch_Tcell);
	normalizeImage("Tcell");
	Ext.CLIJ2_push(Tcell);
	
	imageCalculator("Add create stack", "tumor","Tcell");
	rename("markers");
	run("Grays");
	imageCalculator("Add stack", "markers","BF_variance");
	//close("Tcell");

	run("Enhance Contrast", "saturated=0.35");
	getDimensions(width, height, channels, slices, frames);
	setBatchMode("show");
	
	nrColumns = frames;
	run("Make Montage...", "columns="+nrColumns+" rows=1 scale=1");
	setBatchMode("show");
	
	if(!bool_loadLabelmap) {
		if(oldCellposeWrapper) {
			run("Cellpose Advanced", "diameter="+cellDiameter+" cellproba_threshold=0.0 flow_threshold="+CellposeFlowThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
			print("[INFO] Using old CellPose wrapper. Update Fiji / PT-BIOP Update site to use the new version.")
		}
		else run("Cellpose ...", "env_path="+env_path+" env_type="+env_type+" model="+CellposeModel+" model_path=path\\to\\own_cellpose_model diameter="+cellDiameter+" ch1=1 ch2=0 additional_flags=[--use_gpu, --flow_threshold, "+CellposeFlowThreshold+", --cellprob_threshold, 0.0]");
	else open(labelmapFile);
	}
	
///	if(!bool_loadLabelmap) run("Cellpose Advanced", "diameter="+cellDiameter+" cellproba_threshold=0.0 flow_threshold=0.5 anisotropy=1.0 diam_threshold=12.0 model=cyto2 nuclei_channel=2 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
	if(!bool_loadLabelmap) run("Cellpose ...", "env_path="+env_path+" env_type="+env_type+" model="+CellposeModel+" model_path=path\\to\\own_cellpose_model diameter="+cellDiameter+" ch1=1 ch2=0 additional_flags=[--use_gpu, --flow_threshold, "+CellposeFlowThreshold+", --cellprob_threshold, 0.0]");
	else open(labelmapFile);
}

labelmap_cells = "labelmap_cells";
rename(labelmap_cells);
run("16-bit");
if(!bool_loadLabelmap) saveAs("zip", output + File.separator + getTitle());
run("glasbey_on_dark");
setMinAndMax(0, 255);
run("Montage to Stack...", "columns="+nrColumns+" rows=1 border=0");
close(labelmap_cells);
close("Montage");
selectWindow("Stack");
rename(labelmap_cells);

Ext.CLIJ2_push(labelmap_cells);
close(labelmap_cells);
Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_cells, labelmap_cells_gapsclosed);
//shrink labels
eroded_in = labelmap_cells_gapsclosed;
if(erosionSteps > 0) {
	for (i = 0; i < erosionSteps; i++) {
		eroded_out = "erode_"+i+1;
		Ext.CLIJ2_erodeSphereSliceBySlice(eroded_in, eroded_out);
		eroded_in = eroded_out;
		if(i>0) {
			eroded_out_previous = "erode_"+i;
			Ext.CLIJ2_release(eroded_out_previous);
		}
	}
	Ext.CLIJ2_multiplyImages(labelmap_cells_gapsclosed, eroded_out, labelmap_eroded);
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_eroded, labelmap_eroded_sizeFiltered, minCellSize, maxCellSize);
	Ext.CLIJ2_release(labelmap_eroded);
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_eroded_sizeFiltered, labelmap_final);
	Ext.CLIJ2_release(labelmap_eroded_sizeFiltered);
}
else labelmap_final = labelmap_cells_gapsclosed;
Ext.CLIJ2_getMaximumOfAllPixels(labelmap_final, nrCells);
Ext.CLIJ2_pull(labelmap_final);
Ext.CLIJ2_clear();

//run("16-bit");
Stack.setSlice(1);
run("glasbey_on_dark");
setMinAndMax(0, 255);
rename(labelmap_cells);

selectWindow(labelmap_cells);
setBatchMode("show");
run("Label image to composite ROIs", "rm=[RoiManager[size="+nrCells+", visible=true]]");	//BIOP plugin

roiManager("Measure");
cellArea = Table.getColumn("Area", "Results");

//selectWindow(original);
//Stack.setChannel(1);
//run("Add Image...", "image="+labelmap_cells+" x=0 y=0 opacity=33 zero");

/*
meshBetweenTouchingLabelsDilated = "meshBetweenTouchingLabelsDilated";
Ext.CLIJ2_push(labelmap_cells);
Ext.CLIJ2_drawMeshBetweenTouchingLabels(labelmap_cells, meshBetweenTouchingLabels);
Ext.CLIJ2_dilateBox(meshBetweenTouchingLabels, meshBetweenTouchingLabelsDilated);
Ext.CLIJ2_pull(meshBetweenTouchingLabelsDilated);
//Ext.CLIJ2_drawDistanceMeshBetweenTouchingLabels(Image_input, Image_destination);
selectWindow(original);
run("Add Image...", "image="+meshBetweenTouchingLabelsDilated+" x=0 y=0 opacity=100 zero");

run("Clear Results");
Ext.CLIJ2_statisticsOfLabelledPixels(Tcell, labelmap_cells);
TcellIntensity = Table.getColumn("MEAN_INTENSITY");
run("Clear Results");
Ext.CLIJ2_statisticsOfLabelledPixels(tumor, labelmap_cells);
TumorIntensity = Table.getColumn("MEAN_INTENSITY");
 */

selectWindow(Tcell);
run("Clear Results"); 
roiManager("Measure");
TcellIntensity = Table.getColumn("Mean");
selectWindow(tumor);
run("Clear Results"); 
roiManager("Measure");
TumorIntensity = Table.getColumn("Mean");

ratio = divideArrays(TcellIntensity, TumorIntensity);
Table.create(CELL_INTENSITY_TABLE);
Table.reset(CELL_INTENSITY_TABLE);
Table.setColumn("index", Array.getSequence(nrCells), CELL_INTENSITY_TABLE);
Table.setColumn("T-cell intensity", TcellIntensity, CELL_INTENSITY_TABLE);
Table.setColumn("Tumor intensity", TumorIntensity, CELL_INTENSITY_TABLE);
Table.setColumn("Marker ratio", ratio, CELL_INTENSITY_TABLE);
Table.setColumn("Area", cellArea, CELL_INTENSITY_TABLE);

intRatioOverArea = divideArrays(ratio, cellArea);
Table.setColumn("Marker ratio over Area", intRatioOverArea);

//Ext.CLIJ2_pushArray(intRatioOverAreaImage, intRatioOverArea, nrCells, 1, 1);
//Ext.CLIJ2_getAutomaticThreshold(intRatioOverAreaImage, autoThresholdMethod, classThreshold);
//Ext.CLIJ2_pull(intRatioOverAreaImage);
//Ext.CLIJ2_release(intRatioOverAreaImage);
//print("AutoThreshold to discriminate between Tcells and tumor cells ("+autoThresholdMethod+"): "+classThreshold);

setBatchMode(true);
newImage("forClusteringAnalysis", "32-bit black", nrCells, 1, 3);
for(i=0; i<nrCells; i++) {
	setSlice(1);
	setPixel(i, 0, log(Table.get("T-cell intensity", i, CELL_INTENSITY_TABLE)));
	setSlice(2);
	setPixel(i, 0, log(Table.get("Tumor intensity", i, CELL_INTENSITY_TABLE)));
	setSlice(3);
	setPixel(i, 0, Table.get("Area", i, CELL_INTENSITY_TABLE));
}

// normalize to mean 0 std 1 and find clusters using K-means
for (i=1; i<= nSlices; i++){
	setSlice(i);
	getStatistics(area, mean, min, max, std);
	//print(mean + " " + std);
	run("Subtract...", "value="+mean);
	run("Divide...", "value="+std);
}

run("k-means Clustering ...", "number_of_clusters="+N_CLUSTERS+" cluster_center_tolerance=0.00010000 enable_randomization_seed randomization_seed=48");
for(i=0; i<nrCells; i++) Table.set("Class", i, getPixel(i, 0), CELL_INTENSITY_TABLE);
Table.update;

// so what class is tumor? the one with the lowest average ratio!
ratio_sum0=0;
ratio_sum1=0;
for(i=0; i<nrCells; i++) {
	if (Table.get("Class", i, CELL_INTENSITY_TABLE) == 0){
		ratio_sum0 += Table.get("Marker ratio over Area", i, CELL_INTENSITY_TABLE);
	}else{
		ratio_sum1 += Table.get("Marker ratio over Area", i, CELL_INTENSITY_TABLE);
	}
}
tcell_class = 0;
tumor_class = 1;
if (ratio_sum0 < ratio_sum1){
	tcell_class = 1;
	tumor_class = 0;
}


Array.getStatistics(ratio, ratioMin, ratioMax, ratioMean, ratioStdDev);
Array.getStatistics(cellArea, areaMin, areaMax, areaMean, areaStdDev);

labelmap_classes = "labelmap_classes";
selectWindow(labelmap_cells);
run("Duplicate...", "title="+labelmap_classes+" duplicate");

setBatchMode(true);

selectWindow(labelmap_classes);

Table.create(CELL_COUNT_TABLE);
oldFrame = -1;
x = newArray(1);
y = newArray(1);
Plot.create("cell Classes", "Area", "Marker ratio");
for (i = 0; i < roiManager("Count"); i++) {
	roiManager("select", i);
	frame = parseInt(substring(Roi.getName, 0, 4));
	if(frame != oldFrame) {		//We are at new frame; write values for oldFrame into table
		if(oldFrame != -1) {	//Start after First frame or run this at the last frame
			rowIndex = Table.size;
			Table.set("Frame", rowIndex, oldFrame, CELL_COUNT_TABLE);
			Table.set("T cells", rowIndex, tCellCount, CELL_COUNT_TABLE);
			Table.set("Tumor cells", rowIndex, tumorCount, CELL_COUNT_TABLE);
		}
		tCellCount = 0;
		tumorCount = 0;
		oldFrame = frame;
	}
	if(Table.get("Class", i, CELL_INTENSITY_TABLE) == tcell_class) {
		roiManager("Set Color", "blue");
		roiManager("Rename", Roi.getName + " Tcell");
		tCellCount++;
		changeValues(1, 1/0, 1);

		Plot.setColor("blue");
		x[0] = Table.get("Area", i, CELL_INTENSITY_TABLE);
		y[0] = Table.get("Marker ratio", i, CELL_INTENSITY_TABLE);
		Plot.add("Circle", x, y);
	}
	else if(Table.get("Class", i, CELL_INTENSITY_TABLE) == tumor_class) {
		roiManager("Set Color", "yellow");
		roiManager("Rename", Roi.getName + " Tumor");
		tumorCount++;
		changeValues(1, 1/0, 2);

		Plot.setColor("#ffa000");
		x[0] = Table.get("Area", i, CELL_INTENSITY_TABLE);
		y[0] = Table.get("Marker ratio", i, CELL_INTENSITY_TABLE);
		Plot.add("Circle", x, y);
	}
	else {
		roiManager("Set Color", "green");
		roiManager("Rename", Roi.getName + " class2");
		tumorCount++;
		changeValues(1, 1/0, 3);

		Plot.setColor("green");
		x[0] = Table.get("Area", i, CELL_INTENSITY_TABLE);
		y[0] = Table.get("Marker ratio", i, CELL_INTENSITY_TABLE);
		Plot.add("Circle", x, y);
	}

}
Plot.addLegend("0__immune cells\n1__ tumor cells\n2", "Auto");
setBatchMode("show");

roiManager("deselect");
Plot.setLogScaleX(true);
Plot.setLogScaleY(true);
Plot.setLimitsToFit();
Plot.setFrameSize(600, 600);
Plot.show();
setBatchMode("show");
saveAs("PNG", output+File.separator + getTitle()+".png");


selectWindow(CELL_COUNT_TABLE);
rowIndex = Table.size;
Table.set("Frame", rowIndex, oldFrame, CELL_COUNT_TABLE);
Table.set("T cells", rowIndex, tCellCount, CELL_COUNT_TABLE);
Table.set("Tumor cells", rowIndex, tumorCount, CELL_COUNT_TABLE);
Table.update;

Ext.CLIJ2_push(labelmap_classes);

//Create cell membrane masks
Ext.CLIJ2_labelToMask(labelmap_classes, mask_Tcell, 1);
Ext.CLIJ2_labelToMask(labelmap_classes, mask_tumor, 2);
mask_Tcell_eroded = binaryErodeGPU(mask_Tcell, membraneThickness);
Tcell_membrane = "Tcell_membrane";
tumor_membrane = "Tumor cell membrane";
Ext.CLIJ2_binarySubtract(mask_Tcell, mask_Tcell_eroded, Tcell_membrane);
mask_tumor_eroded = binaryErodeGPU(mask_tumor, membraneThickness);
Ext.CLIJ2_binarySubtract(mask_tumor, mask_tumor_eroded, tumor_membrane);

//Create interface and Not-in-interface masks
Tcell_interface = "T cell interface";
tumor_interface = "Tumor cell interface";
Tcell_membrane_not_in_interface = "T cell membrane not in interface";
tumor_membrane_not_in_interface = "Tumor cell membrane not in interface";
//dilate cell masks recursively
mask_Tcell_dilated = binaryDilateGPU(mask_Tcell, membraneThickness);
mask_tumor_dilated = binaryDilateGPU(mask_tumor, membraneThickness);
//Create interface masks and not-in-interface masks
Ext.CLIJ2_binaryIntersection(mask_Tcell, mask_tumor_dilated, Tcell_interface);
Ext.CLIJ2_binaryIntersection(mask_tumor, mask_Tcell_dilated, tumor_interface);
Ext.CLIJ2_binaryXOr(Tcell_interface, Tcell_membrane, Tcell_membrane_not_in_interface);
Ext.CLIJ2_binaryXOr(tumor_interface, tumor_membrane, tumor_membrane_not_in_interface);

//Pull images from GPU and release memory
setBatchMode(true);
Ext.CLIJ2_pull(Tcell_membrane);
setBatchMode("show");
Ext.CLIJ2_pull(tumor_membrane);
setBatchMode("show");
Ext.CLIJ2_pull(Tcell_interface);
setBatchMode("show");
run("Magenta");
Ext.CLIJ2_pull(tumor_interface);
setBatchMode("show");
run("Green");
Ext.CLIJ2_pull(Tcell_membrane_not_in_interface);
run("biop-Azure");
Ext.CLIJ2_pull(tumor_membrane_not_in_interface);
run("Yellow");
setBatchMode("show");
//Ext.CLIJ2_clear();

//Overlay membranes and interfaces on original image
selectWindow(original);
originalID = getImageID();
setBatchMode("hide");
selectWindow(Tcell_membrane);
Tcell_membraneID = getImageID();
selectWindow(tumor_membrane);
tumor_membraneID = getImageID();
selectWindow(Tcell_interface);
Tcell_interfaceID = getImageID();
selectWindow(tumor_interface);
tumor_interfaceID = getImageID();
selectWindow(Tcell_membrane_not_in_interface);
Tcell_membrane_not_in_interfaceID = getImageID();
selectWindow(tumor_membrane_not_in_interface);
tumor_membrane_not_in_interfaceID = getImageID();

Overlay.remove;
for (f = 1; f <= frames; f++) {
	selectImage(Tcell_interfaceID);
	setSlice(f);
	selectImage(tumor_interface);
	setSlice(f);
	selectImage(Tcell_membrane_not_in_interface);
	setSlice(f);
	selectImage(tumor_membrane_not_in_interface);
	setSlice(f);
	selectImage(originalID);
	Stack.setFrame(f);	
	run("Add Image...", "image=["+Tcell_interface+"] x=0 y=0 opacity="+OPACITY+" zero");
	Overlay.setPosition(0);
	run("Add Image...", "image=["+tumor_interface+"] x=0 y=0 opacity="+OPACITY+" zero");
	Overlay.setPosition(0);
	run("Add Image...", "image=["+Tcell_membrane_not_in_interface+"] x=0 y=0 opacity="+OPACITY+" zero");
	Overlay.setPosition(0);	
	run("Add Image...", "image=["+tumor_membrane_not_in_interface+"] x=0 y=0 opacity="+OPACITY+" zero");
	Overlay.setPosition(0);	
}
setBatchMode("show");

selectWindow(original);
run("To ROI Manager");	//Copy overlay to ROI manager
Overlay.remove;
run("From ROI Manager");
roiManager("Show All without labels");

	
for (f = 1; f <= frames; f++) {
	selectImage(originalID);
	Stack.setFrame(f);
	Overlay.copy;
	x = (f-1)%montageWidth;
	y = floor((f-1)/montageWidth);
	selectImage(montageID);
	Overlay.paste;
//	Overlay.moveTo(width*x, height*y);
	Overlay.setPosition(0);		//N.B. Doesn't work for these overlays, because it is somehow treated as a single item, copied from the ROI manager
	Overlay.moveSelection(4*(f-1), x*width, y*height);
	Overlay.moveSelection(4*(f-1)+1, x*width, y*height);
	Overlay.moveSelection(4*(f-1)+2, x*width, y*height);
	Overlay.moveSelection(4*(f-1)+3, x*width, y*height);
}
selectImage(montageID);
for (f = 1; f <= frames; f++) {
	x = (f-1)%montageWidth;
	y = floor((f-1)/montageWidth);
	Overlay.drawString(f, x*width+width/2, y*height+height);
}

//saveAs("TIF", output+File.separator + getTitle() +".tif");

//Measure intensities
selectWindow(original);
run("Duplicate...", "duplicate title=MOI channels="+ch_MOI);
MOI = "MOI";
Ext.CLIJ2_push(MOI);

for (i = 0; i < frames; i++) {
	showStatus("Measuring intensities... "+i+1+"/"+frames);	//Doesn't work here, probably because of GPU computing
	showProgress(i, frames);
//	print("frame " + i+1 + " of " + frames);
//	selectWindow(MOI);
//	Stack.setFrame(i+1);
//	Ext.CLIJ2_pushCurrentSlice(MOI);
	Ext.CLIJ2_copySlice(MOI, MOI_slice, i);
	
	Ext.CLIJ2_copySlice(Tcell_interface, mask2D, i);
	Ext.CLIJ2_getMeanOfMaskedPixels(MOI_slice, mask2D, mean);
	Table.set("mean Tcell interface", i, mean, CELL_COUNT_TABLE);

	Ext.CLIJ2_copySlice(Tcell_membrane_not_in_interface, mask2D, i);
	Ext.CLIJ2_getMeanOfMaskedPixels(MOI_slice, mask2D, mean);
	Table.set("mean Tcell NOT in interface", i, mean, CELL_COUNT_TABLE);

	Ext.CLIJ2_copySlice(tumor_interface, mask2D, i);
	Ext.CLIJ2_getMeanOfMaskedPixels(MOI_slice, mask2D, mean);
	Table.set("mean tumor interface", i, mean, CELL_COUNT_TABLE);

	Ext.CLIJ2_copySlice(tumor_membrane_not_in_interface, mask2D, i);
	Ext.CLIJ2_getMeanOfMaskedPixels(MOI_slice, mask2D, mean);
	Table.set("mean tumor NOT in interface", i, mean, CELL_COUNT_TABLE);
}
meanTCellInterface = Table.getColumn("mean Tcell interface", CELL_COUNT_TABLE);
meanTCellNotInInterface = Table.getColumn("mean Tcell NOT in interface", CELL_COUNT_TABLE);
TCell_membrane_ratio = divideArrays(meanTCellInterface, meanTCellNotInInterface);
meanTumorInterface = Table.getColumn("mean tumor interface", CELL_COUNT_TABLE);
meanTumorNotInInterface = Table.getColumn("mean tumor NOT in interface", CELL_COUNT_TABLE);
Tumor_membrane_ratio = divideArrays(meanTumorInterface, meanTumorNotInInterface);
Table.setColumn("T cell membrane ratio", TCell_membrane_ratio, CELL_COUNT_TABLE);
Table.setColumn("Tumor membrane ratio", Tumor_membrane_ratio, CELL_COUNT_TABLE);

selectWindow(CELL_COUNT_TABLE);
Table.update;
Table.save(output+File.separator+CELL_COUNT_TABLE + ".tsv")
selectWindow(CELL_INTENSITY_TABLE);
Table.update;
Table.save(output+File.separator+CELL_INTENSITY_TABLE + ".tsv")

duplicate_table(CELL_COUNT_TABLE, DOUBLETS_COUNT_TABLE);
Table.showRowIndexes(true);
//Table.setLocationAndSize(0, 0, 0, 0);
for(i=frames-1; i>=0; i--) {
	if(Table.get("T cells",i) == 1 && Table.get("Tumor cells",i) == 1 && toString(Table.get("mean Tcell interface",i)) != "NaN") continue;
	else Table.deleteRows(i, i);
}
//Table.setLocationAndSize(0, 0, 800, 800);
selectWindow(DOUBLETS_COUNT_TABLE);
Table.update;
Table.save(output+File.separator+DOUBLETS_COUNT_TABLE + ".tsv")

graph_allCells = create_interface_intensity_graph(CELL_COUNT_TABLE);
saveAs("PNG", output+File.separator + getTitle() +".png");
graph_doublets = create_interface_intensity_graph(DOUBLETS_COUNT_TABLE);
saveAs("PNG", output+File.separator + getTitle() +".png");


if(CellOfInterest == "T cell") {
	mask = "Tcell_membrane";
	mask_in_interface = Tcell_interface;
}
else if(markerLoc == "membrane") {
	mask = "Tumor cell membrane";
	mask_in_interface = tumor_interface;
}
else {
	mask = "Tumor cell membrane";
	mask_in_interface = tumor_interface;
}
MOI_image = MOI;


setLineWidth(1);
setBatchMode(true);
roiManager("reset");
run("Clear Results");
print("\\Clear");
if(isOpen("membrane_profiles")) close("membrane_profiles");
if(isOpen("interface_profiles")) close("interface_profiles");

//selectWindow(mask);
//selectWindow(MOI_image);

roiManager("reset");
selectWindow(mask);
setAutoThreshold("Default dark");
run("Analyze Particles...", "add stack");
resetThreshold();
nrcells = roiManager("count");
traceLengths = newArray(nrcells);

selectImage(MOI_image);
metadata = split(getImageInfo(), "\n");
medians = split(metadata[1], ',');
if(medians.length == 1) medians = newArray(channels);
run("Duplicate...", "title=MOI_membrane duplicate channels="+ch_MOI);
run("Subtract...", "value="+medians[ch_MOI-1]+" stack");
membrane_profiles = create_straightened_membrane_MIP_profiles("MOI_membrane", traceLengths);
rename("membrane_profiles");
close("MOI_membrane");
membrane_profiles = "membrane_profiles";
changeValues(0.0, 0.0, -10000);
setThreshold(-9999, 65535);
run("NaN Background");
setMinAndMax(0, 500);

interface_profiles_mask = create_straightened_membrane_MIP_profiles(mask_in_interface, traceLengths);
rename("interface_profiles_mask");
interface_profiles_mask = "interface_profiles_mask";
setThreshold(1E-5, 65535);
run("Convert to Mask");
run("Duplicate...", "title=non-interface_profiles_mask");
run("Invert");
run("Divide...", "value=255.0");
run("32-bit");
setThreshold(0.5, 65535);
run("NaN Background");
setMinAndMax(0, 1);

selectImage(interface_profiles_mask);
run("Divide...", "value=255.0");
run("32-bit");
setThreshold(0.5, 65535);
run("NaN Background");
setMinAndMax(0, 1);

imageCalculator("Multiply create", "membrane_profiles","interface_profiles_mask");
rename("membrane_profiles_interface");
setMinAndMax(0, 500);
setBatchMode("show");
imageCalculator("Multiply create", "membrane_profiles","non-interface_profiles_mask");
rename("membrane_profiles_non-interface");
setBatchMode("show");
setMinAndMax(0, 500);

getDimensions(width, height, channels, slices, frames);
for(i=0; i<nrcells; i++) {
	selectImage("membrane_profiles_interface");
	makeRectangle(0, i, traceLengths[i], 1);
	interfaceMean = getValue("Mean");

	selectImage("membrane_profiles_non-interface");
	makeRectangle(0, i, traceLengths[i], 1);
	nonInterfaceMean = getValue("Mean");
	nonInterfaceStdDev = getValue("StdDev");
	
	profile = getProfile();

	Table.set("interface Mean", i, interfaceMean, "Results");
	Table.set("non-interface Mean", i, nonInterfaceMean, "Results");
	Table.set("non-interface stdDev", i, nonInterfaceStdDev, "Results");
	Table.set("Ratio interface/non-interface", i, interfaceMean/nonInterfaceMean, "Results");
	Table.set("stdDevs above mean", i, (interfaceMean - nonInterfaceMean)/nonInterfaceStdDev, "Results");

	//Internal control: create control interface if there is none at a random location and perform measurements
	controlInterfaceLength = 25;
	if(isNaN(interfaceMean)) {
		position = round(random*traceLengths[i]-1);
		//print(i, position, traceLengths[i]);
		Array.rotate(profile, position);
		controlInterface=Array.slice(profile, 0, controlInterfaceLength);
		controlNonInterface = Array.slice(profile, controlInterfaceLength, profile.length);
		controlInterfaceMean = meanOfArray(controlInterface);
		controlNonInterfaceMean = meanOfArray(controlNonInterface);
		controlNonInterfaceStdDev = stdDevOfArray(controlNonInterface);
		Table.set("control stdDevs above mean", i, (controlInterfaceMean - controlNonInterfaceMean)/controlNonInterfaceStdDev, "Results");
	}
	else Table.set("control stdDevs above mean", i, NaN, "Results");

	//Normalize membrane traces in image
	selectImage("membrane_profiles");
	makeRectangle(0, i, traceLengths[i], 1);
	run("Subtract...", "value="+nonInterfaceMean);
	run("Divide...", "value="+nonInterfaceStdDev);
}
Table.update;
Table.save(output+File.separator+"Membrane intensities.tsv")


setBatchMode(false);
//setLineWidth(3);
run("Distribution...", "parameter=[control stdDevs above mean] or=30 and=-5-10");
setForegroundColor(255, 0, 0);
setLineWidth(3);
//drawLine(109, 138, 109, 10);
//run("Select None");
saveAs("PNG", output + File.separator + getTitle());
run("Distribution...", "parameter=[stdDevs above mean] or=30 and=-5-10");
setForegroundColor(255, 0, 0);
setLineWidth(3);
//drawLine(109, 138, 109, 10);
//run("Select None");
saveAs("PNG", output + File.separator + getTitle());

roiManager("Set Color", "#3300ffff");
setForegroundColor(255, 255, 255);
setLineWidth(1);


//Detect and align interfaces

// Required when only running this part
selectWindow("membrane_profiles");
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);
nrcells = height;
traceLengths = newArray(nrcells);
for(i=0; i<height; i++) {
	x=0;
	while(!isNaN(getPixel(x, i))) x++;
	traceLengths[i] = x;
}

selectWindow("interface_profiles_mask");
shift = newArray(nrcells);
for(i=0; i<nrcells; i++) {
	if (traceLengths[i] == 0) continue;
	makeRectangle(0, i, traceLengths[i], 1);
	profile = getProfile();
	interfaceLength = sumArray(profile);
	if(profile[0]!=1) index = firstIndexOfArray(profile, 1);
	else {
		Array.reverse(profile);
		index=0;
		while(profile[index] == 1) index++;
		index = traceLengths[i]-index;
		Array.reverse(profile);
	}
	shift[i] = -index + 0 + interfaceLength/2;
	Array.rotate(profile, shift[i]);
	for(x=0; x<profile.length; x++) {
		setPixel(x, i, profile[x]);
	}
}
updateDisplay();

selectWindow("membrane_profiles");
for(i=0; i<nrcells; i++) {
	if (traceLengths[i] == 0) continue;
	makeRectangle(0, i, traceLengths[i], 1);
	profile = getProfile();
	Array.rotate(profile, shift[i]);
	for(x=0; x<profile.length; x++) {
		setPixel(x, i, profile[x]);
	}
}
updateDisplay();

profileLength = 200;
makeRectangle(0, 0, profileLength, nrcells);
profile = getProfile();
xValues = multiplyArraywithScalar(Array.getSequence(profileLength), pixelWidth);
Plot.create("Average membrane profile", "membrane perimeter position ("+unit+")", "normalized intensity", xValues, profile);

//Returns the first index at which a value occurs in an array
function firstIndexOfArray(array, value) {
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			break;
		}
	}
	return a;
}
//Returns the last index at which a value occurs in an array
function lastIndexOfArray(array, value) {
	for (a=lengthOf(array)-1; a>0; a--) {
		if (array[a]==value) {
			break;
		}
	}
	return a;
}
//Returns the sum of all elements of an arrays, neglecting NaNs
function sumArray(array) {
	sum=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
	}
	return sum;
}
//Multiplies all elements of an array with a scalar and returns the new array
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]=array[a]*scalar;
	}
	return multiplied_array;
}



restoreSettings;


//// FUNCTIONS ////



function create_straightened_membrane_MIP_profiles(target, traceLengths) {
	selectImage(target);
	setBatchMode("hide");
	nrcells = roiManager("count");
	names = newArray(nrcells);
	concatString = "";
	for (i=0; i<nrcells; i++) {
		showStatus("Straightening membranes... "+i+1+"/"+nrcells);
		showProgress(i, nrcells);
		selectImage(target);
		roiManager("select", i);
		if(Roi.getType != "freeline") {	//Do this only in the first pass
			run("Area to Line");
			//Roi.getSplineAnchors(x, y)
			//Roi.setPolylineSplineAnchors(x, y);
			roiManager("update");
		}
//		run("Interpolate", "interval=0.5");	//Do not smooth, because it will create a smaller straigthened profile!
		run("Straighten...", "title=membrane_"+i+" line="+membraneThicknessForStraightening);
		names[i] = getTitle;
	//	IDs = getImageID();
		getDimensions(width, height, channels, slices, frames);
		traceLengths[i] = width;
		concatString += "image"+i+1+"="+names[i]+" ";
	}
	selectImage(target);
	setSlice(1);
	roiManager("Deselect");
	run("Select None");
	setBatchMode("show");
	
	Array.getStatistics(traceLengths, minTraceLength, maxTraceLength, meanTraceLength, stdDevTraceLength);
	for (i=0; i<nrcells; i++) {
		showStatus("Resizing images... "+i+1+"/"+nrcells);
		showProgress(i, nrcells);
		selectImage(names[i]);
		run("Canvas Size...", "width="+maxTraceLength+" height="+membraneThicknessForStraightening+" position=Center-Left zero");
	}
	run("Concatenate...", concatString);
	straigtened_membranes = "straigtened_membranes";
	straigtened_membranes_max = "straigtened_membranes_max";
	rename(straigtened_membranes);
	setBatchMode("show");
	run("CLIJ2 Macro Extensions", "cl_device=");
	Ext.CLIJ2_clear();
	Ext.CLIJ2_push(straigtened_membranes);
	close(straigtened_membranes);
	if(MEASURE_STRAIGHTENED_MEMBRANE_METHOD == "MIP") Ext.CLIJ2_maximumYProjection(straigtened_membranes, straigtened_membranes_max);
	else if(MEASURE_STRAIGHTENED_MEMBRANE_METHOD == "Mean") Ext.CLIJ2_meanYProjection(straigtened_membranes, straigtened_membranes_max);
	Ext.CLIJ2_pull(straigtened_membranes_max);
	Ext.CLIJ2_release(straigtened_membranes_max);
	setBatchMode("show");
	
	return straigtened_membranes_max;
}


function create_interface_intensity_graph(table) {
	if(CellOfInterest == "T cell") {
		Array.getStatistics(Table.getColumn("mean Tcell interface", table), minY, maxY, meanY, stdevY);
		Array.getStatistics(Table.getColumn("mean Tcell NOT in interface", table), minX, maxX, meanX, stdevX);
		Plot.create("Channel "+ch_MOI+" in T cell membrane plot - "+table, "mean Tcell NOT in interface", "mean Tcell interface");
		Plot.add("Circle", Table.getColumn("mean Tcell NOT in interface", table), Table.getColumn("mean Tcell interface", table));
		Plot.setStyle(0, "blue,#a0a0ff,1.0,Circle");
		Plot.setColor("gray");
		Plot.drawLine(0, 0, maxOf(maxX, maxY), maxOf(maxX, maxY));
		Plot.setFrameSize(600, 600);
		Plot.setLimits(0,maxOf(maxX, maxY),0,maxOf(maxX, maxY));
		Plot.show();
		setBatchMode("show");
	}
	else if (CellOfInterest == "Tumor cell") {
		Array.getStatistics(Table.getColumn("mean tumor interface", table), minY, maxY, meanY, stdevY);
		Array.getStatistics(Table.getColumn("mean tumor NOT in interface", table), minX, maxX, meanX, stdevX);
		Plot.create("Channel "+ch_MOI+" in Tumor cell membrane plot - "+table, "mean tumor NOT in interface", "mean tumor interface");
		Plot.add("Circle", Table.getColumn("mean tumor NOT in interface", table), Table.getColumn("mean tumor interface", table));
		Plot.setStyle(0, "blue,#a0a0ff,1.0,Circle");
		Plot.setColor("gray");
		Plot.drawLine(0, 0, maxOf(maxX, maxY), maxOf(maxX, maxY));
		Plot.setFrameSize(600, 600);
		Plot.setLimits(0,maxOf(maxX, maxY),0,maxOf(maxX, maxY));
		Plot.show();
		setBatchMode("show");
	}
	return getTitle();
}

//Duplicate a table
function duplicate_table(inputTable, outputTable){
	Table.create(outputTable);
	headings = split(Table.headings(inputTable), "\t");
    for (col=0; col<headings.length; col++) {
    	if(col==0 && headings[col]==" ") continue;
    	else {
			col_values = Table.getColumn(headings[col], inputTable);
			Table.setColumn(headings[col], col_values, outputTable);
    	}
    }
	Table.update(outputTable);
}

//Normalize the image with the median value of Otsu thresholds in every frame
function normalizeImage(image) {
	getDimensions(width, height, channels, slices, frames);
	thresholds = newArray(frames);
	for (i = 0; i < frames; i++) {
		Stack.setFrame(i+1);
		setAutoThreshold("Otsu dark");
		getThreshold(thresholds[i], upper);
	}
	resetThreshold;
	median = medianOfArray(thresholds);

	selectWindow(image);
	if(bitDepth() != 32) run("32-bit");
	run("Divide...", "value="+median+" stack");
	run("Enhance Contrast", "saturated=0.35");
}

//Returns a 32-bit binary eroded image with a number of iterations, in GPU memory. The input binaryImage should already be in GPU memory. 
function binaryErodeGPU(binaryImage, iterations) {
	image_in = binaryImage;
	for (i = 0; i < iterations; i++) {
		image_out = binaryImage+i+1;
		Ext.CLIJ2_erodeSphereSliceBySlice(image_in, image_out);
		image_in = image_out;
		if(i>0) {
			image_out_previous = binaryImage+i;
			Ext.CLIJ2_release(image_out_previous);
		}
	}
	return image_out;
}

//Returns a 32-bit binary dilated image with a number of iterations, in GPU memory. The input binaryImage should already be in GPU memory. 
function binaryDilateGPU(binaryImage, iterations) {
	image_in = binaryImage;
	for (i = 0; i < iterations; i++) {
		image_out = binaryImage+i+1;
		Ext.CLIJ2_dilateSphereSliceBySlice(image_in, image_out);
		image_in = image_out;
		if(i>0) {
			image_out_previous = binaryImage+i;
			Ext.CLIJ2_release(image_out_previous);
		}
	}
	return image_out;
}

//Returns the median of the array
function medianOfArray(array) {
	array_sorted = Array.copy(array);
	Array.sort(array_sorted);
	return array_sorted[floor(array_sorted.length/2)];
}

//Divides the elements of two arrays and returns the new array
function divideArrays(array1, array2) {
	divArray=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		divArray[a]=array1[a]/array2[a];
	}
	return divArray;
}

//Returns the average of all elements of an arrays, neglecting NaNs
function averageArray(array) {
	sum=0;
	nans=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
		else nans+=1;
	}
	return sum/(array.length-nans);
}


//Returns the mean of the array
function meanOfArray(array) {
	Array.getStatistics(array, min, max, mean, stdDev);
	return mean;
}


//Returns the stdDev of the array
function stdDevOfArray(array) {
	Array.getStatistics(array, min, max, mean, stdDev);
	return stdDev;
}


//Returns the first index at which a value occurs in an array
function firstIndexOfArray(array, value) {
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			break;
		}
	}
	return a;
}


//Returns the last index at which a value occurs in an array
function lastIndexOfArray(array, value) {
	for (a=lengthOf(array)-1; a>0; a--) {
		if (array[a]==value) {
			break;
		}
	}
	return a;
}


//Returns the sum of all elements of an arrays, neglecting NaNs
function sumArray(array) {
	sum=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
	}
	return sum;
}


//Multiplies all elements of an array with a scalar and returns the new array
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]=array[a]*scalar;
	}
	return multiplied_array;
}
