#@ String (value="<html><p style='font-size:12px; color:navy; font-weight:bold'>Input settings</p></html>", visibility="MESSAGE") input_message
#@ File[] (label = "Input image stack", style = "file") inputImages
#@ File (label = "Output folder", style = "directory") outputFolder
#@ String (label = "Basename of output files", value = "Panel_1", min=1) baseName
#@ Integer (label = "Brightfield channel", value = 1, min=1) ch_BF
#@ Integer (label = "Tumor cell marker channel", value = 3, min=1) ch_Tumor
#@ String (label = "Independent (phenotypic) immune cell markers (comma-separated list)", value = "2,4", min=1) phenotypicImmuneChannelsString
//#@ Integer (label = "number of different cell types (including tumor cells)", value = 3, min=1) nrCellTypes
#@ String (label = "Functional immune cell marker channels (comma-separated list)", value = "2,3") functionalImmuneChannelsString
#@ String (label = "Functional immune cell marker channels used in tumor-immune cell classification", value = "2,3") functionalImmuneChannelsForClusteringString
#@ String (value="<html><p style='font-size:12px; color:navy; font-weight:bold'>Options</p></html>", visibility="MESSAGE") options_message
//#@ Integer (label = "Marker Of Interest (MOI)", value = 2, min=1) ch_MOI
#@ String (label = "Tumor / Immune cell marker localization", choices={"nuclear","membrane"}, style="radioButtonHorizontal", value="nuclear") markerLoc
//#@ String (label = "MOI expressed in cell type", choices={"Immune cell","Tumor cell"}, style="radioButtonHorizontal", value="Immune cell") CellOfInterest
#@ Boolean (label = "Smooth brightfield before cell segmentation", value = false) smoothBrightfield
#@ Double (label = "cell diameter", value = 20, min=1) cellDiameter
#@ Integer (label = "Shrink cells with (pixels)", value = 2, min=0) erosionSteps
#@ Integer (label = "Membrane thickness (pixels)", value = 3, min=1) membraneThickness
#@ String (label = "Thresholding immune cell channels", choices={"automatic","manual"}, style="radioButtonHorizontal", value="automatic") thresholdingMethod
#@ Boolean (label = "Load measurements from disk", value = false) loadMeasurements
#@ File (label = "Measurements file", style = "file", value="-") measurementsFile
#@ Boolean (label = "Load labelmap from disk", value = false) loadLabelmap
#@ File (label = "Labelmap file", style = "file", value="-") labelmapFile
#@ Boolean (label = "Limit number of analyzed cells", value = false) bool_limitNrImages
#@ Integer (label = "...to", value = 100, min=1) limitNrImages
#@ String (value="<html><p style='font-size:12px; color:navy; font-weight:bold'>Cellpose environment settings</p></html>", visibility="MESSAGE") cellpose_message
#@ File (label = "Cellpose environment path", style="Directory", value="D:\\Software\\Python\\Cellpose3\\venv") env_path
#@ String (label = "Cellpose environment type", choices={"conda","venv"}, style="listBox") env_type
#@ String (label = "Cellpose model", default="cyto3") CellposeModel
#@ Boolean (label = "Debug mode (show extra images and info)", value = false) debugMode


/* Macro to analyze cell-cell interactions in Imaging Flow Cytometry data
 * More info on https://github.com/BioImaging-NKI/ImageStreamAnalysis/
 * 
 * ► Requires the following Fiji update sites:
 * - CLIJ
 * - CLIJ2
 * - IJPB-plugins
 * - PTBIOP, with proper settings for the Fiji Cellpose wrapper
 * 
 * ► A working Cellpose Python environment
 * 
 * 
 * Authors: Bram van den Broek & Rolf Harkes, The Netherlands Cancer Institute, b.vd.broek@nki.nl
 */


List.setCommands;
if (List.get("Cellpose ...")!="") oldCellposeWrapper = false;
else if (List.get("Cellpose Advanced")!="") oldCellposeWrapper = true;
else exit("ERROR: Cellpose wrapper not found! Update Fiji and activate the PT-BIOP Update site.");

nrCellTypes = 2;	//K-means Clustering only for Tumor and Immune cells

//Other parameters
useLogIntensityForClusterning = false;	//FOR NOW LEAVE THIS OFF - Log is taken anyway
removeBadDetections = false;
CellposeFlowThreshold = 0.5;
immuneChannelsLUTs = newArray("biop-Azure", "Green", "Red", "biop-Amber", "Magenta");
overlayClassesMembranes_bool = false;

// ---- CONSTANTS ---- //
VERSION = 2.8;
CELL_COUNT_TABLE = "Image statistics";
CELL_INTENSITY_TABLE = "Cell statistics";
DOUBLETS_COUNT_TABLE = "Single Tumor Single Immune cell";
OPACITY = 50;		//Of the membrane overlays
minCellSize = 50;	//In pixels
maxCellSize = 10000;
RESULT_DECIMALS = 3;
MONTAGE_ASPECT_RATIO = 1.6;
MINIMUM_CLASS_FRACTION = 1; //in percentage
//FONTCOLOR = "256,165,0";
FONTCOLOR = "Magenta";
LABELFONTSIZE = 8;

// -----STORE SETTINGS ----- //
List.clear();
MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
List.set("Date",DayNames[dayOfWeek]+" "+dayOfMonth +"-"+MonthNames[month]+"-"+year);
List.set("Time",hour +":"+minute+":"+second);
List.set("Version", VERSION);
List.set("inputImages", arrayToString(inputImages, ","));
List.set("output", outputFolder);
List.set("ch_BF", ch_BF);
List.set("ch_Tumor", ch_Tumor);
List.set("phenotypicImmuneChannelsString", phenotypicImmuneChannelsString);
List.set("nrCellTypes", nrCellTypes);
List.set("functionalImmuneChannelsString", functionalImmuneChannelsString);
List.set("functionalImmuneChannelsForClusteringString", functionalImmuneChannelsForClusteringString);
//List.set("ch_TC",ch_TCR);
List.set("cellDiameter",cellDiameter);
List.set("erosionSteps",erosionSteps);
//List.set("ch_MOI",ch_MOI);
//List.set("CellOfInterest",CellOfInterest);
//List.set("membraneThickness",membraneThickness);
List.set("CellposeFlowThreshold", CellposeFlowThreshold);
List.set("minCellSize", minCellSize);
List.set("maxCellSize", maxCellSize);
List.set("markerLoc", markerLoc);
List.set("bool_limitNrImages", bool_limitNrImages);
List.set("limitNrImages", limitNrImages);

if(!File.exists(outputFolder)) File.makeDirectory(outputFolder);

//interfaceDilationSteps = round(interfaceThickness/2)
phenotypicImmuneChannels = string_to_int_array(phenotypicImmuneChannelsString);
functionalImmuneChannels = string_to_int_array(functionalImmuneChannelsString);
functionalImmuneChannelsForClustering = string_to_int_array(functionalImmuneChannelsForClusteringString);


run("Close All");
if(isOpen(CELL_INTENSITY_TABLE)) close(CELL_INTENSITY_TABLE);
if(isOpen(CELL_COUNT_TABLE)) close(CELL_COUNT_TABLE);
if(isOpen(DOUBLETS_COUNT_TABLE)) close(DOUBLETS_COUNT_TABLE);
if(isOpen("labelmap_cells-Morphometry")) close("labelmap_cells-Morphometry");
roiManager("reset");
run("Clear Results");
print("\\Clear");

run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();
run("Set Measurements...", "area mean min median stack redirect=None decimal="+RESULT_DECIMALS);

//Get max dimension of all input images
setBatchMode(true);
maxWidth = 0;
maxHeight = 0;
for (i = 0; i < inputImages.length; i++) {
	open(inputImages[i]);
	getDimensions(width, height, channels, slices, frames);
	maxWidth = maxOf(maxWidth, width);
	maxHeight = maxOf(maxHeight, height);
	if(i>0) { if(channels != nrChannels) exit("The input images do not have the same number of channels!"); }
	else nrChannels = channels;
}
//Expand images, replacing background with the medians
concatString = "";
for (i = 0; i < inputImages.length; i++) {
	showStatus("Resizing images "+i+1+"/"+inputImages.length);
	selectImage(File.getName(inputImages[i]));
	run("Canvas Size...", "width="+maxWidth+" height="+maxHeight+" position=Center zero");
	metadata = split(getImageInfo(), "\n");
	medians = split(metadata[1], ',');
	if(medians.length == 1) medians = newArray(channels);
	if(debugMode) { print("Backgrounds (medians) of "+File.getName(inputImages[i])+":"); Array.print(medians); }
	getDimensions(width, height, channels, slices, frames);
	for(c=1; c<=channels; c++) {
		Stack.setChannel(c);
		for(f=1; f<=frames; f++) {
			Stack.setFrame(f);
			changeValues(0, 0, medians[c-1]);
		}
	}
	concatString += " image"+i+1+"=["+File.getName(inputImages[i])+"]";
}
if(inputImages.length > 1) run("Concatenate...", "  title=input_stack"+concatString);
setBatchMode("show");

//baseName = File.getNameWithoutExtension(inputImages);
original = getTitle();
getDimensions(width, height, channels, slices, totalFrames);
ch_array = newArray(channels);
for (c = 1; c <= channels; c++) {
	Stack.setPosition(c, 1, 1);
	ch_array[c-1] = getMetadata("Label");
}
Stack.setPosition(1, 1, 1);
//Overlay.remove;
//metadata = split(getImageInfo(), "\n");
//medians = split(metadata[1], ',');
//if(medians.length == 1) medians = newArray(channels);
//if(debugMode) { print("Medians:"); Array.print(medians); }

if(bool_limitNrImages == true && frames > limitNrImages) {
	run("Duplicate...", "duplicate title=[First "+limitNrImages+" cells] frames=1-"+limitNrImages);
	original = getTitle();
	getDimensions(width, height, channels, slices, frames);
	nrImages = frames;
}
else nrImages = totalFrames;

selectWindow(original);
setBatchMode("hide");
Stack.setDisplayMode("grayscale");
run("32-bit");
run("Grays");

//Set active channels to display
Stack.setChannel(ch_BF);
Stack.setDisplayMode("composite");
activeChannels = newArray(channels);
activeChannels[ch_BF-1] = 1;
activeChannels[ch_Tumor-1] = 1;
for(i=0; i<phenotypicImmuneChannels.length; i++) activeChannels[phenotypicImmuneChannels[i]-1] = 1;
Stack.setActiveChannels(arrayToString(activeChannels, ""));
 

//Create montage
montageWidth = floor(sqrt(MONTAGE_ASPECT_RATIO)*floor(sqrt(nrImages))+1);
montageHeight = floor(nrImages/montageWidth + round(0.4999 + (nrImages/montageWidth)%1));

run("Make Montage...", "columns="+montageWidth+" rows="+montageHeight+" scale=1");
rename("Montage of "+original);
montage = getTitle();
montageID = getImageID();
Stack.setDisplayMode("composite");
Stack.setActiveChannels(arrayToString(activeChannels, ""));


//Set display scaling & colors
selectWindow(montage);
setBatchMode("hide");

for(i=0; i<phenotypicImmuneChannels.length; i++) {
	selectWindow(montage);
	Stack.setChannel(phenotypicImmuneChannels[i]);
	run("Biop-Azure");
	run("Enhance Contrast", "saturated=0.35");
	getMinAndMax(min, max);
	setMinAndMax(medians[phenotypicImmuneChannels[i]-1], max);
	
	selectWindow(original);
	Stack.setChannel(phenotypicImmuneChannels[i]);
	run("Biop-Azure");
	setMinAndMax(medians[phenotypicImmuneChannels[i]-1], max);
}

for(i=0; i<functionalImmuneChannels.length; i++) {
	selectWindow(montage);
	Stack.setChannel(functionalImmuneChannels[i]);
	run(immuneChannelsLUTs[i]);
	run("Enhance Contrast", "saturated=0.35");
	getMinAndMax(min, max);
	setMinAndMax(medians[functionalImmuneChannels[i]-1], max);
	
	selectWindow(original);
	Stack.setChannel(functionalImmuneChannels[i]);
	run(immuneChannelsLUTs[i]);
	setMinAndMax(medians[functionalImmuneChannels[i]-1], max);
}

//Set Tumor B&C
selectWindow(montage);
Stack.setChannel(ch_Tumor);
run("Yellow");
run("Enhance Contrast", "saturated=0.1");
getMinAndMax(min, max);
setMinAndMax(medians[ch_Tumor-1], max);
selectWindow(original);
Stack.setChannel(ch_Tumor);
run("Yellow");
setMinAndMax(medians[ch_Tumor-1], max);

//Set brightfield B&C
selectWindow(montage);
Stack.setChannel(ch_BF);
run("Grays");
setMinAndMax(medians[ch_BF-1]-100, medians[ch_BF-1]+400);
selectWindow(original);
Stack.setChannel(ch_BF);
setMinAndMax(medians[ch_BF-1]-100, medians[ch_BF-1]+400);

label = newArray(channels);
selectWindow(original);
for(c=1; c<=channels; c++) {
	Stack.setChannel(c);
	label[c-1] = getInfo("slice.label");
}
Stack.setChannel(1);
setBatchMode("show");

selectWindow(montage);
for(c=1; c<=channels; c++) {
	Stack.setChannel(c);
	run("Set Label...", "label="+label[c-1]);
}
Stack.setChannel(1);
setBatchMode("show");


//Preprocess and segment cells
selectWindow(original);
if(markerLoc == "nuclear") {
	selectWindow(original);
	Stack.setChannel(ch_BF);
	run("Duplicate...", "title=BF_variance duplicate channels="+ch_BF);
	if(smoothBrightfield) run("Mean...", "radius=1");
	run("Variance...", "radius=2 stack");
	
	selectWindow(original);
	Stack.setChannel(ch_Tumor);
	tumor = "tumor";
	run("Duplicate...", "title="+tumor+ " duplicate channels="+ch_Tumor);
	normalizeImage("tumor");
	Ext.CLIJ2_push(tumor);
	
	for(i=0; i<phenotypicImmuneChannels.length; i++) {
		selectWindow(original);
		Stack.setChannel(phenotypicImmuneChannels[i]);
		run("Duplicate...", "title=ImCell_"+i+" duplicate channels="+phenotypicImmuneChannels[i]);
		normalizeImage("ImCell_"+i);
		if(i>0) {
			imageCalculator("Add 32-bit stack", "ImCell_"+i,"ImCell_"+i-1);
			close("ImCell_"+i-1);
			//close("ImCell_"+i);
			//rename("ImCell_"+i);
		}
	}
	rename("ImCell");
	ImCell = "ImCell";
	Ext.CLIJ2_push(ImCell);
	
	imageCalculator("Add create stack", "tumor","ImCell");
	rename("nuclei");
	close("ImCell");

	//Merge "membrane" and nucleus channels
	run("Merge Channels...", "c1=BF_variance c2=nuclei create");
	Stack.setChannel(2);
	run("Green");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(1);
	run("Magenta");
	run("Enhance Contrast", "saturated=0.35");
	getDimensions(width, height, channels, slices, frames);
	rename("for_segmentation");
	setBatchMode("show");
	if(loadLabelmap == false) {
		run("Make Montage...", "columns=1 rows="+nrImages+" scale=1");
		setBatchMode("show");
		setBatchMode(false);
		if(oldCellposeWrapper) {
			run("Cellpose Advanced", "diameter="+cellDiameter+" cellproba_threshold=0.0 flow_threshold="+CellposeFlowThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=2 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
			print("[INFO] Using old CellPose wrapper. Update Fiji / PT-BIOP Update site to use the new version.")
		}
		else run("Cellpose ...", "env_path="+env_path+" env_type="+env_type+" model="+CellposeModel+" model_path=path\\to\\own_cellpose_model diameter="+cellDiameter+" ch1=1 ch2=2 additional_flags=[--use_gpu, --flow_threshold, "+CellposeFlowThreshold+", --cellprob_threshold, 0.0]");
	}
	else {
		open(labelmapFile);
		//run("Make Montage...", "columns=1 rows="+nrImages+" scale=1");
		setBatchMode("show");
	}
}
else if(markerLoc == "membrane") {	//Sum the membrane markers and the variance of the Brightfield image
	selectWindow(original);
	Stack.setChannel(ch_BF);
	run("Duplicate...", "title=BF_variance duplicate channels="+ch_BF);
	if(smoothBrightfield) run("Mean...", "radius=1");
	run("Variance...", "radius=2 stack");
	normalizeImage("BF_variance");
	
	selectWindow(original);
	Stack.setChannel(ch_Tumor);
	tumor = "tumor";
	run("Duplicate...", "title="+tumor+ " duplicate channels="+ch_Tumor);
	normalizeImage("tumor");
	Ext.CLIJ2_push(tumor);
	
//	selectWindow(original);
//	Stack.setChannel(ch_ImCell);
//	ImCell = "ImCell";
//	run("Duplicate...", "title="+ImCell+ " duplicate channels="+ch_ImCell);
//	normalizeImage("ImCell");
//	Ext.CLIJ2_push(ImCell);
	
	for(i=0; i<phenotypicImmuneChannels.length; i++) {
		selectWindow(original);
		Stack.setChannel(phenotypicImmuneChannels[i]);
		run("Duplicate...", "title=ImCell_"+i+" duplicate channels="+phenotypicImmuneChannels[i]);
		normalizeImage("ImCell_"+i);
		if(i>0) {
			imageCalculator("Add 32-bit stack", "ImCell_"+i,"ImCell_"+i-1);
			close("ImCell_"+i-1);
			//close("ImCell_"+i);
			//rename("ImCell_"+i);
		}
	}
	rename("ImCell");
	setBatchMode("show");
	ImCell = "ImCell";
	Ext.CLIJ2_push(ImCell);

	//Add all segmentation markers together
	imageCalculator("Add create stack", "tumor","ImCell");
	rename("markers");
	run("Grays");
	imageCalculator("Add stack", "markers","BF_variance");
	rename("for_segmentation");
	close("ImCell");
	close("BF_variance");
	
	run("Enhance Contrast", "saturated=0.35");
	getDimensions(width, height, channels, slices, frames);
	setBatchMode("show");

	if(loadLabelmap == false) {
		run("Make Montage...", "columns=1 rows="+nrImages+" scale=1");
		setBatchMode("show");
		setBatchMode(false);
		if(oldCellposeWrapper) {
			run("Cellpose Advanced", "diameter="+cellDiameter+" cellproba_threshold=0.0 flow_threshold="+CellposeFlowThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
			print("[INFO] Using old CellPose wrapper. Update Fiji / PT-BIOP Update site to use the new version.")
		}
		else run("Cellpose ...", "env_path="+env_path+" env_type="+env_type+" model="+CellposeModel+" model_path=path\\to\\own_cellpose_model diameter="+cellDiameter+" ch1=1 ch2=0 additional_flags=[--use_gpu, --flow_threshold, "+CellposeFlowThreshold+", --cellprob_threshold, 0.0]");
	}
	else {
		open(labelmapFile);
		//run("Make Montage...", "columns=1 rows="+nrImages+" scale=1");
		setBatchMode("show");
	}
}
setBatchMode(true);
if(loadLabelmap == false) saveAs("zip", outputFolder + File.separator + baseName +"_labelmap_cells");
labelmap_cells_temp = "labelmap_cells_temp";
labelmap_cells = "labelmap_cells";
rename(labelmap_cells_temp);
run("Conversions...", " ");	//ensure no scaling
run("16-bit");				//Convert to 16-bit before pushing, because of this issue: https://forum.image.sc/t/strange-behavior-of-clij2s-closeindexgapsinlabelmap-possible-bug/88761

//Create labelmap stacks from segmentation montage and create a new montage
selectWindow(labelmap_cells_temp);
run("glasbey_on_dark");
setMinAndMax(0, 255);
run("Montage to Stack...", "columns=1 rows="+totalFrames+" border=0");
run("Slice Remover", "first="+nrImages+1+" last="+totalFrames+" increment=1");
close(labelmap_cells_temp);
close("Montage");
selectWindow("Stack");
rename(labelmap_cells_temp);
setBatchMode("show");

Ext.CLIJ2_push(labelmap_cells_temp);
close(labelmap_cells_temp);
Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_cells_temp, labelmap_cells, minCellSize, maxCellSize);
Ext.CLIJ2_pull(labelmap_cells);
run("glasbey_on_dark");
setMinAndMax(0, 255);
setBatchMode("show");

selectImage(labelmap_cells);
run("Make Montage...", "columns="+montageWidth+" rows="+montageHeight+" scale=1");
Stack.setSlice(1);
run("glasbey_on_dark");
setMinAndMax(0, 255);
labelmap_cells_montage_temp = labelmap_cells + "_montage_temp";
labelmap_cells_montage = labelmap_cells + "_montage";
rename(labelmap_cells_montage_temp);
run("32-bit");
setBatchMode("show");

Ext.CLIJ2_push(labelmap_cells_montage_temp);
//shrink labels
eroded_in = labelmap_cells_montage_temp;
if(erosionSteps > 0) {
	for (i = 0; i < erosionSteps; i++) {
		eroded_out = "erode_"+i+1;
		Ext.CLIJ2_erodeSphere(eroded_in, eroded_out);
		eroded_in = eroded_out;
		if(i>0) {
			eroded_out_previous = "erode_"+i;
			Ext.CLIJ2_release(eroded_out_previous);
		}
	}
}
if(erosionSteps > 0) {
	Ext.CLIJ2_multiplyImages(eroded_out, labelmap_cells_montage_temp, labelmap_cells_montage);
	Ext.CLIJ2_release(eroded_out);
	Ext.CLIJ2_pull(labelmap_cells_montage);
}
else labelmap_cells_montage = labelmap_cells_montage_temp;
run("glasbey_on_dark");
setMinAndMax(0, 255);
setBatchMode("show");

//Get interactions
selectWindow(labelmap_cells_montage);
run("Analyze Regions", "area centroid");
x_montage = Table.getColumn("Centroid.X", "labelmap_cells_montage-Morphometry");
y_montage = Table.getColumn("Centroid.Y", "labelmap_cells_montage-Morphometry");
run("Region Adjacency Graph", "show image=["+montage+"]");
print("\\Update:");

//create (inside) edges excluding cell interfaces (2 pixels wide)
Ext.CLIJ2_erodeSphere(labelmap_cells_montage, mask_cells_eroded);	//dilate 1 pixel
Ext.CLIJ2_dilateSphere(labelmap_cells_montage, mask_cells_dilated);	//erode 1 pixel
Ext.CLIJ2_binarySubtract(mask_cells_dilated, mask_cells_eroded, mask_cells_edges);
Ext.CLIJ2_release(mask_cells_eroded);
Ext.CLIJ2_release(mask_cells_dilated);

//create membrane labelmap (better than Ext.CLIJ2_reduceLabelsToLabelEdges, because the membrane gets messes up when dilated afterwards)
Ext.CLIJ2_detectLabelEdges(labelmap_cells_montage, mask_membranes);	//Creates a 2-pixel think mask
Ext.CLIJ2_pull(mask_membranes);
run("Grays");
setBatchMode("show");

//Subtract to create interface masks (2 pixels wide), and dilate if needed
//mask_interfaces = "mask_interfaces";
Ext.CLIJ2_subtractImages(mask_membranes, mask_cells_edges, mask_interfaces_temp);
mask_interfaces = binaryDilateGPU(mask_interfaces_temp, maxOf(membraneThickness-1,1));
Ext.CLIJ2_pull(mask_interfaces);
rename("mask_interfaces");
setBatchMode("show");
Ext.CLIJ2_release(mask_cells_edges);
Ext.CLIJ2_release(mask_interfaces_temp);

//Dilate membranes and labelmap
labelmap_membranes = "labelmap_membranes";
Ext.CLIJ2_dilateLabels(labelmap_cells_montage, labelmap_cells_dilated, maxOf(membraneThickness-1,1));
mask_membranes_dilated = binaryDilateGPU(mask_membranes, maxOf(membraneThickness-1,1));
Ext.CLIJ2_multiplyImages(mask_membranes_dilated, labelmap_cells_dilated, labelmap_membranes);
Ext.CLIJ2_pull(labelmap_membranes);
run("glasbey_on_dark");
setBatchMode("show");
//Ext.CLIJ2_release(mask_membranes_dilated);
Ext.CLIJ2_pull(mask_membranes_dilated);
run("Cyan");

//Create interfaces labelmap
labelmap_interfaces = "labelmap_interfaces";
Ext.CLIJ2_multiplyImages(mask_interfaces, labelmap_cells_dilated, labelmap_interfaces);
Ext.CLIJ2_pull(labelmap_interfaces);
run("glasbey_on_dark");
setBatchMode("show");
Ext.CLIJ2_release(labelmap_cells_dilated);

//Create noninterfaces labelmap
labelmap_noninterfaces = "labelmap_noninterfaces";
Ext.CLIJ2_subtractImages(labelmap_membranes, labelmap_interfaces, labelmap_noninterfaces);
Ext.CLIJ2_pull(labelmap_noninterfaces);
run("glasbey_on_dark");
setBatchMode("show");

// THIS HAS THE PROBLEM THAT SOME MEMBRANES WILL GROW OUTSIDE THE IMAGE, CAUSING A DOUBLE MEASUREMENT AND LATER A CRASH
//selectWindow(labelmap_membranes);
//run("Montage to Stack...", "columns="+montageWidth+" rows="+montageHeight+" border=0");
//run("Slice Remover", "first="+frames+1+" last="+montageWidth*montageHeight+" increment=1");
//rename("labelmap_membranes_separate_images");


setBatchMode(false);
selectWindow(labelmap_cells);
run("Label image to composite ROIs", "rm=[RoiManager[size=, visible=true]]");	//BIOP plugin
nrCells = roiManager("count");
selectImage(original);
setBatchMode(true);

//Measure intensities
//run("Set Measurements...", "mean area standard median stack redirect=None decimal=3");
run("Set Measurements...", "area mean standard stack display redirect=None decimal=3");
Table.create(CELL_INTENSITY_TABLE);
Table.reset(CELL_INTENSITY_TABLE);
cellIDs = addScalarToArray(Array.getSequence(nrCells), 1);
Table.setColumn("cell ID", cellIDs, CELL_INTENSITY_TABLE);
for(i=0; i<nrCells; i++) {
	roiManager("select", i);
	Table.set("image ID", i, parseInt(substring(Roi.getName, 0, 4)), CELL_INTENSITY_TABLE);	//Get the frame number from the ROI name
}
imageIDs = Table.getColumn("image ID");

selectImage(original);
//selectImage(montage);	//For later use/more efficient measurements. Get the imageID from ROIs (see a few lines up)
getDimensions(width, height, channels, slices, frames);
roiManager("Deselect");
for(c=1; c<=channels; c++) {
	run("Clear Results");
	showStatus("Measuring intensities in channel "+c+"/"+channels);
	Stack.setChannel(c);
	
	if(loadMeasurements == false) {
		roiManager("Measure");
//		Ext.CLIJ2_statisticsOfLabelledPixels(montage, labelmap_membranes);
		if(c==1) {
			cellArea = Table.getColumn("Area", "Results");
			Table.setColumn("Area [µm^2]", cellArea, CELL_INTENSITY_TABLE);
		}
		meanList = Table.getColumn("Mean", "Results");
		meanList = subtract_scalar_from_array(meanList, medians[c-1]);	//subtract background - N.B. This uses only the medians of the last input image! 
		Table.setColumn("Mean ch"+c+" ("+label[c-1]+")", meanList, CELL_INTENSITY_TABLE);
		//to do: measure the membranes - MorpholibJ on the montage with the montage labelmaps
	}
}
//Get total immune cell marker intensity - NOT NORMALIZED (YET)!
for(i=0; i<nrCells; i++) {
	sum = 0;
	for(c=0; c<phenotypicImmuneChannels.length; c++) {
		sum += Table.get("Mean ch"+phenotypicImmuneChannels[c]+" ("+label[phenotypicImmuneChannels[c]-1]+")", i, CELL_INTENSITY_TABLE);
	}
	Table.set("Phenotypic immune Cell intensity", i, sum, CELL_INTENSITY_TABLE);
}

if(loadMeasurements) {
	close(CELL_INTENSITY_TABLE);
	open(measurementsFile);
	Table.rename(File.getName(measurementsFile), CELL_INTENSITY_TABLE);
	cellArea = Table.getColumn("Area [µm^2]");
}


ImCellIntensity = Table.getColumn("Phenotypic immune Cell intensity", CELL_INTENSITY_TABLE);
TumorIntensity = Table.getColumn("Mean ch"+ch_Tumor+" ("+label[ch_Tumor-1]+")", CELL_INTENSITY_TABLE);
ratio = divideArrays(ImCellIntensity, TumorIntensity);
Table.setColumn("Immune/tumor marker ratio", ratio, CELL_INTENSITY_TABLE);
Table.update();

//Document interactions
selectWindow(labelmap_cells_montage+"-RAG");
label1 = Table.getColumn("Label 1");
label2 = Table.getColumn("Label 2");
//setOption("ExpandableArrays", true);
//interactions = newArray(nrCells);
nrInteractionsPerImage = newArray(nrImages);

selectWindow(CELL_INTENSITY_TABLE);
for(i=0; i<nrCells; i++) {
	interactionsA = lookup(label1, i+1, label2);	// A->B
	interactionsB = lookup(label2, i+1, label1);	// B->A
	if(interactionsA.length > 0 || interactionsB.length > 0) {
		interactions = Array.concat(interactionsA, interactionsB);
		string = arrayToString(interactions, ",");
		Table.set("Interaction with (IDs)", i, string, CELL_INTENSITY_TABLE);
		nrInteractionsPerImage[imageIDs[i]-1] += interactions.length/2;	//Add the number of interactions (only once)
//		print(i, imageIDs[i], nrInteractionsPerImage[imageIDs[i]-1]);
	}
	else Table.set("Interaction with (IDs)", i, "", CELL_INTENSITY_TABLE);
}
Table.update();
Table.setLocationAndSize(300, 100, 1000, 800);
close(labelmap_cells_montage+"-RAG");

/*
//Get the number of interacting cells per image (quite elaborate at the moment - need to find a better way!)
nrInteractionsPerImage = newArray(nrImages);
xxx = Table.getColumn("Interaction with (IDs)";
for(i=0; i<nrImages; i++) {
	interactionsInCurrentImageArray = lookup(imageIDs, i+1, xxx, CELL_INTENSITY_TABLE);
	for(k=0; k<interactionsInCurrentImageArray.length; k++) {
		split(interactionsInCurrentImageArray[k], ",");
	}
	nrInteractionsPerImage = uniqueValuesInArray(array);
}
*/
//Sum other immune marker channels data - normalize first
functionalImmuneChannelsForClusteringSum = newArray(nrCells);
for(i=0; i<functionalImmuneChannelsForClustering.length; i++) {
	data = Table.getColumn("Mean ch"+functionalImmuneChannelsForClustering[i]+" ("+label[functionalImmuneChannelsForClustering[i]-1]+")", CELL_INTENSITY_TABLE);
	normalizedData = normalizeArrayToMedian0AndStddev1(data, 75);
	Table.setColumn("Mean ch"+functionalImmuneChannelsForClustering[i]+" ("+label[functionalImmuneChannelsForClustering[i]-1]+")", data, CELL_INTENSITY_TABLE);
	functionalImmuneChannelsForClusteringSum = addArrays(functionalImmuneChannelsForClusteringSum, data);
}
//Also create an array of *all* immune cell markers (for visualization of the clusters)
allImmuneCellMarkersIntensity = addArrays(ImCellIntensity, functionalImmuneChannelsForClusteringSum);
Table.setColumn("immune cell markers for clustering", allImmuneCellMarkersIntensity, CELL_INTENSITY_TABLE);

setBatchMode(true);


//Prepare image for clustering - exclude the crap detections
setBatchMode(true);

//Remove crap detections (positive in all other immune channels)
badDetections = newArray(nrCells);
if(removeBadDetections == true) {
	for(i=0; i<functionalImmuneChannels.length; i++) {
		for (k = nrCells-1; k >= 0; k--) {
			if(PositiveInfunctionalImmuneChannels[k] == functionalImmuneChannels.length) {
				showStatus("Detecting crap (positive in all other immune channels)...");
				//print("cell "+k+1+" in image "+imageIDs[k]+" is crap!);
				badDetections[k] = 1;
			}
		}
		Table.update;
		selectWindow(labelmap_membranes_dilated);
		updateDisplay();
		print(sumArray(badDetections)+" objects removed (positive in all other immune channels)");
	}
}
nrCells = roiManager("count");


if(functionalImmuneChannelsForClustering.length > 0) newImage("forClusteringAnalysis", "32-bit black", nrCells, 1, 4);
else newImage("forClusteringAnalysis", "32-bit black", nrCells, 1, 2);
k=0;
for(i=0; i<nrCells; i++) {
	if(badDetections[i] != 1) {
		setSlice(1);
		if(useLogIntensityForClusterning) setPixel(k, 0, log(maxOf(ImCellIntensity[i],0.1)));
		else setPixel(k, 0, ImCellIntensity[i]);
		setSlice(2);
		if(useLogIntensityForClusterning) setPixel(k, 0, log(maxOf(TumorIntensity[i],0.1)));
		else setPixel(k, 0, TumorIntensity[i]);
		if(functionalImmuneChannelsForClustering.length > 0) {
			setSlice(3);
			if(useLogIntensityForClusterning) setPixel(k, 0, log(maxOf(functionalImmuneChannelsForClusteringSum[i],0.1)));
			else setPixel(k, 0, functionalImmuneChannelsForClusteringSum[i]);
		}
		setSlice(4);
		setPixel(k, 0, cellArea[i]);
		k++;
	}
}
nrGoodCells = k;
makeRectangle(0, 0, k-1, 1);
run("Crop");
run("Select None");

if(debugMode) setBatchMode("show");
//run("Histogram", "bins=256 use x_min=0 x_max=1 y_max=Auto");
//setBatchMode("show");
//waitForUser("first");
//selectWindow("forClusteringAnalysis");

//Add the minimum value + 1
for(i=0; i<nSlices; i++) {
	setSlice(i+1);
	min = getValue("Min");
	run("Add...", "value="+(1-min)+" slice");
}
for(i=0; i<nSlices-1; i++) {
	setSlice(i+1);
	run("Log", "slice");
}

setBatchMode("show");


// normalize to mean 0 std 1 and find clusters using K-means
for (i=1; i<= nSlices; i++) {
	setSlice(i);
	getStatistics(area, mean, min, max, std);
	// normalize to mean 0 std 1
	run("Subtract...", "value="+mean);
	run("Divide...", "value="+std);
}
//run("Histogram", "bins=256 use x_min=0 x_max=1 y_max=Auto");
//setBatchMode("show");
//waitForUser("after normalization");
//selectWindow("forClusteringAnalysis");

//Remove extreme values
//	changeValues(-9999, -2, -2);
//	changeValues(5, 9999, 5);

if(debugMode) setBatchMode("show");
run("k-means Clustering ...", "number_of_clusters="+nrCellTypes+" cluster_center_tolerance=0.00010000 enable_randomization_seed randomization_seed=48");

//MANUAL ALTERNATIVE: Imcell/Tumor marker ratio >1 && Imcell marker > 10
//newImage("Clusters", "8-bit black", nrCells, 1, 1);
//for (i=1; i<= nSlices; i++) {
//	setSlice(i);
//	for (k=0; k<nrCells; k++) {
//		if(getPixel(k, 0) >= ratio[k] && ImCellIntensity[k] >= 10) setPixel(k, 0, 1);
//		else if (getPixel(k, 0) < ratio[k] && ImCellIntensity[k] < 10) setPixel(k, 0, 0);
//	}
//}

selectWindow("Clusters");
if(debugMode) setBatchMode("show");
classArray = newArray(nrCells);
k=0;
for(i=0; i<nrCells; i++) {
	if(badDetections[i] != 1) {
		classArray[i] = getPixel(k, 0);
		k++;
	}
	else classArray[i] = 2;
}

// so which class is tumor? the one with the lowest average immune/tumor marker ratio!
ratio_sum0=0;
ratio_sum1=0;
for(i=0; i<nrCells; i++) {
	if (classArray[i] == 0){
		ratio_sum0 += ratio[i];
	}else if (classArray[i] == 1){
		ratio_sum1 += ratio[i];
	}
}
imCell_class = 0;
tumor_class = 1;
invert_classes = false;
if (ratio_sum0 < ratio_sum1){
	imCell_class = 1;
	tumor_class = 0;
	invert_classes = true;
}

selectImage("Clusters");
if(invert_classes == true) run("XOR...", "value=1");
run("Add...", "value=1");
getDimensions(clustersWidth, clustersHeight, clustersChannels, clustersSlices, clustersFrames);
run("Canvas Size...", "width="+clustersWidth+1+" height=1 position=Center-Right zero");	//Add a zero for the background
clusters = "Clusters";
labelmap_classes_montage = "labelmap_classes_montage";
Ext.CLIJ2_push(clusters);
Ext.CLIJ2_generateParametricImage(labelmap_membranes, clusters, labelmap_classes_montage);
Ext.CLIJ2_pull(labelmap_classes_montage);
getLut(reds, greens, blues);
reds[1]=0; greens[1]=127; blues[1]=255;
reds[2]=255; greens[2]=255; blues[2]=0;
setLut(reds, greens, blues);
setMinAndMax(0, 255);
setBatchMode("show");

if(debugMode) print("Tumor cell class:\t"+tumor_class);
if(debugMode) print("Immune cell class:\t"+imCell_class);
//N.B. We assume that the 'crap' class is always 3!

for(i=0; i<nrCells; i++) {
	if(classArray[i] == tumor_class) Table.set("Class", i, "Tumor", CELL_INTENSITY_TABLE);
	else if(classArray[i] == imCell_class) Table.set("Class", i, "Immune", CELL_INTENSITY_TABLE);
	else Table.set("Class", i, "Crap", CELL_INTENSITY_TABLE);
}
/* Class labels
for(i=0; i<nrCells; i++) {
	if(classArray[i] == tumor_class) Table.set("Class label", i, tumor_class+1, CELL_INTENSITY_TABLE);
	else if(classArray[i] == imCell_class) Table.set("Class label", i, imCell_class+1, CELL_INTENSITY_TABLE);
	else Table.set("Class label", i, 3, CELL_INTENSITY_TABLE);
}
*/
Table.update;

Array.getStatistics(ratio, ratioMin, ratioMax, ratioMean, ratioStdDev);
Array.getStatistics(cellArea, areaMin, areaMax, areaMean, areaStdDev);

labelmap_classes = "labelmap_classes";
selectWindow(labelmap_cells);
run("Select None");
run("Duplicate...", "title="+labelmap_classes+" duplicate");


//Rename ROIs, write cell_count_table and create clustering plot
selectWindow(labelmap_classes);
Table.create(CELL_COUNT_TABLE);
numbersArray = addScalarToArray(Array.getSequence(nrImages), 1);
Table.setColumn("ImageID", numbersArray);
oldFrame = 0;
x = newArray(1);
y = newArray(1);
Plot.create("Cell Classes Plot", "Tumor cell marker intensity", "All Immune cell markers used for clustering intensity");
for (i = 0; i < nrCells; i++) {
	roiManager("select", i);
	frame = parseInt(substring(Roi.getName, 0, 4));	//Get the frame number from the ROI name
	if(frame > oldFrame) {		//We are at new frame; write values for oldFrame into table
		if(oldFrame != 0) {		//Start after First frame or run this at the last frame
			rowIndex = frame-1;
//			Table.set("Image ID", rowIndex, oldFrame, CELL_COUNT_TABLE);
			Table.set("Immune cells", rowIndex, ImCellCount, CELL_COUNT_TABLE);
			Table.set("Tumor cells", rowIndex, tumorCount, CELL_COUNT_TABLE);
			if(removeBadDetections == true) Table.set("crap", rowIndex, crapCount, CELL_COUNT_TABLE);
		}
		ImCellCount = 0;
		tumorCount = 0;
		crapCount = 0;
		oldFrame = frame;
	}
	
	if(classArray[i] == imCell_class) {
		roiManager("Set Color", "blue");
//		roiManager("Rename", Roi.getName + " Immune");
		ImCellCount++;
		changeValues(1, 1/0, 1);

		Plot.setColor("#0033ff");
		x[0] = TumorIntensity[i];
		y[0] = allImmuneCellMarkersIntensity[i];
		Plot.add("Circle", x, y);
	}
	else if(classArray[i] == tumor_class) {
		roiManager("Set Color", "yellow");
//		roiManager("Rename", Roi.getName + " Tumor");
		tumorCount++;
		changeValues(1, 1/0, 2);

		Plot.setColor("#ffa000");
		x[0] = TumorIntensity[i];
		y[0] = allImmuneCellMarkersIntensity[i];
		Plot.add("Circle", x, y);
	}
	else {	//third class - currently bad detections (positive in all additional immune cell markers) 
		roiManager("Set Color", "gray");
//		roiManager("Rename", Roi.getName + " crap");
		crapCount++;
		changeValues(1, 1/0, 3);

		Plot.setColor("#ff0000");
		x[0] = TumorIntensity[i];
		y[0] = allImmuneCellMarkersIntensity[i];
		Plot.add("circle", x, y);
	}
	if(frame == 1) {	//Ugly way of fixing frame 1, but it works
		rowIndex = 0;
		Table.set("Immune cells", rowIndex, ImCellCount, CELL_COUNT_TABLE);
		Table.set("Tumor cells", rowIndex, tumorCount, CELL_COUNT_TABLE);
		if(removeBadDetections == true) Table.set("crap", rowIndex, crapCount, CELL_COUNT_TABLE);
	}
}
Plot.setLegend("0__immune cells\n1__tumor cells", "Auto");
Plot.setStyle(1, "#ffa000,none,1.0,Circle");
//write values for the last frame
rowIndex = frame-1;
//Table.set("Image ID", rowIndex, oldFrame, CELL_COUNT_TABLE);
Table.set("Immune cells", rowIndex, ImCellCount, CELL_COUNT_TABLE);
Table.set("Tumor cells", rowIndex, tumorCount, CELL_COUNT_TABLE);
if(removeBadDetections == true) Table.set("crap", rowIndex, crapCount, CELL_COUNT_TABLE);
setBatchMode("show");

Table.setLocationAndSize(1300, 100, 1000, 800);
Table.setColumn("Interactions", nrInteractionsPerImage, CELL_COUNT_TABLE);

roiManager("deselect");
Plot.setLogScaleX(true);
Plot.setLogScaleY(true);
Plot.setLimitsToFit();
Plot.getLimits(xMin, xMax, yMin, yMax);
Plot.setLimits(0.1, xMax, 0.1, yMax);
Plot.setFrameSize(600, 600);
Plot.show();
setBatchMode("show");
saveAs("tif", outputFolder + File.separator + baseName +"_celltype_clustering.tif");


selectWindow(montage);
getLocationAndSize(xpos, ypos, width, height);

//Classify the immune cells subtypes
functionalImmuneChannelsThresholdsLog = newArray(functionalImmuneChannels.length);
functionalImmuneChannelsThresholds = newArray(functionalImmuneChannels.length);
functionalImmuneChannelsColumnNames = newArray(functionalImmuneChannels.length);
PositiveInfunctionalImmuneChannels = newArray(nrCells);
for(i=0; i<functionalImmuneChannels.length; i++) {
	data = Table.getColumn("Mean ch"+functionalImmuneChannels[i]+" ("+label[functionalImmuneChannels[i]-1]+")", CELL_INTENSITY_TABLE);
	dataImageID = arrayTo1DImage(data);
	rename("ch"+functionalImmuneChannels[i]+" : "+label[functionalImmuneChannels[i]-1]);
	min = getValue("Min");
	run("Subtract...", "value="+min+" slice");
	run("Add...", "value=1 slice");
	run("Log");	//natural logarithm - intensity distributions are usually log-normal
	setAutoThreshold("MaxEntropy dark");
	getThreshold(functionalImmuneChannelsThresholdsLog[i], upper);
	resetThreshold();

	run("Histogram", "bins=100 x_min=0 x_max=10 y_max=Auto");
	rename("(log)Histogram of "+label[functionalImmuneChannels[i]-1]);
	histogram = getImageID();

	selectImage(histogram);
	setBatchMode("show");
	setLocation(xpos+width, ypos);

	if(thresholdingMethod == "manual") {
		selectWindow(montage);
		run("Duplicate...", "title=montage_ch"+functionalImmuneChannels[i]+" duplicate channels="+functionalImmuneChannels[i]);
		run("Grays");
		data_and_background = prependToArray(1, data);	//Add background intensity - will become 0 after taking logartihm
		log_data_and_background = logArray(data_and_background);
		Ext.CLIJ2_pushArray(dataImageGPU, log_data_and_background, data_and_background.length, 1, 1);
		Ext.CLIJ2_replaceIntensities(labelmap_cells_montage, dataImageGPU, labelmap_cell_log_intensities);
		Ext.CLIJ2_pull(labelmap_cell_log_intensities);
		rename("Thresholding image for ch"+functionalImmuneChannels[i]+" : "+label[functionalImmuneChannels[i]-1]);
		setBatchMode("show");
		selectWindow("Thresholding image for ch"+functionalImmuneChannels[i]+" : "+label[functionalImmuneChannels[i]-1]);
//		run("Add Image...", "image="+mask_membranes+" x=0 y=0 opacity=30");
		run("Add Image...", "image="+labelmap_classes_montage+" x=0 y=0 opacity=30");
		run("Threshold...");
		setAutoThreshold("MaxEntropy dark no-reset");	//Just to activate the no-reset button
		setMinAndMax(0,40);	//Make it almost invisible, but not too much, because the histogram in the threshold window is based on these values.
		run("Add Image...", "image=[montage_ch"+functionalImmuneChannels[i]+"] x=0 y=0 opacity=50");
		setThreshold(functionalImmuneChannelsThresholdsLog[i], upper);
		updateDisplay();
		setBatchMode("show");

		waitForUser("Automatic log(intensity) threshold for "+label[functionalImmuneChannels[i]-1]+" (ch"+functionalImmuneChannels[i]+") set at "+functionalImmuneChannelsThresholdsLog[i]+". Adjust if necessary, using the threshold window.\nAlso look at the (log)histogram (top-right).");
		getThreshold(functionalImmuneChannelsThresholdsLog[i], upper);
	}
	
	//Draw a red line at the threshold location
	selectImage(histogram);
	setColor("red");
	linePosition = Math.map(functionalImmuneChannelsThresholdsLog[i], 0, 10, 0, 256);
//		Overlay.drawLine(23 + linePosition, 13, 23 + linePosition, 152);
	Overlay.drawLine(20 + linePosition, 11, 20 + linePosition, 138);
	Overlay.show();

	selectImage(dataImageID);
	close();
	//transform the threshold back to linear space
	functionalImmuneChannelsThresholds[i] = exp(functionalImmuneChannelsThresholdsLog[i]) - 1 + min;

	functionalImmuneChannelsColumnNames[i] = "ch"+functionalImmuneChannels[i]+" ("+label[functionalImmuneChannels[i]-1]+")";
	for(k=0; k<nrCells; k++) {
		if(data[k] > functionalImmuneChannelsThresholds[i]) {
			Table.set(functionalImmuneChannelsColumnNames[i], k, 1, CELL_INTENSITY_TABLE);
				PositiveInfunctionalImmuneChannels[k] += 1;
		}
	}
}
Table.update;
print("Intensity thresholds for immune cell markers (natural logarithm):");
Array.print(functionalImmuneChannelsThresholdsLog);
print("Thresholds converted to real intensities:");
Array.print(functionalImmuneChannelsThresholds);
List.set("Immune cell markers thresholds (natural logarithm)", arrayToString(functionalImmuneChannelsThresholdsLog, ","));
List.set("Immune cell markers thresholds", arrayToString(functionalImmuneChannelsThresholds, ","));
//To do: save this in settings/log

//Overlay label classes
selectWindow(montage);
if(overlayClassesMembranes_bool) run("Add Image...", "image="+labelmap_classes_montage+" x=0 y=0 opacity="+OPACITY+" zero");
else run("Add Image...", "image="+mask_membranes+" x=0 y=0 opacity="+OPACITY+" zero");


//Create short labels
label_short = newArray(label.length);
for(i=0; i<label.length; i++) {
	if(indexOf(label[i],"-") != -1) label_short[i] = label[i].substring(0, indexOf(label[i],"-"));
	else label_short[i] = label[i];
}

setFont("SansSerif", LABELFONTSIZE, "antialiased");

setBatchMode(true);
for(i=0; i<nrCells; i++) {
	class = Table.getString("Class", i, CELL_INTENSITY_TABLE);

	subclass = "";
	if(class == "Immune") {
		for(k=0; k<functionalImmuneChannelsColumnNames.length; k++) {
			if(Table.get(functionalImmuneChannelsColumnNames[k], i, CELL_INTENSITY_TABLE) == true) {
				subclass += label_short[functionalImmuneChannels[k]-1] + "+ | ";
			}
		}
		if(subclass != "") subclass = subclass.substring(0, lastIndexOf(subclass, " | "));
//		roiManager("select", i);
//		roiManager("Rename", Roi.getName + " - "+subclass);
		Table.set("Immune subtype", i, subclass, CELL_INTENSITY_TABLE);
		setColor(0,127,255);
	}
	else if(class == "Tumor") {
//		roiManager("select", i);
//		roiManager("Rename", Roi.getName + " - Tumor");
		Table.set("Immune subtype", i, "Tumor", CELL_INTENSITY_TABLE);
		setColor("orange");
	}
	//Draw ROI numbers as overlay in the class color
	Overlay.drawString(i+1, x_montage[i] - LABELFONTSIZE/2, y_montage[i] + LABELFONTSIZE/2);

}
roiManager("deselect");
Table.update;
updateDisplay();


//Update Image statistics table
for (i = 0; i < nrCells; i++) {
	interactions_string = Table.getString("Interaction with (IDs)", i, CELL_INTENSITY_TABLE);
	interactions = string_to_int_array(interactions_string);
	if(interactions.length > 0) Table.set("interacting cell #1", i, Table.getString("Immune subtype", interactions[0]-1, CELL_INTENSITY_TABLE), CELL_INTENSITY_TABLE);
	else Table.set("interacting cell #1", i, "", CELL_INTENSITY_TABLE);
	if(interactions.length > 1) Table.set("interacting cell #2", i, Table.getString("Immune subtype", interactions[1]-1, CELL_INTENSITY_TABLE), CELL_INTENSITY_TABLE);
	else Table.set("interacting cell #2", i, "", CELL_INTENSITY_TABLE);
	if(interactions.length > 2) Table.set("interacting cell #3", i, Table.getString("Immune subtype", interactions[2]-1, CELL_INTENSITY_TABLE), CELL_INTENSITY_TABLE);
	else Table.set("interacting cell #3", i, "", CELL_INTENSITY_TABLE);
	if(interactions.length > 3) Table.set("interacting cell #4", i, Table.getString("Immune subtype", interactions[3]-1, CELL_INTENSITY_TABLE), CELL_INTENSITY_TABLE);
	else Table.set("interacting cell #4", i, "", CELL_INTENSITY_TABLE);
}
selectWindow(CELL_INTENSITY_TABLE);
Table.update;



Table.save(outputFolder + File.separator + baseName + "_"+CELL_INTENSITY_TABLE+".tsv");
roiManager("deselect");

selectWindow(CELL_COUNT_TABLE);
Table.save(outputFolder + File.separator + baseName + "_"+ CELL_COUNT_TABLE +".tsv");

Table.create("Multiplets");
for(i=0; i<=10; i++) {
	Table.set("nr of interactions", i, i, "Multiplets");
	Table.set("count", i, occurencesInArray(nrInteractionsPerImage, i), "Multiplets");
}
Table.update;
Table.save(outputFolder + File.separator + baseName + "_multiplets.tsv");




settings = List.getList;
File.saveString(settings, outputFolder + File.separator + baseName + "_settings.txt");

/*

//Measure interface intensities
selectImage(montage);
run("Select None");
for(i=0; i<functionalImmuneChannelsForClustering.length; i++) {
	Stack.setChannel(functionalImmuneChannelsForClustering[i]);
//	Ext.CLIJ2_pushCurrentSlice(montage);
//	Ext.CLIJ2_statisticsOfLabelledPixels(montage, Image_labelmap);
	run("Duplicate...", "title="+label[functionalImmuneChannelsForClustering[i]-1]+" duplicate channels="+functionalImmuneChannelsForClustering[i]);
	run("Intensity Measurements 2D/3D", "input="+label[functionalImmuneChannelsForClustering[i]-1]+" labels="+labelmap_interfaces+" mean stddev max min median mode skewness kurtosis numberofvoxels volume");
	interfaces_mean = Table.getColumn("Mean", label[functionalImmuneChannelsForClustering[i]-1]+"-intensity-measurements");
	run("Intensity Measurements 2D/3D", "input="+label[functionalImmuneChannelsForClustering[i]-1]+" labels="+labelmap_noninterfaces+" mean stddev max min median mode skewness kurtosis numberofvoxels volume");
	noninterfaces_mean = Table.getColumn("Mean", label[functionalImmuneChannelsForClustering[i]-1]+"-intensity-measurements");
	Table.create(label[functionalImmuneChannelsForClustering[i]-1]);
	Table.setColumn("Interface mean", interfaces_mean, label[functionalImmuneChannelsForClustering[i]-1]);
	Table.setColumn("Non-interface mean", noninterfaces_mean, label[functionalImmuneChannelsForClustering[i]-1]);
//NOTE: doesn't work yet, because the labels have to be matched! (not every cell has an interface)
}

*/

//.CLIJ2_statisticsOfLabelledPixels(montage, labelmap_interfaces);


exit;


Ext.CLIJ2_clear();


Ext.CLIJ2_push(labelmap_classes);

//Create cell membrane masks
Ext.CLIJ2_labelToMask(labelmap_classes, mask_ImCell, 1);
Ext.CLIJ2_labelToMask(labelmap_classes, mask_tumor, 2);
mask_ImCell_eroded = binaryErodeGPU(mask_ImCell, membraneThickness);
Ext.CLIJ2_binarySubtract(mask_ImCell, mask_ImCell_eroded, ImCell_membrane);
mask_tumor_eroded = binaryErodeGPU(mask_tumor, membraneThickness);
Ext.CLIJ2_binarySubtract(mask_tumor, mask_tumor_eroded, tumor_membrane);

//Create interface and Not-in-interface masks
ImCell_interface = "Immune cell interface";
tumor_interface = "Tumor cell interface";
ImCell_membrane_not_in_interface = "Immune cell membrane not in interface";
tumor_membrane_not_in_interface = "Tumor cell membrane not in interface";
//dilate cell masks recursively
mask_ImCell_dilated = binaryDilateGPU(mask_ImCell, membraneThickness);
mask_tumor_dilated = binaryDilateGPU(mask_tumor, membraneThickness);
//Create interface masks and not-in-interface masks
Ext.CLIJ2_binaryIntersection(mask_ImCell, mask_tumor_dilated, ImCell_interface);
Ext.CLIJ2_binaryIntersection(mask_tumor, mask_ImCell_dilated, tumor_interface);
Ext.CLIJ2_binaryXOr(ImCell_interface, ImCell_membrane, ImCell_membrane_not_in_interface);
Ext.CLIJ2_binaryXOr(tumor_interface, tumor_membrane, tumor_membrane_not_in_interface);

//Pull images from GPU and release memory
setBatchMode(true);
Ext.CLIJ2_pull(ImCell_membrane);
Ext.CLIJ2_pull(tumor_membrane);
Ext.CLIJ2_pull(ImCell_interface);
run("Magenta");
Ext.CLIJ2_pull(tumor_interface);
run("Green");
Ext.CLIJ2_pull(ImCell_membrane_not_in_interface);
run("biop-Azure");
Ext.CLIJ2_pull(tumor_membrane_not_in_interface);
run("Yellow");
//Ext.CLIJ2_clear();

//Overlay membranes and interfaces on original image
selectWindow(original);
originalID = getImageID();
//setBatchMode("hide");
selectWindow(ImCell_membrane);
ImCell_membraneID = getImageID();
selectWindow(tumor_membrane);
tumor_membraneID = getImageID();
selectWindow(ImCell_interface);
ImCell_interfaceID = getImageID();
selectWindow(tumor_interface);
tumor_interfaceID = getImageID();
selectWindow(ImCell_membrane_not_in_interface);
ImCell_membrane_not_in_interfaceID = getImageID();
selectWindow(tumor_membrane_not_in_interface);
tumor_membrane_not_in_interfaceID = getImageID();

selectWindow(original);
Overlay.remove;
selectWindow(montage);
Overlay.remove;
setBatchMode("hide");

for (f = 1; f <= frames; f++) {
	selectImage(ImCell_interfaceID);
	setSlice(f);
	selectImage(tumor_interface);
	setSlice(f);
	selectImage(ImCell_membrane_not_in_interface);
	setSlice(f);
	selectImage(tumor_membrane_not_in_interface);
	setSlice(f);

	x = (f-1)%montageWidth;
	y = floor((f-1)/montageWidth);

	addOverlays(ImCell_interface, originalID, f, montageID, x, y);
	addOverlays(tumor_interface, originalID, f, montageID, x, y);
	addOverlays(ImCell_membrane_not_in_interface, originalID, f, montageID, x, y);
	addOverlays(tumor_membrane_not_in_interface, originalID, f, montageID, x, y);
}
selectWindow(montage);
Overlay.show();
setBatchMode("show");

setBatchMode("show");
run("To ROI Manager");	//Copy overlay to ROI manager
Overlay.remove;
run("From ROI Manager");
roiManager("Show All without labels");
/*
for (f = 1; f <= frames; f++) {
	selectImage(originalID);
	Stack.setFrame(f);
	Overlay.copy;
	x = (f-1)%montageWidth;
	y = floor((f-1)/montageWidth);
	selectImage(montageID);
	Overlay.paste;
Overlay.setPosition(1,1,1);
waitForUser(f);
//	Overlay.moveTo(width*x, height*y);
	Overlay.moveSelection(4*(f-1), x*width, y*height);
	Overlay.moveSelection(4*(f-1)+1, x*width, y*height);
	Overlay.moveSelection(4*(f-1)+2, x*width, y*height);
	Overlay.moveSelection(4*(f-1)+3, x*width, y*height);
}
Overlay.show();
*/

selectImage(montageID);
for (f = 1; f <= frames; f++) {
	x = (f-1)%montageWidth;
	y = floor((f-1)/montageWidth);
	Overlay.drawString(f, x*width+width/2, y*height+height);
}



exit;

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
	
	Ext.CLIJ2_copySlice(ImCell_interface, mask2D, i);
	Ext.CLIJ2_getMeanOfMaskedPixels(MOI_slice, mask2D, mean);
	Table.set("mean ImCell interface", i, mean, CELL_COUNT_TABLE);

	Ext.CLIJ2_copySlice(ImCell_membrane_not_in_interface, mask2D, i);
	Ext.CLIJ2_getMeanOfMaskedPixels(MOI_slice, mask2D, mean);
	Table.set("mean ImCell NOT in interface", i, mean, CELL_COUNT_TABLE);

	Ext.CLIJ2_copySlice(tumor_interface, mask2D, i);
	Ext.CLIJ2_getMeanOfMaskedPixels(MOI_slice, mask2D, mean);
	Table.set("mean tumor interface", i, mean, CELL_COUNT_TABLE);

	Ext.CLIJ2_copySlice(tumor_membrane_not_in_interface, mask2D, i);
	Ext.CLIJ2_getMeanOfMaskedPixels(MOI_slice, mask2D, mean);
	Table.set("mean tumor NOT in interface", i, mean, CELL_COUNT_TABLE);
}
meanImCellInterface = Table.getColumn("mean ImCell interface", CELL_COUNT_TABLE);
meanImCellNotInInterface = Table.getColumn("mean ImCell NOT in interface", CELL_COUNT_TABLE);
ImCell_membrane_ratio = divideArrays(meanImCellInterface, meanImCellNotInInterface);
meanTumorInterface = Table.getColumn("mean tumor interface", CELL_COUNT_TABLE);
meanTumorNotInInterface = Table.getColumn("mean tumor NOT in interface", CELL_COUNT_TABLE);
Tumor_membrane_ratio = divideArrays(meanTumorInterface, meanTumorNotInInterface);
Table.setColumn("Immune cell membrane ratio", ImCell_membrane_ratio, CELL_COUNT_TABLE);
Table.setColumn("Tumor membrane ratio", Tumor_membrane_ratio, CELL_COUNT_TABLE);

selectWindow(CELL_COUNT_TABLE);
Table.update;
Table.save(outputFolder+File.separator+CELL_COUNT_TABLE + ".tsv")
selectWindow(CELL_INTENSITY_TABLE);
Table.update;
Table.save(outputFolder+File.separator+CELL_INTENSITY_TABLE + ".tsv")

duplicate_table(CELL_COUNT_TABLE, DOUBLETS_COUNT_TABLE);
Table.showRowIndexes(true);
//Table.setLocationAndSize(0, 0, 0, 0);
for(i=frames-1; i>=0; i--) {
	if(Table.get("Immune cells",i) == 1 && Table.get("Tumor cells",i) == 1 && toString(Table.get("mean ImCell interface",i)) != "NaN") continue;
	else Table.deleteRows(i, i);
}
//Table.setLocationAndSize(0, 0, 800, 800);
graph_allCells = create_interface_intensity_graph(CELL_COUNT_TABLE);
saveAs("PNG", outputFolder+File.separator + getTitle() +".png");
graph_doublets = create_interface_intensity_graph(DOUBLETS_COUNT_TABLE);
saveAs("PNG", outputFolder+File.separator + getTitle() +".png");













//// FUNCTIONS ////

function addOverlays(sourceImageID, targetImageID, frame, montageImageID, x, y) {
	selectImage(targetImageID);
	Stack.setFrame(frame);	
	run("Add Image...", "image=["+sourceImageID+"] x=0 y=0 opacity="+OPACITY+" zero");
	Overlay.show();
	selectImage(montageImageID);
	run("Add Image...", "image=["+sourceImageID+"] x="+x*width+" y="+y*height+" opacity="+OPACITY+" zero");
}


function create_interface_intensity_graph(table) {
	if(CellOfInterest == "Immune cell") {
		Array.getStatistics(Table.getColumn("mean ImCell interface", table), minY, maxY, meanY, stdevY);
		Array.getStatistics(Table.getColumn("mean ImCell NOT in interface", table), minX, maxX, meanX, stdevX);
		Plot.create("Channel "+ch_MOI+" in Immune cell membrane plot - "+table, "mean ImCell NOT in interface", "mean ImCell interface");
		Plot.add("Circle", Table.getColumn("mean ImCell NOT in interface", table), Table.getColumn("mean ImCell interface", table));
		Plot.setStyle(0, "blue,#a0a0ff,1.0,Circle");
		Plot.setColor("red");
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
		Plot.setColor("red");
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
		col_values = Table.getColumn(headings[col], inputTable);
		Table.setColumn(headings[col], col_values, outputTable);
    }
	Table.update(outputTable);
}

//Normalize the image with the median value of Otsu thresholds in every frame
function normalizeImage(image) {
//	getDimensions(width, height, channels, slices, frames);
//	thresholds = newArray(frames);
//	for (i = 0; i < frames; i++) {
//		Stack.setFrame(i+1);
//		setAutoThreshold("Otsu dark");
//		getThreshold(thresholds[i], upper);
//	}
//	resetThreshold;
//	median = medianOfArray(thresholds);
	setAutoThreshold("Otsu dark stack");	//Threshold on the stack - better results in case only a few cells are positive for this staining
	getThreshold(lower, upper);
	
	selectWindow(image);
	if(bitDepth() != 32) run("32-bit");
//	run("Divide...", "value="+median+" stack");
	run("Divide...", "value="+lower+" stack");
	run("Enhance Contrast", "saturated=0.35");
	print("Normalizing "+image+" with threshold: "+lower);
}



//Returns a 32-bit binary eroded image with a number of iterations, in GPU memory. The input binaryImage should already be in GPU memory. 
function binaryErodeGPU(binaryImage, iterations) {
	image_in = binaryImage;
	for (i = 0; i < iterations; i++) {
		image_out = binaryImage+i+1;
		Ext.CLIJ2_erodeSphere(image_in, image_out);
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
		Ext.CLIJ2_dilateSphere(image_in, image_out);
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

//Adds a scalar to all elements of an array
function addScalarToArray(array, scalar) {
	added_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		added_array[a]=array[a] + scalar;
	}
	return added_array;
}

//Returns, as array, all values in the array return_array at the same indices as lookup_value in the ref_array
function lookup(ref_array, lookup_value, return_array) {
	indices = indexOfArray(ref_array, lookup_value);
	return_values = newArray(indices.length);
	for (i=0; i<indices.length; i++) {
		return_values[i] = return_array[indices[i]];
	}
	return return_values;
}

//Returns, as array, the indices at which a value occurs within an array
function indexOfArray(array, value) {
	count=0;
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			count++;
		}
	}
	if (count>0) {
		indices=newArray(count);
		count=0;
		for (a=0; a<lengthOf(array); a++) {
			if (array[a]==value) {
				indices[count]=a;
				count++;
			}
		}
		return indices;
	}
	return newArray(0);
}

//Converts an array into a string, elements separated by 'separator'
function arrayToString(array, separator) {
	outputString = "";
	for (i = 0; i < array.length; i++) {
		outputString += toString(array[i]) + separator;
	}
	return substring(outputString, 0, outputString.length - separator.length);
}

//Converts a string of numbers into an array with integers
function string_to_int_array(string) {
	array = split(string, ",");
	intArray = newArray(array.length);
	for (i = 0; i < array.length; i++) intArray[i] = parseInt(array[i].trim());
	return intArray;
}

//Convert a color into a hexadecimal code
function color_to_hex(color) {
	colorArray = split(color,",,");
	hexcolor = "#" + IJ.pad(toHex(colorArray[0]),2) + IJ.pad(toHex(colorArray[1]),2) + IJ.pad(toHex(colorArray[2]),2);
	return hexcolor;
}

//Returns the sum of all elements of an arrays, neglecting NaNs
function sumArray(array) {
	sum=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
	}
	return sum;
}

//Divides the elements of two arrays and returns the new array
function divideArrays(array1, array2) {
	divArray=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		divArray[a]=array1[a]/array2[a];
	}
	return divArray;
}

//Adds two arrays of equal length element-wise
function addArrays(array1, array2) {
	added_array=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		added_array[a]=array1[a] + array2[a];
	}
	return added_array;
}

//Subtract a scalar from all elements of an array
function subtract_scalar_from_array(array, scalar) {
	subtracted_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		subtracted_array[a]=array[a] - scalar;
	}
	return subtracted_array;
}

//Divides all elements of an array by a scalar
function divideArraybyScalar(array, scalar) {
	divided_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		divided_array[a]=array[a]/scalar;
	}
	return divided_array;
}

//Normalizes the array to mean 0 and stddev 1
function normalizeArrayToMedian0AndStddev1(array, percentile) {
	array_sorted = Array.copy(array);
	Array.sort(array_sorted);
	array_sorted = Array.slice(array_sorted, 0, round(array_sorted.length*percentile/100));
//	Array.show(array, array2);
	median = array_sorted[floor(array_sorted.length/2)];
	Array.getStatistics(array_sorted, min, max, mean, stdDev);
	tempArray = subtract_scalar_from_array(array, median);
	normalizedArray = divideArraybyScalar(tempArray, stdDev);
	return normalizedArray;
}

//Puts the values in a 1D 32-bit image and returns the Image ID
function arrayTo1DImage(array) {
	newImage("1D Image", "32-bit", array.length, 1, 1);
	imageID = getImageID();
	for(i=0; i<array.length; i++) {
		setPixel(i, 0, array[i]);
	}
	updateDisplay();
	return imageID;
}

//Returns the values in a 1D image (m,1,1) or (1,m,1) as an array
function image1DToArray(image) {
	getDimensions(width, height, channels, slices, frames);
	array = newArray(maxOf(width, height));
	if(height == 1) {
		for(x=0; x<width; x++) {
			array[x] = getPixel(x, 0);
		}
	}
	else if (width == 1) {
		for(y=0; y<height; y++) {
			array[y] = getPixel(0, y);
		}
	}
	else exit("Error in function 'image1DToArray': 1D image expected");
	return array;
}

//get the automatic threshold value of the values in an array, as if it was an image
function getThresholdValueOfArray(array, threshold_method, bool_noReset, min, max) {
	imageID = arrayTo1DImage(array);
	selectImage(imageID);
	setMinAndMax(min, max);
	if(bool_noReset) setAutoThreshold(threshold_method+" dark no-reset");
	else setAutoThreshold(threshold_method+" dark");
	getThreshold(lower, upper);
	close();
	return lower;
}

//Returns an array with unique values in the array
function uniqueValuesInArray(array) {
	count=0;
	uniqueArray = newArray(array.length);
	for (a=0; a<lengthOf(array); a++) {
		if(!occursInArray(Array.trim(uniqueArray, count), array[a]) && !matches(array[a],"None")) {
			uniqueArray[count]=array[a];
			count++;
		}
	}
	return Array.trim(uniqueArray, count);
}

//Returns the number of times the value occurs within the array
function occurencesInArray(array, value) {
	count=0;
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			count++;
		}
	}
	return count;
}


//Returns the minimum of the array
function minOfArray(array) {
	max=0;
	for (a=0; a<lengthOf(array); a++) {
		max=maxOf(array[a], max);
	}
	min=max;
	for (a=0; a<lengthOf(array); a++) {
		min=minOf(array[a], min);
	}
	return min;
}


//Returns the maximum of the array
function maxOfArray(array) {
	min=0;
	for (a=0; a<lengthOf(array); a++) {
		min=minOf(array[a], min);
	}
	max=min;
	for (a=0; a<lengthOf(array); a++) {
		max=maxOf(array[a], max);
	}
	return max;
}


//Prepends the value to the array
function prependToArray(value, array) {
	temparray=newArray(lengthOf(array)+1);
	for (i=0; i<lengthOf(array); i++) {
		temparray[i+1]=array[i];
	}
	temparray[0]=value;
	array=temparray;
	return array;
}


//Returns an array containing the 10log values of all elements
function logArray(array) {
	log_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		log_array[a]= Math.log(array[a]);
	}
	return log_array;
}
