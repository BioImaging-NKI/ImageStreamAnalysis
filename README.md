# ImageStream Analysis
Scripts for the analysis of tumor-immune cell interactions in Imaging Flow Cytrometry data.

**Note**: This GitHub site is currently mostly a placeholder. Instructions on how to run them will follow 'soon'.
The two `.ijm` Fiji macros currently have different functionalities, but will at some point be merged into a single script.
`Analyze_ImageStream_stack_2celltypes_+_membrane_intensity_analysis_1.5.ijm` can be used to separate tumor cells from immune cells and to quantify 'marker of interest' enrichment at the interface between these cells and visualize this in various ways.
`ClassifyCellTypes_2_8.ijm` is a more general script that works for more cell types. First, tumor cells and immune cells are determined using K-means clustering. Next, immune cell marker abundance is classified as positive or nagative, compared to a user-set or automatic threshold level. This script currently lacks the membrane enrichment quantification and visualizations.
Finally, the Jupyter Notebook `Summarize_Interactions.ipynb` can be used to quantify interactions between different cells types, for example:

![image](https://github.com/user-attachments/assets/434d8ce4-383f-419c-94bb-7b0dfee70e2d)

For combining exported ImageStream images to create input tiff files for the macros, see our [ImageStreamCombiner repository](https://github.com/BioImaging-NKI/ImageStreamCombiner).
