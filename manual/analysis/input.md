# Input

## IDEAS exports

The input is an export from IDEAS: a zip with one folder per **dataset**, and in it one image per object and channel,
`<object>_Ch<N>.ome.tif`. Every object (an ImageStream *event*) becomes one image in the analysis, in the order of the
object numbers, and the object number is kept with every result (`Object ID`).

Only the channels with a name are read. An object that misses one of them is left out (the log says how many).

Which objects are in the export is decided in IDEAS: the gates and the population that was exported. The analysis
sees nothing else, which matters for the [cell contacts](interactions.md#what-the-selection-in-ideas-means).

??? note "Technical: image size, background and padding"

    IDEAS images differ in size per object. Every image is centred in the size of the largest one, and the border is
    filled with the channel's **background**: the median over the objects of each object's median. The background is
    also what is subtracted from the cell means (`Mean ch<N>`), and the lower end of the display in the dashboard.
    These are the same rules as ImageStreamCombiner (BioImaging-NKI), so an export read directly gives the same
    images as a Combiner stack.

    The extent of every object's own image is kept, so that cells cut by the edge of their image can be found
    ([Segmentation](segmentation.md#cells-at-the-image-edge)).

## ImageStreamCombiner stacks

A `.tif` made by [ImageStreamCombiner](https://github.com/BioImaging-NKI/ImageStreamCombiner) (one frame per event)
can be analysed with a settings file (`input.image`). Channels are then numbered by their position in the stack.

??? note "Technical: what a Combiner stack lacks"

    A Combiner stack orders the events as the files were listed in the zip and drops the object numbers, so its
    results cannot be traced back to IDEAS objects, and "the first N events" are other events than in the export.
    The padding is recognised as the outer rows and columns that hold exactly the background in every channel.

## Channels, roles and pixel size

Every channel that is needed gets a name and a role (brightfield, tumor marker, immune marker(s), subtype markers,
marker of interest): see [The import page](../guide/import-page.md#2-ideas-exports). The pixel size (µm) is set per
export; IDEAS files do not contain it. It converts areas to µm² and positions along the membrane to µm.

## The first N events

*Only the first N events* analyses only the first N objects (by object number) of every dataset, for a quick test of
the settings.
