# Input

## IDEAS exports

The input is an export from IDEAS: a zip file with one folder per **dataset**. Each folder holds one image per object
and channel, named `<object>_Ch<N>.ome.tif`. Every object (an ImageStream *event*) becomes one image in the analysis.
The images are ordered by object number, and the object number is kept with every result in the `Object ID` column.

Only the channels that have been given a name are read. An object that lacks an image for one of these channels is
left out, and the log reports how many objects were left out.

Which objects the export contains is decided in IDEAS, by the gates and by the population that was exported. The
analysis sees no other objects. This matters when you interpret the
[cell contacts](interactions.md#what-the-selection-in-ideas-means).

??? note "Technical: image size, background and padding"

    In IDEAS, the image size differs from object to object. Every image is therefore centered in an image the size of
    the largest one, and the border around it is filled with the channel's **background**. The background is computed
    as follows: first the median of each object's image, and then the median of those values over all objects. The
    same background is subtracted from the mean intensity of each cell (`Mean ch<N>`), and it sets the lower end of the
    display range in the dashboard. These are the same rules that ImageStreamCombiner (BioImaging-NKI) uses, so an
    export that is read directly gives the same images as a Combiner stack.

    The extent of each object's own image is kept, so that cells cut by the edge of their image can be detected (see
    [Segmentation](segmentation.md#cells-at-the-image-edge)).

## ImageStreamCombiner stacks

A `.tif` file made by [ImageStreamCombiner](https://github.com/BioImaging-NKI/ImageStreamCombiner), with one frame per
event, can also be analyzed. It is given in a settings file (`input.image`). The channels are then numbered by their
position in the stack.

??? note "Technical: what a Combiner stack lacks"

    A Combiner stack orders the events in the order in which the files were listed in the zip, and it does not keep the
    object numbers. Its results can therefore not be traced back to IDEAS objects, and "the first N events" of a
    stack are different events from the first N objects of the export. The padding is recognized as the outer rows
    and columns that contain exactly the background value in every channel.

## Channels, roles and pixel size

Every channel that the analysis needs is given a name and a role: brightfield, tumor marker, immune marker(s), subtype
markers or marker of interest. See [The import page](../guide/import-page.md#2-ideas-exports) for how to set them. The
pixel size (in µm) is set for each export, because IDEAS files do not contain it. It is used to convert areas to µm²
and positions along the membrane to µm.

## The first N events

With *Only the first N events*, only the first N objects (by object number) of each dataset are analyzed. Use it for a
quick test of the settings.
