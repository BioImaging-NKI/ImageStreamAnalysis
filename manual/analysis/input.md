# Input

## IDEAS exports

The input is an export from IDEAS: a zip file with one folder per **dataset**. Each folder holds one image per object
and channel, named `<object>_Ch<N>.ome.tif`. If a file-name prefix was set on export, the names read
`<prefix>_<object>_Ch<N>.ome.tif` (the prefix is often the name of the data file); the object number is then the last
number before `_Ch`. Every object (an ImageStream *event*) becomes one image in the analysis. The images are ordered
by object number, and the object number is kept with every result in the `Object ID` column.

If one folder holds files with **different prefixes**, for example the exports of two data files put into the same
folder, each prefix becomes a dataset of its own, named *folder - prefix*. Their object numbers can be the same, so
they are never mixed. A small extra dataset of this kind is often a few files copied into the wrong folder: check it,
and untick it on the import page if it does not belong there.

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

## Several exports of one sample

Sometimes one sample ends up in more than one export: it was acquired in two runs, or exported in two parts. These
parts are not replicates; together they are one dataset. Give them the same name under *Combine as* on the
[import page](../guide/import-page.md#2-ideas-exports) (they may come from different zips), and they are read and
analysed as **one dataset** of that name: one result folder, one dashboard and one row in the overview, with every
number computed on all their cells.

- The parts need the same channel names and the same pixel size, and they belong to one experiment, condition and
  replicate. The import page lists what does not fit, and the run does not start until it is fixed.
- Each part keeps its own **background**: a second acquisition can have a slightly different background, and every
  event is corrected with the background of its own part.
- Object numbers can repeat between parts, so the tables get a column `Source` (the zip and dataset of each event)
  next to `Object ID`. In the dashboard, the *Part* menu shows the events of one part at a time, so you can check
  that the parts look alike.

??? note "Technical: how the parts are joined"

    Each part is read as described above, with its own background per channel. The events of the parts follow each
    other (part after part, each in order of object number). The pixels of each part are shifted by the difference
    between the highest background of the parts and its own, so that all parts share one background. Every
    measurement subtracts the background, so the values above background stay those of the part. Segmentation,
    classification and every other step then see one dataset. `input.max_events` counts the events of all parts
    together, starting with the first part.

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
markers or markers of interest. See [The import page](../guide/import-page.md#2-ideas-exports) for how to set them. The
pixel size (in µm) is set for each export, because IDEAS files do not contain it. It is used to convert areas to µm²
and positions along the membrane to µm.

## The first N events

With *Only the first N events*, only the first N objects (by object number) of each dataset are analyzed. Use it for a
quick test of the settings.
