import numpy as np
import pytest

from conftest import IDEAS_CHANNELS, IDEAS_IDS, make_events, write_ideas_zip
from imagestream_analysis import ideas


def test_scan_finds_datasets_objects_and_channels(ideas_zip):
    info = ideas.scan(ideas_zip["zip"])
    assert list(info) == ["Sample A"]
    ds = info["Sample A"]
    assert ds.channels == IDEAS_CHANNELS
    assert ds.objects() == sorted(IDEAS_IDS)


def test_read_dataset_orders_by_object_and_keeps_the_numbers(ideas_zip):
    stack = ideas.read_dataset(ideas_zip["zip"], "Sample A", [1, 6, 9], ["BF", "T", "Tu"], pixel_size=0.5)
    assert stack.object_ids.tolist() == sorted(IDEAS_IDS)
    assert stack.channel_numbers == [1, 6, 9] and stack.channel_names == ["BF", "T", "Tu"]
    assert stack.pixel_size == 0.5 and stack.n_events == len(IDEAS_IDS)
    order = ideas_zip["order"]
    assert np.array_equal(stack.channel(6), ideas_zip["data"][order, 2].astype(np.float32))
    assert np.array_equal(stack.channel(9), ideas_zip["data"][order, 3].astype(np.float32))
    with pytest.raises(ValueError, match="not in"):
        stack.channel(3)  # not loaded


def test_smaller_images_are_padded_with_the_background(tmp_path):
    data, _ = make_events(4)
    path = write_ideas_zip(tmp_path / "x.zip", data, [1, 2, 3, 4], [1, 2, 3, 4], shrink=(2,))
    stack = ideas.read_dataset(path, "Sample A", [1, 2])
    assert stack.data.shape[2:] == data.shape[2:]
    bf = stack.channel(1)[2]
    assert np.array_equal(bf[5:-5, 5:-5], data[2, 0, 5:-5, 5:-5])
    assert (bf[:5] == round(stack.medians[0])).all()  # the border: the median background of the channel


def test_objects_without_all_channels_are_left_out(tmp_path):
    import zipfile

    data, _ = make_events(3)
    path = write_ideas_zip(tmp_path / "x.zip", data, [10, 20, 30], [1, 2, 3, 4])
    with zipfile.ZipFile(path) as z:  # rewrite without Ch3 of object 20
        files = {n: z.read(n) for n in z.namelist() if n != "Sample A/20_Ch3.ome.tif"}
    with zipfile.ZipFile(path, "w") as z:
        for n, b in files.items():
            z.writestr(n, b)
    assert ideas.read_dataset(path, "Sample A", [1, 3]).object_ids.tolist() == [10, 30]
    assert ideas.read_dataset(path, "Sample A", [1, 2]).object_ids.tolist() == [10, 20, 30]



def test_file_names_with_a_prefix(tmp_path):
    """IDEAS can put a prefix before the object number (<prefix>_<object>_Ch<N>, often the name of the data file).
    One prefix: the folder is the dataset. Several prefixes in one folder (several data files exported together,
    whose object numbers may clash): one dataset per prefix."""
    import zipfile

    data, _ = make_events(3)
    path = write_ideas_zip(tmp_path / "x.zip", data, [10, 20, 30], [1, 2, 3, 4])
    with zipfile.ZipFile(path) as z:
        files = {n: z.read(n) for n in z.namelist()}
    with zipfile.ZipFile(path, "w") as z:
        for n, b in files.items():
            z.writestr(n.replace("Sample A/", "Sample A/40_A1D3_2_"), b)
    info = ideas.scan(path)
    assert list(info) == ["Sample A"] and info["Sample A"].objects() == [10, 20, 30]
    stack = ideas.read_dataset(path, "Sample A", [1, 2])
    assert stack.object_ids.tolist() == [10, 20, 30] and np.array_equal(stack.channel(1), data[:, 0].astype(np.float32))
    with zipfile.ZipFile(path, "w") as z:  # a second export in the same folder, with the same object numbers
        for n, b in files.items():
            z.writestr(n.replace("Sample A/", "Sample A/40_A1D3_2_"), b)
            z.writestr(n.replace("Sample A/", "Sample A/5_B2_4_"), b)
    info = ideas.scan(path)
    assert list(info) == ["Sample A - 40_A1D3_2", "Sample A - 5_B2_4"]
    assert all(d.objects() == [10, 20, 30] for d in info.values())
    assert ideas.read_dataset(path, "Sample A - 5_B2_4", [1]).n_events == 3


def test_max_events_and_errors(ideas_zip):
    stack = ideas.read_dataset(ideas_zip["zip"], "Sample A", [1], max_events=3)
    assert stack.object_ids.tolist() == sorted(IDEAS_IDS)[:3] and stack.n_events_in_file == len(IDEAS_IDS)
    with pytest.raises(ValueError, match="no dataset"):
        ideas.read_dataset(ideas_zip["zip"], "nope", [1])
    with pytest.raises(ValueError, match="no channel"):
        ideas.read_dataset(ideas_zip["zip"], "Sample A", [5])


def test_thumbnails_per_channel(ideas_zip):
    thumbs = ideas.thumbnails(ideas_zip["zip"], "Sample A", n=4)
    assert sorted(thumbs) == IDEAS_CHANNELS
    assert all(url.startswith("data:image/png;base64,") for url in thumbs.values())


@pytest.mark.parametrize("dtype", ["uint8", "uint16", "int16", "float32"])
@pytest.mark.parametrize("byteorder", ["<", ">"])
def test_the_plain_tiff_reader_gives_what_tifffile_gives(dtype, byteorder):
    """IDEAS writes plain TIFFs: read straight from the first IFD (fast), with tifffile for anything else."""
    import io

    import tifffile

    from imagestream_analysis.ideas import _plain_tiff

    rng = np.random.default_rng(1)
    image = (rng.random((37, 23)) * 200).astype(dtype)
    buf = io.BytesIO()
    tifffile.imwrite(buf, image, byteorder=byteorder, rowsperstrip=5)  # several strips
    got = _plain_tiff(buf.getvalue())
    assert got is not None and got.dtype == image.dtype and np.array_equal(got, image)
    packed = io.BytesIO()
    tifffile.imwrite(packed, image, compression="zlib")
    assert _plain_tiff(packed.getvalue()) is None  # compressed: left to tifffile


def test_parts_of_a_sample_keep_their_own_background(tmp_path):
    """ideas.combine: a second part with a higher background and smaller images is shifted to the common
    background, so that the values above background stay its own."""
    data, _ = make_events(3)
    a = write_ideas_zip(tmp_path / "a.zip", data, [1, 2, 3], [1, 2, 3, 4])
    b = write_ideas_zip(tmp_path / "b.zip", data[:2] + np.uint16(40), [1, 2], [1, 2, 3, 4], shrink=(0, 1))
    sa = ideas.read_dataset(a, "Sample A", [1, 2], ["BF", "M"])
    sb = ideas.read_dataset(b, "Sample A", [1, 2], ["BF", "M"])
    both = ideas.combine([sa, sb], "both", ["a", "b"])
    assert both.n_events == 5 and both.name == "both" and both.object_ids.tolist() == [1, 2, 3, 1, 2]
    assert both.sources.tolist() == [0, 0, 0, 1, 1] and both.source_names == ["a", "b"]
    assert both.medians == [max(x, y) for x, y in zip(sa.medians, sb.medians, strict=True)]
    for i in range(2):  # above background: unchanged in each part
        above_a = both.data[:3, i].astype(float) - both.medians[i]
        assert np.allclose(above_a, sa.data[:, i].astype(float) - sa.medians[i], atol=0.5)
        y0, x0, y1, x1 = both.extents[3]
        above_b = both.data[3, i, y0:y1, x0:x1].astype(float) - both.medians[i]
        assert np.allclose(above_b, sb.data[0, i].astype(float) - sb.medians[i], atol=0.5)
    assert (both.extents[3] == [5, 5, data.shape[2] - 5, data.shape[3] - 5]).all()
    with pytest.raises(ValueError, match="different channels"):
        ideas.combine([sa, ideas.read_dataset(b, "Sample A", [1, 3], ["BF", "X"])], "x", ["a", "b"])
