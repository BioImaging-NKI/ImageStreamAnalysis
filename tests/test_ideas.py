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
