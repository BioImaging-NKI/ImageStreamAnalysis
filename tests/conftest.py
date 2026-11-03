"""Synthetic ImageStream data: every event has one tumor cell (large) touching one T cell (small)."""

import numpy as np
import pytest
import tifffile

H, W = 60, 60
CHANNELS = ["BF", "MOI", "Tumor", "T cell"]
MEDIANS = [800, 10, 20, 20]


def disk(cy, cx, r):
    yy, xx = np.mgrid[:H, :W]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r


def make_events(n=10, seed=0, second_marker=False):
    """Returns data (T, C, Y, X) uint16 and the labels (T, Y, X): tumor label 2t+1, T cell label 2t+2.
    second_marker: a fifth channel (SECOND) with a second marker of interest, on the T cells (enriched at their side
    of the interface), where the first one is on the tumor cells."""
    rng = np.random.default_rng(seed)
    data = np.zeros((n, len(CHANNELS) + second_marker, H, W), dtype=np.float64)
    labels = np.zeros((n, H, W), dtype=np.int32)
    xx = np.mgrid[:H, :W][1]
    for t in range(n):
        r_tumor = 12 + t % 3
        tumor = disk(30, 22, r_tumor)
        t_cell = disk(30, 22 + r_tumor + 8, 8) & ~tumor
        labels[t][tumor] = 2 * t + 1
        labels[t][t_cell] = 2 * t + 2
        for c, m in enumerate(MEDIANS + ([SECOND[1]] if second_marker else [])):
            data[t, c] = m + rng.normal(0, 2, (H, W))
        data[t, 0][tumor | t_cell] -= 200  # brightfield: cells darker
        data[t, 2][tumor] += 400 + 20 * t
        data[t, 3][t_cell] += 600 + 10 * t
        contact = tumor & (np.abs(xx - (22 + r_tumor)) <= 2)
        data[t, 1][tumor] += 50
        data[t, 1][contact] += 150  # marker of interest enriched at the tumor side of the interface
        if second_marker:
            data[t, 4][t_cell] += 60
            data[t, 4][t_cell & (np.abs(xx - (22 + r_tumor)) <= 3)] += 120  # at the T-cell side of the interface
    return np.clip(data, 0, 65535).astype(np.uint16), labels


SECOND = ("Second", 15)  # the second marker of interest of make_events(second_marker=True): name, background


def write_stack(path, data, names):
    medians = MEDIANS + [SECOND[1]] * (len(names) - len(MEDIANS))
    tifffile.imwrite(path, data, imagej=True, resolution=(2.0, 2.0),
                     metadata={"axes": "TCYX", "unit": "um", "Labels": names * len(data),
                               "Properties": {"Medians": "\n" + ",".join(str(m) for m in medians)}})
    return path


@pytest.fixture
def synthetic(tmp_path):
    data, labels = make_events()
    image = write_stack(tmp_path / "events.tif", data, CHANNELS)
    labelmap = tmp_path / "labels.tif"
    tifffile.imwrite(labelmap, labels.astype(np.uint16), imagej=True, metadata={"axes": "TYX"})
    return {"image": image, "labelmap": labelmap, "data": data, "labels": labels, "tmp": tmp_path}


@pytest.fixture
def synthetic_two_markers(synthetic):
    """As synthetic, with a fifth channel: a second marker of interest, on the T cells (the first is on the tumor)."""
    data, _ = make_events(second_marker=True)
    image = write_stack(synthetic["tmp"] / "events_two_markers.tif", data, [*CHANNELS, SECOND[0]])
    return synthetic | {"image": image, "data": data}


def write_ideas_zip(path, data, object_ids, channels, dataset="Sample A", shrink=()):
    """An IDEAS-like export: <dataset>/<object>_Ch<N>.ome.tif for every event and channel. data: (T, C, Y, X), the
    channels numbered `channels`; events in `shrink` are written 10 px smaller (IDEAS images differ in size)."""
    import io
    import zipfile

    with zipfile.ZipFile(path, "w") as z:
        for t, obj in enumerate(object_ids):
            for i, ch in enumerate(channels):
                im = data[t, i][5:-5, 5:-5] if t in shrink else data[t, i]
                buf = io.BytesIO()
                tifffile.imwrite(buf, im, ome=True)
                z.writestr(f"{dataset}/{obj}_Ch{ch}.ome.tif", buf.getvalue())
    return path


IDEAS_CHANNELS = [1, 3, 6, 9]  # BF, MOI, T cell, tumor in the IDEAS numbering (arbitrary)
IDEAS_IDS = [512, 7, 90, 33, 1001, 64, 3, 250, 18, 777]  # not in order, as IDEAS object numbers are not


@pytest.fixture
def ideas_zip(tmp_path):
    """The synthetic events as an IDEAS zip, with the labels in object-number order."""
    data, labels = make_events()
    # conftest channel order: BF, MOI, Tumor, T cell -> IDEAS Ch1, Ch3, Ch9, Ch6
    data = data[:, [0, 1, 3, 2]]
    path = write_ideas_zip(tmp_path / "export.zip", data, IDEAS_IDS, IDEAS_CHANNELS)
    order = np.argsort(IDEAS_IDS)
    labelmap = tmp_path / "labels_by_object.tif"
    tifffile.imwrite(labelmap, labels[order].astype(np.uint16), imagej=True, metadata={"axes": "TYX"})
    return {"zip": path, "data": data, "labels": labels, "order": order, "labelmap": labelmap, "tmp": tmp_path}
