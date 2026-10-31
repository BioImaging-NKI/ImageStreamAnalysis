import numpy as np

from imagestream_analysis import imagej


def test_otsu_bin_splits_two_peaks():
    hist = np.zeros(256)
    hist[40:50] = 100
    hist[200:210] = 100
    k = imagej.otsu_bin(hist)
    assert 49 <= k < 200


def test_otsu_threshold_dark_between_classes():
    img = np.full((20, 20), 100.0, dtype=np.float32)
    img[5:15, 5:15] = 1000.0
    thr = imagej.otsu_threshold_dark(img)
    assert 100 < thr <= 1000


def test_circular_kernel_radius_2_is_imagej_shape():
    k = imagej.circular_kernel(2)
    assert k.shape == (5, 5)
    assert k.sum() == 21  # 5x5 without the corners
    assert not k[0, 0] and k[0, 1] and k[2, 0]


def test_variance_filter_constant_and_step():
    assert np.allclose(imagej.variance_filter(np.full((10, 10), 7.0)), 0)
    img = np.zeros((10, 10))
    img[:, 5:] = 10
    v = imagej.variance_filter(img)
    assert v[5, 0] == 0 and v[5, 4] > 0 and v.dtype == np.float32


def test_erode_keeps_image_border_and_uses_cross():
    m = np.ones((5, 5), bool)
    assert imagej.erode(m, 1).all()  # outside the image counts as the edge pixel
    m = np.zeros((5, 5), bool)
    m[1:4, 1:4] = True
    e = imagej.erode(m, 1)
    assert e.sum() == 1 and e[2, 2]
    d = imagej.dilate(e, 1)
    assert d.sum() == 5  # von Neumann neighbourhood


def test_close_index_gaps_keeps_order():
    labels = np.array([[0, 7, 7], [3, 0, 12]])
    out = imagej.close_index_gaps(labels)
    assert out.tolist() == [[0, 2, 2], [1, 0, 3]]


def test_median_upper():
    assert imagej.median_upper([4, 1, 3, 2]) == 3
    assert imagej.median_upper([5, 1, 3]) == 3
