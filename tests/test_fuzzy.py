from geolab.tui.fuzzy import fuzzy_filter, fuzzy_match

NAMES = [
    "Geodetic → ECEF",
    "ECEF → Geodetic",
    "Inverse problem (Vincenty)",
    "Direct problem (Vincenty)",
    "Ellipsoid",
    "2D ellipse",
    "Meridian and parallel arcs",
]


def test_subsequence_required():
    assert fuzzy_match("xyz", "Ellipsoid") is None
    assert fuzzy_match("lpd", "Ellipsoid") is not None


def test_case_and_spaces_ignored():
    assert fuzzy_match("IN v", "inverse") is not None


def test_empty_query_matches_everything_in_order():
    assert [m.index for m in fuzzy_filter("", NAMES)] == list(range(len(NAMES)))


def test_prefix_beats_later_match():
    assert [NAMES[m.index] for m in fuzzy_filter("ell", NAMES)][:2] == ["Ellipsoid", "2D ellipse"]


def test_word_starts_preferred():
    # "pa" should land on "parallel", not the "a" in "Meridian".
    score, positions = fuzzy_match("mpa", "Meridian and parallel arcs")
    assert positions == (0, 13, 14)


def test_consecutive_beats_scattered():
    a = fuzzy_match("ecef", "ECEF → Geodetic")
    b = fuzzy_match("ecef", "Geodetic → ECEF")
    assert a[0] > b[0]
    assert b[1] == (11, 12, 13, 14)


def test_positions_are_the_best_alignment():
    # The greedy first "e" would be at 0; the best match uses the whole word "ellipse".
    assert fuzzy_match("ellipse", "2D ellipse")[1] == (3, 4, 5, 6, 7, 8, 9)
