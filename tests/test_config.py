import pytest

from geolab import Ellipsoid
from geolab.config import config_dir, ellipsoids_path, load_ellipsoids, save_ellipsoids


def test_config_dir_follows_xdg(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    assert config_dir() == tmp_path / "geolab"
    assert ellipsoids_path() == tmp_path / "geolab" / "ellipsoids.toml"


def test_config_dir_defaults_to_dot_config(monkeypatch, tmp_path):
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    assert config_dir() == tmp_path / ".config" / "geolab"


def test_missing_file_is_empty(tmp_path):
    assert load_ellipsoids(tmp_path / "nope.toml") == ([], [])


def test_round_trip(tmp_path):
    path = tmp_path / "sub" / "ellipsoids.toml"
    ellipsoids = [
        Ellipsoid("Mine", 6378000.0, 298.257223563),
        Ellipsoid('Odd "name" ğ', 6377000.5, 300.0),
        Ellipsoid("Ball", 6371000.0, 0.0),
    ]
    save_ellipsoids(ellipsoids, path)
    assert load_ellipsoids(path) == (ellipsoids, [])


def test_saved_file_is_readable(tmp_path):
    path = tmp_path / "ellipsoids.toml"
    save_ellipsoids([Ellipsoid("Mine", 6378000.0, 298.25)], path)
    text = path.read_text()
    assert text.startswith("# Custom ellipsoids saved by geolab")
    assert '[[ellipsoid]]\nname = "Mine"\na = 6378000.0\ninv_f = 298.25\n' in text


def test_save_empty_list(tmp_path):
    path = tmp_path / "ellipsoids.toml"
    save_ellipsoids([], path)
    assert load_ellipsoids(path) == ([], [])


def test_save_leaves_no_temp_file(tmp_path):
    save_ellipsoids([Ellipsoid("Mine", 6378000.0, 298.25)], tmp_path / "e.toml")
    assert [p.name for p in tmp_path.iterdir()] == ["e.toml"]


def test_hand_written_integers_accepted(tmp_path):
    path = tmp_path / "e.toml"
    path.write_text('[[ellipsoid]]\nname = "Int"\na = 6378000\ninv_f = 300\n')
    assert load_ellipsoids(path) == ([Ellipsoid("Int", 6378000.0, 300.0)], [])


def test_malformed_toml(tmp_path):
    path = tmp_path / "e.toml"
    path.write_text("[[ellipsoid]\nname = ")
    ellipsoids, warnings = load_ellipsoids(path)
    assert ellipsoids == []
    assert len(warnings) == 1 and "could not read" in warnings[0]


def test_wrong_top_level_type(tmp_path):
    path = tmp_path / "e.toml"
    path.write_text('ellipsoid = "x"\n')
    assert load_ellipsoids(path)[0] == []
    assert "array of tables" in load_ellipsoids(path)[1][0]


@pytest.mark.parametrize(
    ("entry", "problem"),
    [
        ('a = 6378000.0\ninv_f = 300.0', "missing name"),
        ('name = "X"\ninv_f = 300.0', "a must be a number"),
        ('name = "X"\na = "big"\ninv_f = 300.0', "a must be a number"),
        ('name = "X"\na = true\ninv_f = 300.0', "a must be a number"),
        ('name = "X"\na = 6378000.0\ninv_f = 0.5', "inverse flattening"),
        ('name = "X"\na = -1.0\ninv_f = 300.0', "semi-major"),
    ],
)
def test_invalid_entries_are_skipped(tmp_path, entry, problem):
    path = tmp_path / "e.toml"
    path.write_text(f'[[ellipsoid]]\n{entry}\n\n[[ellipsoid]]\nname = "Good"\na = 6378000.0\ninv_f = 300.0\n')
    ellipsoids, warnings = load_ellipsoids(path)
    assert [e.name for e in ellipsoids] == ["Good"]
    assert len(warnings) == 1
    assert "#1" in warnings[0] and problem in warnings[0]
