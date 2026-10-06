import sys

import pytest


def test_main_without_textual_explains_how_to_install(monkeypatch):
    import geolab.tui

    monkeypatch.delitem(sys.modules, "geolab.tui.app", raising=False)
    monkeypatch.setitem(sys.modules, "textual", None)  # makes "import textual" fail
    monkeypatch.setitem(sys.modules, "textual.app", None)
    with pytest.raises(SystemExit, match=r"geolab\[tui\]"):
        geolab.tui.main()
