def main() -> None:
    try:
        from .app import GeolabApp
    except ModuleNotFoundError as exc:
        if exc.name is None or exc.name.partition(".")[0] != "textual":
            raise
        raise SystemExit("The geolab TUI needs Textual. Install it with: pip install 'geolab[tui]'")
    GeolabApp().run()
