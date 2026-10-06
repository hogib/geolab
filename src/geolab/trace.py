from dataclasses import dataclass, field


@dataclass(frozen=True)
class Step:
    name: str
    value: float
    formula: str = ""


@dataclass(frozen=True)
class Section:
    title: str


@dataclass
class Trace:
    """Records the intermediate values of a calculation ("show working").

    Calculations accept ``trace: Trace | None``; check it with ``is not None``.
    """

    entries: list[Step | Section] = field(default_factory=list)

    def step(self, name: str, value: float, formula: str = "") -> float:
        self.entries.append(Step(name, value, formula))
        return value

    def section(self, title: str) -> None:
        self.entries.append(Section(title))

    @property
    def steps(self) -> list[Step]:
        return [e for e in self.entries if isinstance(e, Step)]

    def __getitem__(self, name: str) -> float:
        """Value of the most recent step called ``name``."""
        for entry in reversed(self.entries):
            if isinstance(entry, Step) and entry.name == name:
                return entry.value
        raise KeyError(name)
