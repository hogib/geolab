"""fzf-style fuzzy matching: the query's characters must appear in order, and
matches score higher when they are consecutive, start words, or come early."""

from dataclasses import dataclass

_SCORE_MATCH = 16
_BONUS_CONSECUTIVE = 12
_BONUS_WORD_START = 10
_BONUS_FIRST_CHAR = 8
_PENALTY_GAP_START = 3
_PENALTY_GAP = 1


@dataclass(frozen=True)
class Match:
    index: int
    """Position of the candidate in the list that was searched."""
    score: int
    positions: tuple[int, ...]
    """Indices of the matched characters in the candidate."""


def _word_start(text: str, i: int) -> bool:
    if i == 0:
        return True
    prev, cur = text[i - 1], text[i]
    return not prev.isalnum() or (prev.islower() and cur.isupper())


def fuzzy_match(query: str, text: str) -> tuple[int, tuple[int, ...]] | None:
    """Best (score, matched positions) of ``query`` in ``text``, or None if it doesn't match.

    Case-insensitive; spaces in the query are ignored. Uses dynamic programming
    over all alignments, so the best-scoring one is found, not just the first.
    """
    q = [c for c in query.lower() if not c.isspace()]
    if not q:
        return 0, ()
    t = text.lower()
    n, m = len(q), len(t)
    neg = float("-inf")
    # best[i][j]: best score with q[:i+1] matched and q[i] matched at t[j]
    best = [[neg] * m for _ in range(n)]
    back: list[list[int]] = [[-1] * m for _ in range(n)]
    for i in range(n):
        running, running_at = neg, -1  # best of best[i-1][k] - gap penalty for k < j-1
        for j in range(m):
            if i > 0 and j >= 2:
                cand = best[i - 1][j - 2]
                # Extend the running maximum by one more gap character.
                running = running - _PENALTY_GAP if running != neg else neg
                if cand != neg and cand - _PENALTY_GAP_START > running:
                    running, running_at = cand - _PENALTY_GAP_START, j - 2
            if t[j] != q[i]:
                continue
            bonus = _SCORE_MATCH
            if _word_start(text, j):
                bonus += _BONUS_WORD_START
            if i == 0:
                score = bonus + (_BONUS_FIRST_CHAR if j == 0 else -min(j, 6))
                best[i][j], back[i][j] = score, -1
                continue
            options = []
            if j >= 1 and best[i - 1][j - 1] != neg:
                options.append((best[i - 1][j - 1] + _BONUS_CONSECUTIVE, j - 1))
            if running != neg:
                options.append((running, running_at))
            if options:
                prev_score, prev_j = max(options)
                best[i][j], back[i][j] = prev_score + bonus, prev_j
    last = best[n - 1]
    end = max(range(m), key=lambda j: last[j], default=-1)
    if end < 0 or last[end] == neg:
        return None
    positions = [end]
    for i in range(n - 1, 0, -1):
        positions.append(back[i][positions[-1]])
    return int(last[end]), tuple(reversed(positions))


def fuzzy_filter(query: str, candidates: list[str]) -> list[Match]:
    """Candidates matching ``query``, best first (ties keep the original order)."""
    matches = []
    for index, text in enumerate(candidates):
        result = fuzzy_match(query, text)
        if result is not None:
            matches.append(Match(index, result[0], result[1]))
    return sorted(matches, key=lambda m: -m.score)
