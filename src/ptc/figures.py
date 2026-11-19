"""README figure: does each reference value follow from the study it cites?

Writes reports/figures/refvalue_consistency_{light,dark}.svg (hand-built SVG, so the
text stays text and both themes are exact). Encoding: a one-hue blue ramp for the
three "explained" outcomes (darkest = exact), orange for "not explained", neutral for
"not checkable". 100 % stacked bars, 2 px surface gaps, 4 px rounded data end.
"""
from __future__ import annotations

from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd

FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"

CATS = [
    ("exact", "Exact"),
    ("rounded", "After rounding"),
    ("absorption_not_recorded", "Oral absorption not recorded"),
    ("not_explained", "Not explained"),
    ("not_checkable", "Not checkable"),
]

THEMES = {
    "light": {"surface": "#fcfcfb", "primary": "#0b0b0b", "secondary": "#52514e",
              "muted": "#898781", "fills": ["#104281", "#2a78d6", "#86b6ef", "#eb6834", "#c3c2b7"]},
    "dark": {"surface": "#1a1a19", "primary": "#f0efec", "secondary": "#c3c2b7",
             "muted": "#898781", "fills": ["#86b6ef", "#3987e5", "#184f95", "#d95926", "#383835"]},
}


def _ink_on(fill: str) -> str:
    r, g, b = (int(fill[i:i + 2], 16) / 255 for i in (1, 3, 5))
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in (r, g, b)]
    lum = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    return "#0b0b0b" if lum > 0.3 else "#ffffff"


def _rows(cons: pd.DataFrame) -> list[tuple[str, int, list[int]]]:
    status = cons["status"].where(~cons["status"].str.startswith("not_checkable"), "not_checkable")
    rows = []
    for label, mask in [("All four", slice(None)), ("ADI", cons["descriptor"] == "ADI"),
                        ("ARfD", cons["descriptor"] == "ARfD"), ("AOEL", cons["descriptor"] == "AOEL"),
                        ("AAOEL", cons["descriptor"] == "AAOEL")]:
        s = status[mask]
        rows.append((label, int(len(s)), [int((s == k).sum()) for k, _ in CATS]))
    return rows


def _segment(x: float, y: float, w: float, h: float, fill: str, round_end: bool) -> str:
    if not round_end or w <= 4:
        return f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{h}" fill="{fill}"/>'
    r = 4
    return (f'<path d="M{x:.1f},{y} H{x + w - r:.1f} Q{x + w:.1f},{y} {x + w:.1f},{y + r} '
            f'V{y + h - r} Q{x + w:.1f},{y + h} {x + w - r:.1f},{y + h} H{x:.1f} Z" fill="{fill}"/>')


def svg(cons: pd.DataFrame, theme: str) -> str:
    t = THEMES[theme]
    rows = _rows(cons)
    W, left, right = 820, 92, 70
    bar_w, bar_h, gap = W - left - right, 22, 2
    top = 112
    row_step = 44
    H = top + row_step * len(rows) + 30
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
           f'font-family="{FONT}" role="img" aria-labelledby="t d">',
           '<title id="t">Do EFSA reference values follow from the study they cite?</title>',
           '<desc id="d">' + escape("; ".join(
               f"{lab} (n={n}): " + ", ".join(f"{name} {c}" for (_, name), c in zip(CATS, counts, strict=True))
               for lab, n, counts in rows)) + "</desc>",
           f'<rect width="{W}" height="{H}" rx="8" fill="{t["surface"]}"/>',
           f'<text x="24" y="36" font-size="17" font-weight="600" fill="{t["primary"]}">'
           "Do EFSA reference values follow from the study they cite?</text>",
           f'<text x="24" y="58" font-size="13" fill="{t["secondary"]}">'
           "Each value recomputed as cited study dose ÷ uncertainty factor · "
           "EFSA OpenFoodTox 3.0, pesticide dossiers</text>"]
    # legend
    lx = 24
    for (_, name), fill in zip(CATS, t["fills"], strict=True):
        out.append(f'<rect x="{lx}" y="76" width="12" height="12" rx="2" fill="{fill}"/>')
        out.append(f'<text x="{lx + 18}" y="86.5" font-size="12.5" fill="{t["secondary"]}">{escape(name)}</text>')
        lx += 18 + len(name) * 6.4 + 26
    for i, (label, n, counts) in enumerate(rows):
        y = top + i * row_step + (0 if i == 0 else 8)
        weight = "600" if i == 0 else "400"
        out.append(f'<text x="{left - 12}" y="{y + 15.5}" font-size="13" font-weight="{weight}" '
                   f'text-anchor="end" fill="{t["primary"]}">{label}</text>')
        nonzero = [j for j, c in enumerate(counts) if c]
        usable = bar_w - gap * (len(nonzero) - 1)
        x = left
        for j, c in enumerate(counts):
            if not c:
                continue
            w = usable * c / n
            fill = t["fills"][j]
            out.append(_segment(x, y, w, bar_h, fill, round_end=(j == nonzero[-1])))
            pct = f"{c / n:.0%}"
            if w >= len(pct) * 7.5 + 10:
                out.append(f'<text x="{x + w / 2:.1f}" y="{y + 15.5}" font-size="12" font-weight="600" '
                           f'text-anchor="middle" fill="{_ink_on(fill)}">{pct}</text>')
            x += w + gap
        out.append(f'<text x="{left + bar_w + 10}" y="{y + 15.5}" font-size="12.5" '
                   f'fill="{t["muted"]}">n = {n:,}</text>')
    out.append(f'<text x="24" y="{H - 12}" font-size="11.5" fill="{t["muted"]}">'
               "Not explained: none of the cited study's doses / UF give the value (a query, not a verdict). "
               "Source: EFSA OpenFoodTox 3.0 (CC BY-ND 4.0).</text>")
    out.append("</svg>")
    return "\n".join(out)


def write(cons: pd.DataFrame, reports: Path) -> list[Path]:
    d = reports / "figures"
    d.mkdir(parents=True, exist_ok=True)
    paths = []
    for theme in THEMES:
        p = d / f"refvalue_consistency_{theme}.svg"
        p.write_text(svg(cons, theme), encoding="utf-8")
        paths.append(p)
    return paths
