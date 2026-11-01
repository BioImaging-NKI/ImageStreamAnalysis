"""The context help of the import page and the dashboard, from the user manual: the guide, analysis, reference and help
pages of manual/ converted to HTML, one section per heading, into src/imagestream_analysis/templates/help.json.

    python tools/manual/make_help.py

The pages show it in a side panel (templates/help.js); their ? buttons point at '<page>#<heading anchor>'.
tests/test_help.py fails when help.json no longer matches the manual, or a ? button points at a heading that is gone.
Screenshots are left out (the help is read next to the real thing); formulas become plain HTML (no MathJax: the
dashboard also works offline); links to pages that are not in the help open the online manual.
"""

import json
import posixpath
import re
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parents[2]
MANUAL = ROOT / "manual"
OUT = ROOT / "src" / "imagestream_analysis" / "templates" / "help.json"
SITE = "https://bioimaging-nki.github.io/ImageStreamAnalysis/"

# the pages in the help, in the order of the manual's contents (mkdocs.yml)
PAGES = {
    "User guide": ["guide/import-page", "guide/dashboard", "guide/overview", "guide/output-folder"],
    "How the analysis works": ["analysis/input", "analysis/segmentation", "analysis/classification",
                               "analysis/interactions", "analysis/membrane",
                               "analysis/comparing-conditions"],
    "Reference": ["reference/settings", "reference/command-line", "reference/fiji"],
    "Help": ["help/faq", "help/glossary"],
}
IN_HELP = {p for pages in PAGES.values() for p in pages}

EXTENSIONS = ["abbr", "admonition", "attr_list", "def_list", "footnotes", "md_in_html", "tables", "toc",
              "pymdownx.details", "pymdownx.keys", "pymdownx.superfences"]

TEX = [(r"\\text\{([^{}]*)\}", r"\1"), (r"\\overline\{([^{}]*)\}", r'<span style="text-decoration:overline">\1</span>'),
       (r"\\frac\{([^{}]*)\}\{([^{}]*)\}", r"(\1) / (\2)"), (r"_\{([^{}]*)\}", r"<sub>\1</sub>"),
       (r"\^\{([^{}]*)\}", r"<sup>\1</sup>"), (r"_([A-Za-z0-9])", r"<sub>\1</sub>"), (r"\\[,;]", " "),
       (r"\\(ln|min|max|exp)\b", r"\1"), (r"\\cdot", "·"), (r"\\times", "×"), (r"\\", "")]


def _tex(formula: str) -> str:
    """A simple TeX formula of the manual as HTML (fractions, sub- and superscripts, \\text, \\overline)."""
    s = formula.strip()
    for _ in range(3):  # nested braces: inner first
        for pattern, repl in TEX[:-1]:
            s = re.sub(pattern, repl, s)
    return re.sub(r"\\", "", s).replace("{", "").replace("}", "")


def _math(text: str) -> str:
    """$$…$$ as a centred line and $…$ inline, as HTML (outside code)."""
    text = re.sub(r"\$\$(.+?)\$\$", lambda m: f'<span class="math display">{_tex(m.group(1))}</span>', text,
                  flags=re.DOTALL)
    return re.sub(r"(?<![\\`$])\$([^$\n`]+?)\$", lambda m: f'<span class="math">{_tex(m.group(1))}</span>', text)


def _source(page: str) -> str:
    text = (MANUAL / f"{page}.md").read_text(encoding="utf-8")
    text = re.sub(r"\A---\n.*?\n---\n", "", text, flags=re.DOTALL)  # front matter
    text = re.sub(r"<figure markdown>.*?</figure>\s*", "", text, flags=re.DOTALL)  # screenshots
    text = re.sub(r"!\[[^\]]*\]\([^)]*\)(\{[^}]*\})?", "", text)  # other images
    text = re.sub(r"^\[[^\]]*\]\([^)]*\)\{ ?\.md-button[^}]*\}\s*$", "", text, flags=re.MULTILINE)  # 'Try it' buttons
    text = re.sub(r"\s*:(material|octicons|fontawesome)-[a-z0-9-]+:", "", text)  # icons
    text = _math(text)
    return text + "\n\n" + (MANUAL / "includes" / "abbreviations.md").read_text(encoding="utf-8")


def _link(page: str, href: str) -> str:
    """An attribute for a link of the page: data-h for a page in the help, else the online manual."""
    if re.match(r"[a-z]+:", href):
        return f'href="{href}" target="_blank" rel="noopener"'
    path, _, anchor = href.partition("#")
    if not path:
        target = page
    else:
        target = posixpath.normpath(posixpath.join(posixpath.dirname(page), path))
        target = re.sub(r"(/index)?\.md$", "", target)
    if target in IN_HELP:
        return f'data-h="{target}#{anchor}"' if anchor else f'data-h="{target}"'
    url = SITE + ("" if target in ("index", ".") else target + "/") + (f"#{anchor}" if anchor else "")
    return f'href="{url}" target="_blank" rel="noopener"'


def page_html(page: str) -> tuple[str, str]:
    """(title, HTML) of a manual page; every h2/h3 with what follows it is a <section data-a="anchor">."""
    md = markdown.Markdown(extensions=EXTENSIONS, extension_configs={"toc": {"permalink": False}})
    html = md.convert(_source(page))
    m = re.match(r'\s*<h1 id="[^"]*">(.*?)</h1>\s*', html, re.DOTALL)
    title, html = (re.sub(r"<[^>]+>", "", m.group(1)), html[m.end():]) if m else (page, html)
    html = re.sub(r"<!--.*?-->\s*", "", html, flags=re.DOTALL)
    html = re.sub(r'<a href="([^"]*)"', lambda a: f"<a {_link(page, a.group(1))}", html)
    parts = re.split(r'(?=<h[23] id=")', html)
    out = [f'<section data-a="">{parts[0].strip()}</section>'] if parts[0].strip() else []
    for p in parts[1:]:
        anchor = re.match(r'<h[23] id="([^"]*)"', p).group(1)
        out.append(f'<section data-a="{anchor}">{p.strip()}</section>')
    return title, "\n".join(out)


def build() -> dict:
    pages = []
    for group, names in PAGES.items():
        for name in names:
            title, html = page_html(name)
            pages.append({"path": name, "group": group, "title": title, "html": html})
    return {"site": SITE, "pages": pages}


def text() -> str:
    return json.dumps(build(), ensure_ascii=False, indent=0) + "\n"


def main() -> None:
    OUT.write_text(text(), encoding="utf-8", newline="\n")
    print(f"{OUT.relative_to(ROOT)}: {OUT.stat().st_size // 1024} kB, "
          f"{sum(len(p) for p in PAGES.values())} pages")


if __name__ == "__main__":
    main()
