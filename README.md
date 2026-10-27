# NE 806 Lecture Notes

LaTeX sources for the NE 806 (Neutronics) lecture notes at Kansas State
University. The notes grew out of MIT 22.106 (Neutron Interactions and
Applications); the `f26` branch carries the Fall 2026 revision.

Open work is tracked in [TODO.md](TODO.md).

## Building

```sh
latexmk                  # draft build -> main.pdf
latexmk -r student.rc    # student build -> build/student/main.pdf
latexmk -c               # remove auxiliary files
```

The draft build shows the transcription markup described below; the
student build hides it. Both need a TeX Live with `memoir`; the
vendored `algorithm2e.sty` covers installs that lack it.

## Layout

| Path | Contents |
|---|---|
| `main.tex` | document skeleton: parts and `\include` list |
| `preamble.tex` | packages and macros |
| `format106.sty` | chapter style, title page, TOC formatting |
| `introduction.tex` | front-matter introduction |
| `nuclear_data/` | Part I: nuclear data and reaction theory |
| `deterministic/` | Part II: deterministic transport methods |
| `monte_carlo/` | Part III: Monte Carlo methods |
| `topics/` | Part IV: adjoint, variational, geometry, criticality safety |
| `appendices/` | nomenclature and other appendices |
| `images/<part>/` | figures, grouped like the chapters |
| `code/<part>/` | listings included with `\lstinputlisting` |
| `biblio106.bib` | bibliography |

## Conventions

- One lecture per file, one `\chapter` per lecture. Label each chapter
  `\label{lec:<file name>}` so `\LECTURE{<file name>}` can refer to it.
- Keep images and code under the folder for the part that uses them,
  and refer to them by full path (e.g. `images/deterministic/sndisc`).
- When a new symbol enters the notes, add it to
  `appendices/nomenclature.tex`.
- `images/topics/shortest_path.pdf` is generated: run
  `pdflatex -shell-escape shortest_path.tex` in `images/topics/`. This
  needs `gnuplot` on the path to draw the curves, and uses the
  `gnuplot-lua-tikz.sty` next to it.

## Transcription markup

The nuclear data lectures were transcribed from handwritten notes and
carry markup tying them to the source scans:

| Macro | Purpose |
|---|---|
| `\sourcepage{n}` | marks where page *n* of the scan begins |
| `\editorial{...}` | note on a transcription choice |
| `\uncertain{...}` | flags a reading that needs checking |
| `\marginnoteinline{...}` | a margin note from the scan |

All four are defined in `preamble.tex` and print nothing in the student
build. Remove them as each chapter is revised.
