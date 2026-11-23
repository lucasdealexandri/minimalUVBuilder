# Reserved-symbol lists

Language-level names that generated FeynRules/UFO symbols must not reuse, read by
`tools/fr_names.language_lists`.

- `mathematica_<version>.txt` — ``Names["System`*"]``
- `feynrules_<version>.txt` — ``Names["FeynRules`*"]``

One symbol per line; lines starting with `#` are headers (Mathematica/FeynRules version, date).
All files of a family are **unioned**, so lists from different versions can coexist (reserving a
few extra names only costs a rename). If a family has no file, a short hand-picked fallback list in
`tools/fr_names.py` is used instead.

Regenerate (from the repository root):

    wolframscript -file scripts/dump_reserved_symbols.wl /path/to/feynrules-current

or open `scripts/dump_reserved_symbols.wl` in a notebook, set `frPath`/`outDir`, and evaluate.

Python keywords and builtins are taken from the running interpreter, so they need no file.
