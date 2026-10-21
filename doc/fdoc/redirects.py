"""
Redirect pages for the URLs of the old Sphinx documentation.

The Sphinx documentation (FLINT 3.0 to 3.6) had one page per source file,
<name>.html in a single directory, and anchors "#c.<symbol>" for
functions and types. `fdoc build --redirects` writes a small page for each
old name that sends the browser to the corresponding new page; for module
pages it also maps the old function anchors to the new function pages.
"""

import json
import os

from html import escape as esc

# Old pages that did not document a module, and where their content is now.
SPECIAL = {
    "index": "index.html",
    "genindex": "index.html",
    "search": "index.html",
    "references": "references.html",
    "introduction": "guide/getting-started/introduction.html",
    "building": "guide/getting-started/building.html",
    "memory": "guide/conventions/memory.html",
    "threading": "guide/conventions/threading.html",
    "portability": "guide/conventions/portability.html",
    "contributing": "guide/development/contributing.html",
    "contributors": "guide/development/contributors.html",
    "bug_reporting": "guide/development/bug-reporting.html",
    "python_flint": "guide/interfaces/python.html",
    "history": "history/index.html",
    "overview": "guide/types/ball-arithmetic.html",
    "using": "guide/types/ball-arithmetic.html",
    "issues": "guide/types/ball-arithmetic.html#conventions",
    "introduction_calcium": "guide/types/choosing.html#exact",
    "examples": "tutorials/getting-started/example-programs.html",
    "examples_arb": "tutorials/numerics/index.html",
    "examples_calcium": "tutorials/exact/index.html",
    "formulas": "algorithms/ball-arithmetic/error-propagation.html",
    "constants": "algorithms/special-functions/constants.html",
    "gamma": "algorithms/special-functions/gamma.html",
    "hypergeometric": "algorithms/special-functions/hypergeometric.html",
    "hurwitz": "algorithms/special-functions/zeta.html",
    "polylogarithms": "algorithms/special-functions/zeta.html#polylog",
    "agm": "algorithms/special-functions/elliptic-modular.html#agm",
    "gr_implementing": "gr/implementing.html",
    "index_integers": "index.html",
    "index_integers_mod": "index.html",
    "index_rationals": "index.html",
    "index_arb": "index.html",
    "index_generic": "index.html",
}

# Old pages that documented a module (one per header file).
MODULES = """
acb acb_calc acb_dft acb_dirichlet acb_elliptic acb_hypgeom acb_mat acb_modular
acb_poly acb_theta acf aprcl arb arb_calc arb_fmpz_poly arb_fpwrap arb_hypgeom
arb_mat arb_poly arf arith bernoulli bool_mat ca ca_ext ca_field ca_mat ca_poly
ca_vec calcium d_mat d_vec decimal dfloat dirichlet dlog double_extras
double_interval ecpp fexpr fexpr_builtin fft fft_small flint fmpq fmpq_mat
fmpq_mpoly fmpq_mpoly_factor fmpq_poly fmpq_vec fmpz fmpz_extras fmpz_factor
fmpz_lll fmpz_mat fmpz_mod fmpz_mod_mat fmpz_mod_mpoly fmpz_mod_mpoly_factor
fmpz_mod_mpoly_q fmpz_mod_poly fmpz_mod_poly_factor fmpz_mod_vec fmpz_mpoly
fmpz_mpoly_factor fmpz_mpoly_q fmpz_poly fmpz_poly_factor fmpz_poly_mat
fmpz_poly_q fmpz_vec fmpzi fq fq_default fq_default_mat fq_default_poly
fq_default_poly_factor fq_embed fq_mat fq_nmod fq_nmod_embed fq_nmod_mat
fq_nmod_mpoly fq_nmod_mpoly_factor fq_nmod_poly fq_nmod_poly_factor fq_nmod_vec
fq_poly fq_poly_factor fq_vec fq_zech fq_zech_embed fq_zech_mat fq_zech_poly
fq_zech_poly_factor fq_zech_vec gr gr_dft gr_domains gr_generic gr_mat gr_mpoly
gr_ore_poly gr_poly gr_series gr_special gr_vec long_extras longlong
machine_vectors mag mp_real mpn_extras mpn_mod mpoly nf nf_elem nfloat nmod
nmod_mat nmod_mpoly nmod_mpoly_factor nmod_poly nmod_poly_factor nmod_poly_mat
nmod_vec padic padic_mat padic_poly padic_radix partitions perm profiler qadic
qfb qqbar qsieve radix thread_pool ulong_extras
""".split()

TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Moved: {title}</title>
<link rel="canonical" href="{target}">
<meta name="robots" content="noindex">
<noscript><meta http-equiv="refresh" content="0; url={target}"></noscript>
<script>
(function () {{
  var map = {anchors};
  var h = decodeURIComponent(location.hash.slice(1));
  var t = "{target_js}";
  if (h) {{
    var k = h.indexOf("c.") === 0 ? h.slice(2) : h;
    if (map[k]) t = map[k];
    else if (t.indexOf("#") < 0) t += "#" + h;
  }}
  location.replace(t);
}})();
</script>
</head>
<body>
<p>This page has moved to <a href="{target}">{target}</a>.</p>
</body>
</html>
"""


def write_redirects(site, outdir):
    """Write <outdir>/<old>.html for every page of the old documentation.
    Returns the number of pages written. Existing files (index.html,
    references.html, which keep their names) are left alone."""
    n = 0
    targets = dict(SPECIAL)
    for name in MODULES:
        if name in site.modules:
            targets[name] = site.modules[name].url
        else:
            targets.setdefault(name, "index.html")
    # old anchors: functions and types of a module ("#c.fmpz_add"), and the
    # release anchors of the version history ("#flint-3-1-0")
    by_module = {}
    for sym_name, sym in site.symbols.items():
        by_module.setdefault(sym.block.module, {})[sym_name] = sym.url
    history = {}
    for page in site.order:
        if page.collection == "history":
            for node in page.body:
                if node.kind == "heading" and node["id"]:
                    history[node["id"]] = "%s#%s" % (getattr(page, "canonical_url", page.url).split("#")[0],
                                                     node["id"])
    for old, target in sorted(targets.items()):
        path = os.path.join(outdir, old + ".html")
        if os.path.exists(path):
            continue
        anchors = by_module.get(old, {}) if old in MODULES else (history if old == "history" else {})
        text = TEMPLATE.format(title=esc(old), target=esc(target), target_js=target.replace('"', '\\"'),
                               anchors=json.dumps(anchors, separators=(",", ":")))
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
        n += 1
    return n
