latex_test_cases = [
    ("Mul(2, Pi, NumberI)", r"2 \pi i"),
    ("Mul(-2, Pi, NumberI)", r"-2 \pi i"),
    ("Mul(-1, -2, -3)", r"-1 \cdot \left(-2\right) \cdot \left(-3\right)"),
    ("Div(-1, 3)", r"-\frac{1}{3}"),
    ("Pow(2, n)", r"{2}^{n}"),
    ("Pow(-1, n)", r"{\left(-1\right)}^{n}"),
    ("Mul(3, Pow(2, n))", r"3 \cdot {2}^{n}"),
    ("Mul(3, Pow(-1, n))", r"3 \cdot {\left(-1\right)}^{n}"),
    ("Mul(-3, Pow(-1, n))", r"-3 \cdot {\left(-1\right)}^{n}"),
    ("Pow(10, Pow(10, -10))", r"{10}^{{10}^{-10}}"),
    ("Pow(2, Div(-1, 3))", r"{2}^{-1 / 3}"),
    ("Pow(2, Div(-1, Mul(3, n)))", r"{2}^{-1 / \left(3 n\right)}"),
    ("Equal(Add(Pow(Sin(x), 2), Pow(Cos(x), 2)), 1)", r"\sin^{2}\!\left(x\right) + \cos^{2}\!\left(x\right) = 1"),
    ("Sum(f(n) + g(n), For(n, a, b))", r"\sum_{n=a}^{b} \left(f(n) + g(n)\right)"),
    ("Sum(f(n), For(n, S), NotEqual(n, 0))", r"\sum_{\textstyle{n  \in S \atop n \ne 0}} f(n)"),
    ("Sum(f(n), For(n, S))", r"\sum_{n  \in S} f(n)"),
    ("Sum(f(n), For(n, a, b), NotEqual(n, 0))", r"\sum_{\textstyle{n=a \atop n \ne 0}}^{b} f(n)"),
    ("Sum(f(n), For(n, a, b))", r"\sum_{n=a}^{b} f(n)"),
    ("Product(f(n) + g(n), For(n, a, b))", r"\prod_{n=a}^{b} \left(f(n) + g(n)\right)"),
    ("Product(f(n), For(n, S), NotEqual(n, 0))", r"\prod_{\textstyle{n  \in S \atop n \ne 0}} f(n)"),
    ("Product(f(n), For(n, S))", r"\prod_{n  \in S} f(n)"),
    ("Product(f(n), For(n, a, b), NotEqual(n, 0))", r"\prod_{\textstyle{n=a \atop n \ne 0}}^{b} f(n)"),
    ("Product(f(n), For(n, a, b))", r"\prod_{n=a}^{b} f(n)"),
    ("Integral(f(x), For(x, -Infinity, Infinity))", r"\int_{-\infty}^{\infty} f(x) \, dx"),
    ("Integral(f(x), For(x, RR))", r"\int_{x \in \mathbb{R}} f(x) \, dx"),
    ("Integral(f(x) + g(x) / h(x), For(x, a, b))", r"\int_{a}^{b} \left(f(x) + \frac{g(x)}{h(x)}\right) \, dx"),
    ("(Hypergeometric2F1Regularized(Div(-1,4),Div(1,4),1/2, (x-1)/2)**2)", r"{\left(\,{}_2{\textbf F}_1\!\left(-\frac{1}{4}, \frac{1}{4}, \frac{1}{2}, \frac{x - 1}{2}\right)\right)}^{2}"),
]

def test_latex(fexpr):
    print("test latex!")
    fexpr.inject(vars=True)
    for formula, expected in latex_test_cases:
        expr = eval(formula)
        latex = expr.latex()
        if latex != expected:
            raise AssertionError("%s:  got '%s', expected '%s'" % (formula, latex, expected))

def latex_report(fexpr):
    fexpr.inject(vars=True)
    formulas = [eval(formula) for formula, expected in latex_test_cases]
    from os.path import expanduser
    from time import clock
    fp = open(expanduser("~/Desktop/latex_report.html"), "w")
    fp.write("""
<!DOCTYPE html>
<html>
<head>
<title>Fredrik Johansson's website</title>
<meta http-equiv="Content-Type" content="text/html;charset=utf-8" >
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
tt { padding: 0.2em; background-color: #f8f8f8; border:1px solid #eee; }
table { border-collapse:collapse; margin: 1em; }
table, th, td { border: 1px solid #aaa; }
th, td { padding:0.3em; padding-top: 0; margin-top: 0; }
table { width: 95%; }
.katex { font-size: 1em !important; } 
</style>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.css" integrity="sha384-AfEj0r4/OFrOo5t7NnNe46zW/tFgW6x/bCJG8FqQCEo3+Aro6EYUG4+cU+KJWu/X" crossorigin="anonymous">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.js" integrity="sha384-g7c+Jr9ZivxKLnZTDUhnkOnsh30B4H0rpLUpJ4jAIKs4fnJI+sEnkvrMWph2EDg4" crossorigin="anonymous"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/contrib/auto-render.min.js" integrity="sha384-mll67QQFJfxn0IYznZYonOWZ644AWYC+Pt2cHqMaRhXVrursRwvLnLaebdGIlYNa" crossorigin="anonymous"
    onload="renderMathInElement(document.body);"></script>
<script>
  document.addEventListener("DOMContentLoaded", function() {
      renderMathInElement(document.body, {
          delimiters: [
            {left: "$$", right: "$$", display: true},
            {left: "$", right: "$", display: false}
          ]
      });
  });
</script>
</head>
<body>
""")
    t1 = clock()
    output = [formula.latex() for formula in formulas]
    t2 = clock()
    fp.write("""<h1>fexpr to LaTeX test sheet</h1>""")
    fp.write("""<p>Converted %i formulas to LaTeX in %f seconds.</p>""" % (len(formulas), (t2-t1)))
    fp.write("""<table>""")
    fp.write("""<tr><th>fexpr (Python input)</th> <th>Generated LaTeX</th> <th>KaTeX display</th>""")
    for formula, latex in zip(formulas, output):
        fp.write("""<tr>""")
        fp.write("""<td><tt>%s</tt></td>""" % formula)
        fp.write("""<td><tt>%s</tt></td>""" % latex)
        fp.write("""<td>$$%s$$</td>""" % latex)
        fp.write("""</tr>""")
    fp.write("""</table>""")

    #for formula, latex in zip(formulas, output):
    #    fp.write("""<table><tr><td><b>fexpr:</b> <tt>%s</tt></td></tr> <tr><td><b>LaTeX:</b> <tt>%s</tt></td></tr> <tr><td>$$%s$$</td></tr></table>""" % (formula, latex, latex))

    fp.write("""</body></html>""")
    fp.close()
