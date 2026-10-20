---
title: Ball semantics, nonfinite values and domains
summary: How the dfloat balls represent the whole real line, which status each function returns at the boundary of its domain, and how predicates and conversions treat uncertain cases.
order: 2
---

The ball formats `dNb` and `dNcb` implement the real and complex numbers
rigorously in the sense of interval arithmetic: every operation returns
a ball containing the image of the input balls under the exact
operation. Since the representation is built on doubles, which have
infinities and NaN, the module fixes how such values are interpreted and
what happens when an input ball meets the boundary of the domain of a
function. These rules are the same for the `dNb` functions and for the
generic rings created with `gr_ctx_init_dfloat`, and they follow those of
`arb_t` and the `gr` interface [@Joh2025gr].

## The whole real line {#whole-line}

A ball with a nonfinite radius or a nonfinite midpoint component
($\pm\infty$ or NaN) represents the whole real line $W$: `[0 +/- inf]`,
`[5 +/- inf]`, `[inf +/- 0]`, `[nan +/- 1]` and a midpoint `(1, nan)`
with radius 0 all denote $W$. Every function accepts any of these forms.
The functions return $W$ only in the normal form $[0 \pm \infty]$ (all
midpoint components $+0$, radius $+\infty$); any other result has finite
midpoint components and a finite nonnegative radius. The exceptions are
`dNb_set` and the data movement methods of the generic rings (`swap`,
`vec_set`, `vec_gather`, `vec_scatter`), which copy their input as it
is. $W$ results from overflow, from nonfinite input, and from functions
that are unbounded on the input ball. NaN never occurs as a value
otherwise, and every form of $W$ prints as `[+/- inf]`.

## Functions defined everywhere {#total}

Operations defined on all of $\mathbb{R}$ always succeed: `add`, `sub`, `mul`,
`sqr`, `neg`, `abs`, `mul_2exp_si`, `exp`, `expm1`, `sin`, `cos`,
`sin_cos`, `atan`, `sinh`, `cosh`, `tanh`, `atan2`, `sgn`, the dot
products, integer powers with a nonnegative exponent, and division by a
ball not containing zero. On $W$ they give $W$, or a ball containing the
range where the function is bounded: $[0 \pm 1]$ for `sin`, `cos`,
`tanh` and `sgn`, $[0 \pm \pi/2]$ for `atan`, and $[0 \pm \pi]$ for
`atan2` (which, like `arb_atan2`, has $\operatorname{atan2}(0, 0) = 0$).
An exact zero times anything, including $W$, is exactly zero, and
$x^0 = 1$ exactly for every $x$.

## Partial functions {#domains}

An operation defined on a subset $D$ of $\mathbb{R}$ or $\mathbb{R}^2$ returns
`GR_SUCCESS` only if the input balls are contained in $D$, `GR_DOMAIN`
only if they are disjoint from $D$, and `GR_UNABLE` otherwise; the output
is then unspecified. Since $W$ contains points on both sides of every
boundary, an operation with a nonfinite argument returns `GR_UNABLE`
unless the other argument alone decides. The domains are:

| Function | Domain $D$ |
|---|---|
| `div`, `inv` | nonzero divisor |
| `sqrt` | $x \ge 0$ |
| `rsqrt`, `log` | $x > 0$ |
| `log1p` | $x > -1$ |
| `asin`, `acos` | $\lvert x \rvert \le 1$ |
| `tan` | $x \notin \pi/2 + \pi\Z$ (never `GR_DOMAIN`, since the poles are irrational) |
| `pow` | $x > 0$; or $x = 0$ and $y \ge 0$; or $x < 0$ and $y \in \Z$ |

For example:

| Operation | Status and result |
|---|---|
| $[0 \pm \infty] / 3$ | `GR_SUCCESS`, $[0 \pm \infty]$ |
| $1/0$, $[0 \pm \infty]/0$ | `GR_DOMAIN` |
| $[0 \pm \infty] / [3 \pm 2]$ | `GR_SUCCESS`, $[0 \pm \infty]$ |
| $1/[1 \pm 2]$, $1/[0 \pm \infty]$ | `GR_UNABLE` |
| $\exp([0 \pm \infty])$ | `GR_SUCCESS`, $[0 \pm \infty]$ |
| $\sin([\infty \pm 0])$ | `GR_SUCCESS`, $[0 \pm 1]$ |
| $\sqrt{-1}$ | `GR_DOMAIN` |
| $\sqrt{[0 \pm 1]}$ | `GR_UNABLE` |
| $\log 0$, $\log [-2 \pm 1]$ | `GR_DOMAIN` |
| $\log [0 \pm \infty]$ | `GR_UNABLE` |
| $\operatorname{asin}(1)$ | `GR_SUCCESS`, $\pi/2$ |
| $\operatorname{asin}([2 \pm 0.5])$ | `GR_DOMAIN` |
| $0^{-1/2}$, $(-2)^{1/2}$ | `GR_DOMAIN` |
| $(-2)^{[0 \pm \infty]}$ | `GR_UNABLE` |
| $[0 \pm \infty]^3$ | `GR_SUCCESS`, $[0 \pm \infty]$ |

## Predicates and conversions {#predicates}

The predicates of the generic rings (`gr_is_zero`, `gr_equal`,
`gr_cmp`, ...) answer `T_TRUE` or `T_FALSE` (or `GR_SUCCESS`) only when
every point of the balls decides the question, and `T_UNKNOWN`
(`GR_UNABLE`) otherwise, in particular for $W$; they may answer
`T_UNKNOWN` in borderline cases. The `dNb_` predicates return 0 or 1
with the same meaning: see `dNb_contains` for which answers are
certain.

For conversions:

- `gr_set_d` of $\pm\infty$ or NaN returns `GR_DOMAIN`, since these are
  not real numbers (`dNb_set_d` gives $W$);
- `gr_set_other` from a nonfinite `arb_t` gives $W$;
- `gr_get_fmpz` (and `gr_get_si`, `gr_get_ui`) succeeds only for an exact
  integer, and returns `GR_DOMAIN` for a ball that contains no integer
  (or only integers out of range) and `GR_UNABLE` otherwise;
- `gr_get_d` and `dNb_get_d` return 0 for $W$;
- `gr_get_interval_mid_rad` returns `GR_UNABLE` for $W$, whose radius
  is not a real number.

Between `dfloat` rings, `gr_set_other` and `_gr_vec_set_other` copy
when the target has at least as many components, truncate with the tail
added to the radius when the source context is strong (canonical), and
renormalise otherwise.

## Plain formats {#plain}

The plain formats `dN` behave like `double`: they contain infinities and
NaN, overflow gives infinity, and invalid operations give NaN. In the
generic rings they are approximate floating-point rings, and their
operations return `GR_DOMAIN` only for an exact singularity (a zero
divisor, $1/\sqrt{0}$, $0^y$ with $\operatorname{Re}(y) \le 0$ and
$y \ne 0$); $\log 0$ is $-\infty$.

## Complex balls {#complex}

Each part of a complex ball follows the conventions of the real balls,
so that $W + 2i$, $3 + Wi$ and $W + Wi$ are the sets $\mathbb{R} + 2i$,
$3 + \mathbb{R} i$ and $\mathbb{C}$, and results use the normal form $[0 \pm \infty]$
for such parts. Operations defined on all of $\mathbb{C}$ always succeed: `add`,
`sub`, `mul`, `sqr`, `neg`, `conj`, `mul_2exp_si`, `exp`, `sin`,
`cos`, `sin_cos`, `sinh`, `cosh`, `sqrt`, `abs`, `arg`, `sgn`, the dot
products and integer powers with a nonnegative exponent. The others
follow the status rules above, with these domains:

| Function | Domain |
|---|---|
| `div`, `inv` | nonzero divisor |
| `rsqrt`, `log` | $z \ne 0$ |
| `tan` | $z \notin \pi/2 + \pi\Z$ (never `GR_DOMAIN`) |
| `tanh` | $z \notin i(\pi/2 + \pi\Z)$ (never `GR_DOMAIN`) |
| `pow` | $x \ne 0$; or $x = 0$ and $\operatorname{Re}(y) > 0$; or $y \in \{0, 1, 2, \ldots\}$ |

The functions with branch cuts (`log`, `sqrt`, `rsqrt`, `pow`, `arg`)
are the principal branches, with the cut $(-\infty, 0]$ continuous from
above ($\log(-1) = \pi i$, $\sqrt{-4} = 2i$), and a ball crossing the
cut gets an enclosure of the values on both sides. For example:

| Operation | Status and result |
|---|---|
| $W/3$ | `GR_SUCCESS`, $W$ |
| $1/0$, $W/0$ | `GR_DOMAIN` |
| $1/([0 \pm 1] + [0 \pm 1]i)$ | `GR_UNABLE` |
| $1/(W + 2i)$ | `GR_SUCCESS`, $[0 \pm 1/2] + [0 \pm 1/2]i$ |
| $\exp(W)$ | `GR_SUCCESS`, $W + 0i$ |
| $\exp(Wi)$ | `GR_SUCCESS`, $[0 \pm 1] + [0 \pm 1]i$ |
| $\sqrt{-4}$ | `GR_SUCCESS`, $2i$ |
| $\sqrt{[0 \pm 1] + [0 \pm 1]i}$ | `GR_SUCCESS` (an enclosure) |
| $\log 0$ | `GR_DOMAIN` |
| $\log(-1 + [0 \pm 1/4]i)$ | `GR_SUCCESS`, imaginary part $[0 \pm \pi]$ |
| $\log(W + 2i)$ | `GR_SUCCESS`, real part $W$ |
| $\tan(\pi/2 + 10^{-3}i)$ | `GR_SUCCESS` |
| $\tan([\pi/2 \pm 0.01])$ | `GR_UNABLE` |
| $0^{2+i}$, $0^0$ | `GR_SUCCESS`, 0 and 1 |
| $0^{-1}$, $0^i$ | `GR_DOMAIN` |
| $[0 \pm 1/2]^{1/2}$ | `GR_SUCCESS` (an enclosure) |

With the crude enclosures used where the precise formulas fail (wide
rectangles, parts equal to $W$), a complex ball function returns
`GR_UNABLE` for an input inside its domain only rarely: for balls whose
lower bound on $|z|$ is below $2^{-1000}$, or tangents at the edge of
what the double range can decide.

Conversions of the complex rings: from the real types (`arb_t`, `fmpz`,
`fmpq`, the real `dfloat` rings) through the real ring of the parts,
with a zero imaginary part; from `acb_t` part by part (a nonfinite part
gives $W$); between complex formats part by part; from any other ring
through `acb_t`. To a real type (`gr_get_fmpz`, `gr_get_d`,
`gr_set_other` into a real ring or `arb_t`), the imaginary part must be
exactly zero: the result is `GR_DOMAIN` if it is certainly not zero and
`GR_UNABLE` if a ball merely contains zero.
