.. _gr-tower:

**gr_tower.h** -- towers of algebraic and transcendental extensions
===============================================================================

(The lazy fields of the second half of this chapter are declared in
**gr_tower_lazy.h**, which includes ``gr_tower.h``.)

A :type:`gr_tower_t` represents a tower of fields

.. math::

    F_0 \subset F_1 \subset \cdots \subset F_n, \quad F_k = F_{k-1}[x] / (m_k)

where each step adjoins an algebraic element `a_k` with monic minimal
polynomial `m_k \in F_{k-1}[x]`, and the base field
`F_0 = K(t_1, \ldots, t_r)` is a rational function field over the
constant field `K` (currently `\mathbb{Q}`) in transcendental generators
`t_j \in \{\pi, \exp(u), \log(u)\}` (where `u` is an element of the
tower defined by the generators preceding `t_j`). Every generator is
embedded in `\mathbb{C}` by a numerical enclosure, which fixes the choice
of conjugate or branch, so that the tower represents a definite subfield
`F_n \subset \mathbb{C}`.

The generators are ordered by their *definition order* (the order in
which they were adjoined, which is compatible with dependencies); the
first `p` generators in this order form a *prefix* of the tower, which
is itself a tower. The algebraic generators have a numbering
`a_1, \ldots, a_n` of their own and the transcendental generators
`t_1, \ldots, t_r` likewise.

Transcendental generators come with a status: `\pi` is transcendental,
while `\exp(u)` and `\log(u)` are transcendental over the field below
only conjecturally (Schanuel's conjecture) in general. An element which
is nonzero as a rational function of such generators is a nonzero
number if the generators are algebraically independent; deciding this
is the subject of Richardson's algorithm, described below.

Elements of `F_n` are represented recursively: an element of `F_k` is a
polynomial over `F_{k-1}` of degree less than `\deg m_k`
(see :func:`gr_ctx_init_gr_poly_quotient`). This representation is
efficient for towers of many small extensions, such as
`\mathbb{Q}(\sqrt{2}, \sqrt{3}, \ldots, \sqrt{47})` of degree `2^{15}`,
where a primitive element would be impractical, as well as for a single
extension of high degree.

Dynamic evaluation
-------------------------------------------------------------------------------

The minimal polynomials `m_k` are not required to be known irreducible
over `F_{k-1}`: a step may be adjoined with a polynomial that merely vanishes
at the given root. The quotient ring `F_{k-1}[x]/(m_k)` then pretends to
be a field (see :func:`gr_ctx_set_is_pretend_field`); it reports itself
as a field via :func:`gr_ctx_is_field` only when `m_k` has been proven
irreducible (and all the steps below it are fields). Whenever an inversion
exposes a factorization of `m_k` (failing with ``GR_UNABLE``), the tower
is *refined* by replacing `m_k` with the factor vanishing at the
enclosure of `a_k`, and the operation is retried. The tower contexts
(lazy and fixed) are genuine fields: the refinement is transparent there.

Since the new modulus divides the old one, refinement never invalidates
existing elements: any polynomial representing an element of the old
quotient ring represents the same number in the refined tower and is
reduced automatically by subsequent operations. Refinement is therefore
transparent to the user, apart from the fact that the value of
:func:`gr_tower_degree` can decrease.

This makes the zero test complete without irreducibility proofs and
without relying on numerical evaluation: an element whose representative
is nonzero is either invertible (hence provably nonzero, since the
evaluation map to `\mathbb{C}` is a ring homomorphism) or exposes a
factorization. Numerical evaluation is used only as a fast path, and
to choose between factors.

Richardson's algorithm
-------------------------------------------------------------------------------

With `\exp` and `\log` generators, an element which is a nonzero rational
function of the generators may still be zero as a number, due to a
relation between the generators. Richardson's algorithm decides this:
if the element is not separated from zero numerically at the current
precision, a `\mathbb{Q}`-linear relation is sought (by LLL) among the
numbers on the *log side* of the tower -- the arguments `a` of the
generators `\exp(a)`, the generators `\log(u)` themselves, and `2 \pi i`
when `\pi` and `i` are present. A candidate relation `\sum c_j l_j = 0`
is verified by a zero test in the tower below its highest generator
(recursively): for a top generator `\exp(a)` the relation itself lives
there, and for a top generator `\log(u)` the equivalent multiplicative
relation `\prod e^{c_j l_j} = 1` is tested, with the branch pinned down by
`|\sum c_j l_j| < 2`. A verified relation is used to eliminate the top
generator, which becomes *algebraic* over the generators preceding it:
of degree 1 (`\log(u) = \ldots`, or `\exp(a) = u` with a unit
coefficient) or a radical `x^c = u` (a dynamic step, since `x^c - u` may
factor; when `u` is rational, the irreducible factor over `\mathbb{Q}`
vanishing at the generator is used, so that e.g. `\exp(i \pi / 3)` gets
the modulus `x^2 - x + 1`). If no relation is found, the precision is
doubled. Under
Schanuel's conjecture a zero element always leads to a relation, so
the algorithm terminates; every answer it gives is proved, the
conjecture being needed only for termination. The precision limit
(the option :macro:`GR_TOWER_OPT_CERTIFY_PREC_LIMIT`) bounds the effort, after which the answer
is unknown.

A relation `c_g u_g + c_j u_j = 0` between the arguments of exactly two
exponential generators with `|c_g| > 1` (as arises with `\exp(x/2)`
adjoined over `\exp(x)`) is not resolved by a radical: the exponential
of the primitive argument `u^* = u_j / (c_g / \gcd(c_g, c_j))` is
adjoined as a new generator, inserted before both, of which both become
powers (linear moduli), so that the degree of the tower does not grow.

*Torsion.* Arguments with rational parts in `\pi i`, `u_j = y_j + r_j \pi
i` (`r_j` the constant term of the coefficient of `\pi` divided by `i`),
lead to relations `\sum c_j u_j + c_k 2 \pi i = 0` which are resolved
modulo the roots of unity: when `\sum c_j r_j + 2 c_k = 0`, the
exponentials become `\exp(y_j) \exp(r_j \pi i)` with a root of unity of
order `\operatorname{lcm}(2 \operatorname{den}(r_j))` (first in the tower,
the root of unity generators already present becoming its powers): for
two exponentials directly through the primitive argument `y^*` as above,
for several by new exponentials `\exp(y_j)`, related by the next search
without roots of unity. Thus `\exp((16 + 30 \pi i)/225)` and `\exp((5 + 42
\pi i)/60)` are expressed through `\exp(1/900)` and `\zeta_{60}` (degree
16), not through a radical of degree 375 over a root of unity of order
125. This is done when the root of unity is present, or when its degree
is smaller than that of the radical it replaces. When the rational part
is hidden (`\pi` in a denominator of the argument, as in the
transformation factors of the theta functions), a relation
`c_g u_g + \sum c_j u_j + k \cdot 2\pi i = 0` in which `c_g` divides the
`c_j` still gives a modulus of degree one when `\exp(-2\pi i k/c_g)` is
`\pm 1` or `\pm i` with `i` present: `\exp(u + \pi i/2) = i \exp(u)`
rather than a root of the reducible `X^2 + \exp(u)^2`.

A generator for a definition (a function and its arguments) which is not
in the tower of the arguments is first looked for in the other towers of
the context (by enclosures, then exact comparison of the arguments), and
takes over the definition id of one found there, so that merges identify
the two (a value cached with an element of an older tower, or adjoined to
a copy of a prefix, would otherwise give a second, unrelated generator
for the same value).

The trigonometric generators enter through their *angles*: `A = 2u`
for `\tan(u)` and `A = 2 \arctan(u)`, whose points on the unit circle
`(\cos A, \sin A) = ((1 - q^2)/(1 + q^2), 2q/(1 + q^2))` (with
`q = \tan(u)`, the generator, or `q = u`) are rational functions of the
generators. When `i` is present, the angles are on the log side as
`i A`, with `e^{i A} = (1 + i q)/(1 - i q)`, so that the relations with
the complex exponentials are found (`e^{2i} = (1 + i \tan 1)/(1 - i
\tan 1)`). Without `i`, a separate *real angular search* runs on the
angles and `\pi` (included numerically when the tower does not contain
it: a relation with it makes `\pi` adjoined): a relation `\sum c_j A_j
= 0` is verified directly below its top generator when that is a
tangent (the angles `2u` live there), and otherwise through the product
of the points, `\prod (\cos A_j + i \sin A_j)^{c_j} = 1` computed with
pairs of real elements, plus `|\sum c_j A_j| < 1`. The elimination is
linear for an arctangent (`\arctan(u) = -\sum_{j \ne g} c_j A_j / (2 c_g)`)
and for a tangent with a unit coefficient (`\tan(u) = \sin A /
(1 + \cos A)` for the point of `A = 2u` given by the others). A tangent
with a larger coefficient is first swapped with another tangent of unit
coefficient when the definition order allows it (so that `\tan(2x)`
becomes a rational function of `\tan(x)`, not `\tan(x)` a root of a
polynomial over `\tan(2x)`); otherwise the lattice spanned by the angles
involved gets a basis (Hermite normal form), new generators `\tan(B_r/2)`
are inserted before the tangents involved, and these become rational
functions of the new ones through the addition formulas, so that the
tower degree does not grow; as a last resort, `\tan(u)` becomes a root
of the real polynomial `C \operatorname{Im}((1 + ix)^{2c}) - S
\operatorname{Re}((1 + ix)^{2c})` (a dynamic step). Thus `4 \arctan(1/5)
- \arctan(1/239) = \pi/4`, `\arctan(\tan 2) = 2 - \pi` and the addition
formulas are decided without complex numbers.

Lambert W values `w = W_k(z)` are on the log side as well, with
`e^w = z / w` in place of the exponential generator, so that their
relations with exponentials and logarithms are found in the same way:
`W(2 \log 2) = \log 2`, `W(3 e^3) = 3`, `e^{W(z)} = z / W(z)`,
`\log W(z) = \log z - W(z)`. When `w` is the top generator of a relation
(with a unit coefficient), the candidate value `w^*` in the tower below
is verified exactly as a solution, `w^* e^{w^*} = z`, and identified with
the branch of the generator by a Krawczyk test: a ball containing the
enclosures of both on which `w e^w - z` has a unique zero. The
two branches `W_0(-\log(2)/2) = -\log 2` and
`W_{-1}(-\log(2)/2) = -2 \log 2` are thus distinguished.

Since eliminated generators are kept (as algebraic generators), elements
in the flat representation remain valid across eliminations, exactly as
across refinements: the reduction modulo the new modulus performs the
substitution. This is why the lazy field can decide equalities in place,
learning relations as it goes. In the nested representation, an
elimination rebuilds the chain of contexts, so :func:`gr_tower_is_zero`
runs the search on a copy of the tower.

Special functions
-------------------------------------------------------------------------------

Values of special functions (the gamma function, the error function, the
Lambert W function, polygamma functions, polylogarithms, the zeta
function, complete elliptic integrals, and constants such as Euler's
`\gamma`) are transcendental generators of their own kinds. Unlike
`\exp` and `\log`, most of them are not covered by Schanuel's conjecture:
their algebraic relations come from functional equations (shifts,
reflection, multiplication theorems), special values and transformations.
The principle is the same as for Richardson's algorithm: every relation
used is an identity, applied exactly, and numerical evaluation only
proves nonvanishing and chooses representatives. Generators of these
kinds have the status ``GR_TOWER_STATUS_CONJECTURAL``: distinct
generators are assumed algebraically independent unless a relation is
known, so that the zero test is sound (an element which is numerically
indistinguishable from zero without an explaining relation gives
``T_UNKNOWN``) but complete only under the corresponding conjectures
(for instance the Rohrlich-Lang conjecture for the gamma function at
rational arguments).

The lazy field applies the identities before creating a generator, so
that values related by them are expressed through the same generator
(*canonical arguments*). The canonical arguments are intrinsic -- they
depend on the value of the argument, not on its representation -- so
that different towers create the same generators, which merges then
identify:

* `\Gamma(z)`, `\psi^{(m)}(z)`: the argument is shifted by an integer
  into the strip `0 \le \operatorname{Re}(z) < 1` (Pochhammer products,
  respectively sums of `(z + k)^{-m-1}`) and reflected (`z \to 1 - z`,
  with `\pi / \sin(\pi z)`, respectively derivatives of
  `\pi \cot(\pi z)`) into `0 < \operatorname{Re}(z) < 1/2`, or
  `\operatorname{Re}(z) \in \{0, 1/2\}` with `\operatorname{Im}(z) > 0`.
  Integers and half-integers give factorials and `\sqrt{\pi}`, the
  digamma function at rationals is given by Gauss's digamma theorem
  (through `\gamma`, logarithms and `\pi`), and `\psi^{(m)}` at `1` and
  `1/2` by zeta values.
* `\Gamma(p/q)` at rationals with denominator `q \le 36`: a normal form
  with respect to all the relations which follow from the reflection
  formula and Gauss's multiplication theorem (the distribution relations,
  which by the Rohrlich-Lang conjecture are all the algebraic relations
  between these values and `\pi`). A reduced row echelon form of the
  relations at level `q` (in the logarithms of the values, of `\pi`, of
  integers and of the numbers `\sin(\pi k/q)`), with the values ordered
  by preference (smallest denominator, then smallest numerator), selects
  `\varphi(q)/2` basis values at level `q`, which contain the basis values
  at every level dividing `q`; every other value is a monomial in the
  basis values times rational powers of `\pi`, of integers and of products
  of sines. For example `\Gamma(1/6) = 2^{-1/3} \sqrt{3/\pi}\, \Gamma(1/3)^2`
  and `\Gamma(2/3) = 2\pi / (\sqrt{3}\, \Gamma(1/3))`. The products of sines
  are written in terms of the root of unity `e^{\pi i/q}` (square roots of
  rationals for `q = 3, 4, 6`), so that values at one level share a
  cyclotomic field rather than introducing independent algebraic
  generators. Values with larger denominators (or when the normal form
  cannot be built) use the reflection formula only and become
  generators, possibly dependent on each other: the zero test then
  computes the distribution relations between the gamma generators at
  rational arguments of each level `N \le 240` exactly (no numerical
  relation search is needed, since they are theorems) and eliminates
  the latest generator involved, which becomes algebraic of degree one
  over the others, `\pi`, roots of integers (of prime power orders, as
  in the lazy layer, or Gauss sums for square roots) and of `\pi`, and
  roots of products of sines (written through `e^{\pi i/N}`). The
  relation is made an integral combination of the distribution relations
  through a Hermite normal form at the level, cached per thread.
* `\psi^{(m)}(p/q)` (`m \ge 1`) and `\zeta(s, p/q) = (-1)^s \psi^{(s-1)}(p/q) / (s-1)!`
  at rationals with denominator `q \le 240`: the additive analogue of
  the gamma normal form. The values `x_k = \zeta(s, k/q)`,
  `0 < k \le q` (`x_q = \zeta(s)`), satisfy the reflection relations
  `x_k + (-1)^s x_{q-k} = (-1)^{s-1} \pi^s P_{s-1}(\cot(\pi k/q)) / (s-1)!`
  (with `P_j(\cot x) = (d/dx)^j \cot x`) and the distribution relations
  `\sum_{j<m} x_{b + jq/m} = m^s x_{bm}` (`m \mid q`), which by a
  conjecture of the Chowla-Milnor type are all the linear relations
  between these values over the algebraic numbers. A reduced row echelon
  form with the same preference order gives a basis of values at each
  level (with `\zeta(s)` for odd `s`, and Catalan's constant in place of
  `\psi'(1/4) = \pi^2 + 8G`); every other value is a rational combination
  of the basis values plus `\pi^s` times an element of a cyclotomic field
  (the cotangents written through roots of unity). For example
  `\psi'(3/4) = \pi^2 - 8G`, `\psi''(1/3) = -4\pi^3/(3\sqrt 3) - 26\zeta(3)`,
  `\psi'(1/6) = 5\psi'(1/3) - 4\pi^2/3`. The row echelon forms are
  cached per thread by weight and level. Values with larger denominators
  are generators (after reflection), and the zero test computes the
  relations between the Hurwitz generators (polygamma values at
  rationals, `\zeta(s)` for odd `s`, Catalan's constant) at each weight
  and level `N \le 240` exactly, through the normal form of the cached
  level; every relation whose pivot is a
  generator makes that generator a rational combination of earlier ones
  plus `\pi^s` times an element of `\mathbb{Q}(\zeta_{\operatorname{lcm}(N,4)})`
  (computed modulo the cyclotomic polynomial, then written through the
  roots of unity of the tower, with no algebraic denominator). On a
  common rational line `a w + b` with irrational `w` (as for
  `\Gamma`), the zero test uses the shift, reflection and multiplication
  relations along the line: they are exact (no `2\pi i` ambiguity, no
  exponential factors), with constants `(a w + c)^{-s}` and the
  reflection terms, whose cotangents are written through
  `\exp(\pi i w)` and roots of unity as
  `i (Y + 1)(1 + Y + \cdots + Y^{r-1}) / (Y^r - 1)` with
  `Y = e^{2\pi i (a w + c)}`, so that
  `\psi'(\sqrt 2) + \psi'(\sqrt 2 + 1/2) = 4 \psi'(2 \sqrt 2)` is decided.
  The digamma function takes part with `s = 1` (`-\psi` in place of
  `\zeta(1, z)`, `P_0(x) = x`), its multiplication theorem
  `\sum_{j<n} \psi(z + j/n) = n \psi(nz) - n \log n` contributing the
  logarithms of the primes: `\psi(2\sqrt 2) = (\psi(\sqrt 2) +
  \psi(\sqrt 2 + 1/2))/2 + \log 2`.
* `\log \Gamma(z)` (``lgamma``, the analytic continuation from the
  positive reals with the cut on the negative reals, as
  :func:`acb_lgamma`) is `\log(\Gamma(z)) + 2\pi i k`, the integer `k`
  determined numerically: the relations of `\log \Gamma` follow from
  those of `\Gamma` and the structure of the logarithms
  (`\log\Gamma(\pi + 1) = \log\Gamma(\pi) + \log \pi`).
* `\Gamma` at irrational arguments on a common rational line
  `a w + b` (`a, b` rational: `\Gamma(\sqrt 2)`, `\Gamma(\sqrt 2 + 1/2)`,
  `\Gamma(2 \sqrt 2)`): the zero test finds the line (an integer relation
  between the arguments, verified exactly) and the shift, reflection and
  multiplication relations along it, and eliminates a generator in terms
  of the others, `\pi`, `\exp(\pi i w)`, `\exp(w \log p)` (from
  `n^{n z}`), roots and linear factors `a w + c`. The result is checked
  numerically, since the logarithmic relations hold modulo `2 \pi i`.
  Since `\Gamma(a w + b)` is, by these relations, a constant times the
  product of the `\Gamma(w + r)` for `r` in `\{(b + i)/a \bmod 1\}`
  (likewise for the Hurwitz zeta function, with a sum), values each of
  whose sets of residues has an element of its own (as two values from
  unrelated computations which happen to lie on a line) have no
  relation, and the linear algebra is skipped.
* `\operatorname{erf}`: odd, so the argument has `\operatorname{Re}(z) > 0`;
  on the imaginary axis `\operatorname{erf}(iy) = i \operatorname{erfi}(y)`
  with a real generator; `\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)`
  and `\operatorname{erfi}(z) = -i \operatorname{erf}(iz)`.
* `\zeta(n)` at integers: Bernoulli numbers (rational multiples of
  `\pi^n` for even `n`); `\zeta(3), \zeta(5), \ldots` are generators.
  At other `s`, the functional equation
  `\zeta(s) = \pi^{s - 1/2} \Gamma((1-s)/2) / \Gamma(s/2)\, \zeta(1-s)`
  gives the canonical half `\operatorname{Re}(s) > 1/2` (or
  `\operatorname{Re}(s) = 1/2`, `\operatorname{Im}(s) \ge 0`): `\zeta(1/3)`
  is expressed through `\zeta(2/3)`.
* Dirichlet `L`-functions (:func:`gr_dirichlet_l`): `L(s, \chi) =
  L(s, \chi^*) \prod_{p \mid q} (1 - \chi^*(p) p^{-s})` for the primitive
  character `\chi^*` of conductor `f` inducing `\chi` mod `q`. At
  integers, `L(s, \chi^*) = f^{-s} \sum_{a \le f} \chi^*(a) \zeta(s, a/f)` in
  the Hurwitz normal form (`-(1/f) \sum \chi^*(a) \psi(a/f)` at `s = 1`,
  Bernoulli polynomials at `s \le 0`), so that `L(1, \chi_{-4}) = \pi/4`
  and `L(2, \chi_{-4}) = G`. Otherwise the functional equation
  `L(s, \chi) = \varepsilon(\chi) (f/\pi)^{1/2 - s}
  \frac{\Gamma((1 - s + a)/2)}{\Gamma((s + a)/2)} L(1 - s, \bar\chi)`
  (`\varepsilon(\chi) = \tau(\chi)/(i^a \sqrt f)` with the Gauss sum
  `\tau(\chi)` written through roots of unity, `a` the parity; conductors
  `f \le 240`, otherwise ``GR_UNABLE`` at the other half) gives the
  canonical half as for `\zeta` (at `s = 1/2`, the smaller Conrey label
  of `\chi`, `\bar\chi`), where `L(s, \chi^*)` is a generator of kind
  ``GR_TOWER_DIRICHLET_L`` (printed as ``dirichlet_l(chi(f, k), s)``);
  its conjugate is `L(\bar s, \bar\chi)`.
* `\zeta(s, a)` for `s` not an integer `\ge 2` (those are polygamma
  values): `-B_{1-s}(a)/(1-s)` at integers `s \le 0` (the continuation
  in `a`, also at the nonpositive integers `a`, as
  :func:`acb_hurwitz_zeta` at `a = 0`; ``GR_DOMAIN`` there for the other
  `s`, a term `(a + k)^{-s}` being infinite); at rational `a`
  (shifted into `(0, 1]`, denominators `q \le 240`) the character sum
  `\zeta(s, p/q) = \frac{q^s}{\varphi(q)} \sum_{\chi \bmod q} \bar\chi(p) L(s, \chi)`,
  so that the distribution relations
  `\sum_{j<m} \zeta(s, (a + j)/m) = m^s \zeta(s, a)` and the relations
  with `L`-functions and `\zeta(s)` hold by construction; at irrational
  `a` a generator of kind ``GR_TOWER_HURWITZ_ZETA`` (two arguments) at
  the shift with `0 < \operatorname{Re}(a) \le 1` (with the terms
  `(a + k)^{-s}`, principal powers). The multiplication relations on
  rational lines `a w + b` are not used at non-integer `s`.
* The Lerch transcendent `\Phi(z, s, a) = \sum_n z^n (n + a)^{-s}`
  (:func:`gr_lerch_phi`): `\zeta(s, a)` at `z = 1`,
  `v^{-s} \sum_{k<v} z^k \zeta(s, (a + k)/v)` at the roots of unity `z`
  of order `v \le 240`, `\operatorname{Li}_s(z)/z` at `a = 1` with an
  integer `s`; and `\operatorname{Li}_s(z) = z \Phi(z, s, 1)` at the roots
  of unity for non-integer `s`.
* `\operatorname{Li}_s(z)`: rational functions for `s \le 0`,
  `-\log(1 - z)` for `s = 1`, zeta values at `z = \pm 1`,
  `\operatorname{Li}_2(1/2) = \pi^2/12 - \log^2(2)/2`; at a root of unity
  `w = e^{2\pi i a/N}` (`3 \le N \le 240`, recognized numerically and
  verified exactly),
  `\operatorname{Li}_s(w) = N^{-s} \sum_{k=1}^{N} w^k \zeta(s, k/N)` in terms
  of Hurwitz zeta values (in their normal form), so that
  `\operatorname{Li}_2(i) = -\pi^2/48 + iG` and
  `\operatorname{Li}_3(i) = -3\zeta(3)/32 + i\pi^3/32`, and the
  distribution relations of polylogarithms at roots of unity follow.
  `\operatorname{Li}_2` at other arguments is reduced under the
  anharmonic group `\{z, 1-z, 1/z, 1/(1-z), 1-1/z, z/(z-1)\}` by the
  reflection and inversion formulas (with logarithms and `\pi^2`, and
  `-i\pi\log z` on the cut `z > 1`, as for :func:`acb_polylog`) to a
  canonical argument: in `[0, 1/2]` for real `z`, otherwise the orbit
  point in the region `\operatorname{Re}(w) \le 1/2`, `|w| \le 1`,
  `|w - 1| \le 1` (with `\operatorname{Im}(w) > 0` on its boundary).
  The zero test finds the distribution relations
  `\operatorname{Li}_2(y^n) = n \sum_{k<n} \operatorname{Li}_2(\zeta_n^k y)`
  (`|y| \le 1`, `n = 2, 3, 4`) between the generators present (matching
  the arguments through the anharmonic orbits, with the special values
  at `0, \pm 1, 1/2, 2`), and eliminates the latest generator occurring
  linearly in the relation after the eliminated ones are substituted. A
  rational value without a generator is adjoined at its canonical point
  when its height is smaller than that of the generator to eliminate
  (which bounds the process). Abel's five-term relation
  `\operatorname{Li}_2(x) + \operatorname{Li}_2(y) - \operatorname{Li}_2(xy)
  - \operatorname{Li}_2(\frac{x(1-y)}{1-xy}) - \operatorname{Li}_2(\frac{y(1-x)}{1-xy})
  = \log\frac{1-x}{1-xy} \log\frac{1-y}{1-xy}` is used likewise for
  real `x, y \in (0, 1)` among the orbit points of the generators'
  arguments, when the five values are found in the orbits.
  `\operatorname{Li}_s` for `s \ge 3` (other than at roots of unity) is
  reduced by the inversion formula
  `\operatorname{Li}_s(z) + (-1)^s \operatorname{Li}_s(1/z) =
  -\frac{(2\pi i)^s}{s!} B_s(\frac12 + \frac{\log(-z)}{2\pi i})`
  (with `\log(-z) = \log z + \pi i` on the cut `z > 1`) to `|z| < 1`,
  or `|z| = 1` with `\operatorname{Im}(z) > 0`; `\operatorname{Li}_3(1/2)
  = 7\zeta(3)/8 - \pi^2\log(2)/12 + \log^3(2)/6`. The zero test finds
  the distribution relations
  `\operatorname{Li}_s(y^n) = n^{s-1} \sum_{k<n} \operatorname{Li}_s(\zeta_n^k y)`
  between the generators (`n = 2, 3, 4`), which are linear with no
  elementary part at these arguments. Landen's relations for
  `\operatorname{Li}_3` (and the other functional equations of the
  higher polylogarithms) are not used.
* `K(m) = (\pi/2)\, {}_2F_1(1/2, 1/2; 1; m)`,
  `E(m) = (\pi/2)\, {}_2F_1(-1/2, 1/2; 1; m)` are evaluated through the
  hypergeometric normal form below: `K(0) = E(0) = \pi/2`, `E(1) = 1`
  (Gauss's sum), `K(1/2) = \Gamma(1/4)^2 / (4 \sqrt{\pi})` and `E(1/2)`
  (Gauss's second summation theorem, with contiguity), and the
  imaginary-modulus transformation `K(m) = K(m/(m-1))/\sqrt{1-m}`,
  `E(m) = \sqrt{1-m}\, E(m/(m-1))` to a canonical argument with
  `|m - 1| \le 1` (`\operatorname{Im}(m) \ge 0` on its boundary), which is
  Pfaff's transformation. The generators `K(m)`, `E(m)` are the closed
  forms of the basis values of the coset `(1/2, 1/2; 1) + \mathbb{Z}^3`;
  the singular values `m = k_r^2` for `r = 2, 3, 4` and their complements
  (Chowla-Selberg: `K` through `\Gamma` at rationals and `\pi`, `E`
  through the elliptic alpha function) are applied when such a generator
  would be created. The zero test finds Legendre's relation
  `E(m) K(1-m) + E(1-m) K(m) - K(m) K(1-m) = \pi/2` and Landen's
  transformation `K(m) = (1 + k_1) K(k_1^2)`,
  `E(m) = (1 + s) E(k_1^2) - s K(m)` (`s = \sqrt{1-m}`,
  `k_1 = (1-s)/(1+s)`) between the generators present, also through the
  transformation (the partner of a generator at `a` may sit at `1 - a`,
  `(a-1)/a`, `1/(1-a)` or `1/a`), and eliminates the latest generator
  involved (degree one, with square roots as root generators).
* `{}_pF_q(a_1, \ldots, a_p; b_1, \ldots, b_q; z)` (the generators, of
  kind ``GR_TOWER_HYPGEOM``, are not regularized; see *Hypergeometric
  functions* below).
* Lambert W: `W_0(0) = 0`, `W_{0,-1}(-1/e) = -1`; the other relations are
  found by Richardson's algorithm (above).

Hypergeometric functions
-------------------------------------------------------------------------------

The values `F(P; z) = {}_pF_q(a; b; z)`, `P = (a_1, \ldots, a_p; b_1,
\ldots, b_q)`, are brought to a normal form in three layers
(``lazy_hypgeom.c``):

1. Degenerate parameters: `z = 0`, cancellation of `a_i = b_j`,
   terminating series (`a_i` a nonpositive integer), poles (`b_j` a
   nonpositive integer: ``GR_DOMAIN``; the regularized function is
   finite there, `\tilde F = \prod_i (a_i)_{n+1} z^{n+1} \tilde F(a + n
   + 1; n + 2, b' + n + 1; z)` for `b_1 = -n`), divergence (`p > q +
   1`), Gauss's sum at `z = 1`, `{}_0F_0 = e^z`, `{}_1F_0 = (1-z)^{-a}`,
   and the reduction of the order when `a_i - b_j = m` is a positive
   integer: `F = (\theta + b_j)_m F_{\mathrm{red}} / (b_j)_m`,
   `\theta = z\, d/dz`, with `F_{\mathrm{red}}` the function without
   `a_i`, `b_j` (expanded exactly through `\theta^k = \sum_l S(k, l) z^l
   D^l` and `D^l F = (a)_l / (b)_l\, F(P + l e)`, `e = (1, \ldots, 1)`).

2. Transformations of one term, choosing a canonical member of the orbit
   of `(P, z)`: for `{}_2F_1`, Pfaff's transformation `F(a, b; c; z) =
   (1-z)^{-a} F(a, c-b; c; z/(z-1))` when `|z - 1| > 1` (or `|z - 1| = 1`
   and `\operatorname{Im}(z) < 0`), then the smaller (in a fixed order of
   the cells) of `P` and its Euler transform `(c-a, c-b; c)` (factor
   `(1-z)^{c-a-b}`), off the cut `z \ge 1`; for `{}_1F_1`, Kummer's
   transformation `F(a; b; z) = e^z F(b-a; b; -z)` to
   `\operatorname{Re}(z) \le 0`. Powers with an irrational exponent are
   written `x^{e_0} x^n`, `e = e_0 + n`, `0 < \operatorname{Re}(e_0) \le
   1`, so that contiguous parameters share one generator.

3. The contiguity module of the coset `P + \mathbb{Z}^{p+q}`. With `r =
   \max(p, q+1)` and `v(P) = (F, \theta F, \ldots, \theta^{r-1} F)`, a
   unit shift of a parameter acts as `v(P') = M v(P)` with `M = I +
   C/a_i` (raising `a_i`) or `I + C/(b_j - 1)` (lowering `b_j`), where
   `C` is the companion matrix of `\theta` modulo the hypergeometric
   operator `\theta \prod_j (\theta + b_j - 1) - z \prod_i (\theta +
   a_i)`; `\det(I + C/a_i)` is a multiple of `\prod_j (b_j - 1 - a_i)`, so
   the shifts are invertible unless the operator is reducible. Every
   value in the coset is thus a rational combination of the `r` basis
   values `G_k = F(P_0 + k e)` at the canonical cell `P_0` (parameters
   shifted to `0 < \operatorname{Re} \le 1`). A reducible pair (`b_j -
   a_i` a positive integer) keeps `b_j - a_i = 1` in the cell, and the
   paths stay in the chamber `b_j - a_i \ge 1`; it gives the relation
   `(\theta + a_i) F(P_0) = a_i F_{\mathrm{red}}(P_0)` between the basis
   values. A numerator `a_i = 1` (a positive integer, in its cell) gives
   another: the operator factors as `\theta L_1` with `L_1 = \prod_j
   (\theta + b_j - 1) - z \prod_{k \ne i} (\theta + a_k)`, so that `L_1 F
   = \prod_j (b_j - 1)`.

Known values give further linear relations between the basis values:
the entries of a declarative table, each a pattern (parameters affine
in symbols, the argument free or a rational constant) with a value in a
small expression language and optional conditions, for example::

    "3 2 | a, b, c | 1+a-b, 1+a-c | 1 | gamma(1+a/2)*gamma(1+a-b)*... | re: 1+a/2-b-c"
    "0 1 | | 1/2 | z | cosh(2*sqrt(z))"
    "2 1 | 1/2, 1/2 | 1 | z | 2*elliptic_k_gen(z)/pi"

Every entry is matched against every point of the coset (in a box of
shifts of its symbols) and of the images of the coset under the
transformations of layer 2; the relations are reduced in order
(reducibility relations, then entries with a constant argument, then by
distance from the cell), each value being computed only when its
relation is independent of the previous ones, until `F(P)` is
determined or the relations are exhausted. The basis values which
remain undetermined become generators. Thus the two entries for
`{}_0F_1(; 1/2; z)` and `{}_0F_1(; 3/2; z)` give every `{}_0F_1(; n +
1/2; z)` in terms of `\cosh(2\sqrt z)` and `\sinh(2 \sqrt z)`, the
entries for `K` and `E` give every `{}_2F_1` with parameters in `(1/2,
1/2, 1) + \mathbb{Z}^3`, the error function entry gives every
`{}_1F_1(n + 1/2; m + 3/2; z)` (`m \ge n`, with the reducibility
relation `(\theta + 1/2) F = e^z/2`), and the summation theorems at `z =
\pm 1, 1/2` apply to every contiguous function of a summable one (the
"contiguous" variants of Kummer's, Gauss's and Bailey's theorems). The
current table contains the elementary cases of `{}_0F_1`, `{}_1F_1`
and `{}_2F_1`, the error function, Kummer's second formula `{}_1F_1(a;
2a; z) = e^{z/2} {}_0F_1(; a + 1/2; z^2/16)`, `K` and `E`, the sums of
Gauss (`z = 1`), Kummer (`z = -1`), Gauss and Bailey (`z = 1/2`), and
the sums of Dixon, Watson and Whipple for `{}_3F_2` at `z = 1` (the
latter only at the point itself: the contiguity module at `z = 1`,
where the rank drops to `p - 1`, is not implemented). Each entry is checked
numerically by the tests.

*Links.* Some entries relate the coset to the values of another
hypergeometric function at another argument; they are used only when a
generator of the context is at that argument (an ``anchor:`` condition,
tested numerically up to Pfaff's transformation, before the exact
evaluation), and one at a time (not inside the evaluation of another
link), so that the generators of the context are reused rather than
duplicated, and no chain of transformations is explored. These are the
quadratic transformations of `{}_2F_1`,

.. math ::

    F(a, b; 2b; z) = (1 - z/2)^{-a} F(\tfrac{a}{2}, \tfrac{a+1}{2}; b + \tfrac12; \tfrac{z^2}{(2-z)^2}) \quad (z \notin [1, \infty)),

    F(a, b; 1 + a - b; z) = (1 + z)^{-a} F(\tfrac{a}{2}, \tfrac{a+1}{2}; 1 + a - b; \tfrac{4z}{(1+z)^2}) \quad (|z| < 1),

    F(a, b; \tfrac{a+b+1}{2}; z) = F(\tfrac{a}{2}, \tfrac{b}{2}; \tfrac{a+b+1}{2}; 4z(1-z)) \quad (\operatorname{Re}(z) < 1/2),

and their inverses (with the principal roots `\sqrt{z}`, `\sqrt{1 - z}`,
valid off the cut); for instance `F(1/4, 3/4; 1; 4z/(1+z)^2) = \sqrt{1 +
z}\, F(1/2, 1/2; 1; z) = 2\sqrt{1 + z} K(z)/\pi` once `K(z)` is a
generator. The complete elliptic integrals are linked by Landen's
transformation in the same way (the zero test also finds it between two
generators already present, see *Special functions*): with `k' = \sqrt{1
- m}` and `k_1 = (1 - k')/(1 + k')`,

.. math ::

    K(m) = (1 + k_1) K(k_1^2), \quad E(m) = (1 + k') E(k_1^2) - k' K(m)
    \quad (m \notin [1, \infty)),

and a new `K(m)` or `E(m)` is expressed through a generator `K`, `E` at a
point of the Landen chain of `m`, down (`m \to k_1^2`) or up (`m \to
4s/(1+s)^2`, `s = \pm \sqrt{m}`, `|s| < 1`), within three steps, rather
than becoming a generator.

The generators of kind ``GR_TOWER_HYPGEOM`` have the argument `z` and
the additional arguments `a_1, \ldots, a_p, b_1, \ldots, b_q` (see
:func:`gr_tower_adjoin_special_multi_flat`), and print as
``hypgeom_2f1(a, b, c, z)`` (``hypgeom_0f1(b, z)``,
``hypgeom_1f1(a, b, z)``, ``hypgeom_pfq([a...], [b...], z)``).

Modular forms and theta functions
-------------------------------------------------------------------------------

The modular functions and forms of a point `\tau` of the upper half-plane
(``lazy_modular.c``) are reduced exactly to the fundamental domain

.. math ::

    \mathcal{F} = \{ -1/2 < \operatorname{Re}(\tau) \le 1/2, |\tau| > 1 \}
        \cup \{ |\tau| = 1, 0 \le \operatorname{Re}(\tau) \le 1/2 \},

`\tau_0 = g \tau` with `g \in PSL_2(\mathbb{Z})` found numerically
(:func:`acb_modular_fundamental_domain_approx`) and then verified and
corrected by exact comparisons, and the values at `\tau` follow from the
transformation laws with exact multipliers (the eighth roots of unity of
:func:`acb_modular_theta_transform`, the 24th roots of unity of
:func:`acb_modular_epsilon_arg`, powers of `c\tau + d`, and `E_2(g\tau) =
(c\tau+d)^2 E_2(\tau) + 6c(c\tau+d)/(\pi i)`). At `\tau_0`, every value is
expressed through `\Lambda = \lambda(\tau_0)` and the complete elliptic
integrals at `\Lambda` (canonical arguments: for `\tau_0 \in \mathcal{F}`,
`|\Lambda - 1| \le 1`, with `\operatorname{Im}(\Lambda) > 0` on the
boundary `\operatorname{Re}(\tau_0) = 1/2`):

.. math ::

    \theta_3^2 = 2K(\Lambda)/\pi, \quad
    \theta_2 = \theta_3 \Lambda^{1/4}, \quad
    \theta_4 = \theta_3 (1 - \Lambda)^{1/4}, \quad
    \eta = (\theta_2 \theta_3 \theta_4 / 2)^{1/3},

    j = 256 \frac{(1 - \Lambda + \Lambda^2)^3}{\Lambda^2 (1 - \Lambda)^2}, \quad
    E_4 = \theta_3^8 (1 - \Lambda + \Lambda^2), \quad
    E_6 = \theta_3^{12} \frac{(1 + \Lambda)(2 - \Lambda)(1 - 2\Lambda)}{2}, \quad
    E_2 = \theta_3^4 \left(\frac{3 E(\Lambda)}{K(\Lambda)} - 2 + \Lambda\right),

with principal roots (the right branches on `\mathcal{F}`), and `E_{2k}`,
`k \ge 4`, by the recurrence of the Weierstrass coefficients. The
relations between these values (Jacobi's identity, `j = E_4^3/\Delta`,
`E_4^3 - E_6^2 = 1728 \Delta`, `E_8 = E_4^2`, ...) thus hold by
construction, and the relations of `K` and `E` known to the zero test
(Legendre's relation, which is the transformation `\tau \to -1/\tau`,
and Landen's transformation, `\tau \to 2\tau`) apply to the modular
values as well. The generator `\lambda(\tau_0)` has the kind
``GR_TOWER_MODULAR_LAMBDA`` and prints as ``modular_lambda(tau0)``.

At a CM point (`\operatorname{Re}(\tau_0)` and `|\tau_0|^2` rational,
`\tau_0` a root of `A\tau^2 + B\tau + C` of discriminant `D`, class
number `h \le 16`, `|D| \le 100000`), `\Lambda` is the algebraic number
selected numerically among the roots of the factor of `\sum_k H_k N^k
M^{h-k}` (`H` the Hilbert class polynomial, `N = 256(1 - x + x^2)^3`, `M
= x^2 (1 - x)^2`) which vanishes at it, so that `j(\tau)` and
`\lambda(\tau)` are algebraic (`j(\sqrt{-5}) = 632000 + 282880 \sqrt 5`,
`j((1 + \sqrt{-163})/2) = -640320^3`). For the fundamental
discriminants of class number one with `|D|` at most the option
``GR_TOWER_OPT_GAMMA_LATTICE_LIMIT`` (36 by default), the formula of
Chowla and Selberg gives the periods,

.. math ::

    |\eta(\tau_0)|^2 = \frac{1}{\sqrt{2\pi |D|}} \prod_{j=1}^{|D|-1}
        \Gamma(j/|D|)^{w \chi_D(j)/4},

with `\eta(\tau_0)^2 = e^{\pi i \operatorname{Re}(\tau_0)/6} |\eta(\tau_0)|^2`
and `\theta_3^2 = 2^{2/3} \eta^2 (\Lambda(1 - \Lambda))^{-1/6}` (times
the sixth root of unity selected numerically), so that `\eta(i) =
\Gamma(1/4)/(2\pi^{3/4})`, `E_2(i) = 3/\pi`, and the elliptic integrals
at the singular moduli of these discriminants are expressed through
gamma values (generalizing `K(1/2)`). Elsewhere `K(\Lambda)` and
`E(\Lambda)` are generators.

*Weierstrass functions.* For the lattice `\mathbb{Z} + \tau \mathbb{Z}`
(the conventions of :func:`acb_elliptic_p`), with `\theta_k =
\theta_k(0, \tau)`,

.. math ::

    e_1 = \tfrac{\pi^2}{3} (\theta_3^4 + \theta_4^4), \quad
    e_2 = \tfrac{\pi^2}{3} (\theta_2^4 - \theta_4^4), \quad
    e_3 = -\tfrac{\pi^2}{3} (\theta_2^4 + \theta_3^4), \quad
    g_2 = \tfrac{4 \pi^4}{3} E_4, \quad g_3 = \tfrac{8 \pi^6}{27} E_6,

    \wp(z) = e_1 + \left(\frac{\pi \theta_3 \theta_4 \theta_2(z)}{\theta_1(z)}\right)^2, \quad
    \wp'(z) = -2\pi^3 (\theta_2 \theta_3 \theta_4)^2 \frac{\theta_2(z) \theta_3(z) \theta_4(z)}{\theta_1(z)^3}, \quad
    \sigma(z) = \frac{\theta_1(z)}{\pi \theta_2 \theta_3 \theta_4} e^{\pi^2 E_2 z^2 / 6},

so that the differential equation `\wp'^2 = 4\wp^3 - g_2 \wp - g_3 =
4(\wp - e_1)(\wp - e_2)(\wp - e_3)`, the periodicity and the values at
the half periods hold by construction, and the addition and
multiplication theorems through those of the theta functions (below).
The Weierstrass zeta function (which needs `\theta_1'(z)/\theta_1(z)`)
and the inverse of `\wp` are not implemented.

The theta functions of a variable `z` (``jacobi_theta(z, tau)``) are
reduced in both arguments: `\tau` to `\tau_0` as above (`z \to
-z/(c\tau + d)` with the factor `\sqrt{i/(c\tau+d)}\, e^{-\pi i c
z^2/(c\tau + d)}`), and then `z` modulo the half-lattice, `z = z_1 + (p +
q\tau_0)/2` with `z_1 = x + y\tau_0`, `x, y \in [-1/4, 1/4)` (exact
floors; `p`, `q` adjusted on the boundary of that box, see below), with
the characteristics: in terms of `\theta[a, b](z) = \sum_n
\exp(\pi i (n + a/2)^2 \tau + 2\pi i (n + a/2)(z + b/2))` (`\theta_1 =
-\theta[1,1]`, `\theta_2 = \theta[1,0]`, `\theta_3 = \theta[0,0]`,
`\theta_4 = \theta[0,1]`),

.. math ::

    \theta[a,b](z_1 + (p + q\tau)/2) = e^{-\pi i q^2 \tau/4 - \pi i q (z_1
        + (b+p)/2)} \theta[a + q, b + p](z_1),

`\theta[a + 2, b] = \theta[a, b]`, `\theta[a, b + 2] = (-1)^a \theta[a, b]`,
and the parity `\theta[a, b](-z) = (-1)^{ab} \theta[a, b](z)`, to the
reduced point `z_0 = \pm z_1` in the fundamental domain of `z \to \pm z +
(\mathbb{Z} + \tau_0\mathbb{Z})/2` given by `0 < y < 1/4`, `-1/4 \le x <
1/4`, or `y \in \{0, 1/4\}`, `0 \le x \le 1/4` (so that `z` and `-z`, or a
point and its reduced point, have the same reduced point). At `\tau_0 = i` and `\tau_0 = e^{\pi i/3}`, the
stabilizer of `\tau_0` acts on `z` (by `z \to iz`, respectively
`e^{\pm \pi i/3} z`), and the least reduced point of the orbit is taken.
At the reduced point `z_0 \ne 0`, `\theta_1(z_0, \tau_0)` and
`\theta_4(z_0, \tau_0)` are generators of the kind
``GR_TOWER_JACOBI_THETA`` (two arguments, printed as
``jacobi_theta_1(z0, tau0)``, ``jacobi_theta_4(z0, tau0)``), and
`\theta_2(z_0, \tau_0)`, `\theta_3(z_0, \tau_0)` are algebraic
generators over them, with the same definition kind and Jacobi's
relations

.. math ::

    \theta_2(z)^2 \sqrt{1 - \Lambda} = \theta_4(z)^2 \sqrt{\Lambda} - \theta_1(z)^2, \quad
    \theta_3(z)^2 \sqrt{1 - \Lambda} = \theta_4(z)^2 - \theta_1(z)^2 \sqrt{\Lambda}

as moduli (adjoined once per point: later evaluations at an equal point
find them by their definition, whatever the representation of the
point). The quasi-periodicity, the transformation laws and the quartic
relations between the theta functions of `z` thus hold by construction.

*Related points.* Every reduced point `z_0` met at `\tau_0` is recorded
in the context with the method of its evaluation, chosen when it is
first met from the points recorded before (relations modulo the
half-lattice, found numerically and verified exactly), so that the
values at related points are not independent generators:

* a torsion point `z_0 = (a + b\tau_0)/N` with `3 \le N \le 6` (the
  reduced point: `z = 1/3` has `z_0 = 1/6`, while `(1 + 2\tau)/5` with
  `z_0 = -1/5 + \tau/10` is not used): the values are algebraic over the theta constants,
  through the root `X = \tilde\wp(z_0)` of the division polynomial
  `\psi_N` (in the units `\tilde\wp = \wp/(\pi^2\theta_3^4)`), the
  ratios `\theta_k(z_0)/\theta_1(z_0)` (square roots of `X - e_j` up to
  constants) and a root of `\theta_1(z_0)^{N^2}` given by the
  multiplication formula at the lattice point `N z_0`;
* `d z_0 = \pm n u` with `u` recorded, `n \le 6`, `d \le 3`: the
  multiplication formulas `\theta_k(nu) = \theta_3 (\theta_1(u)/\theta_3)^{n^2}
  S_k(n)` with `S_k(n)` rational in the ratios at `u` (the recursion of
  the division polynomials), then for `d > 1` the division as for torsion
  points (with `\tilde\wp(nu) = \tilde\wp(u) - \psi_{n-1}\psi_{n+1}/\psi_n^2`);
* `z_0 = \pm(u + v)` with `u`, `v` and `u - v` recorded: rational, by
  the addition formulas `\theta_k(u+v)\theta_k(u-v)\theta_4^2 =
  \theta_k(u)^2\theta_4(v)^2 - \theta_{5-k}(u)^2\theta_1(v)^2`;
* `u`, `u + z_0`, `u - z_0` recorded: `\theta_1(z_0)^2` and
  `\theta_4(z_0)^2` from two of these addition formulas (linear
  equations), the other values from Jacobi's relations;
* `z_0 = \pm(nu \pm v)` otherwise (`n \le 3`): the ratios at `z_0`
  rational in the values at `u` and `v` (Jacobi's addition formulas for
  `\mathrm{sn}, \mathrm{cn}, \mathrm{dn}` from the multiple `nu`), and
  `\theta_1(z_0)` a generator;
* `2 z_0 = \pm(u \pm v)` otherwise: with `P = z + w`, `Q = z - w` and
  `S_{ab} = \theta_a(P)\theta_b(Q) + \theta_b(P)\theta_a(Q)`, Jacobi's
  formulas for `\theta_a(z+w)\theta_b(z-w)` give `\theta_1\theta_4(z) /
  \theta_2\theta_3(z) = S_{14}/S_{23}` and its images under the
  permutations of 2, 3, 4, so that the ratios `r_k = \theta_k/\theta_1`
  at `z = (P + Q)/2` are monomials in the square roots of the `S_{ab}`
  (`r_2 = \sqrt{S_{23} S_{24} / (S_{14} S_{13})}`, and so on, with
  `r_2 r_3 r_4` rational; the signs of the square roots chosen
  numerically, unique by the enclosures, and checked for consistency);
  `\theta_1(z_0)` is a generator. The point `w = (P - Q)/2` then has
  `r_k(w)/r_k(z)` rational in the `S_{ab}` and the differences `D_{ab}`
  (`r_4(w)/r_4(z) = \theta_3 S_{12} / (\theta_2 D_{13})`, say), and
  `\theta_1(z)\theta_1(w)` a monomial in the same square roots and
  `\sqrt{D_{12} D_{13}}`;
* otherwise the generators above.

*Rebasing.* Where the methods above would make the values at a new
point algebraic over those at points recorded with generators of their
own (roots of polynomials, square roots), the roles are exchanged: the
new point (or a related one) gets the generators, the earlier points
become rational over it, and the context records the value of each
earlier generator as a rational function of the new ones. When
elements involving such a generator and elements involving one of the
new ones meet in an operation, the generator becomes algebraic in
their tower with the linear modulus `X - \text{value}` (the generators
of the values are moved before it: one rebuild of the tower for all the
records that apply), so that the elements built on it, before or after,
meet those built on the new generators. (Not when a tower merely
contains both: the towers of a long session hold generators of
unrelated computations.) The cases:

* `d z_0 = \pm n u` with `d > 1` and `u` of generators: with `an + cd
  = 1`, the point `b = a z_0 \pm c u` (so that `n b = z_0`, `d b = \pm
  u`) gets the generators, `u` and `z_0` are its multiples (the division
  by 3 of `\theta(3z)` would be a root of degree 9 with large
  coefficients: `\theta(3z)`, `\theta(2z)`, `\theta(z)` cost the same
  in any order);
* `u`, `P = u + x`, `Q = u - x` with `P`, `Q` of generators (`z = (P + Q)/2`
  with `P = z + w`, `Q = z - w` first, and `x = w`; or `x = z_0` in the
  case of the linear equations above): `x` gets the generators, `P`
  the ratios of Jacobi's formulas from `u`, `x` (its `\theta_1` stays a
  generator) and `Q` the addition formulas;
* `z_0 = \pm(2u \pm v)` with `v` of generators: as the previous case with
  `x = -(u \pm v)`, `P = v`, `Q = z_0`.

The values at the points are then those of the order in which the
points with generators came first, whatever the actual order (the
addition formulas with the values at `z`, `w`, `z \pm w` take 0.06 to
0.25 seconds in all 24 orders, against up to seven seconds without the
rebasing). Elements built on the earlier
generators other than the values themselves (`\wp(z + w)` as a rational
function of `\theta_1(z + w)` and `\theta_4(z + w)`, say) become
rational functions of the new generators, which may be large.

The values are cached with the record (those of the 32 most recently
used points; the ratios `\theta_k/\theta_3` at the 8 most recent points
`\tau_0` and the theta constants themselves at the last one, together
with the anchor chosen for each `\tau_0`, so that values computed later
come from the same anchor), so that repeated evaluations at a point cost
a lookup; the bound keeps the cache from holding the towers of a long
session alive. At the CM points where `\lambda` has degree more than 2
(from a class polynomial: `\lambda(1/2 + i)` is a root of a quartic
with coefficients of five digits), the torsion points and the divisions
are not used (the roots of the division polynomials over that field are
costly: tens of seconds for `\psi_3`), nor the other relations except
the doubling `z_0 = \pm 2u` (the inverses of the denominators of their
formulas in that field take minutes, hours for `n = 6`): the values
there are generators, and identities between the values at related
points are not decided.
Thus the duplication formula `\theta_1(2z)\theta_2\theta_3\theta_4 =
2\theta_1(z)\theta_2(z)\theta_3(z)\theta_4(z)`, the addition formulas and
the addition theorem of `\wp` are decided whatever the order of
evaluation of `\theta(z)`, `\theta(w)`, `\theta(z \pm w)`, as rational
identities in the values at the points with generators (the forms
`S_{ab}`, `D_{ab}` above remain for the cases the rebasing leaves out).
Relations among more than three points (`z_0 = u + v + w`, say) are not
used.

*Commensurable points.* A point `\tau_0 \in \mathcal{F}` which is not a
CM point is compared with the points `\tau_1` whose `\lambda(\tau_1)` is
already a generator of the context (the *anchors*): when `\tau_0 =
\gamma \tau_1` with `\gamma \in GL_2^+(\mathbb{Q})` of determinant
`2^k 3^l` (an integer relation between `1, \tau_1, \tau_0, \tau_0 \tau_1`
found by LLL and verified exactly; `k \le 4` when `l = 0`, otherwise
`l \le 2` and `k \le 1`), the values at `\tau_0` are computed from those
at `\tau_1` rather than from a new generator. Two points commensurable
beyond these limits have generators of their own, and the points
commensurable with both (within the limits) could come from either: of
the anchors found, those whose values come from the most recently used
generator `\lambda` are preferred (the *root*: the values of an
evaluation in progress, `j(\tau)`, `\eta(\tau)` and `\eta(3\tau)` say,
then come from one root, rather than `\eta(3\tau)` from the root of an
earlier evaluation at `(3\tau + 1)/4`, with a cheaper chain to `3\tau`,
which would leave the identity between them undecided), and of these the
one with the shortest chain of steps. A chain with a tripling and a
further step (a halving, another tripling) gives values in towers of
large degree, where everything at `\tau_0` is costly (the theta
functions of `z` at `4i/\pi` through two thirdings from `1/3 + \pi i/4`:
3.5 seconds, against 0.01 at a generator of its own; `j(\tau)`,
`\eta(3\tau)` after `\lambda((2\tau + 1)/3)`: minutes). There
`\lambda(\tau_0)` is a new generator *linked* to the anchor: the values
at `\tau_0` are those of a new generator (the fourth roots of
`\lambda` and `1 - \lambda`, `K(\lambda)`, `\theta_3 = \sqrt{2K/\pi}`,
`E(\lambda)`), and the context records their values through the chain,
computed when first needed, as rebase records (see the theta functions
of `z` above): where values from the two meet, those at `\tau_0` are the
values of the chain, and the identities between them hold as with the
chain (Jacobi's modular equation of degree 3 between `\lambda(\tau)`
and `\lambda(3\tau)` with `3\tau` from a linked `(3\tau + 1)/2`);
elsewhere they cost what they cost in a new field. With `\gamma = g M`, `g
\in SL_2(\mathbb{Z})`, `M = \begin{pmatrix} a & b \\ 0 & d
\end{pmatrix}`, the steps are the duplication formulas

.. math ::

    \theta_3(2w)^2 = \frac{\theta_3^2 + \theta_4^2}{2}, \quad
    \theta_4(2w)^2 = \theta_3 \theta_4, \quad
    \theta_2(2w)^2 = \frac{\theta_3^2 - \theta_4^2}{2}, \quad
    2 E_2(2w) - E_2(w) = \theta_3(2w)^4 + \theta_2(2w)^4

(and their inverses for `w \to w/2`; the square roots are chosen
numerically), Jacobi's modular equation of degree 3 for `w \to 3w` and
`w \to w/3`: with `u = (\theta_2(w)/\theta_3(w))^{1/2}` and `v` the
same at `3w`,

.. math ::

    u^4 - v^4 - 2uv(1 - u^2 v^2) = 0, \quad
    \frac{\theta_3(w)^2}{\theta_3(3w)^2} = 1 + \frac{2v^3}{u}, \quad
    3 E_2(3w) - E_2(w) = \sum_{k=2}^4 \theta_k(w)^2 \theta_k(3w)^2

(the root of the quartic equals `\pm v` numerically and is adjoined as an
algebraic generator), the translations and the transformation laws of
`g`. The modular equations of these levels thus hold by construction:
Landen's `\lambda(2\tau) = ((1 - k')/(1 + k'))^2`, the modular polynomial
`\Phi_2(j(\tau), j(2\tau)) = 0`, eta quotients such as `\theta_4(2\tau) =
\eta(\tau)^2/\eta(2\tau)`, the Hauptmodul `t = (\eta(\tau)/\eta(3\tau))^{12}`
of `\Gamma_0(3)` with `j(\tau) = (t + 27)(t + 243)^3/t^3` and `j(3\tau) =
(t + 27)(t + 3)^3/t`, and so on. When only `\lambda` is needed, the
computation uses the ratios `\theta_k/\theta_3` and no elliptic integral.
Matrices `M` with `3 \mid a` and `3 \mid d` (`\tau_1 + 1/3`, say), other
levels (primes `\ge 5`), and levels with a factor 3 inside a conjugation
are not used: values at points related by such transformations are
independent generators, and the corresponding identities are not
decided (``T_UNKNOWN``).

The factors `\exp(\pi i \alpha)` of the transformation laws are written
`\exp(\pi i \operatorname{Re}(\alpha)) \exp(-\pi \operatorname{Im}(\alpha))`
when `\operatorname{Re}(\alpha)` is rational (a root of unity times a real
exponential), which keeps the multiplicative relations between them of
small degree.

Conjugation maps `F(u)` to `F(\bar u)` (`W_k(u)` to `W_{-k}(\bar u)`) off
the branch cuts, and the values at real points are real generators,
which the real fields accept.

Further relations -- Landen's relations of the trilogarithm,
polylogarithms at other algebraic arguments, isogenies of elliptic
integrals beyond Landen's -- are not applied yet; see the design notes.

Types and macros
-------------------------------------------------------------------------------

.. type:: gr_tower_gen_struct

    A generator of a tower: its current kind (algebraic, or one of
    ``GR_TOWER_EXP``, ``GR_TOWER_LOG``, ``GR_TOWER_PI``, ``GR_TOWER_TAN``,
    ``GR_TOWER_ATAN``), the kind of its
    definition (which differs from the current kind after an elimination;
    for algebraic generators also ``GR_TOWER_ROOT``,
    ``GR_TOWER_ROOT_OF_UNITY`` and :macro:`GR_TOWER_TAN_PI`),
    a status flag, its index (`k` or `j`), the context of `F_k` if
    algebraic, the argument of `\exp` or `\log` as a flat element, an
    enclosure, a name, a definition id and its definition order.
    ``GR_TOWER_GEN(T, d)`` is the generator with definition order *d*,
    ``GR_TOWER_STEP(T, k)`` the algebraic generator `a_{k+1}` and
    ``GR_TOWER_TRANS(T, j)`` the transcendental generator `t_{j+1}`.

.. type:: gr_tower_struct

.. type:: gr_tower_t

    The tower object.

.. macro:: GR_TOWER_TAN_PI

    The definition kind of the real algebraic generator
    `\tan(\pi / M)` (``def_param`` `M`) of the tangent normal form of the
    real trigonometric constants (see the lazy real field).

.. macro:: GR_TOWER_STATUS_PROVEN
           GR_TOWER_STATUS_DYNAMIC
           GR_TOWER_STATUS_INDEPENDENT
           GR_TOWER_STATUS_SCHANUEL
           GR_TOWER_STATUS_CONJECTURAL

    Status of a generator. For a step: the minimal polynomial is known
    to be irreducible over the field below (``PROVEN``), or merely
    assumed to be so (``DYNAMIC``). For a transcendental generator:
    ``PROVEN`` (its alias ``INDEPENDENT``) means that the generator is
    known to be algebraically independent of all the generators
    preceding it in the definition order, which is what the structural
    arguments of the module use (currently only `\pi` at the head of a
    tower has this status: a proof of transcendence alone does not give
    it when another transcendental generator precedes); ``SCHANUEL``
    and ``CONJECTURAL`` mean that the independence is conjectural (under
    Schanuel's conjecture, or beyond it for the values of special
    functions).

Memory management and basic access
-------------------------------------------------------------------------------

.. function:: void gr_tower_init(gr_tower_t T, gr_ctx_t base)

    Initializes *T* to the trivial tower over the field *base*.
    The base context must remain valid during the lifetime of *T*.

.. function:: void gr_tower_clear(gr_tower_t T)

.. function:: gr_tower_struct * gr_tower_heap_init(gr_ctx_t base)
              void gr_tower_heap_clear(gr_tower_struct * T)
              int gr_tower_get_str(char ** s, const gr_tower_t T)
              int gr_tower_gen_get(gr_ptr res, const gr_tower_t T, slong d)
              gr_ctx_struct * gr_tower_field_ptr(const gr_tower_t T)
              slong gr_tower_length_si(const gr_tower_t T)
              slong gr_tower_num_gens_si(const gr_tower_t T)
              slong gr_tower_num_trans_si(const gr_tower_t T)
              const char * gr_tower_gen_name(const gr_tower_t T, slong d)
              int gr_tower_gen_status(const gr_tower_t T, slong d)
              int gr_tower_gen_kind(const gr_tower_t T, slong d)

    Heap-allocated towers, the string written by :func:`gr_tower_write`,
    the generator with definition order *d* as an element of the top
    field, and non-inline accessors (for language bindings; the Python
    wrapper's ``gr_tower`` class is built on them).

.. function:: void gr_tower_set(gr_tower_t res, const gr_tower_t T)

    Sets *res* to a copy of *T*. The copy has its own contexts, so
    refinements of one tower do not affect the other.

.. function:: slong gr_tower_length(const gr_tower_t T)

    The number of algebraic steps `n`. The number of transcendental
    generators is ``T->num_trans`` and the total number of generators
    is :func:`gr_tower_num_gens`.

.. function:: slong gr_tower_num_gens(const gr_tower_t T)
              slong gr_tower_prefix_length(const gr_tower_t T, slong p)
              slong gr_tower_prefix_num_trans(const gr_tower_t T, slong p)
              slong gr_tower_order_code(const gr_tower_t T, slong d)
              slong gr_tower_gid_order(const gr_tower_t T, slong gid)
              slong gr_tower_find_def_order(const gr_tower_t T, ulong def_id)

    The number of generators, the number of algebraic, respectively
    transcendental, generators among the first *p* generators in
    definition order, the code of the generator with definition order *d*
    (`k > 0` for the algebraic generator `a_k`, `-j` for the
    transcendental generator `t_j`), and the definition order of the
    generator with the given gid (a tower-local id which is stable under
    reorderings) or definition id (-1 if absent). The definition order
    of a generator can change when the tower is reordered (by the
    relation search); ``T->structure_version`` is incremented then.

.. function:: void gr_tower_gen_set_def_id(gr_tower_t T, slong d, ulong def_id)

    Sets the definition id of the generator with definition order *d*.
    Definition ids identify generators across towers: the merges
    (:func:`gr_tower_absorb`) identify generators with equal nonzero
    ids. The lazy fields assign them to the generators they create; a
    user of the tower API may assign its own.

.. function:: void gr_tower_set_prefix(gr_tower_t res, const gr_tower_t T, slong p)
              void gr_tower_set_subset(gr_tower_t res, const gr_tower_t T, const int * mark)

    Sets *res* to a copy of the first *p* generators of *T* in definition
    order, respectively of the generators with nonzero ``mark[d]``
    (indexed by definition order, in that order). The marked set must be
    closed under the dependencies of the definitions: the generators in
    the moduli and arguments of marked generators must be marked. The
    statuses carry over (a modulus irreducible over the field of all
    preceding generators is irreducible over a subfield).

.. function:: gr_ctx_struct * gr_tower_base(const gr_tower_t T)
              gr_ctx_struct * gr_tower_field_at(const gr_tower_t T, slong k)
              gr_ctx_struct * gr_tower_field(const gr_tower_t T)

    Contexts of `F_0`, `F_k` and `F_n`. Elements of the tower are
    elements of these contexts. Note that the context of `F_k` is a
    generic quotient ring: its :func:`gr_is_zero` and :func:`gr_inv`
    are complete only when the tower has been fully refined; use the
    ``gr_tower`` versions of these operations, or the context created by
    :func:`gr_ctx_init_tower_field`, for complete versions.

.. function:: slong gr_tower_step_degree(const gr_tower_t T, slong k)
              slong gr_tower_degree(const gr_tower_t T)
              const gr_poly_struct * gr_tower_step_minpoly(const gr_tower_t T, slong k)

    The degree `[F_k : F_{k-1}]`, the degree `[F_n : F_0]`, and the
    polynomial `m_k`. Here `1 \le k \le n`.

.. function:: int gr_tower_write(gr_stream_t out, const gr_tower_t T)
              void gr_tower_print(const gr_tower_t T)

Adjoining generators
-------------------------------------------------------------------------------

.. function:: int gr_tower_adjoin_algebraic(gr_tower_t T, const gr_poly_t m, const acb_t z, int status, const char * name)

    Adjoins the root of the monic polynomial *m* over the current top
    field `F_n` nearest to the midpoint of *z*. The root is certified,
    by an interval Newton step, to be the only root of *m* in a disk
    around the midpoint of *z* whose radius is that of *z*, inflated
    by up to `|z|/4` when *z* is tighter than the roots can be
    separated: *z* may be a rigorous enclosure or an approximation to a
    few digits (``1.4142`` selects `\sqrt 2` among the roots of
    `x^2 - 2`). Returns ``GR_UNABLE`` if no such disk is found (several
    roots near *z*, or insufficient precision), or ``GR_DOMAIN`` if *m*
    is not monic of positive degree. If *name* is ``NULL``, a default
    name is generated. *status* should be ``GR_TOWER_STATUS_PROVEN``
    only if *m* is known to be irreducible over `F_n`.

.. function:: int gr_tower_adjoin_qqbar(gr_tower_t T, const qqbar_t x, const char * name)

    Adjoins the algebraic number *x* using its minimal polynomial over
    `\mathbb{Q}`, which may factor over a nontrivial tower (the step is
    then dynamic). The base field must be `\mathbb{Q}`.

.. function:: int gr_tower_adjoin_root_ui(gr_tower_t T, gr_srcptr x, ulong n, const char * name)
              int gr_tower_adjoin_sqrt(gr_tower_t T, gr_srcptr x, const char * name)

    Adjoins the principal *n*-th root (square root) of the element *x* of
    the top field. The generator records its definition (``def_kind``
    ``GR_TOWER_ROOT``, ``def_param`` `n`, and *x* as its argument), which
    the lazy field uses to identify repeated roots.

.. function:: int gr_tower_adjoin_root_of_unity(gr_tower_t T, ulong n, const char * name)
              int gr_tower_adjoin_root_fmpz(gr_tower_t T, const fmpz_t p, ulong n, const char * name)

    Adjoins the root of unity `\exp(2 \pi i / n)` (``def_kind``
    ``GR_TOWER_ROOT_OF_UNITY``, ``def_param`` `n`) or the principal
    *n*-th root of the positive integer *p* (``GR_TOWER_ROOT`` with a
    constant argument; *p* may be negative, as for `\sqrt{-163}`, the
    root of `x^n - p`), using the minimal polynomial over `\mathbb{Q}`.
    Such *structured* generators are identified by their definitions
    when towers are merged: a root of prime power order `\ell^e` is a
    power of a generator of the same kind of order `\ell^f`, `f \ge e`,
    and roots of composite orders are products of roots of prime power
    orders (`\zeta_n = \prod \zeta_{q_i}^{c_i}` with
    `c_i = (n/q_i)^{-1} \bmod q_i`, and likewise for radicals up to a
    power of *p*). Roots of prime power orders of different primes
    generate linearly disjoint fields, so the representation stays
    canonical (irreducible moduli) and sparse. When a merge meets a
    generator of order `\ell^f` with `f < e`, a generator of order
    `\ell^e` is adjoined in front of the tower and the old one is
    re-expressed as its power (a linear modulus). Generators which end
    up after a generator inserted in front of them are defined over a
    larger field, so their proofs of irreducibility are dropped unless
    a structural rule covers the new prefix: a root of unity of prime
    power order `\ell^e` stays irreducible over roots of unity of orders
    prime to `\ell`, and principal roots of positive integers with
    pairwise coprime radicands (or equal radicands with coprime orders),
    each irreducible over `\mathbb{Q}`, generate a field of degree the
    product of the orders, also in the presence of `i` (Mordell's theorem
    on real radicals, which generalizes Besicovitch's theorem on roots of
    distinct primes); an `n`-th root of a prime `p` stays irreducible
    over roots of unity of orders prime to `p` and roots of integers
    `a` of orders `m` with `p \nmid a m` (the field is unramified at `p`,
    and `x^n - p` is Eisenstein there). The last rule also lets merges
    skip the search for such a root in the target tower whatever the
    statuses of its steps (`\sqrt{11}` over
    `\mathbb{Q}(\zeta_{16}, \sqrt{2}, \ldots, \sqrt{7})` in the
    ``examples/dft.c -tower -input 1`` benchmark, where this skip made
    the computation an order of magnitude faster). A generator whose proof was dropped is then looked
    for in the field generated by the new prefix (lattice reduction, for
    prefixes of degree at most ``GR_TOWER_EXPRESS_DEGREE_LIMIT``): for
    instance `\sqrt{2}` is irreducible over `\mathbb{Q}(i)`, but once
    `\zeta_8` has been moved in front it becomes
    `\zeta_8 - \zeta_8^3` (a linear modulus), keeping the representation
    canonical; failing that, dynamic evaluation discovers the
    factorization when a zero test needs it.
    Square roots of integers are also found through Gauss sums in the
    roots of unity of a tower: `\sqrt{p^*} = \sum_a (a/p)\, \zeta_p^a`
    for odd primes *p* (`p^* = \pm p \equiv 1 \bmod 4`),
    `\sqrt{2} = \zeta_8 + \zeta_8^{-1}` and `\sqrt{-1} = \zeta_4`. A
    merge maps `\sqrt{c}` to such an expression when the roots of
    unity are present (except in real views, which keep real
    generators), and the exact zero test turns a dynamic step
    `x^2 - c` which comes after them into a step of degree one (for
    instance `\sqrt{3}` after `\zeta_3` and `i`), instead of inverting
    in a ring which is not a field. Before the modular proofs, the zero
    test also applies two structural rules which only look at the
    values of the generators before a step, whatever their statuses:
    `x^n - p` (*p* prime, unsplit) over roots of unity and roots of
    integers generating a field unramified at *p* (Eisenstein), and
    `\Phi_{\ell^e}` over a field unramified at `\ell` (linear
    disjointness from the totally ramified `\mathbb{Q}(\zeta_{\ell^e})`;
    `i` over `\mathbb{Q}(\zeta_{111}, 3^{1/37})`, which modular proofs
    miss since primes `\equiv 1 \bmod 4` split `x^2 + 1`).

.. function:: int gr_tower_adjoin_pi(gr_tower_t T, const char * name)
              int gr_tower_adjoin_exp(gr_tower_t T, gr_srcptr u, const char * name)
              int gr_tower_adjoin_log(gr_tower_t T, gr_srcptr u, const char * name)
              int gr_tower_adjoin_exp_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
              int gr_tower_adjoin_log_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)

    Adjoins a transcendental generator `\pi`, `\exp(u)` or `\log(u)`,
    where *u* is an element of the top field, given in the nested
    representation or as a flat element in the given (current or
    superseded) context of the tower's flat machinery. The base field
    changes, so the nested contexts of the tower are rebuilt: elements
    in the nested representation held by the caller (*u* among them)
    are no longer elements of :func:`gr_tower_field`, but they remain
    valid elements of the superseded contexts, which the tower keeps
    until :func:`gr_tower_clear` (so they can be used with those
    contexts and cleared after the call), while flat elements remain
    valid elements of the tower. Returns ``GR_DOMAIN`` for
    `\exp(0)`, `\log(0)` and `\log(1)`, and ``GR_UNABLE`` if it cannot
    be decided whether *u* is such a value. No other relation is searched
    for at this point: the generator is adjoined as conjecturally
    transcendental, and relations are discovered by zero tests.

.. function:: int gr_tower_adjoin_free(gr_tower_t T, const char * name)

    Adjoins a free transcendental generator (kind ``GR_TOWER_FREE``): a
    formal variable, algebraically independent of all the generators by
    definition (status ``GR_TOWER_STATUS_INDEPENDENT``), without a
    numerical value (its enclosure is indeterminate, and
    :func:`gr_tower_get_acb` fails on elements involving it) and taking
    part in no relation search. The top field is then a rational
    function field over the previous one, with the complete arithmetic
    and zero test of the flat representation; algebraic steps over it
    need an enclosure to select their root, so they can be adjoined
    only through generators with values.

.. function:: int gr_tower_adjoin_tan(gr_tower_t T, gr_srcptr u, const char * name)
              int gr_tower_adjoin_atan(gr_tower_t T, gr_srcptr u, const char * name)
              int gr_tower_adjoin_tan_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
              int gr_tower_adjoin_atan_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)

    Adjoins the real transcendental generator `\tan(u)` or `\arctan(u)`
    (kinds ``GR_TOWER_TAN`` and ``GR_TOWER_ATAN``) for a real *u*, which
    the caller must ensure (`\tan(u)` also requires `u \not\equiv \pi/2
    \pmod \pi`); `u = 0` gives ``GR_DOMAIN``, as does a *u* which is
    numerically not real. Their enclosures are real (with an exactly zero
    imaginary part), so that the generators count as real elsewhere (the
    real forms of the lazy real field). Their relations with each other,
    with `\pi` and with the exponentials and logarithms are found by
    Richardson's algorithm on their *angles* (below).

.. function:: int gr_tower_adjoin_special(gr_tower_t T, int kind, slong param, gr_srcptr u, const char * name)
              int gr_tower_adjoin_special_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)

    Adjoins the special function value of the given kind at *u* as a
    transcendental generator: ``GR_TOWER_GAMMA`` (`\Gamma(u)`),
    ``GR_TOWER_ERF`` (`\operatorname{erf}(u)` for *param* = 0,
    `\operatorname{erfi}(u)` for *param* = 1), ``GR_TOWER_LAMBERTW``
    (`W_k(u)`, `k` = *param*), ``GR_TOWER_POLYGAMMA`` (`\psi^{(m)}(u)`,
    `m` = *param*), ``GR_TOWER_POLYLOG`` (`\operatorname{Li}_s(u)`,
    `s` = *param*), ``GR_TOWER_ZETA`` (`\zeta(u)`), ``GR_TOWER_ELLIPTIC_K``
    and ``GR_TOWER_ELLIPTIC_E`` (`K(u)`, `E(u)` with the parameter
    `m = u`), or ``GR_TOWER_CONSTANT`` (no argument; *param* =
    ``GR_TOWER_CONST_EULER`` or ``GR_TOWER_CONST_CATALAN``). The status is
    ``GR_TOWER_STATUS_CONJECTURAL`` (``GR_TOWER_STATUS_SCHANUEL`` for
    Lambert W). ``GR_TOWER_HYPGEOM`` (`{}_pF_q`) has several arguments
    and is adjoined with :func:`gr_tower_adjoin_special_multi_flat`. The caller is responsible for canonical arguments (see
    above): a value which is algebraic or expressible through other
    generators must not be adjoined. Returns ``GR_DOMAIN`` for `u = 0`
    where the value is trivial or a pole (`\operatorname{erf}`, `W`,
    `\operatorname{Li}_s`, `\Gamma`, `\psi^{(m)}`), and ``GR_UNABLE`` if
    the value cannot be evaluated numerically. The generators are
    identified across towers (merges, maps) by their kind, parameter and
    argument.

.. function:: int gr_tower_adjoin_special_multi_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_struct * u, slong nargs, const fmpz_mpoly_ctx_t mctx, const char * name)
              slong _gr_tower_special_num_args(int kind, slong param)

    Adjoins the value of a special function of several arguments: *u*
    is a vector of *nargs* flat elements, of which ``u[0]`` becomes the
    argument of the generator and the others its additional arguments.
    For ``GR_TOWER_HYPGEOM``, *param* is ``GR_TOWER_HYPGEOM_PARAM(p, q)``
    and `u = (z, a_1, \ldots, a_p, b_1, \ldots, b_q)` for the value
    `{}_pF_q(a; b; z)` (not regularized; numerically through
    :func:`acb_hypgeom_2f1` and friends, and for `p = q + 1 \ge 3` only
    inside the unit disk). *nargs* must be
    ``_gr_tower_special_num_args(kind, param)``. Generators are
    identified across towers by their kind, parameter and all their
    arguments.

Numerical evaluation
-------------------------------------------------------------------------------

.. function:: int gr_tower_trans_get_acb(acb_t res, gr_tower_t T, slong j, slong prec)

    Enclosure of the transcendental generator `t_j`.

.. function:: int gr_tower_step_get_acb(acb_t res, gr_tower_t T, slong k, slong prec)

    Sets *res* to an enclosure of the generator `a_k`, refined to
    *prec* bits of relative accuracy using interval Newton iteration.
    The enclosure is cached in the tower.

.. function:: int gr_tower_get_acb_at(acb_t res, gr_srcptr x, slong k, slong prec, gr_tower_t T)
              int gr_tower_get_acb(acb_t res, gr_srcptr x, slong prec, gr_tower_t T)

    Sets *res* to an enclosure of the element *x* of `F_k` (of `F_n`).

Dynamic refinement and complete operations
-------------------------------------------------------------------------------

.. function:: int gr_tower_refine(gr_tower_t T)

    Processes all zero divisors recorded by the quotient ring contexts,
    refining the corresponding minimal polynomials. Returns nonzero if the
    tower changed.

.. function:: int gr_tower_refine_step(gr_tower_t T, slong k, const gr_poly_t g)

    Given a monic proper factor *g* of `m_k`, replaces `m_k` by *g* or
    `m_k / g`, whichever vanishes at the enclosure of `a_k`.

.. function:: truth_t gr_tower_is_zero(gr_srcptr x, gr_tower_t T)
              truth_t gr_tower_equal(gr_srcptr x, gr_srcptr y, gr_tower_t T)
              int gr_tower_inv(gr_ptr res, gr_srcptr x, gr_tower_t T)
              int gr_tower_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_tower_t T)

    Zero test, equality test, inversion and division of elements of the
    top field. For algebraic towers these are complete: they refine the
    tower as needed and always return a definite answer (``GR_DOMAIN``
    for division by zero). An element which is nonzero as a rational
    function of conjecturally transcendental generators is decided by
    Richardson's algorithm on a copy of the tower (up to the precision
    limit, after which the answer is ``T_UNKNOWN`` / ``GR_UNABLE``).

Modular irreducibility proofs
-------------------------------------------------------------------------------

A step adjoined without a proof of irreducibility (a *dynamic* step,
:macro:`GR_TOWER_STATUS_DYNAMIC`) is correct by dynamic evaluation, but
a proof (:macro:`GR_TOWER_STATUS_PROVEN`) makes some decisions cheaper or
definite: a nonzero element of a tower whose steps are all proven is
invertible without computing the inverse, a non-rational element of such
a tower is definitely non-rational, and structural rules apply to proven
generators. For towers over `\mathbb{Q}`, proofs are attempted at
*places* of the tower: a chain of simple roots `\alpha_1, \ldots,
\alpha_k` of the moduli in a finite field `\mathbb{F}_{\ell^m}`
(each `\alpha_j` a root of `m_j(\alpha_1, \ldots, \alpha_{j-1}, x)`)
defines a homomorphism from the ring `\mathbb{Z}[a_1, \ldots, a_k]`
(with the leading coefficients inverted) to `\mathbb{F}_{\ell^m}`,
whose local ring at the kernel is a discrete valuation ring since the
Jacobian of the triangular system is nonsingular there. A monic
factorization of `m_{k+1}` over the tower would therefore reduce to a
factorization over the residue field, so `m_{k+1}` is irreducible over
the tower if its image at the place is irreducible; and a `p`-th root
in the tower of an element `a` would reduce to a `p`-th root of the
image of `a`, so `a` is not a `p`-th power in the tower if its image is
not one in the residue field. By Capelli's theorem, `X^n - a` is
irreducible over the tower if `a` is not a `p`-th power for any prime
`p \mid n` and, when `4 \mid n`, `-a/4` is not a fourth power; this
proves the steps adjoined by :func:`gr_tower_adjoin_root_ui` (the
principal `n`-th roots of the lazy field) in almost all cases where
they are irreducible. The proofs are relative to the moduli of the
lower steps, as the status of a step is.

With transcendental generators, the place also gives them values in
`\mathbb{F}_{\ell}`. A modulus then has a root at the place only for
some of the values: a radical `X^n - a(t)` for about one value of `t`
in `n`, so that over the radicals of `t` and `1 - t` of orders 4 (the
theta constants `\lambda^{1/4}`, `(1-\lambda)^{1/4}`, say) one random
choice gives a place of degree one about once in sixteen. When a
modulus has no simple root, the values of the generators which first
occur from the earliest step sharing a generator with it are drawn
again (the chain of roots being taken back to that step), up to 48
times per place. In the merges of the towers of the modular functions,
whose lattice searches for roots of radicals not in the target tower
otherwise took a third of the time, the absence of a root is then
usually shown at a place.

.. function:: int gr_tower_prove_step_modular(gr_tower_t T, slong k, slong tries)

    Tries to prove that `m_k` is irreducible over `F_{k-1}` at places of
    the tower, with an effort proportional to *tries* (a few times
    *tries* primes with places over the prime fields, then a few primes
    with places over small extension fields). Sets the status of the
    step to proven and returns 1 on success. Reducible moduli are never
    proven; irreducible ones may fail to be (a modulus which is
    reducible modulo every prime, like `x^4 + 1`, or a chain of roots
    which exists only at primes of low density). The steps below `k`
    must be proven (except for a modulus which is Eisenstein in a
    transcendental generator): over a quotient ring which is not a field,
    a place sees only one of its factors, and the irreducibility over
    that factor says nothing about the one the tower is refined to later.

.. function:: int gr_tower_prove_modular(gr_tower_t T, slong tries)

    Attempts a proof for every dynamic step of the tower which has not
    been tried with the current moduli (recorded in ``proof_version`` of
    the generator; the attempt is repeated after a refinement), in
    order, the steps above a step which remains dynamic being left for
    a later attempt. Returns 1
    if all generators are proven afterwards. The zero test calls this
    before resorting to an inversion, and the lazy field before leaving
    a rationality test undecided; the option :macro:`GR_TOWER_OPT_MODULAR_TRIES`
    is the number of primes. After a failed modular attempt on a step
    of a tower over `\mathbb{Q}`, Trager's method (below) is tried,
    which decides the question for small degrees: an irreducible `x^4 + 1` is proven, and a modulus
    which is reducible over the tower is refined.

.. function:: int gr_tower_prove_step_trager(gr_tower_t T, slong k, slong degree_limit)

    Decides the irreducibility of `m_k` over `F_{k-1}` by Trager's
    method (see the factorization below), for a tower whose steps below
    `k` are proven, over `\mathbb{Q}` or over a rational function field
    `\mathbb{Q}(t_1, \ldots, t_r)`, provided that `\deg(m_k) [F_{k-1} :
    F_0] \le` *degree_limit*. If `m_k` is irreducible, the step is marked
    proven; otherwise the tower is refined with its irreducible factor
    vanishing at the generator (the factors which do not vanish being
    divided out one at a time), and the refined step is marked proven.
    Linear steps are proven directly. Returns 1 if the step is proven
    (after a refinement, if any), 0 if the method does not apply, the
    degree limit is exceeded, the modulus is not squarefree, or no
    shift making the norm squarefree was found. Unlike the modular
    proofs this is a decision procedure, but the norm has degree
    `\deg(m_k) [F_{k-1} : F_0]` and is factored over `\mathbb{Z}` (or
    as a multivariate polynomial over `\mathbb{Z}`), so it is limited to
    small absolute degrees: the automatic proofs of
    :func:`gr_tower_prove_modular` use it for towers over `\mathbb{Q}`
    only, with the limit of the option
    :macro:`GR_TOWER_OPT_TRAGER_DEGREE_LIMIT` (48 by default).

.. function:: int gr_tower_binomial_irreducible_modular(gr_tower_t T, const fmpz_mpoly_q_t a, const fmpz_mpoly_ctx_t actx, ulong n, slong tries)

    Proves `X^n - a` irreducible over the field of the tower by the
    criterion above, for the element *a* in the flat representation
    (in the context *actx* of the tower's flat machinery). Returns 1 on
    success; the steps of the tower must be proven.

    (Internally, the same test reports a prime `p \mid n` for which `a`
    was a `p`-th power at many informative places, the evidence used
    by :macro:`GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT`.)

.. function:: int gr_tower_poly_no_roots_modular(const gr_poly_t g, gr_tower_t T, slong tries)

    Returns 1 if the polynomial *g* over the top field of *T* is shown to
    have no root in that field, and 0 if this is not known. About *tries*
    places over prime fields are tried (the transcendental generators
    being given random values): if *g* had a root `r`, then at a place
    where the coefficients of *g* are integral and its leading
    coefficient does not vanish, `\phi(r)` would be a root of the
    reduction of *g* (the local ring at a place where the moduli have
    simple roots is regular, hence integrally closed). The statement is
    relative to the moduli of the steps being irreducible; the function
    serves to skip lattice searches for roots (when towers are merged and
    in root finding), where a wrong answer only costs a dynamic step.

The effort of the automatic uses is set by the options
:macro:`GR_TOWER_OPT_NO_ROOTS_TRIES` (places tried by
:func:`gr_tower_poly_no_roots_modular`), :macro:`GR_TOWER_OPT_MODULAR_TRIES`
(the budget of a modular proof: a schedule of places over prime fields,
then over extension fields, of a size proportional to it) and
:macro:`GR_TOWER_OPT_TRAGER_DEGREE_LIMIT` (the largest absolute degree
`\deg(m_k) [F_{k-1} : \mathbb{Q}]` for which the automatic proofs of
towers over `\mathbb{Q}` fall back to Trager's method).

.. macro:: GR_TOWER_DEFERRED_IDEAL_SIZE

    Number of base coefficients above which the flat conversion of a
    modulus is deferred (see the flat representation below).

Factorization of polynomials
-------------------------------------------------------------------------------

Polynomials over a field `F_k` of a tower are factored by Trager's
method. The *norm* `N(h) \in F_0[x]` of `h \in F_k[x]` is the product of
the conjugates of `h` over `F_0`, computed as the iterated resultant of
`h` with the moduli `m_k, \ldots, m_1` (the coefficients of `h` viewed
as polynomials in `a_k` over `F_{k-1}[x]`, and so on down); it has
degree `\deg(h) [F_k : F_0]`. For squarefree `h`, the norm of
`h(x - \theta)` is squarefree for all but finitely many shifts
`\theta = \sum_j s_j a_j`; for such a shift, the irreducible factors of
`h` are `\gcd(h(x - \theta), N_i)(x + \theta)` for the irreducible
factors `N_i` of the norm over `F_0`. Over `F_0 = \mathbb{Q}` the norm
is factored by :func:`fmpz_poly_factor`; over `F_0 = \mathbb{Q}(t_1,
\ldots, t_r)` its denominators are cleared and it is factored as a
polynomial in `t_1, \ldots, t_r, x` by :func:`fmpz_mpoly_factor`.
Squarefree decomposition is by Yun's algorithm, and the gcds are
computed by the subresultant PRS. The factorization is exact: it
requires the steps up to `k` to be proven (which is attempted first:
modular proofs, then Trager's method for the steps themselves), so
that `F_k` is a field.

The option :macro:`GR_TOWER_OPT_FACTOR_DEGREE_LIMIT` bounds `\deg(h)
[F_k : F_0]` for the squarefree parts `h` factored by
:func:`gr_tower_poly_factor`.

.. function:: int gr_tower_poly_factor(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, slong k, gr_tower_t T)
              int gr_tower_poly_factor_limit(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, slong k, slong degree_limit, gr_tower_t T)

    Factors the nonzero polynomial *f* over `F_k` as `c \prod_i
    f_i^{e_i}` with `c \in F_k` (the leading coefficient) and distinct
    monic irreducible `f_i`, stored in *fac* (a vector of polynomials
    over `F_k`) and *mult*. Returns ``GR_UNABLE`` if a step up to `k`
    cannot be proven irreducible, if a squarefree part exceeds the
    degree limit, if no suitable shift is found, or if the base field is
    not `\mathbb{Q}` or `\mathbb{Q}(t_1, \ldots, t_r)`; ``GR_DOMAIN``
    for the zero polynomial. The field contexts
    (:func:`gr_ctx_init_tower_field`) implement :func:`gr_factor` for
    polynomials (``GR_METHOD_POLY_FACTOR``) by this function.

.. function:: int gr_tower_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t f, slong k, gr_tower_t T)

    The roots of *f* in `F_k` with their multiplicities (the linear
    factors of the factorization). The field contexts implement
    :func:`gr_poly_roots` by this function.

Absolute representation
-------------------------------------------------------------------------------

.. function:: int gr_tower_get_coeffs_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T)
              int gr_tower_set_coeffs_at(gr_ptr res, gr_srcptr coeffs, slong k, gr_tower_t T)

    Conversion between an element of `F_k` and its vector of `[F_k : F_0]`
    coordinates in the monomial basis
    `a_1^{i_1} \cdots a_k^{i_k}`, ordered with `i_1` varying fastest.

.. function:: int gr_tower_multiplication_matrix(gr_mat_t res, gr_srcptr x, gr_tower_t T)
              int gr_tower_charpoly(gr_poly_t res, gr_srcptr x, gr_tower_t T)

    The matrix of multiplication by *x* with respect to the monomial
    basis of `F_n` over `F_0`, and its characteristic polynomial.

.. function:: int gr_tower_get_fmpz_poly_minpoly(fmpz_poly_t res, gr_srcptr x, gr_tower_t T)
              int gr_tower_get_qqbar(qqbar_t res, gr_srcptr x, gr_tower_t T)

    The minimal polynomial of the number represented by *x* over
    `\mathbb{Q}`, respectively the number itself as a :type:`qqbar_t`.
    The base field must be `\mathbb{Q}`. The annihilating polynomial of
    least degree is computed by incremental elimination on the powers of
    *x*; if the tower is not yet fully refined this polynomial may be
    reducible, in which case the correct factor is selected numerically.

Expressing numbers in a tower
-------------------------------------------------------------------------------

.. function:: int gr_tower_sqrt_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T)
              int gr_tower_sqrt(gr_ptr res, gr_srcptr x, gr_tower_t T)

    Sets *res* to the principal square root of *x* if it lies in `F_k`
    (`F_n`): with `F_k = F_{k-1}(a)`, `a^2 = r` after completing the
    square, `x = u + v a` is a square iff `u^2 - r v^2` is a square `s` in
    `F_{k-1}` and `(u \pm s)/2` is a square `p^2` there, and then
    `\sqrt{x} = p + v a / (2p)`; this recurses to the base field, where
    rational numbers and rational functions are tested directly. No
    lattice reduction is involved. Returns ``GR_DOMAIN`` if *x* is not a
    square in `F_k`, and ``GR_UNABLE`` if a step below `k` has degree
    other than 2.

Above the degree :macro:`GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT` (an option),
lattice-based expression searches are skipped when adjoining roots and
merging towers; roots are then adjoined without searching, dynamic
evaluation keeping zero tests complete.

.. function:: int gr_tower_express(gr_ptr res, const gr_poly_t q, const acb_t z, gr_tower_t T)
              int gr_tower_express_limit(gr_ptr res, const gr_poly_t q, const acb_t z, slong prec_limit, gr_tower_t T)

    Given a polynomial *q* over the top field and an enclosure *z*
    isolating a single root of *q*, attempts to find the element of the
    top field which is that root. A candidate is found by an integer
    relation search (LLL) over the monomial basis, then verified exactly
    (`q(x) = 0`, using the complete zero test) and identified numerically
    with the root in *z*. Returns ``GR_SUCCESS`` if found and ``GR_UNABLE``
    if no representation was found within the precision limit (which
    does not prove that the number is not in the field).

.. function:: int gr_tower_express_qqbar(gr_ptr res, const qqbar_t x, gr_tower_t T)

    Expresses the algebraic number *x* in the tower (base field `\mathbb{Q}`).
    Returns ``GR_DOMAIN`` if the degree of *x* does not divide the degree
    of the tower, which proves that *x* is not in the field.

Maps between towers
-------------------------------------------------------------------------------

.. type:: gr_tower_map_struct

.. type:: gr_tower_map_t

    A homomorphism from a source tower to a target tower over the same
    constant field, given by the images of the first *length* generators
    of the source (in definition order) as flat elements of the target.
    Flat elements survive extensions of the target, so a map remains valid
    when the target grows (its images may need conversion to the current
    polynomial context, see :func:`gr_tower_map_sync`).

.. function:: void gr_tower_map_init(gr_tower_map_t map, gr_tower_t source, gr_tower_t target)
              void gr_tower_map_clear(gr_tower_map_t map)
              void gr_tower_map_set(gr_tower_map_t res, const gr_tower_map_t map)
              int gr_tower_map_fit_length(gr_tower_map_t map, slong len)
              void gr_tower_map_sync(gr_tower_map_t map)
              fmpz_mpoly_q_struct * gr_tower_map_image(gr_tower_map_t map, slong d)
              int gr_tower_map_set_image(gr_tower_map_t map, slong d, gr_srcptr x, slong k)
              int gr_tower_map_set_image_flat(gr_tower_map_t map, slong d, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t mctx)

    A newly initialized map is defined on no generators; images are set
    by the functions below or by *set_image* (from an element of `F_k`
    of the target) and *set_image_flat*. The image of the generator with
    definition order *d* is an element of the context ``map->mctx``,
    which *sync* brings to the current context of the target.

.. function:: int gr_tower_map_set_inclusion(gr_tower_map_t map)

    Sets the map to the inclusion of the source into a target whose first
    generators (in definition order) are copies of the source's generators.

.. function:: int gr_tower_map_apply_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_map_t map)
              int gr_tower_map_apply(gr_ptr res, gr_srcptr x, gr_tower_map_t map)
              int gr_tower_map_apply_flat(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, gr_tower_map_t map)
              int gr_tower_map_apply_poly(gr_poly_t res, const gr_poly_t f, gr_tower_map_t map)

    Applies the map to an element of `F_k` of the source, of its top
    field, to a flat element of the source, or to a polynomial over the
    top field of the source. The result is an element of (a polynomial
    over) the top field of the target, respectively a flat element in
    the context ``map->mctx`` (after the synchronization performed by
    this function). Returns ``GR_DOMAIN`` if *x* involves generators on
    which the map is not defined.

.. function:: int gr_tower_absorb(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, int flags)
              int gr_tower_absorb_prefix(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, slong p, int flags)

    Extends *U* (the target of *map*) by the generators of *B* with
    definition order from *length* (of the map) to *p* - 1 (or all), and
    extends the map accordingly. A generator with the same definition id
    as a generator of *U* is identified with it. Otherwise, an algebraic
    generator whose minimal polynomial has become linear is substituted;
    a minimal polynomial whose image has coefficients in the base field
    of *U* is factored there (the factor vanishing at the generator
    replaces it); if *flags* contains ``GR_TOWER_MERGE_EXPRESS``, a
    low-effort expression search is tried, unless a place of *U* shows
    that the polynomial has no root in *U*
    (:func:`gr_tower_poly_no_roots_modular`); otherwise the generator is
    adjoined (as a dynamic step). A transcendental generator is adjoined
    with the image of its argument.

.. function:: int gr_tower_merge(gr_tower_t U, gr_tower_map_t mapA, gr_tower_map_t mapB, gr_tower_t A, gr_tower_t B, int flags)

    Sets *U* (which must be initialized) to a tower containing *A* and *B*,
    with maps (initialized by this function) from *A* and *B* into *U*.

.. function:: int gr_tower_eliminate(gr_tower_t U, gr_tower_map_t map, gr_tower_t T)

    Sets *U* to a tower isomorphic to *T* with the steps of degree 1
    removed, with the isomorphism as *map*.

Tower field context
-------------------------------------------------------------------------------

A fixed tower is built step by step and its top field used as a ``gr``
context; a typical use is a field given in advance, in which many
computations take place (linear algebra over `\mathbb{Q}(\sqrt 2,
\sqrt 3)`, say), where the lazy field's bookkeeping of towers is not
wanted:

.. code-block:: c

    gr_ctx_t QQ, K;
    gr_tower_t T;
    gr_ptr x;

    gr_ctx_init_fmpq(QQ);
    gr_tower_init(T, QQ);

    /* adjoin sqrt(2) and sqrt(3): the arguments are elements of the
       current top field */
    {
        gr_ctx_struct * top = gr_tower_field(T);
        GR_TMP_INIT(x, top);
        GR_MUST_SUCCEED(gr_set_ui(x, 2, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
        GR_TMP_CLEAR(x, top);
    }
    {
        gr_ctx_struct * top = gr_tower_field(T);
        GR_TMP_INIT(x, top);
        GR_MUST_SUCCEED(gr_set_ui(x, 3, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
        GR_TMP_CLEAR(x, top);
    }

    /* the field Q(sqrt 2, sqrt 3) as a gr context: complete zero tests,
       division, generators a1, a2 */
    gr_ctx_init_tower_field(K, T);
    {
        gr_vec_t gens;
        gr_vec_init(gens, 0, K);
        GR_MUST_SUCCEED(gr_gens(gens, K));
        GR_TMP_INIT(x, K);
        GR_MUST_SUCCEED(gr_add(x, gr_vec_entry_ptr(gens, 0, K), gr_vec_entry_ptr(gens, 1, K), K));
        GR_MUST_SUCCEED(gr_inv(x, x, K));     /* 1/(sqrt 2 + sqrt 3) = sqrt 3 - sqrt 2 */
        gr_println(x, K);
        GR_TMP_CLEAR(x, K);
        gr_vec_clear(gens, K);
    }
    gr_ctx_clear(K);
    gr_tower_clear(T);

Adjoining a step whose polynomial is only known to vanish at the root
(:func:`gr_tower_adjoin_algebraic` with ``GR_TOWER_STATUS_DYNAMIC``) is
allowed: the field context stays a field, and the tower is refined when
an operation exposes a factorization. Transcendental generators
(:func:`gr_tower_adjoin_pi`, :func:`gr_tower_adjoin_exp`,
:func:`gr_tower_adjoin_log`) rebuild the tower over a rational function
field, which invalidates the contexts created before (see below). The
Python wrapper's ``gr_tower`` class follows this pattern with
``adjoin_*`` methods and ``field()``.

.. function:: void gr_ctx_init_tower_field(gr_ctx_t ctx, gr_tower_t T)
              int gr_tower_field_get_acb(acb_t res, gr_srcptr x, slong prec, gr_ctx_t ctx)

    Initializes *ctx* to a ``gr`` context for the field `F_n`, `n` being
    the number of algebraic steps of *T* at that moment. Elements
    have the same representation as elements of :func:`gr_tower_field`,
    and arithmetic is delegated to that context, but zero tests, equality
    tests and inversions are complete (dynamically refining *T*).
    The context refers to *T*, which must outlive it. It remains valid
    when further algebraic steps are appended to *T* afterwards: its
    elements then lie in the subfield `F_n` of the new top field, and
    :func:`gr_set_other` converts them into a context made for the
    larger field (and, generally, between contexts of the same tower
    from a lower to a higher level). Adjoining a transcendental
    generator rebuilds the tower over a new base field, which makes the
    contexts created before *stale*: their operations fail with
    ``GR_UNABLE``, while clearing their elements remains valid (the
    tower keeps the superseded nested contexts until it is cleared). Conversion to an ``acb``
    context (:func:`gr_set_other`) is supported.

Flat representation
-------------------------------------------------------------------------------

The nested representation is the structural one (it carries the
triangular set, dynamic refinement and the certificates), but for
arithmetic it is often better to represent elements as multivariate
rational functions in the generators, reduced modulo the triangular set
with denominators cleared. Since the leading monomials `c_k x_k^{d_k}` of
the triangular set are pairwise coprime, it is a Gröbner basis for the
lex order with `x_n > \cdots > x_1`, so reduced numerators are canonical
up to scaling. In this representation inversion is free (swap numerator
and denominator), sparse expressions stay sparse and multiplication uses
:type:`fmpz_mpoly` arithmetic, at the cost of non-canonical denominators
(equality tests cross-multiply). Both representations coexist; conversion
in either direction goes through the coordinate vector.

The flat *ideal* (the moduli converted to polynomials in all the
generators) is built from the nested moduli. Over `\mathbb{Q}`, the
conversion of a modulus with non-constant coefficients and more than
:macro:`GR_TOWER_DEFERRED_IDEAL_SIZE` base coefficients is deferred until
a reduction actually needs it: the reduction of an element then proceeds
by one division per step from the top down (a reduction by the modulus
of step `k` raises only degrees in the variables of lower steps, so the
steps needed are known as the reduction goes), converting the moduli on
first use. The moduli of a splitting tower have binomial size in the
middle, and symmetric functions of the roots -- which reduce through the
last, linear steps -- never touch them. A polynomial whose degree in
every variable is already below that of the corresponding modulus is
reduced (the reduction by a modulus of step `k` raises only degrees in
the variables of lower steps, never that of `a_k` above `\deg m_k`),
so :func:`gr_tower_flat_reduce` returns at once for such input; this
is the common case after a multiplication of elements which are linear
in the generators, and avoids a division per step.

.. type:: gr_tower_flat_struct

.. type:: gr_tower_flat_t

    The flat machinery for one tower: a polynomial context in which the
    variable ``GR_TOWER_FLAT_VAR_D(F, d)`` = ``cap - 1 - d`` represents the
    generator with definition order *d* (``GR_TOWER_FLAT_VAR(F, k)`` and
    ``GR_TOWER_FLAT_TVAR(F, j)`` are the variables of `a_k` and `t_j`),
    so that later generators are more significant in lex order and the
    triangular set is a Groebner basis whatever the kinds of the
    generators are; the tower can grow up to the capacity without
    renumbering. The reduced triangular set, and superseded polynomial
    contexts (kept alive, each with the list of generators its variables
    stood for, so that elements created with them remain convertible even
    after reorderings). Every tower has an embedded flat machinery,
    :func:`gr_tower_flat`. Reduction is by pseudo-division from the top
    generator down; the leading coefficients of the moduli are polynomials
    in the transcendental variables only, so the fraction is multiplied
    through by them.

.. function:: void gr_tower_flat_init(gr_tower_flat_t F, gr_tower_t T, slong cap)
              void gr_tower_flat_clear(gr_tower_flat_t F)
              int gr_tower_flat_ensure(gr_tower_flat_t F)
              gr_tower_flat_struct * gr_tower_flat(gr_tower_t T)

    Initialization and maintenance: *ensure* replaces the polynomial
    context if the tower has grown beyond the capacity or has been
    reordered (returning 1) and rebuilds the triangular set if the tower
    has been refined or restructured.

.. function:: void gr_tower_flat_convert(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t old_mctx, gr_tower_flat_t F)
              int gr_tower_flat_reduce(fmpz_mpoly_q_t x, gr_tower_flat_t F)
              int gr_tower_flat_set_nested_at(fmpz_mpoly_q_t res, gr_srcptr x, slong k, gr_tower_flat_t F)
              int gr_tower_flat_poly_get_nested_at(gr_ptr res, const fmpz_mpoly_t f, slong k, gr_tower_flat_t F)
              int gr_tower_flat_get_nested_at(gr_ptr res, const fmpz_mpoly_q_t x, slong k, gr_tower_flat_t F)
              int gr_tower_flat_get_acb(acb_t res, const fmpz_mpoly_q_t x, slong prec, gr_tower_flat_t F)
              slong gr_tower_flat_level(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
              slong gr_tower_flat_alg_level(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
              slong gr_tower_flat_def_order(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
              int gr_tower_flat_compose(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, fmpz_mpoly_q_struct ** imgs, gr_tower_flat_t F)
              truth_t gr_tower_flat_num_is_zero(fmpz_mpoly_q_t x, gr_tower_flat_t F)

    Conversion from a superseded context, reduction, conversions to and
    from the nested representation of `F_k`, numerical evaluation, the
    length of the shortest prefix containing all generators occurring in
    *x* (respectively the highest algebraic generator, and the highest
    definition order), substitution of elements of *F*'s tower for the
    variables of another context, and the complete zero test of the
    numerator (numerics, then the nested test, which may refine the tower).

.. function:: void gr_ctx_init_tower_field_flat(gr_ctx_t ctx, gr_tower_t T)
              int gr_tower_flat_set_nested(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_tower_flat_get_nested(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_tower_field_flat_get_acb(acb_t res, gr_srcptr x, slong prec, gr_ctx_t ctx)

    A ``gr`` context for the top field of a fixed tower in the flat
    representation, presented as an extension of the constant field with
    generators `t_1, \ldots, t_r, a_1, \ldots, a_n`, conversions to and
    from :func:`gr_ctx_init_tower_field`, and numerical evaluation.

Tuning options
-------------------------------------------------------------------------------

The limits and policies of the algorithms are read from a table of
options (indexed by the constants below) that a tower refers to through
its member ``options``: :var:`gr_tower_default_options` for a tower made
by :func:`gr_tower_init`, the table of the lazy field for its towers
(see :func:`gr_tower_lazy_ctx_set_option`). Towers derived from a tower
(prefixes, copies, merges) share its table. Changing an option affects
the computations that follow; the elements created before remain valid.

.. var:: const slong gr_tower_default_options[GR_TOWER_OPT_NUM_OPTIONS]

.. function:: void gr_tower_set_options(gr_tower_t T, const slong * options)

    Makes *T* read its options from the table *options* (not copied), or
    from the defaults if *options* is ``NULL``.

.. macro:: GR_TOWER_OPTION(T, k)

    The value of the option *k* for the tower *T*.

.. function:: const char * gr_tower_option_name(slong option)
              slong gr_tower_option_find(const char * name)
              slong gr_tower_option_default(slong option)
              int gr_tower_option_valid(slong option, slong value)
              void gr_tower_options_init(slong * options)
              int gr_tower_options_set(slong * options, slong option, slong value)

    The name of an option (``"prec_limit"`` for
    ``GR_TOWER_OPT_PREC_LIMIT``; ``NULL`` for an invalid index), the
    index of the option with the given name (-1 if none), the default
    value, and whether a value is in the range of the option. A table
    for :func:`gr_tower_set_options` (``GR_TOWER_OPT_NUM_OPTIONS``
    entries) is filled with the defaults by :func:`gr_tower_options_init`
    and changed by :func:`gr_tower_options_set`, which returns
    ``GR_DOMAIN`` for a value out of range. (The Python interface takes
    the options by their names.)

.. macro:: GR_TOWER_OPT_VERBOSE
           GR_TOWER_OPT_PRINT_FLAGS
           GR_TOWER_OPT_PRINT_DIGITS
           GR_TOWER_OPT_PREC_LIMIT
           GR_TOWER_OPT_CERTIFY_PREC_LIMIT
           GR_TOWER_OPT_NUMERIC_PREC_LIMIT
           GR_TOWER_OPT_SMOOTH_LIMIT
           GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT
           GR_TOWER_OPT_EXPRESS_PREC
           GR_TOWER_OPT_TRAGER_DEGREE_LIMIT
           GR_TOWER_OPT_FACTOR_DEGREE_LIMIT
           GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT
           GR_TOWER_OPT_MODULAR_TRIES
           GR_TOWER_OPT_MODULAR_TERMS_LIMIT
           GR_TOWER_OPT_NO_ROOTS_TRIES
           GR_TOWER_OPT_RATIONALIZE_LIMIT
           GR_TOWER_OPT_RELATION_COST_LIMIT
           GR_TOWER_OPT_SQRT_BUDGET
           GR_TOWER_OPT_GAUSS_SUM_LIMIT
           GR_TOWER_OPT_DENSE_LIMIT
           GR_TOWER_OPT_INV_DENSE_DEGREE_LIMIT
           GR_TOWER_OPT_INV_LINEAR_DEGREE_LIMIT
           GR_TOWER_OPT_DEFERRED_IDEAL_SIZE
           GR_TOWER_OPT_GROW_THRESHOLD
           GR_TOWER_OPT_ROOT_OF_UNITY_ORDER_LIMIT
           GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT
           GR_TOWER_OPT_REALIFY_ORDER_LIMIT
           GR_TOWER_OPT_TRIG_PI_LIMIT
           GR_TOWER_OPT_TRIG_ALGEBRAIC_LIMIT
           GR_TOWER_OPT_TRIG_FORM
           GR_TOWER_OPT_SPLIT_IMAGINARY
           GR_TOWER_OPT_COMPOSITE_RADICALS
           GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT
           GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT
           GR_TOWER_OPT_ANNIHILATING_FLAT_LIMIT
           GR_TOWER_OPT_GAMMA_LATTICE_LIMIT
           GR_TOWER_OPT_GAMMA_LINE_LIMIT
           GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT
           GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT
           GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT
           GR_TOWER_OPT_GAUSS_DIGAMMA_LIMIT
           GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT
           GR_TOWER_OPT_DENSE_FORM_SPARSITY
           GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT
           GR_TOWER_OPT_SPLIT_DEGREE_LIMIT
           GR_TOWER_OPT_MINPOLY_DEGREE_LIMIT
           GR_TOWER_OPT_INV_DENSE_ALG
           GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT

The options (default values in parentheses; precisions are in bits):

* :macro:`GR_TOWER_OPT_VERBOSE` (0): with 1, the relations found
  between generators (generators becoming of degree one) are printed;
  with 2, also statistics.
* :macro:`GR_TOWER_OPT_PRINT_FLAGS` (``GR_TOWER_PRINT_SYMBOLIC |
  GR_TOWER_PRINT_DEFS``), :macro:`GR_TOWER_OPT_PRINT_DIGITS` (6): printing
  of lazy field elements (:func:`gr_tower_lazy_ctx_set_print`).
* :macro:`GR_TOWER_OPT_PREC_LIMIT` (256): precision up to which a cheap
  numerical separation from zero is attempted before an exact test
  (zero tests, comparisons, signs).
* :macro:`GR_TOWER_OPT_CERTIFY_PREC_LIMIT` (4096): precision limit of
  the certifications (roots, signs and real parts of generators, the
  selection among candidates by their enclosures) and of the relation
  searches of the zero test, after which the answer is unknown.
* :macro:`GR_TOWER_OPT_NUMERIC_PREC_LIMIT` (65536): when the exact
  methods leave a zero test or a sign undecided (no structure theorem
  covers the element), the numerical evaluation is pushed to this
  precision before giving up: `1 - \operatorname{erf}(100)` is
  separated from zero at 14500 bits.
* :macro:`GR_TOWER_OPT_SMOOTH_LIMIT` (3512): number of primes for the
  trial division of the radicands of roots in lazy fields, and of the
  Gaussian integers whose logarithms are decomposed over the prime
  logarithms; a cofactor beyond it is kept as one radicand (or one
  opaque logarithm).
* :macro:`GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT` (64),
  :macro:`GR_TOWER_OPT_EXPRESS_PREC` (256): lattice-based expression
  searches in fields up to this degree, at a precision of this many
  bits plus four per degree of the field
  (:macro:`GR_TOWER_MERGE_EXPRESS_PREC`).
* :macro:`GR_TOWER_OPT_TRAGER_DEGREE_LIMIT` (48),
  :macro:`GR_TOWER_OPT_MODULAR_TRIES` (6),
  :macro:`GR_TOWER_OPT_MODULAR_TERMS_LIMIT` (50000),
  :macro:`GR_TOWER_OPT_NO_ROOTS_TRIES` (24): effort of the automatic
  irreducibility proofs (Trager's method up to this degree; the number
  of places tried by the modular proofs, which are not attempted for
  moduli of more than this many terms; the places tried by the modular
  proofs of the absence of roots).
* :macro:`GR_TOWER_OPT_FACTOR_DEGREE_LIMIT` (128): factorization over
  towers.
* :macro:`GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT` (64): in lazy fields,
  the roots of a polynomial are found by factoring it over the tower
  of its coefficients when the degree of the tower times that of the
  polynomial is at most this (3/8 of it with transcendental generators,
  over which factoring is costlier).
* :macro:`GR_TOWER_OPT_RATIONALIZE_LIMIT` (2000): size of the
  denominators rationalized.
* :macro:`GR_TOWER_OPT_RELATION_COST_LIMIT` (200000): multiplicative
  relations costlier to verify are not verified.
* :macro:`GR_TOWER_OPT_SQRT_BUDGET` (2000): recursive calls of the
  square root search in towers of quadratic steps (the search branches,
  up to three calls one level down; 200 was too few for the ten
  quadratic steps behind cos(pi/257) as nested square roots).
* :macro:`GR_TOWER_OPT_GAUSS_SUM_LIMIT` (10000): square roots of primes
  up to this bound expressed as Gauss sums in the cyclotomic fields
  present.
* :macro:`GR_TOWER_OPT_DENSE_LIMIT` (100000): size of the Kronecker arrays
  of dense products.
* :macro:`GR_TOWER_OPT_INV_DENSE_DEGREE_LIMIT` (512),
  :macro:`GR_TOWER_OPT_INV_LINEAR_DEGREE_LIMIT` (64): inverses by the
  dense multivariate method in fields up to the first degree, and (in
  lazy fields) by a linear system up to the second.
  :macro:`GR_TOWER_OPT_INV_DENSE_ALG` (0) selects the dense inverse
  algorithm in fields of several generators for testing: 0 the modular
  one, then linear algebra; 1 the modular one only; 2 linear algebra
  only.
* :macro:`GR_TOWER_OPT_DEFERRED_IDEAL_SIZE` (12): moduli with more
  coefficients than this over the base field are converted to the flat
  representation only when needed.
* :macro:`GR_TOWER_OPT_GROW_THRESHOLD` (6): number of generators beyond
  which a lazy field extends towers in place instead of merging them.
* :macro:`GR_TOWER_OPT_ROOT_OF_UNITY_ORDER_LIMIT` (360): logarithms of
  roots of unity of orders up to this are rational multiples of
  `\pi i`.
* :macro:`GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT` (100000): roots of unity
  of orders up to this may become generators (one generator for the
  orders present, see :macro:`GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT`).
* :macro:`GR_TOWER_OPT_REALIFY_ORDER_LIMIT` (480): real elements written
  with roots of unity of orders up to this (and whose common level is
  at most this) are rewritten in real generators in real fields.
* :macro:`GR_TOWER_OPT_TRIG_PI_LIMIT` (240): sine, cosine and tangent of
  `r \pi` for denominators of `r` up to this in the tangent normal form.
* :macro:`GR_TOWER_OPT_TRIG_ALGEBRAIC_LIMIT` (1000): beyond the
  tangent normal form, the trigonometric values at `r \pi` are
  algebraic numbers (roots of unity, or real algebraic numbers) for
  denominators of `r` up to this; beyond it, they are built from
  exponentials.
* :macro:`GR_TOWER_OPT_TRIG_FORM` (:macro:`GR_TOWER_TRIG_EXPONENTIAL`):
  with :macro:`GR_TOWER_TRIG_TANGENT`, sine, cosine and tangent of real
  arguments in complex lazy fields take the real forms (through
  `\tan(x/2)`, and the tangent normal form at rational multiples of
  `\pi`), as in real fields, instead of exponentials.
* :macro:`GR_TOWER_OPT_SPLIT_IMAGINARY` (0),
  :macro:`GR_TOWER_OPT_COMPOSITE_RADICALS` (0): the choice of generators
  for square roots (:func:`gr_tower_lazy_ctx_set_gen_flags`).
* :macro:`GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT` (0): the largest degree
  `\varphi(N)` of a root of unity `\zeta_N` of composite order used as
  one generator for the roots of unity whose orders divide `N`; beyond
  it (and with 0), roots of unity of prime power orders are used. One
  field is at best about twice as fast for dense arithmetic within it
  (sums of products of random elements of `\mathbb{Q}(\zeta_{60})` and
  `\mathbb{Q}(\zeta_{120})`), and slower for orders with several odd
  prime factors (`\mathbb{Q}(\zeta_{105})`, `\mathbb{Q}(\zeta_{255})`:
  twice as slow), and computations involving many orders build large
  fields (the DFT benchmark with `x_k = (k+2)^{1000}`, `N = 105` is
  an order of magnitude slower with a cap of 64 or 128), whence the
  default.
* :macro:`GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT` (120): when towers
  are merged, conjugate roots of one integer polynomial of degree up
  to this are recognized by it (and build a splitting tower with the
  Cauchy moduli).
* :macro:`GR_TOWER_OPT_ANNIHILATING_FLAT_LIMIT` (1024): annihilating
  polynomials of elements of towers of degree up to this are found by
  linear algebra over the flat representation.
* :macro:`GR_TOWER_OPT_GAMMA_LATTICE_LIMIT` (36),
  :macro:`GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT` (240): normal forms of
  `\Gamma` and Hurwitz zeta values at rationals up to these denominators.
* :macro:`GR_TOWER_OPT_GAMMA_LINE_LIMIT` (12): the relations between
  `\Gamma` values on a rational line `a w + b` are searched for
  multipliers `|a|` up to this.
* :macro:`GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT` (64): Hurwitz zeta values
  and polygamma functions up to this weight take part in the relation
  searches and the normal forms.
* :macro:`GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT` (240): levels of the
  relations between `\Gamma` and Hurwitz zeta generators searched by the
  zero test.
* :macro:`GR_TOWER_OPT_GAUSS_DIGAMMA_LIMIT` (30): Gauss's digamma theorem
  for denominators up to this.
* :macro:`GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT` (1048576): the dense
  forms of elements of number fields of one generator up to this degree
  (0: none), see the lazy field below.
* :macro:`GR_TOWER_OPT_DENSE_FORM_SPARSITY` (8): an element given by a
  polynomial in the generator longer than 256 takes the dense form only
  if at least one coefficient in this many is nonzero (0: whatever its
  sparsity), see the lazy field below.
* :macro:`GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT` (0): in lazy fields, a
  merged tower of several proven steps over `\mathbb{Q}` with monic
  integral univariate moduli (radicals, generic algebraic numbers; no
  roots of unity, tangents or transcendental generators) of degree at
  most this is replaced by `\mathbb{Q}(\theta)` for a primitive element
  `\theta = \ell \sum c_i a_i` (small integers `c_i` making
  `1, \theta, \ldots, \theta^{D-1}` independent, `\ell` the least
  integer making `\theta` integral), the definitions of the tower being
  recorded as polynomials in `\theta`; its elements then have the dense
  form. It loses the structure of the steps (sparse elements become
  dense with large coefficients), so it pays only for small degrees:
  sums of products of random elements are 1.7 times faster in
  `\mathbb{Q}(\sqrt{2}, \sqrt{3})`, but 7 times slower in
  `\mathbb{Q}(\sqrt{2}, \sqrt{3}, \sqrt{5})`, 30 times slower with
  `\sqrt{7}` added, and fields of degree 20 later merged with
  transcendental generators make zero tests by inversion very slow.
* :macro:`GR_TOWER_OPT_SPLIT_DEGREE_LIMIT` (4096),
  :macro:`GR_TOWER_OPT_MINPOLY_DEGREE_LIMIT` (512): in lazy fields, the
  zero test of an algebraic number `f + g` in a tower of nominal degree
  beyond the first limit, not all proven, with `f` and `g` involving
  disjoint sets of generators (the difference of two numbers built
  independently, whose towers a merge has stacked), decides `f = -g`
  in towers of their own generators: by `P(-g) = 0` for the minimal
  polynomial `P` of `f` (or the other way round; the smaller tower, of
  degree up to the second limit), a zero test in the tower of `-g`
  alone, and the numerical isolation of the roots of `P`. With 0, never.
* :macro:`GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT` (64): in lazy fields,
  a root `\sqrt[n]{x}` which none of the cheaper rules finds in the
  field is looked for there when `x` is likely a `p`-th power for a
  prime `p \mid n` and `p` times the degree of the tower is at most
  this (with transcendental generators, 3/8 of it, halved for each
  generator beyond the first: the norms of the factorization are
  multivariate): the steps of the tower are proven first if needed (the
  modular proofs, and Trager's method within its limit), and the
  places of the modular proof of irreducibility of `X^n - x` at which
  `x` is a unit and `p \mid q - 1` must all have shown `x` to be a
  `p`-th power, at least eight of them (more places are tried in that
  case only; an element which is not a `p`-th power is one at about
  one such place in `p`). The roots of `X^p - x` are then found by
  factoring over the tower, the principal one is identified by its
  enclosure, and `\sqrt[n]{x} = \sqrt[n/p]{\sqrt[p]{x}}` (principal
  roots). Otherwise, as with 0, a new generator is adjoined whose
  modulus may be reducible, which the zero tests reduce when they meet
  it (`\sqrt{1/(\pi^{1/6} + 2\pi + 1)^2}` is the inverse of `\pi^{1/6}
  + 2\pi + 1` at once, rather than a generator of degree 2 over
  `\mathbb{Q}(\pi^{1/6})`). Over the profile catalog, the check costs
  about 0.2 seconds in all (most of it the modular test, which the
  adjunction makes anyway).

Lazy field (gr_tower_lazy.h)
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_tower_lazy(gr_ctx_t ctx, gr_ctx_t base, int flags)

    Initializes *ctx* to a field of numbers over *base*
    (currently `\mathbb{Q}`) in which each element carries a pointer to
    a tower and a prefix length in that tower, and is stored in the flat
    representation. *flags* selects the field: the complex field of
    all the numbers the towers can represent by default, its real
    subfield with :macro:`GR_TOWER_LAZY_REAL`, and the algebraic
    numbers with :macro:`GR_TOWER_LAZY_ALGEBRAIC` (the four contexts
    print as "Complex field (lazy towers)", "Real field (lazy towers)",
    "Complex algebraic field (lazy towers)" and "Real algebraic field
    (lazy towers)"); :macro:`GR_TOWER_MERGE_EXPRESS` enables the
    low-effort expression searches when towers are merged (see below).
    In a real context, an operation whose result would not be real fails
    with ``GR_DOMAIN`` (`i`, roots and logarithms of negative numbers,
    inverse trigonometric functions outside their real domains,
    conversions of non-real numbers, non-real polynomial roots are
    dropped by :func:`gr_poly_roots`, and :func:`gr_factor` gives
    quadratic factors for the pairs of nonreal roots); since all
    elements are real, the tests are sign tests. The real and algebraic
    fields can also be created as views of a complex field sharing its
    elements (:func:`gr_ctx_init_tower_lazy_view`). In an algebraic context, `\pi` is not
    available and the elementary functions only at the arguments where
    their values are algebraic (`\exp(0)`, `\log(1)`, `\cos(0)`, `1^y`,
    `x^{p/q}`, and so on -- decided by the Lindemann--Weierstrass and
    Gelfond--Schneider theorems), the others failing with
    ``GR_DOMAIN``. The restrictions apply to the operations called by
    the user, not to the intermediate values of the implementations
    (`\arccos(1/2) = \pi/3` in the real field goes through complex
    logarithms). Towers are created and extended in place as needed
    (for instance by :func:`gr_sqrt`, which first tries to express the
    root in the tower of its argument, and by :func:`gr_exp`,
    :func:`gr_log` and :func:`gr_pi`, which adjoin transcendental
    generators; generators with equal definitions are shared), and
    arithmetic between elements of different towers goes through a tower
    containing both, found by the definition ids of the generators
    involved (or constructed); the operands are mapped by substitution.
    Results are stored with the shortest prefix of their tower in which
    they can be represented. Zero and equality tests are complete for
    algebraic numbers, and for numbers involving `\exp` and `\log`
    they are decided by Richardson's algorithm, the relations found
    being recorded in the towers (shared by all elements referring to
    them). A nonzero polynomial in one generator `\exp(u)`, `\log(u)`
    or `\pi` over the number field of proven algebraic steps below it
    (`u` algebraic) is a nonzero number by the Hermite--Lindemann
    theorem, which decides `\exp(10^{-10000}) - 1` without a
    33000-bit evaluation. The relation search behind a zero test looks at the
    generators the element involves (with those their definitions
    involve), not at the others in the tower: `\arctan(\tan(260515))
    - 260515 + 82924 \pi` is decided the same way whether or not
    `\exp(i)` is in the tower (where the relation `\exp(260515 i) =
    \exp(i)^{260515}`, irrelevant and costly, would otherwise be
    found), and a nonzero element of the field is a nonzero number
    when the transcendental generators it involves are proven
    transcendental, whatever the others. A relation search at the
    default precision is also run when `\exp` or `\log` of an element
    is adjoined, so that for instance `\exp(\log(x))` is `x` at once.
    An exponential `\exp(\sum c_j \log u_j + r \pi i)` with rational
    `c_j`, `r` and algebraic `u_j` is the product of principal powers
    `u_j^{c_j}` (powers of principal roots) and a root of unity, without
    a new generator (`\cos(\arctan(1/3))`, which goes through
    `\exp((\pi i + 2 \log 5 - 4 \log(2+i))/4)`, is algebraic at once).

    Towers are shared by all elements built in them, so a tower can grow
    long (many generators accumulated by earlier, unrelated computations).
    To keep such history from slowing down new computations, a
    transcendental generator (`\exp`, `\log`, `\tan`, `\arctan`) for an
    argument stored with a short prefix of a long tower (at least four
    generators beyond the prefix) is not appended to the long tower but
    adjoined to a copy of the prefix, keeping the definition id and name
    of an equal generator already defined elsewhere; the relation search
    for it then runs over the generators of the prefix only, and its
    values are expressed without the generators of the rest of the tower
    (for instance, after sixteen identities for `\sum_{n=1}^N \cos(n a)`
    have filled a tower with the tangents of multiples of `a/2`, the
    identity `\tan a \tan b = (\cos(a-b) - \cos(a+b)) / (\cos(a-b) +
    \cos(a+b))` for other angles is then found without a relation
    search). Radicals (:func:`gr_sqrt`, `n`-th roots) are still
    adjoined to the tower of their argument, since the expressions and
    factorizations found for them there are too costly to rediscover
    (polynomial roots are sought over the generators their coefficients
    involve, see below), unless the tower continues with at least
    twelve generators beyond the prefix (Watson's sum after Dixon's in
    one context took four minutes, the roots of its sines appended
    after the gamma values of Dixon's sum, which the relation searches
    of its merges then visited; a second). For the same reason the tower of
    one operand of an operation, containing the definitions of the
    other, is used for the result only if no generators of other
    computations lie between them (`i \pi` with `i` in a tower of its own grown in
    place by an earlier computation, say: the result would carry that
    computation's generators in its prefix), and the cached generators
    (the roots of integers and of unity, algebraic numbers) whose
    tower is collected move to a tower where they come first.
    Relations between generators of different copies are found when
    the copies meet: a merge which adds
    transcendental generators to a tower which already has some is
    followed by a relation search at the default precision. Likewise,
    `\exp(r \pi i + y)` with rational `r` is split as `\exp(r \pi i)
    \exp(y)` when the
    root of unity `\exp(r \pi i)` is already in the tower of the
    argument, rather than defining a new exponential. When a relation
    makes an exponential a radical `x^c = u` with `u = r w^c`
    (`r` rational, `w` a monomial in the generators, or a square up to
    the content for `c = 2`), its modulus is the factor of `x^c - u`
    vanishing at it, `w^e f(x/w)` with `f` the factor of `y^c - r` over
    `\mathbb{Q}` (degree 1 when `x/w = \pm 1`).

    Every definition (a generator adjoined to some tower, or found
    expressible in one) gets a name unique in the context, ``a1``,
    ``a2``, ... for algebraic definitions and ``t1``, ``t2``, ... for
    exponentials and logarithms in order of creation (`\pi` and `i`
    keep their names), which follows the generator into every tower
    containing it, so that elements print consistently across towers
    and the symbols can be read back; :func:`gr_gens` lists the
    definitions made so far which are still present in some tower, in
    order. Generators created inside a tower by the relation search (a
    primitive exponential `\exp(u/m)` inserted before its powers, say)
    are named when they are printed or listed at the latest.

    Towers in which no element lives any more are collected: every
    element counts as a reference to its tower, and the towers whose
    count has dropped to zero are freed at the end of the outermost
    operation on the context (in batches). Towers holding the canonical
    definitions of the registries (`\pi` and other constant
    transcendentals, roots of unity, radicals of integers) are kept; the
    cache of algebraic numbers and the recorded expressions of
    definitions are moved to other towers containing them, or dropped.
    The definitions of a collected tower which no other tower contains
    are forgotten: they no longer appear in :func:`gr_gens`, and their
    names can no longer be parsed (printed forms with definitions remain
    readable). Memory therefore stays bounded when a long session keeps
    few elements alive (solving hundreds of polynomials over
    `\mathbb{Q}(\sqrt{2})` leaves a handful of towers).

    The context is thread-safe when FLINT is built with pthreads: every
    operation on it (including the ``gr_tower_lazy_*`` functions below)
    holds a recursive mutex of the context, since operations may refine
    or merge the towers of the shared registry and update the cached
    representation of elements. Operations on one context are thus
    serialized, except for the arithmetic of the two kinds of elements
    below; distinct contexts are independent. The lower-level
    :func:`gr_ctx_init_tower_field` context is not thread-safe.

    Representations: an element holds exactly one of three
    representations (a tagged union): a rational number (``fmpq``), a
    flat element (``fmpz_mpoly_q`` over the flat context of its tower,
    the general case), or a dense form. Rational numbers live in the
    trivial tower, which never changes and is not reference counted, so
    that arithmetic, comparisons and conversions between rational
    elements are done directly on their numerators and denominators,
    without the lock (rational results of operations in other towers are
    moved there). An element of a tower whose first generator `a` is
    algebraic, proven, of degree `2 \le d \le`
    :macro:`GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT` (`2^{20}` by default;
    0 disables dense forms) with a monic integral modulus, and which
    involves only `a` with an integer denominator (`\mathbb{Q}(\sqrt{2})`,
    `\mathbb{Q}(\zeta_5)`, `\mathbb{Q}(2^{1/3})`, `\mathbb{Q}(\zeta_{9973})`),
    can be stored in dense form: an ``fmpq_poly`` in `a` of length less
    than `d` (normalized, without padding, with a positive denominator
    coprime to the content of the numerator), referring to an immutable
    descriptor of `a` and its modulus kept by the tower. Products are
    computed with FLINT's polynomial multiplication followed by a
    reduction by the modulus (folding for prime cyclotomic moduli,
    division with remainder for dense moduli of large degree), inverses
    by the resultant, sums and differences by ``fmpq_poly`` additions.
    Sums, differences, products, inverses and quotients, operations with
    rational scalars, comparisons and the predicates of elements with
    the same descriptor and tower, or with rational elements, are
    computed in this form without the lock when the result already lives
    in that tower (a fresh result goes through the locked operation
    once). A polynomial longer than 256 takes the dense form only when
    at least one coefficient in
    :macro:`GR_TOWER_OPT_DENSE_FORM_SPARSITY` (8 by default) is nonzero,
    so that sparse elements of large cyclotomic fields (a cosine
    `(\zeta^k + \zeta^{-k})/2` in `\mathbb{Q}(\zeta_{9973})`) stay flat.
    The dense form is made at the end of the locked operations on the
    outermost level, for the operands and the result, and only when
    every flat operand can take it: if one stays flat, nothing is
    converted (a running sum of sparse terms would otherwise change form
    at every step). Operations needing the flat form of a rational or
    dense operand work on a temporary flat copy (released when the
    outermost locked operation returns), so that an operand is never
    rewritten while another thread may be reading it without the lock;
    a flat operand converted to the dense form is read by other threads
    only under the lock. The cost of a *Horner step* `x \gets x y + c`
    (a multiplication followed by the addition of a constant, the inner
    loop of polynomial evaluation by Horner's rule) in these
    representations is comparable to that of ``nf_elem`` in the number
    field for small degrees and several times lower for large degrees,
    and comparable to ``fmpq`` for rational numbers (see the design
    notes for timings).

    Number fields: when the tower of the operands has only algebraic
    generators whose moduli are monic univariate integer polynomials,
    each in its own variable (`\mathbb{Q}(\zeta_N)` with the roots of
    unity of prime power orders, `\mathbb{Q}(\sqrt{-163})`, radicals of
    integers, and products of such fields), elements with integer
    denominators are multiplied densely, like ``nf_elem`` in Calcium: the
    generators are replaced by powers of one variable (Kronecker
    substitution, with room for the products of reduced elements), the
    integer polynomials are multiplied with :func:`fmpz_poly_mul`, and
    the product is reduced digit by digit by the moduli (only the nonzero
    coefficients of a sparse modulus such as `x^{64} + 1` being used, and
    `1 + x + \cdots + x^{q-1}` for a prime `q` applied as `x^q = 1`
    followed by one subtraction, while a dense modulus of degree at
    least 64 is applied by :func:`fmpz_poly_rem`). In a compositum
    `\mathbb{Q}(\zeta_{q_1}, \ldots, \zeta_{q_m})` of cyclotomic fields
    of pairwise coprime orders (the decomposition of
    `\mathbb{Q}(\zeta_n)`, `n = q_1 \cdots q_m`), products are instead
    one cyclic convolution of length `n` (the Good-Thomas index map
    `x_k \mapsto x^{c_k}`, `c_k` the idempotents of the Chinese remainder
    theorem, identifies the tensor product of the
    `\mathbb{Z}[x_k]/(x_k^{q_k} - 1)` with `\mathbb{Z}[x]/(x^n - 1)`)
    followed by the reductions by the `\Phi_{q_k}`, when this is shorter
    than the Kronecker array (1155 against 4389 entries for
    `\mathbb{Q}(\zeta_{1155})`); polynomial and matrix products use the
    same layout. A factor which is a monomial, or a power of a root of
    unity such as `\zeta_n^j` (whose reduced form is a tensor product of
    powers `\zeta_{q_k}^{e_k}`, each a monomial or, for `e_k \ge d_k`,
    a short sum), is applied fibre by fibre to the compact array of the
    other factor. The conversions between the
    sparse and the packed representations read and write the packed
    exponents directly. Inverses are computed by the extended
    Euclidean algorithm over `\mathbb{Q}` with a single generator (by
    the closed formula for a quadratic one). With several generators,
    the inverse of an integer numerator `X` is computed modulo primes
    `p` at which all moduli split into distinct linear factors (primes
    `p \equiv 1 \pmod n` for roots of unity of order `n`, and
    `p \equiv 1 \pmod{8|c|}` for `\sqrt{c}`): modulo such a prime the
    algebra is `\mathbb{F}_p^D`, reached by one Vandermonde matrix per
    generator applied along its axis of the coefficient array, where
    inversion is pointwise. The adjugate `N(X)/X` (`N` the norm) and
    `N(X)` are reconstructed by the Chinese remainder theorem from
    batches of primes of doubling size until the product of the primes
    exceeds twice `\|X\|_1 \|N(X)/X\|_\infty G + |N(X)|`, where `G`
    bounds the growth of coefficients in the reduction by the moduli;
    as `X \cdot N(X)/X - N(X)` vanishes modulo every prime used, this
    proves the result. The cost is `O(D (d_1 + \cdots + d_m))` per prime
    against `O(D^3)` for a linear system with the multiplication matrix,
    which remains the fallback (up to degree 512) when primes of this
    kind are too rare; the primes found are kept with the tower. Towers
    whose moduli are not univariate (triangular steps such as
    `\sqrt{1+\sqrt{2}}`, or roots of polynomials with algebraic
    coefficients) of degree up to 64 over `\mathbb{Q}` invert and divide
    through the multiplication matrix as well, its columns `x b_j` on the
    monomial basis obtained from one another by multiplications by
    single generators (2.3–3.2 times faster than the inverse of the
    denominator in the nested representation, which remains the
    fallback). The last inverse is also cached in the tower (fraction-free elimination
    divides by the same pivot many times). Products of sparse elements
    in towers of several generators (fewer than `L/8` term products for
    a Kronecker array of length `L = \prod (2 d_k - 1)`) are computed
    sparsely. Sums with integer denominators only cancel the gcd of the
    denominators, and arithmetic with rational scalars
    (:func:`gr_add_si`, :func:`gr_mul_fmpq`, :func:`gr_div_ui`, ...)
    works on the element directly. Polynomial products
    (:func:`gr_poly_mullow`, from length 4) and matrix products
    (:func:`gr_mat_mul`, from 3 x 3) over such fields pack all the
    coefficients or entries into integer polynomials with common
    denominators (by rows and columns for matrices) and use one
    :func:`fmpz_poly_mul` or :func:`fmpz_poly_mat_mul`. Products with
    rational numbers are scalar multiplications. Determinants use
    ``fmpq_mat`` for rational matrices, the cofactor expansion up to
    4 x 4, LU decomposition (recursive, with the fast matrix products)
    over number fields of degree `D < n` and the division-free Berkowitz
    algorithm otherwise; linear systems use fraction-free LU over number
    fields of degree `D \ge n/2` and LU otherwise (see
    ``examples/number_field_bench.c``: over `\mathbb{Q}(\zeta_{64})`, a
    16 x 16 determinant is more than an order of magnitude faster by
    Berkowitz than by LU, while over `\mathbb{Q}(\sqrt{-163})` LU is
    fastest).

    Generators of number fields: the square root of a negative rational
    number `-A B^2 / D^2` with `A` squarefree is `B \sqrt{-A} / D`
    with `\sqrt{-A}` a single generator (a root of `x^2 + A`), except
    for `\sqrt{-1} = i` and `\sqrt{-3} = 2 \zeta_3 + 1`, which are
    roots of unity: `\mathbb{Q}(\sqrt{-163})` has degree 2, like the
    quadratic fields of Calcium. Otherwise roots of rational numbers are
    canonical products of roots of primes and roots of unity (so
    `\sqrt{-5} \sqrt{-7} = -\sqrt{35}` is found by merging the towers
    `\mathbb{Q}(\sqrt{-5})` and `\mathbb{Q}(\sqrt{-7})`, whose
    compositum contains `\sqrt{35}`, and roots of unity of composite
    orders are products of roots of prime power orders, like
    `\zeta_{255} = \zeta_3^a \zeta_5^b \zeta_{17}^c`).

    Polynomials over the field: :func:`gr_poly_gcd`, :func:`gr_poly_xgcd`
    and :func:`gr_poly_resultant` use the subresultant PRS (the
    Euclidean remainder sequence over a function field with algebraic
    generators has exponential coefficient growth), and
    :func:`gr_poly_roots` finds the roots of a quadratic by the formula
    (with a structured square root of the discriminant), of a binomial
    `a x^n + b` as a principal root times the roots of unity, and
    otherwise, after the roots which are algebraic numbers of small
    degree and height have been recognized numerically and verified
    exactly (a quadratic one as a square root; only guesses of degree
    below `n D` are tried, `D` being the degree of the field of the
    coefficients, since a root of degree `n D` over `\mathbb{Q}` is a
    generic root which the factorization below finds more cheaply),
    factors the polynomial
    exactly over the common tower of its coefficients
    (:func:`gr_tower_poly_factor`, when the tower's steps are proven and
    the norm degree is at most 64, or 24 with transcendental
    generators): the linear factors give roots in the tower, the other
    factors are solved separately (by the formulas for quadratics and
    binomials), and a root of an irreducible factor of degree `\ge 3`
    is adjoined as a proven step, over which the cofactor is factored
    again (a splitting tower). Beyond the limits, each root is
    expressed in the tower when it lies there (`\pi` and `e` for
    `(x - \pi)(x - e)`; the lattice search is skipped when a place of
    the tower shows that there is no such root) and adjoined as a
    dynamic step otherwise, a remaining quadratic factor being solved
    by the formula. The common tower of the coefficients is reduced to
    the generators the coefficients involve (with those their
    definitions involve), in a new tower when the shared tower has
    others: roots of `x^5 - \sqrt{2} x - 1` are sought over
    `\mathbb{Q}(\sqrt{2})` even when the tower of `\sqrt{2}` has been
    extended by the roots of other polynomials or by unrelated radicals
    (a splitting tower of degree 240 instead of one of degree 5760).
    :func:`gr_factor` of a polynomial gives the linear factors
    `x - r` for the roots, and over the real fields the quadratic
    factors `(x - r)(x - \overline{r})` for the pairs of nonreal
    roots.

    Square roots of elements whose flat numerator and denominator are
    squares (or minus squares) as polynomials are read off
    (`\sqrt{(1-\pi)^2} = \pi - 1`, `\sqrt{-(\pi+1)^2} = i (\pi + 1)`);
    otherwise square roots are found exactly through quadratic steps
    (:func:`gr_tower_sqrt`) before any lattice search, repeated roots
    of the same element are identified, and roots of likely powers are
    found by factoring (see
    :macro:`GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT`). A positive rational content is
    taken out of an `n`-th root (`\sqrt{4 \pi} = 2 \sqrt{\pi}`), and the
    root of a monomial in roots of unity and roots of integers is a
    monomial in such roots of higher orders (`\sqrt{4 \sqrt 2} = 2 \cdot
    2^{1/4}`). Roots of rational numbers are
    factored into roots of primes (`\sqrt{6} = \sqrt{2} \sqrt{3}`,
    `\sqrt{-2} = i \sqrt{2}`), and `\exp(r \pi i)` with rational `r`
    is a root of unity: both are products of powers of canonical
    generators of prime power orders, each the generator of a tower of
    its own registered in the context (see
    :func:`gr_tower_adjoin_root_of_unity`), `i` being
    `\exp(2 \pi i / 4)`. A root of a smaller order is not a power of a
    root of a larger order created earlier (`\sqrt{2}` computed after
    `2^{1/128}` lives in `\mathbb{Q}(\sqrt{2})`, not in a field of
    degree 128) until elements of the two towers meet in an operation,
    when the towers are merged; `\exp(\sum_j c_j \log p_j + r \pi i)` with
    rational `c_j`, `r` and positive rational `p_j` is the product
    `\prod_j p_j^{c_j} \exp(r \pi i)` of structured radicals and a root
    of unity (`\exp((\log(2i) - \pi i/2)/2)` is `\sqrt{2}`).
    Logarithms of rational and Gaussian rational
    numbers are decomposed over the logarithms of primes and of Gaussian
    primes, plus `\pi i`; `\log(c \exp(u)) = \log c + u` is taken
    from the definition when `\operatorname{Im} u \in (-\pi, \pi]`.

    :func:`gr_conj` fixes real generators, sends `i` to `-i` and roots
    of unity to their inverses, conjugates `\exp`,
    `\log` and root generators through their definitions (away from the
    branch cut of `\log`, or with the cut handled when the argument is
    a negative real number), and conjugates any other algebraic
    generator as an algebraic number (the conjugate, a root of the same
    polynomial, is adjoined next to the generator; the operation is
    restarted when this changes the tower); :func:`gr_re`, :func:`gr_im`,
    :func:`gr_arg`, :func:`gr_abs` (as `\sqrt{x \bar{x}}`) and
    :func:`gr_sgn` follow. The trigonometric and hyperbolic functions
    and their inverses are expressed through `\exp`, `\log` and
    `\sqrt` (the inverse functions being checked against the principal
    branch numerically, with ``GR_UNABLE`` on the branch cuts), and
    :func:`gr_floor`, :func:`gr_ceil`, :func:`gr_trunc`,
    :func:`gr_nint` and :func:`gr_cmp` are decided for real elements by
    combining enclosures with exact zero tests. Relations found in one
    tower are shared: when towers are merged, generators are identified
    by their definitions (arguments compared exactly), so that
    `\exp(i a)` computed twice is one generator.

    :func:`gr_pow` handles rational exponents through roots and other
    exponents as `\exp(y \log x)`; :func:`gr_set_str` parses expressions
    with ``^`` (``examples/huge_expr.c -tower`` evaluates two algebraic
    numbers given by huge expressions in nested radicals and proves them
    equal, much faster than the Calcium and ``qqbar`` versions of the
    same example).

    :func:`gr_poly_roots` is supported: a polynomial with rational
    coefficients is factored over `\mathbb{Q}`; quadratic factors are
    solved by radicals, and the roots of an irreducible factor of degree
    `d \ge 3` are algebraic numbers (``qqbar``), each in its own
    tower. Every algebraic generator adjoined from an algebraic number
    remembers its *origin* polynomial: when two conjugate roots of the
    same polynomial meet (in an arithmetic operation, or through
    :func:`gr_conj`), the second is adjoined as a root of the Cauchy
    modulus `p(x) / \prod (x - \alpha_i)` over the first, without a
    certification of irreducibility, so that a *splitting tower*
    `\alpha_1, \ldots, \alpha_{d-1}` emerges as roots are combined, the
    last root being given by Vieta's formula; symmetric functions of the
    roots then reduce to the coefficients (``examples/hilbert_matrix_ca.c
    -tower`` verifies the trace and determinant of the Hilbert matrix
    against its eigenvalues up to `n = 22`, where Calcium
    reaches `n = 10` with the Vieta option). Otherwise each squarefree
    factor is
    handled over a
    common tower of its coefficients, its roots being expressed in that
    tower when possible and adjoined as dynamic steps otherwise (the
    linear factor being divided out before the next root is sought).
    Multiplicities are returned as usual.

    Elements are reduced with respect to the triangular set of moduli;
    a monic univariate integer modulus (a cyclotomic polynomial, say) is
    applied by univariate division, sums of reduced elements whose
    denominators involve only the transcendental generators are not
    reduced again, and the denominator of an inverse is rationalized so
    that denominators involve only the transcendental generators: a
    monomial in the algebraic generators through their moduli
    (`1/a = -(c_1 + c_2 a + \cdots + c_d a^{d-1})/c_0`), a polynomial in
    the highest algebraic generator through the adjugate of its
    multiplication matrix modulo the modulus (the norm being free of
    that generator), within a size limit (:macro:`GR_TOWER_OPT_RATIONALIZE_LIMIT`);
    a larger denominator is left in place, as a valid fraction, rather
    than inverted at enormous cost. ``examples/dft.c -tower`` runs the
    exact DFT benchmark of Calcium: with roots of unity as inputs
    (``-input 3``) lengths up to 128 are proved (Calcium fails to prove
    the identity beyond length 32), with `1/(1 + k \pi)` (``-input 4``)
    length 64 and with `1/(1 + \sqrt{k} \pi)` (``-input 5``) length 32,
    where Calcium is much slower. ``examples/machin.c -tower`` proves
    all Machin-like formulas tabulated in ``mp_real`` (up to 48 terms),
    where Calcium proves only those with few terms. (Timings are kept in
    the design notes rather than here.)

    *flags* may contain ``GR_TOWER_MERGE_EXPRESS`` to enable low-effort
    expression searches when towers are merged.

    Linear algebra over this field (and over :func:`gr_ctx_init_tower_field`)
    uses division-free algorithms by default, since inverses of tower
    elements are large.

.. macro:: GR_TOWER_LAZY_REAL
           GR_TOWER_LAZY_ALGEBRAIC

    Flags of :func:`gr_ctx_init_tower_lazy` selecting the real and the
    algebraic subfield.

.. function:: void gr_ctx_init_tower_lazy_view(gr_ctx_t ctx, gr_ctx_t parent, int field_flags)
              int gr_tower_lazy_ctx_same_state(gr_ctx_t ctx1, gr_ctx_t ctx2)

    Initializes *ctx* to a *view* of the lazy field *parent*: a context
    sharing its towers, registry, generator names and lock (and hence
    its elements), restricted further by *field_flags*
    (:macro:`GR_TOWER_LAZY_REAL`, :macro:`GR_TOWER_LAZY_ALGEBRAIC`, added
    to those of *parent*). The shared state is reference counted:
    *parent* and its views may be cleared in any order.
    :func:`gr_set_other` between contexts of the same state
    (:func:`gr_tower_lazy_ctx_same_state`) copies the element and checks
    that it lies in the target field; between different states it
    transfers the element structurally (the generators it involves are
    recreated in the target state from their definitions, with the
    roots matched by their enclosures, so that even roots closer than
    their printed approximations keep their identity). The real field
    is thus the real subfield of a complex field whose computations it
    shares: its values are computed through the complex numbers where
    that is the natural way, and only the values entering it are
    checked. As everywhere in ``gr``, the operations of a view assume
    that their operands are elements of the view: an element of the
    complex field that is used directly in its real view without
    going through :func:`gr_set_other` (which is possible since the
    elements have the same representation) is not checked for realness,
    and the sign tests, comparisons and rounding functions of the view
    then act on its real part.

    Membership of the subfield is enforced where values enter it from
    outside the operations: conversions, :func:`gr_set_str` (whose
    definitions may introduce auxiliary nonreal generators, as long as
    the value is real) and :func:`gr_gens`, which lists only the
    definitions lying in the subfield (in a real field, not `i` or
    `\exp(i)`). The realness test is numerical when the imaginary part
    is separated from zero, and otherwise the exact test `\overline{x}
    = x`. The results of the elementary functions, of root finding and of
    factorization in a real field are rewritten in real terms when their
    representation involves roots of unity `\zeta_j = e^{2 \pi i / n_j}`
    and otherwise only real generators: `x = \operatorname{Re}(x) =
    \sum A \cos(2 \pi \sum_j e_j / n_j)` over the terms `A \prod_j
    \zeta_j^{e_j}` of its numerator, with the cosines in the *tangent
    normal form*. The real trigonometric constants at angles in
    `(\pi/D)\mathbb{Z}` (sines, cosines, tangents, and the real elements
    of the cyclotomic fields they generate) lie in
    `\mathbb{Q}(\zeta_c)^+` with `c = 4D` (`D` odd) or `c = 2D` (`D`
    even), which is generated by `\tau = \tan(\pi/M)` with `M = D`
    (`D` odd), `D/2` (`D \equiv 2 \bmod 4`) or `2D` (`4 \mid D`), the
    smallest `M` with `\operatorname{lcm}(M, 4) = c`. With
    `w = e^{2\pi i/M} = (1 + i\tau)^2/(1 + \tau^2)`, every root of unity of
    order dividing `c` is `i^a w^b`, so a real element of the cyclotomic
    field is a polynomial in `\tau`, computed modulo the minimal
    polynomial of `\tau` with Gaussian arithmetic; `\tau` is a generator
    with ``def_kind`` :macro:`GR_TOWER_TAN_PI` (``def_param`` `M`,
    printed as ``tan(pi/M)``), or its square root form when it is
    quadratic (`\tan(\pi/3) = \sqrt 3`, `\tan(\pi/8) = \sqrt 2 - 1`,
    `\tan(\pi/12) = 2 - \sqrt 3`). This is the algebraic counterpart of
    the half-angle substitution used for the trigonometric functions at
    other real arguments (below): the values at one level share one
    generator, also those in a quadratic subfield. Thus
    `\sin(2\pi/5) = (11\tau - \tau^3)/8` and
    `\cos(\pi/5) = (7 - \tau^2)/8` with `\tau = \tan(\pi/5)`,
    `\sin(\pi/7) = (\tau^5 - 20\tau^3 + 19\tau)/16` with
    `\tau = \tan(\pi/7)`, and the real factors of `x^4 + 1` are
    `x^2 \pm \sqrt 2 x + 1`. The trigonometric functions of a real field
    at rational multiples of `\pi` (denominators up to 240), the
    elementary parts of the gamma and polygamma values at rationals (in
    the reflection formula, Gauss's digamma theorem, and the cotangent
    parts of the Hurwitz normal form) and the real forms above all use
    this form: `\Gamma(3/5) = \pi (15\tau - \tau^3)/(10\,\Gamma(2/5))`,
    `\psi'(4/5) = (15 - \tau^2)\pi^2/5 - \psi'(1/5)`, and
    `\psi(1/7)` through `\tan(\pi/7)` and the logarithms of its
    polynomials. (Tangents are computed at the level of the doubled angle,
    `\tan(r\pi) = \sin(2r\pi)/(1 + \cos(2r\pi))`, whose field is
    smaller.) So that the real generators so introduced are
    not re-expressed through the roots of unity when towers are merged,
    a state with a real view merges with
    :macro:`GR_TOWER_MERGE_REAL_FIRST`: a real algebraic generator which
    is expressible in the target tower only through nonreal algebraic
    generators is adjoined and moved before them instead (the steps after
    it become dynamic, and are proven or refined by Trager's method when
    the norms are small: a fifth root of unity gets a quadratic modulus
    over `\mathbb{Q}(\tan(\pi/5))`), and a root of unity placed in front
    of a tower by its definition is placed after the tangent generators
    already there; thus the real quadratic factors of `x^5 - 3` have
    their coefficients in `\mathbb{Q}(3^{1/5}, \tan(\pi/5))`.

    The trigonometric functions of real arguments in a real field use
    the real generators: `\sin x` and `\cos x` are rational functions of
    `\tan(x/2)` (`\cos 1 = (1 - t^2)/(1 + t^2)` with `t = \tan(1/2)`),
    `\tan x` is `\tan(x)`, and `\arctan`, `\arcsin` (`\arctan(y /
    \sqrt{1 - y^2})`) and `\arccos` (`\pi/2 - \arcsin`) use `\arctan`
    generators; a rational multiple `r \pi` in the argument is split off
    and handled by the addition formulas with the algebraic values at
    `r \pi` (so that no generator has `\pi` in its argument), odd
    functions get generators of positive arguments, and `\arctan(y)` is
    recognized as `r \pi` for a rational `r` with a small denominator
    when `\tan(r \pi) = y` (verified exactly). A complex field computing
    with the same state (a view) keeps its exponential forms; relations
    between the two forms are found when both occur in a tower (which
    then contains `i`).

.. function:: void gr_tower_lazy_ctx_set_gen_flags(gr_ctx_t ctx, int flags)
              int gr_tower_lazy_ctx_gen_flags(gr_ctx_t ctx)

    How new generators are chosen. With
    :macro:`GR_TOWER_GENS_SPLIT_IMAGINARY`, the square root of a negative
    rational number is `i \sqrt{A}` (the generators `i` and
    `\sqrt{p}` for the primes `p` dividing `A`, shared with the real
    square roots) rather than the single generator `\sqrt{-A}`, which
    is the default (0). By default, roots of unity of composite orders
    are products of roots of unity of prime power orders (linearly
    disjoint fields: canonical, sparse representations, monomials for
    the powers of roots of unity), and square roots of positive
    rationals are products of square roots of primes. With
    :macro:`GR_TOWER_GENS_COMPOSITE_ROOTS` (which sets the option
    :macro:`GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT` to
    :macro:`GR_TOWER_COMPOSITE_ROOTS_DEGREE`, 128, when it is 0;
    clearing it sets the option to 0), a root of unity is a power of
    one generator `\zeta_N` (a root of `\Phi_N`) for `N` the least common
    multiple of the orders requested, as in Calcium: arithmetic within
    one cyclotomic field is dense univariate (3–4 times faster for
    `\mathbb{Q}(\zeta_{60})` and `\mathbb{Q}(\zeta_{255})`), but
    computations involving many orders build large fields (the DFT
    benchmark is several times slower); towers merged in this mode adjoin
    a missing root of unity as `\zeta_L` for `L` the least common multiple
    of its order and those present, which become powers of it, up to
    the degree limit and :macro:`GR_TOWER_CYCLOTOMIC_ORDER_LIMIT`, beyond
    which roots of prime power orders are used. The flags are a view of
    the options :macro:`GR_TOWER_OPT_SPLIT_IMAGINARY`,
    :macro:`GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT` and
    :macro:`GR_TOWER_OPT_COMPOSITE_RADICALS`. With
    :macro:`GR_TOWER_GENS_COMPOSITE_RADICALS`, the square root of a
    positive rational `A B^2 / D^2` (`A` squarefree) is `B \sqrt{A} / D`
    with `\sqrt{A}` one generator (relations such as
    `\sqrt{6} = \sqrt{2}\sqrt{3}` are then found by the merges rather
    than by construction). Elements created before a change remain
    valid and compatible with the new ones.

.. function:: int gr_tower_lazy_ctx_set_option(gr_ctx_t ctx, slong option, slong value)
              slong gr_tower_lazy_ctx_get_option(gr_ctx_t ctx, slong option)

    Sets or gets a tuning option (see the table above) of the lazy
    field, shared by its views and its towers; a change takes effect
    at once, also for the existing elements. Returns ``GR_DOMAIN`` for
    a value out of the range of the option. The Python interface takes
    the names (:func:`gr_tower_option_name`), also as keyword arguments
    of the constructors (``ComplexField_tower(cyclotomic_degree_limit=64)``).

.. function:: void gr_tower_lazy_ctx_set_print(gr_ctx_t ctx, int flags, slong digits)
              int gr_tower_lazy_ctx_print_flags(gr_ctx_t ctx)
              int gr_tower_lazy_ctx_field_flags(gr_ctx_t ctx)

    Printing of elements, by any combination of
    :macro:`GR_TOWER_PRINT_NUMERIC` (a numerical value with *digits*
    digits), :macro:`GR_TOWER_PRINT_SYMBOLIC` (the expression in the
    generators) and :macro:`GR_TOWER_PRINT_DEFS` (the definitions of the
    generators the expression involves, transitively, in definition
    order); the default is the expression with the definitions. For
    instance `2 + e^2 + 3 e^{e^2}` prints as
    ``3*t2+t1+2 {t1 = exp(2); t2 = exp(t1)}``, as ``3*t2+t1+2``
    without the definitions, as ``4863.92`` numerically, and as
    ``4863.92 {3*t2+t1+2 where t1 = exp(2); t2 = exp(t1)}`` with
    everything. The definitions are written as ``sqrt(x)``,
    ``root(x, n)``, ``exp(x)``, ``log(x)``, ``exp(2*pi*i/n)`` and, for an
    algebraic generator without a closed form,
    ``root(m(a), approximation)`` with the minimal polynomial *m* over the
    field below, in the order of creation of the generators (a definition
    only involves earlier ones). A rational number prints as such. The
    field flags are those of :func:`gr_ctx_init_tower_lazy`.

    The numerical value omits a real or imaginary part which is zero
    to within the working precision (as for `\sqrt{-163}`, printed as
    ``12.7671*i``).

    The printed form with the definitions is read back by
    :func:`gr_set_str`: the definitions are evaluated in order, each in
    the names defined before it -- the names are local to the string,
    so a string printed in one context or session is read in another --
    then the expression in all of them (a ``root(m, approximation)``
    definition takes the unique root of *m* within the precision of the
    approximation -- the digits printed for a root are enough to
    separate it from the other roots of its polynomial -- and the string
    is rejected with ``GR_UNABLE`` when several roots agree to the given
    digits). A string without definitions is read in the generators of
    the context by their names, with ``pi``, ``i`` and the elementary
    functions (``sqrt(2)``, ``exp(a1)``). Conversions between lazy
    contexts (:func:`gr_set_other`), for instance from the complex field
    into the real or algebraic subfield, do not go through strings but
    transfer the generators structurally; the conversion fails if the
    number does not belong to the target field.

.. macro:: GR_TOWER_PRINT_NUMERIC
           GR_TOWER_PRINT_SYMBOLIC
           GR_TOWER_PRINT_DEFS

The symbolic expression of an element (:func:`gr_get_fexpr`) is
``Where(expr, Def(a_1, Sqrt(2)), Def(t_1, Exp(a_1)), ...)`` with the
generators as the symbols ``a_1``, ``t_1``, ... (``Pi`` and
``NumberI`` for `\pi` and `i`), the definitions as ``Exp``,
``Log``, ``Sqrt``, ``Pow(x, Div(1, n))``, ``Exp(Div(Mul(2, Pi,
NumberI), n))`` and ``PolynomialRootNearest(List(...),
Decimal(...))`` (with the origin polynomial over `\mathbb{Q}` when
the generator has one, otherwise the coefficients of the minimal
polynomial over the field below as expressions). The LaTeX form
follows from :func:`fexpr_get_str_latex`.

.. function:: int gr_tower_lazy_get_qqbar(qqbar_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_get_acb(acb_t res, const gr_tower_lazy_elem_t x, slong prec, gr_ctx_t ctx)
              int gr_tower_lazy_root_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong n, gr_ctx_t ctx)
              gr_tower_struct * gr_tower_lazy_get_tower(slong * level, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              const fmpz_mpoly_q_struct * gr_tower_lazy_get_data(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              void gr_tower_lazy_ctx_stats(gr_ctx_t ctx)
              slong gr_tower_lazy_ctx_num_towers(gr_ctx_t ctx)

    Conversions and access to the underlying tower representation
    of an element (the tower stays valid while *x* lives in it: towers
    without elements are collected). :func:`gr_tower_lazy_ctx_num_towers`
    returns the number of towers of the context after collecting the
    unused ones. The ``qqbar`` conversion rebuilds the algebraic part
    of the tower over `\mathbb{Q}` when the tower has transcendental
    generators which neither the element nor the moduli involve.
    :func:`gr_tower_lazy_get_data` returns the element as a
    ``fmpz_mpoly_q`` over the flat context of its tower, converting
    the element to that representation in place if it is currently
    stored as a rational number or in dense form (see
    :type:`gr_tower_lazy_elem_struct`); the pointer stays valid until the
    next operation involving *x*.
    :func:`gr_tower_lazy_ctx_stats` prints the number of towers, their
    largest degree and the number of steps; with the environment
    variable ``GR_TOWER_STATS_VERBOSE`` set, also the sizes of the
    moduli of every tower with at least two steps (flat terms and
    coefficient bits, or "deferred", and nested coefficients).

.. type:: gr_tower_lazy_elem_struct

.. type:: gr_tower_lazy_elem_t

    The elements of the lazy fields (``elem_size`` of the context); a
    :type:`gr_tower_lazy_elem_t` is an array of length one of
    :type:`gr_tower_lazy_elem_struct`, permitting it to be passed by
    reference. The structure is public so that its data can be read
    without copies, but it must not be written: an element refers to a
    tower shared with other elements, and the context converts, reduces
    and moves the data of elements lazily, under its lock. The fields:

    * ``slong repr``: which member of the union ``elem`` holds the
      element, one of the constants below.
    * ``fmpq elem.q`` (:macro:`GR_TOWER_LAZY_REPR_RATIONAL`): a rational
      number.
    * ``fmpz_mpoly_q_struct elem.flat.data`` and
      ``fmpz_mpoly_ctx_struct * elem.flat.mctx``
      (:macro:`GR_TOWER_LAZY_REPR_FLAT`): a fraction over the flat
      context of the tower. The context and the reduction of the data
      may lag behind the tower (they are brought up to date by the next
      operation or accessor), so that the data should be read through
      :func:`gr_tower_lazy_get_data`.
    * ``fmpq_poly_struct elem.dense.poly``
      (:macro:`GR_TOWER_LAZY_REPR_DENSE`): a canonical polynomial of
      degree less than `d` in the first generator `a` of the tower, of
      degree `d` (the modulus is obtained with
      :func:`gr_tower_lazy_get_fmpq_poly`). The other member of this
      branch, ``elem.dense.nf``, is internal.
    * ``F``, ``level``, ``structure_version``, ``reduced_version``,
      ``shallow``: the tower (through ``F->T``) and the bookkeeping of
      the lazy updates (internal).

    The layout may change between versions of FLINT; the accessors
    below are the stable interface. A rational or dense element may be
    read while other threads read it; reading a flat element races with
    the lock-holding thread, which may convert it to the dense form, so
    that it should be read through the locked accessors.

.. macro:: GR_TOWER_LAZY_REPR_FLAT
           GR_TOWER_LAZY_REPR_RATIONAL
           GR_TOWER_LAZY_REPR_DENSE

    The representations of :type:`gr_tower_lazy_elem_struct`.

.. function:: int gr_tower_lazy_repr(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)

    The current representation of *x* (one of the constants above). It
    is a property of the storage, not of the value: the same number may
    be flat in one element and dense in another, and an element may
    change representation after an operation that involves it.

.. function:: int gr_tower_lazy_get_fmpq_poly(fmpq_poly_t res, fmpz_poly_t modulus, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)

    Sets *res* to *x* as a polynomial in the first generator `a` of
    its tower, of degree less than the degree of `a`, and *modulus*
    (unless ``NULL``) to the minimal polynomial of `a` (monic, with
    integer coefficients), whatever the representation of *x*. A rational
    number is a constant polynomial with modulus `x`. Returns
    ``GR_DOMAIN`` if *x* involves another generator or `a` is not
    algebraic, and ``GR_UNABLE`` if the modulus of `a` is not monic
    integral or *x* has an algebraic denominator. The polynomial is
    canonical when the step of `a` is proven: equal elements give equal
    polynomials. For example, `(1 + \sqrt{2})^3` gives `7 + 5x` with
    modulus `x^2 - 2`.

.. function:: void gr_tower_lazy_init(gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              void gr_tower_lazy_clear(gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              void gr_tower_lazy_swap(gr_tower_lazy_elem_t x, gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              void gr_tower_lazy_set_shallow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_randtest(gr_tower_lazy_elem_t res, flint_rand_t state, gr_ctx_t ctx)
              int gr_tower_lazy_write(gr_stream_t out, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_set_str(gr_tower_lazy_elem_t res, const char * s, gr_ctx_t ctx)
              int gr_tower_lazy_get_fexpr(fexpr_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_zero(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
              int gr_tower_lazy_one(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
              int gr_tower_lazy_pi(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
              int gr_tower_lazy_i(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
              int gr_tower_lazy_gens(gr_vec_t vec, gr_ctx_t ctx)
              int gr_tower_lazy_set(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_set_si(gr_tower_lazy_elem_t res, slong c, gr_ctx_t ctx)
              int gr_tower_lazy_set_ui(gr_tower_lazy_elem_t res, ulong c, gr_ctx_t ctx)
              int gr_tower_lazy_set_fmpz(gr_tower_lazy_elem_t res, const fmpz_t c, gr_ctx_t ctx)
              int gr_tower_lazy_set_fmpq(gr_tower_lazy_elem_t res, const fmpq_t c, gr_ctx_t ctx)
              int gr_tower_lazy_set_other(gr_tower_lazy_elem_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int gr_tower_lazy_get_si(slong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_get_ui(ulong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_get_fmpz(fmpz_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_get_fmpq(fmpq_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              truth_t gr_tower_lazy_is_zero(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              truth_t gr_tower_lazy_is_one(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              truth_t gr_tower_lazy_equal(const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_neg(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_add(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_sub(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_mul(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_div(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_inv(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_add_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx)
              int gr_tower_lazy_add_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx)
              int gr_tower_lazy_add_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx)
              int gr_tower_lazy_add_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx)
              int gr_tower_lazy_sub_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx)
              int gr_tower_lazy_sub_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx)
              int gr_tower_lazy_sub_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx)
              int gr_tower_lazy_sub_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx)
              int gr_tower_lazy_mul_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx)
              int gr_tower_lazy_mul_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx)
              int gr_tower_lazy_mul_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx)
              int gr_tower_lazy_mul_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx)
              int gr_tower_lazy_div_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx)
              int gr_tower_lazy_div_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx)
              int gr_tower_lazy_div_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx)
              int gr_tower_lazy_div_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx)
              int gr_tower_lazy_pow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_sqrt(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_abs(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_conj(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_exp(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_log(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_poly_mullow(gr_tower_lazy_elem_struct * res, const gr_tower_lazy_elem_struct * p1, slong len1, const gr_tower_lazy_elem_struct * p2, slong len2, slong n, gr_ctx_t ctx)
              int gr_tower_lazy_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
              int gr_tower_lazy_poly_factor(gr_poly_t c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
              int gr_tower_lazy_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
              int gr_tower_lazy_mat_det(gr_tower_lazy_elem_t res, const gr_mat_t A, gr_ctx_t ctx)
              int gr_tower_lazy_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

    The methods of the ``gr`` interface of the lazy fields, which the
    generic element functions of ``gr.h`` (:func:`gr_add` and so on)
    dispatch to, and which may also be called directly; the semantics
    are those of the generic functions. Operations on rational elements,
    and on dense elements of the same field (with rational scalars),
    take lock-free paths; the others hold the context's lock.
    :func:`gr_tower_lazy_set_other` converts an element *x* of the
    context *x_ctx*. In the real and algebraic views (see
    :func:`gr_ctx_init_tower_lazy_view`), operations whose result
    leaves the subfield return ``GR_DOMAIN``.

.. function:: int gr_tower_lazy_sin(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_cos(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_tan(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_sinh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_cosh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_tanh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_asin(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_acos(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_atan(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_asinh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_acosh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_atanh(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_re(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_im(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_arg(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_sgn(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_csgn(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              truth_t gr_tower_lazy_is_real(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_cmp(int * res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_cmpabs(int * res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_floor(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_ceil(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_trunc(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_nint(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)

    The elementary functions and real-number operations of the lazy
    field, as installed in its method table. The inverse functions
    return ``GR_UNABLE`` on their branch cuts (where the logarithmic
    formulas may differ from the principal branch). As for
    :type:`qqbar_t`, the rounding functions act on the real part
    (``floor(2 + i sqrt(2)) = 2``), *csgn* is the sign of the real part,
    or of the imaginary part when the real part is zero, and *cmpabs*
    compares `x \bar x` and `y \bar y`; *cmp* returns ``GR_DOMAIN`` for
    elements known to be nonreal (``GR_UNABLE`` when realness cannot be
    decided). All of these are exact: enclosures decide when they
    separate, and the exact zero test decides otherwise.

.. function:: int gr_tower_lazy_gamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_rgamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_lgamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_beta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
              int gr_tower_lazy_digamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_polygamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_erf(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_erfc(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_erfi(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_lambertw(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_lambertw_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t k, gr_ctx_t ctx)
              int gr_tower_lazy_zeta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_hurwitz_zeta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t a, gr_ctx_t ctx)
              int gr_tower_lazy_polylog(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_lerch_phi(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t a, gr_ctx_t ctx)
              int gr_tower_lazy_dirichlet_l(gr_tower_lazy_elem_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_tower_lazy_elem_t s, gr_ctx_t ctx)
              int gr_tower_lazy_dilog(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_elliptic_k(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_elliptic_e(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
              int gr_tower_lazy_euler(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
              int gr_tower_lazy_catalan(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
              int gr_tower_lazy_hypgeom_0f1(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t b, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
              int gr_tower_lazy_hypgeom_1f1(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t a, const gr_tower_lazy_elem_t b, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
              int gr_tower_lazy_hypgeom_2f1(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t a, const gr_tower_lazy_elem_t b, const gr_tower_lazy_elem_t c, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
              int gr_tower_lazy_hypgeom_pfq(gr_tower_lazy_elem_t res, const gr_vec_t a, const gr_vec_t b, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)
              int gr_tower_lazy_hypgeom_pfq_vec(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_struct * a, slong p, const gr_tower_lazy_elem_struct * b, slong q, const gr_tower_lazy_elem_t z, int flags, gr_ctx_t ctx)

    The special functions of the lazy field, as installed in its method
    table, with canonical arguments (see *Special functions* above).
    Poles give ``GR_DOMAIN``. The order of ``polygamma`` must be an
    integer; ``polylog`` with a non-integer order is evaluated at the
    roots of unity (through the Lerch transcendent), and ``lerch_phi`` at
    `z = 1`, at the roots of unity, and at `a = 1` with an integer `s`
    (``GR_UNABLE`` otherwise). :func:`gr_dirichlet_l` dispatches to
    ``gr_tower_lazy_dirichlet_l`` for the lazy fields. In the real fields, a value
    which is not real gives ``GR_DOMAIN``; in the algebraic fields, only
    algebraic values (such as `\Gamma(n)` and `\zeta(-n)`) are returned,
    and ``GR_UNABLE`` otherwise, since the transcendence of special
    function values is not known in general. The printed forms
    (``gamma(...)``, ``erf(...)``, ``erfi(...)``, ``lambertw(x, k)``,
    ``digamma(...)``, ``polygamma(m, x)``, ``polylog(s, x)``,
    ``zeta(...)``, ``elliptic_k(...)``, ``elliptic_e(...)``, ``euler``,
    ``catalan``) are read back by the parser (the hypergeometric
    functions, of several arguments, are not). The hypergeometric
    functions take the flag 1 for the regularized functions (as
    :func:`gr_hypgeom_pfq`); their normal form is described in
    *Hypergeometric functions* above.

.. function:: int gr_tower_lazy_modular_lambda(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_modular_j(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_modular_delta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_dedekind_eta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_eisenstein_e(gr_tower_lazy_elem_t res, ulong k, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_eisenstein_g(gr_tower_lazy_elem_t res, ulong k, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_jacobi_theta(gr_tower_lazy_elem_t res1, gr_tower_lazy_elem_t res2, gr_tower_lazy_elem_t res3, gr_tower_lazy_elem_t res4, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_jacobi_theta_1(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_jacobi_theta_2(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_jacobi_theta_3(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_jacobi_theta_4(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_weierstrass_p(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_weierstrass_p_prime(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_weierstrass_sigma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_elliptic_roots(gr_tower_lazy_elem_t e1, gr_tower_lazy_elem_t e2, gr_tower_lazy_elem_t e3, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)
              int gr_tower_lazy_elliptic_invariants(gr_tower_lazy_elem_t g2, gr_tower_lazy_elem_t g3, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx)

    The modular lambda function, the `j`-invariant, the discriminant
    `\Delta = \eta^{24}`, the Dedekind eta function, the Eisenstein series
    `E_k` (normalized to constant term 1) and `G_k = 2\zeta(k) E_k` for even
    `k \ge 2` (``GR_DOMAIN`` for other `k`, as for the complex numbers;
    ``GR_UNABLE`` for `k > 1000`), the Jacobi theta functions (with the conventions of
    :func:`acb_modular_theta`, `q = e^{\pi i \tau}`), and the Weierstrass
    functions `\wp`, `\wp'`, `\sigma`, the roots `e_1, e_2, e_3` and the
    invariants `g_2, g_3` of the lattice `\mathbb{Z} + \tau \mathbb{Z}`
    (as :func:`acb_elliptic_p`; `\wp` at a lattice point gives
    ``GR_DOMAIN``), in the normal form
    described in *Modular forms and theta functions* above. They are
    installed in the method table of the lazy field (``gr_modular_j`` and
    so on). A point with `\operatorname{Im}(\tau) \le 0` gives
    ``GR_DOMAIN``; ``GR_UNABLE`` is returned when `\tau` cannot be reduced
    exactly (for instance when the sign of `\operatorname{Im}(\tau)` is
    undecided). In the algebraic fields, only algebraic values (`j` and
    `\lambda` at CM points) are returned. ``modular_lambda``,
    ``modular_j``, ``modular_delta`` and ``dedekind_eta`` are read back by
    the parser.

The Python binding in ``src/python/flint_ctypes.py`` offers the four
lazy fields as ``ComplexField_tower`` (``CC_tower``), ``RealField_tower``
(``RR_tower``), ``ComplexAlgebraicField_tower`` (``QQbar_tower``) and
``RealAlgebraicField_tower`` (``AA_tower``), with ``set_print`` for the
display, ``set_gens`` for the choice of generators, ``gens()``, ``fexpr()`` and ``latex()`` on elements and
``tower()`` to inspect the tower an element lives in, as well as fixed
towers: ``gr_tower`` (``adjoin_sqrt``, ``adjoin_root``,
``adjoin_qqbar``, ``adjoin_algebraic``, ``adjoin_root_of_unity``,
``adjoin_pi``, ``adjoin_exp``, ``adjoin_log``, each returning the new
generator as an element of the new top field) and ``TowerField``
(``tower.field()``), documented by doctests. The binding runs the
Calcium test cases of that file
against the lazy field (``test_tower``, ``test_tower_trigonometric``),
including identities which Calcium cannot decide, such as
`\arccos(\cos(\sqrt{2} - 1)) = \sqrt{2} - 1`. Two further test
functions collect the open issues of the Calcium issue tracker which
the tower field settles (``test_tower_calcium_issues``: cancellation of
`a^{50} - a^{51} a^{-1}` for `a = 417/(962 \pi + 80808)`,
`\exp((\log(2i) - \pi i/2)/2) = \sqrt{2}`, the algebraicity of
`\cos(\arccos(\sqrt{2} - \sqrt{3})/3)`, `\exp(\log(M)) = M` for a
matrix, and so on) and examples from the documentation of Sage's
``QQbar`` (``test_tower_sage_examples``: nested radicals, the roots of
`x^5 - x - 1`, the discriminant of a cubic against its roots, the
regular 34-gon in radicals, Lehmer's polynomial and the ARPREC identity
for `\alpha^{630} - 1`), each running quickly on the tower field
where Sage's ``QQbar`` is reported to take minutes or not to finish.

Two randomized tests exercise the lazy field against independent
oracles: ``t-lazy_qqbar`` evaluates random algebraic expressions
(rationals, `i`, roots of unity, radicals, field operations, principal
roots, conjugates, real parts) both in the lazy field and in ``qqbar``,
and ``t-lazy_numeric`` evaluates random expressions with `\exp`,
`\log`, `\pi`, roots, trigonometric functions and their inverses in
parallel in ball arithmetic (the enclosures must overlap, a value whose
ball excludes zero must not be found zero) and checks identities holding
by construction (`\exp(\log x) = x`, `x = \operatorname{Re} x + i
\operatorname{Im} x`, `\sin^2 x + \cos^2 x = 1`, `\arctan(\tan x) = x`,
and so on). They run with a larger ``FLINT_TEST_MULTIPLIER`` for longer
searches. ``t-modular`` checks that the modular proofs never prove a
reducible modulus (`X^n - b^n` over towers containing `b`, `X^2 - 2`
over `\mathbb{Q}(\sqrt 2)`, `X^4 + 4`) and prove most irreducible
radical steps, and ``t-threads`` runs several threads on one lazy
context with shared elements.
