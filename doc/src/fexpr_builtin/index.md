---
title: fexpr_builtin.h
header: fexpr_builtin.h
description: The builtin symbols of symbolic expressions, with their meaning and LaTeX rendering, and functions to look them up.
category: Exact real and complex numbers
order: 9
---

This module defines the symbols with a predefined meaning that can be
used in symbolic expressions of type `fexpr_t` (see [](fexpr)). The
builtin symbols cover logic, sets, arithmetic, calculus, linear
algebra and several hundred mathematical functions, following the
conventions of the Fungrim formula database [@Fungrim]. Every builtin
symbol has a LaTeX rendering, used by `fexpr_write_latex`; the symbols
that [](ca) and [](qqbar) understand also have a numerical meaning.

By convention, builtin symbol names start with an uppercase letter and
have at least two characters. Single letters and names starting with a
lowercase letter are left for variables, with the exception of the
Greek letter names in [](#greek), which are builtin so that they render
as Greek letters.

For each builtin symbol `Symbol`, the header defines the C constant
`FEXPR_Symbol`, its index in the table of builtin symbols, to be used
with `fexpr_set_symbol_builtin`, `fexpr_is_builtin_symbol`,
`fexpr_call_builtin1` and related functions. Below, the symbols are
documented by name, grouped by topic, with example expressions and their
LaTeX rendering. Expressions are written in the notation of
`fexpr_print`.

## C interface {#c-interface}

::: type
fexpr_builtin_symbol

The enumeration of builtin symbols: `FEXPR_AGM`, `FEXPR_AGMSequence`,
..., `FEXPR_zeta_`, in strictly increasing alphabetical (ASCII) order of
the names, followed by `FEXPR_BUILTIN_LENGTH`, the number of builtin
symbols. The alphabetical order allows lookup by binary search.
:::

::: type
fexpr_symbol_info
_fexpr_latex_writer

An entry of the table of builtin symbols: the enumerator `symbol`, the
name `string`, the LaTeX code `latex_string` used for the symbol itself
(empty if it has a writer function), and `latex_writer`, a pointer of
type `_fexpr_latex_writer` (that is, `void (*)(calcium_stream_t, const
fexpr_t, ulong)`) to the function that renders calls of the symbol, or
`NULL` for the default rendering.
:::

::: variable
fexpr_builtin_table

The table `const fexpr_symbol_info fexpr_builtin_table[FEXPR_BUILTIN_LENGTH]`,
indexed by `fexpr_builtin_symbol`.
:::

::: function
slong fexpr_builtin_lookup(const char * s)

Returns the index of the builtin symbol with the name `s`, or $-1$ if
`s` is not the name of a builtin symbol. The lookup is a binary search
in `fexpr_builtin_table`.
:::

::: function
const char * fexpr_builtin_name(slong n)

Returns a pointer to the name of the builtin symbol with index `n`, a
string owned by the table that must not be modified or freed.
:::

::: function
slong fexpr_builtin_length(void)

Returns the number of builtin symbols, `FEXPR_BUILTIN_LENGTH`.
:::

::: macro
FEXPR_SYMBOL_Pos
FEXPR_SYMBOL_Neg
FEXPR_SYMBOL_Add
FEXPR_SYMBOL_Sub
FEXPR_SYMBOL_Mul
FEXPR_SYMBOL_Div
FEXPR_SYMBOL_Pow

The header words of the expressions consisting of the symbols `Pos`,
`Neg`, `Add`, `Sub`, `Mul`, `Div` and `Pow`, for fast comparisons when
recognising arithmetic operations.
:::

## Variables and iteration {#variables}

The following symbols bind variables locally. They are syntactic
constructs that do not denote mathematical objects on their own.

::: macro
For

`For(x, ...)` binds the symbol `x` as a local variable in its parent
expression; the further arguments give a range, set or point whose
meaning depends on the parent operator. `For(x, S)` iterates over a set,
optionally followed by a filter predicate `P(x)`, as in `Set(f(x),
For(x, S), P(x))` ($\left\{ f(x) : x \in
S\,\mathbin{\operatorname{and}}\, P(x) \right\}$) or `Sum(f(x), For(x,
S))` ($\sum_{x \in S} f(x)$). `For(x, a, b)` lets `x` run between
endpoints, as in `Sum(f(n), For(n, a, b))` ($\sum_{n=a}^{b} f(n)$) or
`Integral(f(x), For(x, a, b))` ($\int_{a}^{b} f(x) \, dx$). `For(x, a)`
lets `x` approach a point or selects a point of differentiation, as in
`Limit(f(x), For(x, a))` ($\lim_{x \to a} f(x)$) or `Derivative(f(x),
For(x, a))` ($f'\!\left(a\right)$), and `For(x, a, n)` gives the order
of a derivative, as in `Derivative(f(x), For(x, a, n))`
(${f}^{(n)}\!\left(a\right)$).

::: details
A sum `Sum(f(n), For(n, a, b))` is empty if $b < a$. An integral
`Integral(f(x), For(x, a, b))` follows the straight-line path from $a$
to $b$; swapping the endpoints negates it.
:::
:::

::: macro
Where
Def

`Where(f(x), Def(x, a))` defines `x` as an alias for `a` within `f(x)`,
so that it is equivalent to `f(a)`; it renders as `Where(f(x), Def(x,
a))` ($f(x)\; \text{ where } x = a$). `Def(f(t), a)` defines `f` as a
function of the dummy variable `t`, and `Def(Tuple(a, b), T)`
destructures a tuple. A `Where` expression may contain several
definitions. `Def` occurs only inside `Where`.
:::

::: macro
Fun

`Fun(x, expr)` is the anonymous function mapping `x` to `expr`, rendered
`Fun(x, Pow(x, 2))` ($x \mapsto {x}^{2}$); `x` is bound within the
expression.
:::

::: macro
Step
Repeat

Notation for sequences: `Step(x, For(x, a, b))` ($a, a + 1, \ldots, b$)
and `Repeat(x, n)` ($\underbrace{x, \ldots, x}_{n \text{ times}}$).
:::

## Booleans and logic {#logic}

Logical values, predicates and connectives. Logical connectives are
rendered as words unless `FEXPR_LATEX_LOGIC` is given or the expression
is wrapped in `Logic`.

::: macro
Equal
NotEqual

`Equal(a, b)` ($a = b$) is `True` if `a` and `b` represent the same
object and `False` otherwise; with more arguments it states that all are
equal. `NotEqual(a, b)` ($a \ne b$) is `Not(Equal(a, b))`.
:::

::: macro
Same

`Same(a, b)` is `a` (equivalently `b`) if `a` and `b` represent the same
object, and `Undefined` otherwise. It asserts or emphasises that
expressions have the same value within a formula; it accepts any number
of arguments.
:::

::: macro
True
False

The logical constants `True` ($\operatorname{True}$) and `False`
($\operatorname{False}$).
:::

::: macro
Not
And
Or
Equivalent
Implies

Logical negation `Not(x)` ($\operatorname{not} x$), conjunction `And(x,
y)` ($x \;\mathbin{\operatorname{and}}\; y$) and disjunction `Or(x, y)`
($x \;\mathbin{\operatorname{or}}\; y$), both with any number of
arguments, equivalence `Equivalent(x, y)` ($x \iff y$), which is `Equal`
applied to logical values, and implication `Implies(x, y)` ($x
\;\implies\; y$).
:::

::: macro
Exists
All

Quantifiers: `Exists(f(x), For(x, S))` ($f(x) \;\text{ for some } x \in
S$) and `All(f(x), For(x, S))` ($f(x) \;\text{ for all } x \in S$),
optionally with a filter predicate as third argument, as in `All(f(x),
For(x, S), P(x))` ($f(x) \;\text{ for all } x \in S \text{ with }
P(x)$).
:::

::: macro
Cases
Case
Otherwise

Piecewise definitions: `Cases(Case(f(x), P(x)), Case(g(x), Otherwise))`
($\begin{cases} f(x), & P(x)\\g(x), & \text{otherwise}\\ \end{cases}$).
Any number of cases can be given. If several predicates hold
simultaneously, the corresponding values must agree; if no predicate
holds and there is no `Otherwise` case, the value is undefined.
:::

## Tuples, lists and sets {#collections}

Collections and set operations.

::: macro
Tuple
List
Set

A tuple `Tuple(a, b, c)` ($\left(a, b, c\right)$), a list `List(a, b,
c)` ($\left[a, b, c\right]$) and a finite set `Set(a, b, c)` ($\left\{a,
b, c\right\}$). A set can also be given by comprehension, as in
`Set(f(x), For(x, S), P(x))` ($\left\{ f(x) : x \in
S\,\mathbin{\operatorname{and}}\, P(x) \right\}$).
:::

::: macro
Item
Length
Cardinality
Concatenation

The element with index `n` of a tuple or list `Item(L, n)` (${L}_{n}$),
the length of a tuple or list `Length(L)` ($\# L$), the cardinality of a
set `Cardinality(S)` ($\# S$), and the concatenation of lists
`Concatenation(L, M)` ($L \,^\frown M$).
:::

::: macro
Element
NotElement
EqualAndElement

Membership `Element(x, S)` ($x \in S$) and `NotElement(x, S)` ($x \notin
S$); `EqualAndElement(a, b, S)` states that `a` equals `b` and is an
element of `S`, rendered `EqualAndElement(a, b, S)` ($a = b \in S$).
:::

::: macro
Union
Intersection
SetMinus
Subset
SubsetEqual

Set operations `Union(A, B)` ($A \cup B$), `Intersection(A, B)` ($A \cap
B$), `SetMinus(A, B)` ($A \setminus B$) and the relations `Subset(A, B)`
($A \subset B$) and `SubsetEqual(A, B)` ($A \subseteq B$).
:::

::: macro
CartesianProduct
CartesianPower
Subsets
Sets
Tuples

The Cartesian product `CartesianProduct(A, B)` ($A \times B$) and power
`CartesianPower(S, n)` (${S}^{n}$), the power set `Subsets(S)`
($\mathscr{P}(S)$) of all subsets of `S`, the class `Sets`
($\operatorname{Sets}$) of all sets, and `Tuples` (the class of all
tuples), `Tuples(S)` (tuples with elements in `S`) and `Tuples(S, n)`
(tuples of length `n` with elements in `S`).
:::

## Numbers and arithmetic {#numbers}

Constants, arithmetic, sets of numbers, intervals and infinities.

::: macro
Undefined
Unknown

`Undefined` ($\mathfrak{u}$) is the special value $\mathfrak{u}$ for
undefined results, for example of $0/0$. `Unknown`
($\operatorname{Unknown}$) is the meta-value used by [](ca) for values
that could not be computed.
:::

::: macro
Pi
NumberI
NumberE
GoldenRatio
Euler
CatalanConstant
KhinchinConstant
GlaisherConstant

The constants `Pi` ($\pi$), the imaginary unit `NumberI` ($i$), the base
of the natural logarithm `NumberE` ($e$), the golden ratio `GoldenRatio`
($\varphi$), Euler's constant `Euler` ($\gamma$), Catalan's constant
`CatalanConstant` ($G$), Khinchin's constant `KhinchinConstant` ($K$)
and Glaisher's constant `GlaisherConstant` ($A$). The verbose names
`NumberI` and `NumberE` leave `i`, `I`, `e` and `E` free as variable
names.
:::

::: macro
RootOfUnity

`RootOfUnity(n)` ($\zeta_{n}$) is the principal root of unity $e^{2 \pi
i / n}$, and `RootOfUnity(n, k)` ($\zeta_{n, k}$) its $k$-th power.
:::

::: macro
Decimal

`Decimal(str)` is the rational number given in decimal floating-point
notation by the string `str`, for example `Decimal("-3.25e-725")`,
rendered `Decimal("-3.25e-725")` ($-3.25 \cdot 10^{-725}$).
:::

::: macro
AlgebraicNumberSerialized
PolynomialRootIndexed
PolynomialRootNearest

Algebraic numbers: `AlgebraicNumberSerialized(List(c_0, ..., c_d),
enclosure)` gives the root of the polynomial with integer coefficients
$c_0, \ldots, c_d$ isolated by the ball `enclosure` (a `RealBall`, or
`Add(RealBall(...), Mul(RealBall(...), NumberI))` for a complex ball);
`PolynomialRootIndexed( List(c_0, ..., c_d), k)` is the $k$-th root
(from 1) in the order of `qqbar_cmp_root_order`; and
`PolynomialRootNearest(List(c_0, ..., c_d), x)` is the root closest to
the number `x`. These are produced and read by the `qqbar_get_fexpr`
functions and `qqbar_set_fexpr` (see [](qqbar)).
:::

::: macro
Enclosure
Approximation
Guess
EqualNearestDecimal

Notation for numerical information: a value enclosed by `x`,
`Enclosure(x)` ($\left(\in x \right)$); an approximation
`Approximation(x, n)`; a conjectured value `Guess(x)` ($\left(?\; x
\right)$); and the statement that `x` rounds to a given decimal string
with a given number of digits, `EqualNearestDecimal(x, "1.25", 3)` ($x =
\text{``1.25''} \;\, {\scriptstyle (\text{nearest } 3 \text{
digits})}$).
:::

::: macro
Pos
Neg
Add
Sub
Mul
Div
Pow

Arithmetic: `Pos(x)` ($+x$), `Neg(x)` ($-x$), `Add(x, y)` ($x + y$),
`Sub(x, y)` ($x - y$), `Mul(x, y)` ($x y$), `Div(x, y)` ($\frac{x}{y}$)
and `Pow(x, y)` (${x}^{y}$). `Add` and `Mul` accept any number of
arguments.
:::

::: macro
Sqrt
Root

The principal square root `Sqrt(x)` ($\sqrt{x}$) and the principal
$n$-th root `Root(x, n)`.
:::

::: macro
Less
LessEqual
Greater
GreaterEqual

Order relations `Less(x, y)` ($x < y$), `LessEqual(x, y)` ($x \le y$),
`Greater(x, y)` ($x > y$) and `GreaterEqual(x, y)` ($x \ge y$). With
more arguments they denote chains, such as `Less(x, y, z)` ($x < y <
z$).
:::

::: macro
NN
ZZ
QQ
RR
CC
Primes
AlgebraicNumbers
RealAlgebraicNumbers

Sets of numbers: the natural numbers including 0 `NN` ($\mathbb{N}$),
the integers `ZZ` ($\mathbb{Z}$), the rationals `QQ` ($\mathbb{Q}$), the
reals `RR` ($\mathbb{R}$), the complex numbers `CC` ($\mathbb{C}$), the
positive primes `Primes` ($\mathbb{P}$), the algebraic numbers
`AlgebraicNumbers` ($\overline{\mathbb{Q}}$) and the real algebraic
numbers `RealAlgebraicNumbers` ($\overline{\mathbb{Q}}_{\mathbb{R}}$).
:::

::: macro
IntegersGreaterEqual
IntegersLessEqual
Range

Integer ranges: `IntegersGreaterEqual(x)` ($\mathbb{Z}_{\ge x}$) and
`IntegersLessEqual(x)` ($\mathbb{Z}_{\le x}$) for an extended real
number $x$, and `Range(a, b)` ($\{a, a + 1, \ldots, b\}$), the integers
from $a$ to $b$ (empty if $a > b$).
:::

::: macro
Interval
OpenInterval
ClosedOpenInterval
OpenClosedInterval
RealBall
OpenRealBall

Real intervals with extended real endpoints: `Interval(a, b)` ($\left[a,
b\right]$), `OpenInterval(a, b)` ($\left(a, b\right)$),
`ClosedOpenInterval(a, b)` ($\left[a, b\right)$), `OpenClosedInterval(a,
b)` ($\left(a, b\right]$); and real balls with real center $m$ and
extended real radius $r$, closed `RealBall(m, r)` ($\left[m \pm
r\right]$) or open `OpenRealBall(m, r)` ($\left(m \pm r\right)$).
:::

::: macro
OpenComplexDisk
ClosedComplexDisk
UpperHalfPlane
UnitCircle
BernsteinEllipse
Lattice

Regions of the complex plane: the open disk `OpenComplexDisk(m, r)`
($D(m, r)$), the closed disk `ClosedComplexDisk(m, r)` ($\overline{D}(m,
r)$) with center $m$ and radius $r$, the upper half-plane
`UpperHalfPlane` ($\mathbb{H}$) of numbers with positive imaginary part,
the unit circle `UnitCircle` ($\mathbb{T}$), the Bernstein ellipse
`BernsteinEllipse(rho)` ($\mathcal{E}_{\rho}$) with foci $\pm 1$ and sum
of semi-axes $\rho$, and the lattice `Lattice(a, b)` ($\Lambda_{(a,
b)}$) generated by `a` and `b`.
:::

::: macro
Infinity
UnsignedInfinity

The positive signed infinity `Infinity` ($\infty$) and the unsigned
infinity `UnsignedInfinity` ($\hat{\infty}$). Signed infinities in other
directions are written as products such as `Mul(NumberI, Infinity)`.
:::

::: macro
RealSignedInfinities
ComplexSignedInfinities
RealInfinities
ComplexInfinities

Sets of infinities: `RealSignedInfinities` ($\{\pm \infty\}$), equal to
$\{+\infty, -\infty\}$, `ComplexSignedInfinities` ($\{[e^{i \theta}]
\infty\}$), equal to $\{e^{i \theta} \infty : \theta \in \mathbb{R}\}$,
and the same sets together with unsigned infinity, `RealInfinities`
($\{\hat{\infty}, \pm \infty\}$) and `ComplexInfinities`
($\{\hat{\infty}, [e^{i \theta}] \infty\}$).
:::

::: macro
ExtendedRealNumbers
ProjectiveRealNumbers
SignExtendedComplexNumbers
ProjectiveComplexNumbers

Extended number systems: `ExtendedRealNumbers`
($\overline{\mathbb{R}}$), equal to $\mathbb{R} \cup \{\pm \infty\}$,
`ProjectiveRealNumbers` ($\hat{\mathbb{R}}_{\infty}$), equal to
$\mathbb{R} \cup \{\tilde \infty\}$, `SignExtendedComplexNumbers`
($\overline{\mathbb{C}}_{[e^{i \theta}] \infty}$), the complex numbers
with all signed infinities, and `ProjectiveComplexNumbers`
($\hat{\mathbb{C}}_{\infty}$) $= \mathbb{C} \cup \{\tilde \infty\}$ (the
Riemann sphere).
:::

::: macro
RealSingularityClosure
ComplexSingularityClosure

The value sets of [](ca): `RealSingularityClosure`
($\overline{\mathbb{R}}_{\text{Sing}}$), equal to $\mathbb{R} \cup \{\pm
\infty, \tilde \infty, \mathfrak{u}\}$ and `ComplexSingularityClosure`
($\overline{\mathbb{C}}_{\text{Sing}}$), the complex numbers together
with all signed infinities, unsigned infinity and `Undefined`.
:::

## Operators and calculus {#operators}

Operators that aggregate values over a range or set, limits, derivatives
and integrals. The variable is bound with `For` (see [](#variables)).

::: macro
Sum
Product

Sums and products: `Sum(f(n), For(n, a, b))` ($\sum_{n=a}^{b} f(n)$),
`Product(f(n), For(n, a, b))` ($\prod_{n=a}^{b} f(n)$), and over sets
with optional filter, `Sum(f(x), For(x, S), P(x))` ($\sum_{\textstyle{x
\in S \atop P(x)}} f(x)$).
:::

::: macro
PrimeSum
PrimeProduct
DivisorSum
DivisorProduct

Sums and products over primes and divisors: `PrimeSum(f(p), For(p))`
($\sum_{p} f(p)$), `PrimeProduct(f(p), For(p))` ($\prod_{p} f(p)$),
`DivisorSum(f(d), For(d, n))` ($\sum_{d \mid n} f(d)$) and
`DivisorProduct(f(d), For(d, n))` ($\prod_{d \mid n} f(d)$); a filter
predicate may follow.
:::

::: macro
Zeros
UniqueZero
Solutions
UniqueSolution
Poles

Sets of zeros, solutions and poles: `Zeros(f(x), For(x, S))`
($\mathop{\operatorname{zeros}\,}\limits_{x \in S} f(x)$),
`Solutions(Equal(f(x), 0), For(x, S))`
($\mathop{\operatorname{solutions}\,}\limits_{x \in S} f(x) = 0$),
`Poles(f(z), For(z, S))` ($\mathop{\operatorname{poles}\,}\limits_{z \in
S} f(z)$); `UniqueZero` and `UniqueSolution` denote the unique element
of such a set.
:::

::: macro
Supremum
Infimum
Minimum
Maximum
ArgMin
ArgMax
ArgMinUnique
ArgMaxUnique

Extrema over sets: `Supremum(f(x), For(x, S))`
($\mathop{\operatorname{sup}\,}\limits_{x \in S} f(x)$), `Infimum(f(x),
For(x, S))` ($\mathop{\operatorname{inf}\,}\limits_{x \in S} f(x)$),
`Minimum(f(x), For(x, S))` ($\mathop{\min\,}\limits_{x \in S} f(x)$),
`Maximum(f(x), For(x, S))` ($\mathop{\max\,}\limits_{x \in S} f(x)$),
the sets of minimizers `ArgMin(f(x), For(x, S))`
($\mathop{\operatorname{arg\,min}\,}\limits_{x \in S} f(x)$) and
maximizers `ArgMax`, and the unique minimizer or maximizer
`ArgMinUnique`, `ArgMaxUnique`.
:::

::: macro
Limit
SequenceLimit
RealLimit
LeftLimit
RightLimit
ComplexLimit
MeromorphicLimit
SequenceLimitInferior
SequenceLimitSuperior

Limits: `Limit(f(x), For(x, a))` ($\lim_{x \to a} f(x)$), the limit of a
sequence `SequenceLimit(a(n), For(n, Infinity))` ($\lim_{n \to \infty}
a(n)$), limits through real arguments `RealLimit`, from the left
`LeftLimit(f(x), For(x, a))` ($\lim_{x \to {a}^{-}} f(x)$) or the right
`RightLimit`, through complex arguments `ComplexLimit`, the limit in the
sense of meromorphic functions `MeromorphicLimit`, and the limit
inferior and superior of sequences `SequenceLimitInferior`,
`SequenceLimitSuperior`.
:::

::: macro
AsymptoticTo

Asymptotic equivalence: `AsymptoticTo(f(n), g(n), n, Infinity)` ($f(n)
\sim g(n), \; n \to \infty$).
:::

::: macro
Derivative
RealDerivative
ComplexDerivative
ComplexBranchDerivative
MeromorphicDerivative

Derivatives: `Derivative(f(x), For(x, a))` ($f'\!\left(a\right)$) and of
higher order `Derivative(f(x), For(x, a, n))`
(${f}^{(n)}\!\left(a\right)$); the variants specify real
differentiation, complex differentiation, complex differentiation along
a branch, and differentiation of meromorphic functions.
:::

::: macro
Integral
Path
CurvePath

Integrals `Integral(f(x), For(x, a, b))` ($\int_{a}^{b} f(x) \, dx$)
along the straight segment, or along a path given as a polygon `Path(a,
b, c)` ($a \rightsquigarrow b \rightsquigarrow c$) or a parametrised
curve `CurvePath(f(t), For(t, a, b))` ($\left(f(t),\, t : a
\rightsquigarrow b\right)$) in place of the endpoints.
:::

::: macro
Residue
ComplexZeroMultiplicity

The residue `Residue(f(z), For(z, a))`
($\mathop{\operatorname{res}}\limits_{z=a} f(z)$) and the order of a
zero `ComplexZeroMultiplicity(f(z), For(z, a))`
($\mathop{\operatorname{ord}}\limits_{z=a} f(z)$) at a point.
:::

::: macro
IsHolomorphicOn
IsMeromorphicOn
AnalyticContinuation

Predicates `IsHolomorphicOn(f(z), For(z, S))` ($f(z) \text{ is
holomorphic on } z \in S$) and `IsMeromorphicOn`, and the analytic
continuation of $f$ along a path, `AnalyticContinuation(f(z), For(z, a,
b))` ($\mathop{\text{Continuation}}\limits_{\displaystyle{z: a
\rightsquigarrow b}} \, f(z)$).
:::

## Matrices and linear algebra {#matrices}

Matrices and matrix groups.

::: macro
Matrix
Row
Column
RowMatrix
ColumnMatrix
DiagonalMatrix
Matrix2x2

Matrices: `Matrix(Row(a, b), Row(c, d))`
($\displaystyle{\begin{pmatrix}a & b \\c & d\end{pmatrix}}$), also
written `Matrix(List(Row(a, b), Row(c, d)))`, the row and column vectors
`RowMatrix(a, b, c)` ($\displaystyle{\begin{pmatrix}a & b &
c\end{pmatrix}}$) and `ColumnMatrix(a, b, c)`
($\displaystyle{\begin{pmatrix}a \\ b \\ c\end{pmatrix}}$), the diagonal
matrix `DiagonalMatrix(a, b, c)` ($\displaystyle{\begin{pmatrix}a & & \\
& b & \\ & & c\end{pmatrix}}$) and `Matrix2x2(a, b, c, d)`
($\displaystyle{\begin{pmatrix}a & b \\ c & d\end{pmatrix}}$). `Column`
is reserved for columns as `Row` is used for rows.
:::

::: macro
ZeroMatrix
IdentityMatrix
HilbertMatrix

Special $n \times n$ matrices: `ZeroMatrix(n)` ($0_{n}$),
`IdentityMatrix(n)` ($I_{n}$) and the Hilbert matrix `HilbertMatrix(n)`
($H_{n}$).
:::

::: macro
Det
Spectrum
SingularValues

The determinant `Det(A)` ($\operatorname{det}(A)$), the set of
eigenvalues `Spectrum(A)` ($\operatorname{spec}(A)$) and the singular
values `SingularValues(A)` ($\sigma(A)$).
:::

::: macro
Matrices
SL2Z
PSL2Z
SpecialLinearGroup
GeneralLinearGroup

Sets of matrices: `Matrices(CC, n)` ($\operatorname{M}_{n}(\mathbb{C})$)
and `Matrices(CC, m, n)` ($\operatorname{M}_{m \times n}(\mathbb{C})$),
the modular groups `SL2Z` ($\operatorname{SL}_2(\mathbb{Z})$) and
`PSL2Z` ($\operatorname{PSL}_2(\mathbb{Z})$), and the groups
`SpecialLinearGroup(n, ZZ)`
($\operatorname{SL}_{n}\!\left(\mathbb{Z}\right)$) and
`GeneralLinearGroup(n, ZZ)`
($\operatorname{GL}_{n}\!\left(\mathbb{Z}\right)$).
:::

## Polynomials, series and rings {#rings}

Polynomials, power series and algebraic structures.

::: macro
Polynomial
Coefficient
PolynomialDegree

`Polynomial(List(c_0, ..., c_n), x)` is the polynomial $c_0 + c_1 x +
\cdots + c_n x^n$ in `x` (the variable may be omitted, as in the output
of `ca_poly_get_fexpr`); `Coefficient(f, x, n)` ($[{x}^{n}] f$) is the
coefficient of $x^n$ and `PolynomialDegree(f)` ($\deg(f)$) the degree.
:::

::: macro
Pol
Ser

Placeholders for polynomials and series in formulas: `Pol(x)`
($\operatorname{Pol}(x)$), `Ser(x)` ($\operatorname{Ser}(x)$).
:::

::: macro
Polynomials
PolynomialFractions
FormalPowerSeries
FormalLaurentSeries
FormalPuiseuxSeries

Rings: `Polynomials(ZZ, x)` ($\mathbb{Z}[x]$), the field of fractions
`PolynomialFractions(QQ, x)` ($\mathbb{Q}(x)$), formal power series
`FormalPowerSeries(ZZ, x)` ($\mathbb{Z}[[x]]$), formal Laurent series
`FormalLaurentSeries(ZZ, x)` ($\mathbb{Z}(\!(x)\!)$) and formal Puiseux
series `FormalPuiseuxSeries(ZZ, x)`
($\mathbb{Z}\!\left\langle\!\left\langle x
\right\rangle\!\right\rangle$).
:::

::: macro
Zero
One
Characteristic

The zero `Zero(R)` ($0_{R}$) and one `One(R)` ($1_{R}$) of a ring and
its characteristic `Characteristic(R)` ($\operatorname{char}(R)$).
:::

::: macro
Rings
CommutativeRings
Fields
QuotientRing
FiniteField

The classes of rings `Rings`, commutative rings `CommutativeRings` and
fields `Fields`; the quotient ring `QuotientRing(R, I)`; and the finite
field `FiniteField(q)` with $q$ elements.
:::

::: macro
QSeriesCoefficient
EqualQSeriesEllipsis
IndefiniteIntegralEqual

Notation from Fungrim formulas: `QSeriesCoefficient(f, q, n)`, the
coefficient of $q^n$ in a $q$-series; `EqualQSeriesEllipsis`, which
states that a $q$-series begins with given terms; and
`IndefiniteIntegralEqual(f(x), g(x), x)` ($\int f(x) \, dx = g(x) +
\mathcal{C}$).
:::

::: macro
Call
CallIndeterminate

Function application: `Call(f, x)` ($f(x)$) is $f$ applied to $x$, and
`CallIndeterminate(f, x, v)` applies $f$ to the indeterminate $v$,
rendered `CallIndeterminate(f, x, v)` ($f(v)$).
:::

## Special functions {#special-functions}

Mathematical functions, following the conventions of the Fungrim formula
database [@Fungrim]. Several functions accept an optional last argument
$r$ for the $r$-th derivative, such as `AiryAi(z, r)`. The functions
supported by numerical evaluation in [](ca) and `fexpr_get_decimal_str`
are a subset.

::: macro
Abs
Sign
Csgn
RealAbs
Re
Im
Arg
Conjugate

Complex parts: `Abs(z)` ($\left|z\right|$), the complex sign `Sign(z)`
($\operatorname{sgn}(z)$), equal to $z/|z|$ for $z \ne 0$; the
real-valued sign extension `Csgn(z)` ($\operatorname{csgn}(z)$); the
absolute value of a real number `RealAbs(x)` ($\left|x\right|$), `Re(z)`
($\operatorname{Re}(z)$), `Im(z)` ($\operatorname{Im}(z)$), the argument
`Arg(z)` ($\arg(z)$) and the complex conjugate `Conjugate(z)`
($\overline{z}$).
:::

::: macro
Max
Min
Floor
Ceil

`Max(a, b)` ($\max(a, b)$) and `Min(a, b)` ($\min(a, b)$), with any
number of arguments, and `Floor(x)` ($\left\lfloor x \right\rfloor$),
`Ceil(x)` ($\left\lceil x\right\rceil$).
:::

::: macro
KroneckerDelta
IsOdd
IsEven
CongruentMod
Divides
Mod

Integer predicates and operations: `KroneckerDelta(i, j)`
($\delta_{(i,j)}$), `IsOdd(n)` ($n \text{ odd}$), `IsEven(n)` ($n \text{
even}$), `CongruentMod(a, b, m)` ($a \equiv b \pmod {m }$), `Divides(d,
n)` ($d \mid n$) and the remainder `Mod(a, n)` ($a \bmod n$).
:::

::: macro
GCD
LCM
XGCD

`GCD(a, b)` ($\gcd(a, b)$) and `LCM(a, b)` ($\operatorname{lcm}(a, b)$),
with any number of arguments, and the extended GCD `XGCD(a, b)`
($\operatorname{xgcd}(a, b)$).
:::

::: macro
IsPrime
Prime
PrimePi

`IsPrime(n)` ($n \text{ prime}$), the $n$-th prime `Prime(n)` ($p_{n}$)
and the prime-counting function `PrimePi(x)` ($\pi(x)$).
:::

::: macro
DivisorSigma
MoebiusMu
EulerPhi
LiouvilleLambda
SquaresR

Arithmetic functions: the divisor function `DivisorSigma(k, n)`
($\sigma_{k}\!\left(n\right)$), the Möbius function `MoebiusMu(n)`
($\mu(n)$), Euler's totient function `EulerPhi(n)` ($\varphi(n)$), the
Liouville function `LiouvilleLambda(n)` ($\lambda(n)$) and the number of
representations as a sum of squares `SquaresR(k, n)`
($r_{k}\!\left(n\right)$).
:::

::: macro
LegendreSymbol
JacobiSymbol
KroneckerSymbol
DiscreteLog

`LegendreSymbol(a, p)` ($\left(\frac{a}{p}\right)$), `JacobiSymbol(a,
n)` ($\left(\frac{a}{n}\right)$), `KroneckerSymbol(a, n)`
($\left(\frac{a}{n}\right)$), and the discrete logarithm `DiscreteLog(n,
b, p)` ($(\epsilon : {b}^{\epsilon} \equiv n \text{ mod }p)$), the least
$\epsilon \ge 0$ with $b^\epsilon \equiv n \pmod p$.
:::

::: macro
Exp
Log

The exponential function `Exp(z)` ($e^{z}$) and the principal branch of
the natural logarithm `Log(z)` ($\log(z)$).
:::

::: macro
Sin
Cos
Tan
Cot
Sec
Csc
Sinc

Trigonometric functions `Sin(z)` ($\sin(z)$), `Cos(z)` ($\cos(z)$),
`Tan(z)` ($\tan(z)$), `Cot(z)` ($\cot(z)$), `Sec(z)` ($\sec(z)$),
`Csc(z)` ($\csc(z)$), and the sinc function `Sinc(z)`
($\operatorname{sinc}(z)$), equal to $\sin(z)/z$.
:::

::: macro
Sinh
Cosh
Tanh
Coth
Sech
Csch

Hyperbolic functions `Sinh(z)` ($\sinh(z)$), `Cosh(z)` ($\cosh(z)$),
`Tanh(z)` ($\tanh(z)$), `Coth(z)` ($\coth(z)$), `Sech(z)`
($\operatorname{sech}(z)$) and `Csch(z)` ($\operatorname{csch}(z)$).
:::

::: macro
Asin
Acos
Atan
Acot
Asec
Acsc
Atan2

Principal branches of the inverse trigonometric functions `Asin(z)`
($\operatorname{asin}(z)$), `Acos(z)` ($\operatorname{acos}(z)$),
`Atan(z)` ($\operatorname{atan}(z)$), `Acot(z)`
($\operatorname{acot}(z)$), `Asec(z)` ($\operatorname{asec}(z)$),
`Acsc(z)` ($\operatorname{acsc}(z)$), and the two-argument arctangent
`Atan2(y, x)` ($\operatorname{atan2}(y, x)$).
:::

::: macro
Asinh
Acosh
Atanh
Acoth
Asech
Acsch

Principal branches of the inverse hyperbolic functions `Asinh(z)`
($\operatorname{asinh}(z)$), `Acosh(z)` ($\operatorname{acosh}(z)$),
`Atanh(z)` ($\operatorname{atanh}(z)$), `Acoth(z)`
($\operatorname{acoth}(z)$), `Asech(z)` ($\operatorname{asech}(z)$) and
`Acsch(z)` ($\operatorname{acsch}(z)$).
:::

::: macro
LambertW

The Lambert W function: the principal branch `LambertW(z)` ($W(z)$) and
the branch with index $k$, `LambertW(z, k)` ($W_{k}(z)$).
:::

::: macro
Factorial
DoubleFactorial
Binomial
RisingFactorial
FallingFactorial

`Factorial(n)` ($n!$), `DoubleFactorial(n)` ($n!!$), the binomial
coefficient `Binomial(n, k)` (${n \choose k}$), and the rising and
falling factorials `RisingFactorial(z, k)` ($\left(z\right)_{k}$) and
`FallingFactorial(z, k)` ($\left(z\right)^{\underline{k}}$).
:::

::: macro
HarmonicNumber
Fibonacci
BellNumber
PartitionsP
LandauG

The harmonic number `HarmonicNumber(n)` ($H_{n}$), the Fibonacci number
`Fibonacci(n)` ($F_{n}$), the Bell number `BellNumber(n)`
($\operatorname{B}_{n}$), the partition function `PartitionsP(n)`
($p(n)$) and Landau's function `LandauG(n)` ($g(n)$).
:::

::: macro
BernoulliB
BernoulliPolynomial
EulerE
EulerPolynomial

Bernoulli numbers `BernoulliB(n)` ($B_{n}$) and polynomials
`BernoulliPolynomial(n, x)` ($B_{n}\!\left(x\right)$), and Euler numbers
`EulerE(n)` ($E_{n}$) and polynomials `EulerPolynomial(n, x)`
($E_{n}\!\left(x\right)$).
:::

::: macro
StirlingCycle
StirlingS1
StirlingS2

Stirling numbers: unsigned of the first kind `StirlingCycle(n, k)`
($\left[{n \atop k}\right]$), signed of the first kind `StirlingS1(n,
k)` ($s\!\left(n, k\right)$) and of the second kind `StirlingS2(n, k)`
($\left\{{n \atop k}\right\}$).
:::

::: macro
SymmetricPolynomial
Cyclotomic
SloaneA

The elementary symmetric polynomial `SymmetricPolynomial(k, List(x, y,
z))` ($e_{k}\!\left(\left[x, y, z\right]\right)$), the cyclotomic
polynomial `Cyclotomic(n, x)` ($\Phi_{n}\!\left(x\right)$), and the
$n$-th term of an OEIS sequence `SloaneA("A000001", n)`
($\text{A000001}\!\left(n\right)$).
:::

::: macro
Gamma
LogGamma
DigammaFunction
DigammaFunctionZero

The gamma function `Gamma(z)` ($\Gamma(z)$), the log-gamma function
`LogGamma(z)` ($\log \Gamma(z)$) (with the branch cut of the principal
logarithm of $\Gamma$ continued analytically from the positive real
axis), the digamma function `DigammaFunction(z)` ($\psi(z)$) and its
derivatives (polygamma functions) `DigammaFunction(z, n)`
(${\psi}^{(n)}\!\left(z\right)$), and the zeros of the digamma function
`DigammaFunctionZero(n)` ($x_{n}$).
:::

::: macro
BetaFunction
BarnesG
LogBarnesG

The beta function `BetaFunction(a, b)` ($\mathrm{B}(a, b)$), the Barnes
G-function `BarnesG(z)` ($G(z)$) and its logarithm `LogBarnesG(z)`
($\log G(z)$).
:::

::: macro
StirlingSeriesRemainder
LogBarnesGRemainder

Remainders of asymptotic series after $n$ terms: of the Stirling series
for $\log \Gamma(z)$, `StirlingSeriesRemainder(n, z)`
($R_{n}\!\left(z\right)$), and of the asymptotic series for $\log G(z)$,
`LogBarnesGRemainder(n, z)`.
:::

::: macro
ChebyshevT
ChebyshevU
LegendreP
JacobiP
HermiteH
LaguerreL
GegenbauerC
SphericalHarmonicY

Orthogonal polynomials and related functions: Chebyshev polynomials
`ChebyshevT(n, x)` ($T_{n}\!\left(x\right)$) and `ChebyshevU(n, x)`
($U_{n}\!\left(x\right)$), Legendre polynomials (and functions)
`LegendreP(n, x)` ($P_{n}\!\left(x\right)$), Jacobi polynomials
`JacobiP(n, a, b, x)`, Hermite polynomials `HermiteH(n, x)`
($H_{n}\!\left(x\right)$), Laguerre polynomials `LaguerreL(n, x)` and
`LaguerreL(n, a, x)`, Gegenbauer polynomials `GegenbauerC(n, lamda, x)`,
and spherical harmonics `SphericalHarmonicY(ell, m, theta, phi)`.
:::

::: macro
LegendrePolynomialZero
GaussLegendreWeight

The $k$-th zero `LegendrePolynomialZero(n, k)` ($x_{n, k}$) of the
Legendre polynomial $P_n$, and the corresponding Gauss–Legendre
quadrature weight `GaussLegendreWeight(n, k)` ($w_{n, k}$).
:::

::: macro
Erf
Erfc
Erfi

The error function `Erf(z)` ($\operatorname{erf}(z)$), the complementary
error function `Erfc(z)` ($\operatorname{erfc}(z)$) and the imaginary
error function `Erfi(z)` ($\operatorname{erfi}(z)$).
:::

::: macro
UpperGamma
LowerGamma
IncompleteBeta
IncompleteBetaRegularized

The upper and lower incomplete gamma functions `UpperGamma(s, z)`
($\Gamma(s, z)$) and `LowerGamma(s, z)` ($\gamma(s, z)$), and the
incomplete beta function `IncompleteBeta(x, a, b)`
($\mathrm{B}_{x}\!\left(a, b\right)$) and its regularised version
`IncompleteBetaRegularized(x, a, b)` ($I_{x}\!\left(a, b\right)$).
:::

::: macro
LogIntegral
ExpIntegralE
ExpIntegralEi
SinIntegral
SinhIntegral
CosIntegral
CoshIntegral
FresnelC
FresnelS

Exponential and trigonometric integrals: the logarithmic integral
`LogIntegral(z)` ($\operatorname{li}(z)$), the generalised exponential
integral `ExpIntegralE(s, z)` ($E_{s}\!\left(z\right)$), the exponential
integral `ExpIntegralEi(z)` ($\operatorname{Ei}(z)$), the sine and
cosine integrals `SinIntegral(z)` ($\operatorname{Si}(z)$),
`CosIntegral(z)` ($\operatorname{Ci}(z)$), their hyperbolic counterparts
`SinhIntegral(z)` ($\operatorname{Shi}(z)$), `CoshIntegral(z)`
($\operatorname{Chi}(z)$), and the Fresnel integrals `FresnelC(z)`
($C(z)$) and `FresnelS(z)` ($S(z)$).
:::

::: macro
AiryAi
AiryBi
AiryAiZero
AiryBiZero

The Airy functions `AiryAi(z)` ($\operatorname{Ai}(z)$) and `AiryBi(z)`
($\operatorname{Bi}(z)$), with derivatives `AiryAi(z, 1)`
($\operatorname{Ai}'\!\left(z\right)$), and their $n$-th zeros on the
negative real axis `AiryAiZero(n)` ($a_{n}$) and `AiryBiZero(n)`
($b_{n}$); a second argument 1 selects the zeros of the derivatives.
:::

::: macro
BesselJ
BesselY
BesselI
BesselK
HankelH1
HankelH2

Bessel functions of the first and second kind `BesselJ(nu, z)`
($J_{\nu}\!\left(z\right)$) and `BesselY(nu, z)`
($Y_{\nu}\!\left(z\right)$), modified Bessel functions `BesselI(nu, z)`
($I_{\nu}\!\left(z\right)$) and `BesselK(nu, z)`
($K_{\nu}\!\left(z\right)$), with derivatives `BesselJ(nu, z, r)`
(${J}^{(r)}_{\nu}\!\left(z\right)$), and the Hankel functions
`HankelH1(nu, z)` ($H^{(1)}_{\nu}\!\left(z\right)$) and `HankelH2(nu,
z)` ($H^{(2)}_{\nu}\!\left(z\right)$).
:::

::: macro
BesselJZero
BesselYZero

The $n$-th positive zeros `BesselJZero(nu, n)` ($j_{\nu, n}$) and
`BesselYZero(nu, n)` ($y_{\nu, n}$) of $J_\nu$ and $Y_\nu$.
:::

::: macro
CoulombF
CoulombG
CoulombH
CoulombC
CoulombSigma

Coulomb wave functions: the regular and irregular functions
`CoulombF(ell, eta, z)` ($F_{\ell,\eta}(z)$) and `CoulombG(ell, eta, z)`
($G_{\ell,\eta}(z)$), the combination `CoulombH(omega, ell, eta, z)`
($H^{\omega}_{\ell,\eta}(z)$), equal to $G \pm i F$ for $\omega = \pm
1$, the normalising constant `CoulombC(ell, eta)`
($C_{\ell}\!\left(\eta\right)$) and the Coulomb phase shift
`CoulombSigma(ell, eta)` ($\sigma_{\ell}\!\left(\eta\right)$).
:::

::: macro
Hypergeometric0F1
Hypergeometric1F1
Hypergeometric1F2
Hypergeometric2F0
Hypergeometric2F1
Hypergeometric2F2
Hypergeometric3F2

Generalised hypergeometric functions ${}_pF_q$ with the parameters
followed by the argument, such as `Hypergeometric0F1(a, z)`
($\,{}_0F_1(a, z)$), `Hypergeometric1F1(a, b, z)` ($\,{}_1F_1(a, b,
z)$), `Hypergeometric2F1(a, b, c, z)` ($\,{}_2F_1(a, b, c, z)$) and
`Hypergeometric3F2(a, b, c, d, e, z)` ($\,{}_3F_2(a, b, c, d, e, z)$);
${}_2F_0$ is understood in the sense of an asymptotic series
(equivalently via $U$).
:::

::: macro
Hypergeometric0F1Regularized
Hypergeometric1F1Regularized
Hypergeometric1F2Regularized
Hypergeometric2F1Regularized
Hypergeometric2F2Regularized
Hypergeometric3F2Regularized

The regularised hypergeometric functions, with each lower parameter $b$
contributing a factor $1/\Gamma(b)$, such as
`Hypergeometric2F1Regularized(a, b, c, z)` ($\,{}_2{\textbf F}_1(a, b,
c, z)$).
:::

::: macro
HypergeometricU
HypergeometricUStar
HypergeometricUStarRemainder

The confluent hypergeometric function of the second kind
`HypergeometricU(a, b, z)` ($U(a, b, z)$), the scaled function
`HypergeometricUStar(a, b, z)` ($U^{*}(a, b, z)$), equal to $z^a U(a, b,
z)$, and the remainder of its asymptotic series after $n$ terms
`HypergeometricUStarRemainder(n, a, b, z)` ($R_{n}\!\left(a, b,
z\right)$).
:::

::: macro
RiemannZeta
RiemannZetaZero
RiemannXi
RiemannHypothesis

The Riemann zeta function `RiemannZeta(s)` ($\zeta(s)$) and its
derivatives `RiemannZeta(s, r)` (${\zeta}^{(r)}\!\left(s\right)$), its
$n$-th nontrivial zero in the upper half-plane `RiemannZetaZero(n)`
($\rho_{n}$), the Riemann xi function `RiemannXi(s)` ($\xi(s)$), and the
proposition `RiemannHypothesis` ($\operatorname{RH}$).
:::

::: macro
HurwitzZeta
LerchPhi
PolyLog
MultiZetaValue

The Hurwitz zeta function `HurwitzZeta(s, a)` ($\zeta(s, a)$), the Lerch
transcendent `LerchPhi(z, s, a)` ($\Phi(z, s, a)$), the polylogarithm
`PolyLog(s, z)` ($\operatorname{Li}_{s}\!\left(z\right)$) and multiple
zeta values `MultiZetaValue(s, t)` ($\zeta(s, t)$).
:::

::: macro
StieltjesGamma
KeiperLiLambda

The Stieltjes constants `StieltjesGamma(n)` ($\gamma_{n}$), the
generalised constants `StieltjesGamma(n, a)`
($\gamma_{n}\!\left(a\right)$) of the Hurwitz zeta function, and the
Keiper–Li coefficients `KeiperLiLambda(n)` ($\lambda_{n}$).
:::

::: macro
DirichletCharacter
DirichletGroup
PrimitiveDirichletCharacters
ConreyGenerator
GaussSum

Dirichlet characters in Conrey's labelling: the character
`DirichletCharacter(q, ell)` ($\chi_{q \, . \, \ell}$) and its value
`DirichletCharacter(q, ell, n)` ($\chi_{q \, . \, \ell}(n)$), the group
of characters modulo $q$ `DirichletGroup(q)` ($G_{q}$), the set of
primitive characters `PrimitiveDirichletCharacters(q)`
($G^{\text{Primitive}}_{q}$), the generator used by Conrey's labelling
`ConreyGenerator(q)` ($g_{q}$), and the Gauss sum `GaussSum(q, chi)`
($G_{q}\!\left(\chi\right)$).
:::

::: macro
DirichletL
DirichletLZero
DirichletLambda
GeneralizedBernoulliB
GeneralizedRiemannHypothesis

The Dirichlet L-function `DirichletL(s, chi)` ($L(s, \chi)$), its zeros
`DirichletLZero(n, chi)` ($\rho_{n, \chi}$), the completed L-function
`DirichletLambda(s, chi)` ($\Lambda(s, \chi)$), the generalised
Bernoulli numbers `GeneralizedBernoulliB(n, chi)` ($B_{n, \chi}$), and
the proposition `GeneralizedRiemannHypothesis` ($\operatorname{GRH}$).
:::

::: macro
AGM
AGMSequence

The arithmetic-geometric mean `AGM(a, b)` ($\operatorname{agm}(a, b)$)
and the $n$-th terms of the AGM iteration `AGMSequence(n, a, b)`
($\operatorname{agm}_{n}\!\left(a, b\right)$).
:::

::: macro
EllipticK
EllipticE
EllipticPi
IncompleteEllipticF
IncompleteEllipticE
IncompleteEllipticPi

Legendre's complete elliptic integrals with parameter $m = k^2$,
`EllipticK(m)` ($K(m)$), `EllipticE(m)` ($E(m)$), `EllipticPi(n, m)`
($\Pi(n, m)$), and the incomplete integrals `IncompleteEllipticF(phi,
m)` ($F(\phi, m)$), `IncompleteEllipticE(phi, m)` ($E(\phi, m)$) and
`IncompleteEllipticPi(n, phi, m)` ($\Pi(n, \phi, m)$).
:::

::: macro
CarlsonRF
CarlsonRG
CarlsonRJ
CarlsonRD
CarlsonRC
CarlsonHypergeometricR
CarlsonHypergeometricT

Carlson's symmetric elliptic integrals `CarlsonRF(x, y, z)` ($R_F(x, y,
z)$), `CarlsonRG(x, y, z)` ($R_G(x, y, z)$), `CarlsonRJ(x, y, z, p)`
($R_J(x, y, z, p)$), `CarlsonRD(x, y, z)` ($R_D(x, y, z)$) and
`CarlsonRC(x, y)` ($R_C(x, y)$), and Carlson's multivariate
hypergeometric functions `CarlsonHypergeometricR` and
`CarlsonHypergeometricT`.
:::

::: macro
JacobiTheta
JacobiThetaQ

Jacobi theta functions $\theta_j$, $j = 1, 2, 3, 4$, in terms of the
lattice parameter $\tau$, `JacobiTheta(j, z, tau)`
($\theta_{j}\!\left(z, \tau\right)$), with derivatives `JacobiTheta(j,
z, tau, r)` (${\theta}^{(r)}_{j}\!\left(z, \tau\right)$), or in terms of
the nome $q$, `JacobiThetaQ(j, z, q)` ($\theta_{j}\!\left(z, q\right)$).
:::

::: macro
DedekindEta
ModularJ
ModularLambda
EisensteinG
EisensteinE
EulerQSeries

Modular forms and functions: the Dedekind eta function
`DedekindEta(tau)` ($\eta(\tau)$), the $j$-invariant `ModularJ(tau)`
($j(\tau)$), the modular lambda function `ModularLambda(tau)`
($\lambda(\tau)$), the Eisenstein series `EisensteinG(k, tau)`
($G_{k}\!\left(\tau\right)$) and their normalised versions
`EisensteinE(k, tau)` ($E_{k}\!\left(\tau\right)$), and the Euler
function `EulerQSeries(q)` ($\phi(q)$), equal to $ \prod_{k \ge 1} (1 -
q^k)$.
:::

::: macro
DedekindSum
DedekindEtaEpsilon
JacobiThetaEpsilon
JacobiThetaPermutation

Quantities in modular transformation formulas: the Dedekind sum
`DedekindSum(h, k)` ($s(h, k)$), the root of unity
`DedekindEtaEpsilon(a, b, c, d)` ($\varepsilon(a, b, c, d)$) in the
transformation of $\eta$ under a matrix $\left(\begin{smallmatrix} a & b
\\ c & d \end{smallmatrix}\right)$, and the root of unity
`JacobiThetaEpsilon(j, a, b, c, d)` ($\varepsilon_{j}\!\left(a, b, c,
d\right)$) and index permutation `JacobiThetaPermutation(j, a, b, c, d)`
($S_{j}\!\left(a, b, c, d\right)$) in the transformation of $\theta_j$.
:::

::: macro
ModularGroupAction
ModularGroupFundamentalDomain
ModularLambdaFundamentalDomain

The action of a matrix on the upper half-plane
`ModularGroupAction(gamma, tau)` ($\gamma \circ \tau$), and the
fundamental domains `ModularGroupFundamentalDomain` ($\mathcal{F}$) of
the modular group and `ModularLambdaFundamentalDomain`
($\mathcal{F}_{\lambda}$) of the modular lambda function.
:::

::: macro
WeierstrassP
WeierstrassZeta
WeierstrassSigma
EllipticRootE

The Weierstrass elliptic functions for the lattice with periods 1 and
$\tau$, `WeierstrassP(z, tau)` ($\wp(z, \tau)$), `WeierstrassZeta(z,
tau)` ($\zeta(z, \tau)$) and `WeierstrassSigma(z, tau)` ($\sigma(z,
\tau)$), and the roots `EllipticRootE(j, tau)`
($e_{j}\!\left(\tau\right)$) of the cubic polynomial associated with
$\wp$.
:::

::: macro
HilbertClassPolynomial
PrimitiveReducedPositiveIntegralBinaryQuadraticForms

The Hilbert class polynomial `HilbertClassPolynomial(D, x)`
($H_{D}\!\left(x\right)$) and the set of primitive reduced positive
definite binary quadratic forms of discriminant $D$,
`PrimitiveReducedPositiveIntegralBinaryQuadraticForms(D)`
($\mathcal{Q}^{*}_{D}$).
:::

## Display forms {#display}

Symbols that only affect LaTeX rendering. They carry no meaning beyond
their argument.

::: macro
Parentheses
Brackets
Braces
AngleBrackets

`Parentheses(x)`, `Brackets(x)`, `Braces(x)` and `AngleBrackets(x)`
denote `x` but render it inside delimiters: `Parentheses(x)`
($\left(x\right)$), `Brackets(x)` ($\left[x\right]$), `Braces(x)`
($\left\{x\right\}$), `AngleBrackets(x)` ($\left\langle
x\right\rangle$).
:::

::: macro
Ellipsis

`Ellipsis` ($\ldots$) renders as an ellipsis, for example to indicate
omitted arguments; it has no builtin meaning.
:::

::: macro
Subscript

`Subscript(x, n)` (${x}_{n}$) renders `x` with a subscript.
:::

::: macro
Logic

`Logic(x)` denotes `x`, but logical connectives within `x` are rendered
with symbols instead of words, as with `FEXPR_LATEX_LOGIC`:
`Logic(And(P, Not(Q)))` ($P \,\land\, \neg Q$).
:::

::: macro
ShowExpandedNormalForm

`ShowExpandedNormalForm(x)` denotes `x`, but renders its expanded normal
form (see `fexpr_expanded_normal_form`) instead of `x` itself. This
triggers a computation that can be expensive.
:::

## Variables {#greek}

Builtin symbols for variables with special LaTeX rendering. Single
letters and other symbols starting with a lowercase letter are free for
use as variables; the following lowercase names are nevertheless builtin
so that they render as Greek letters. A trailing underscore gives a
subscripted variable: `alpha_(n)` renders as $\alpha_n$.

::: macro
alpha
beta
gamma
delta
epsilon
zeta
eta
theta
vartheta
iota
kappa
lamda
mu
nu
xi
pi
rho
sigma
tau
phi
varphi
chi
omega
ell

Lowercase Greek letters and $\ell$, such as `alpha` ($\alpha$),
`epsilon` ($\varepsilon$), `lamda` ($\lambda$), spelled without b, and
`ell` ($\ell$).
:::

::: macro
alpha_
beta_
gamma_
delta_
epsilon_
zeta_
eta_
theta_
vartheta_
iota_
kappa_
lamda_
mu_
nu_
xi_
pi_
rho_
sigma_
tau_
phi_
varphi_
chi_
omega_
ell_

Subscripted lowercase Greek letters, such as `alpha_(n)` ($\alpha_{n}$)
or `ell_(i, j)` ($\ell_{i, j}$).
:::

::: macro
GreekGamma
Delta
Theta
Lamda
Xi
GreekPi
Sigma
Phi
Psi
Omega

Uppercase Greek letters, such as `GreekGamma` ($\Gamma$), `Delta`
($\Delta$) and `Omega` ($\Omega$). `GreekGamma` and `GreekPi` are named
so as not to clash with the functions `Gamma` and the constant `Pi`.
:::

::: macro
GreekGamma_
Delta_
Theta_
Lamda_
Xi_
GreekPi_
Sigma_
Phi_
Psi_
Omega_

Subscripted uppercase Greek letters, such as `Delta_(n)` ($\Delta_{n}$).
:::

::: macro
XX

`XX(n)` ($\text{x}_{n}$) renders as a subscripted upright x, for indexed
variables $\text{x}_1, \text{x}_2, \ldots$.
:::
