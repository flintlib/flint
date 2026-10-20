---
title: Zero testing
summary: Why deciding x = 0 is hard for exact real and complex numbers, what is known about decidability, and how ca combines numerical evaluation in ball arithmetic (which proves x ≠ 0) with exact methods (formal reduction, evaluation as an algebraic number, rewriting in a normal form, factoring) that prove x = 0.
complexity: Proving $x \ne 0$ costs evaluations at precisions $64, 128, \ldots$ up to about $-\log_2 |x|$ bits (bounded by `CA_OPT_PREC_LIMIT`); proving $x = 0$ costs at least one evaluation at the maximum precision plus the exact methods, of which the algebraic evaluation can be expensive for numbers of high degree.
order: 40
---

The zero test is the core of exact real and complex arithmetic. Every
decision in [](ca), from comparisons and divisions to degrees of
polynomials and ranks of matrices, reduces to deciding whether some
number $x$ equals zero. The function page [](ca/check-is-zero) gives a
user-level summary and examples.

## Decidability {#decidability}

For algebraic numbers the problem is decidable: a number given by an
expression with rational numbers, arithmetic operations and roots can
be evaluated to its minimal polynomial and an isolating ball (see
[](algorithms/exact/qqbar-arithmetic)), which is zero if and only if
the minimal polynomial is $t$. The cost can, however, be exponential in
the size of the expression, since degrees multiply.

With transcendental functions the situation is different. Richardson
[@Ric1968] showed that it is undecidable whether an expression built
from a variable $x$, rational numbers, $\pi$, $\log 2$, $+$, $\times$,
$\sin$ and $|\cdot|$ vanishes identically on $\mathbb{R}$. For
*constants*, the *constant problem* is open: it is not known whether
equality of numbers built from rational numbers with arithmetic, $\exp$
and $\log$ is decidable. Richardson [@Ric1997; @Ric2009] gave an
algorithm for such exp-log constants that is correct whenever it
terminates and that terminates if Schanuel's conjecture holds (see
[](algorithms/exact/field-towers#completeness)); see also van der
Hoeven [@vdH2006]. Even for specific numbers such as $e + \pi$, it is
not known whether they are rational.

In practice, the numbers that arise in a computation usually have a
structure that can be exploited, and a semi-decision procedure that
proves inequality numerically and equality by exact methods, and that
answers "unknown" when neither succeeds, is useful. This is the approach
of [](ca) [@Joh2021calcium].

## Proving inequality {#numerical}

If $x \ne 0$, then evaluating $x$ in ball arithmetic at a sufficiently
high precision gives a ball that excludes zero: the radius of the
computed ball tends to zero as the precision grows, for the operations
that occur in `ca_t` numbers, as long as no function is evaluated on a
branch cut or at a singularity. A ball that excludes zero is a proof.
`ca_get_acb_raw` evaluates the representation of $x$ (a rational function
of the generators) from enclosures of the generators, which are computed
recursively and cached in the extension numbers.

The required precision is about $-\log_2 |x|$ plus the precision lost in
the evaluation, so numbers that are extremely close to zero without
being zero are expensive; the search is bounded by
`CA_OPT_PREC_LIMIT`. Ramanujan's constant $e^{\pi \sqrt{163}} -
640320^3 - 744 \approx -7.5 \cdot 10^{-13}$ is separated from zero at
128 bits.

Numerical evaluation alone can never prove $x = 0$: a ball containing
zero is compatible with both answers. A number that cannot be separated
from zero is either zero or smaller than the precision limit allows to
detect.

## Proving equality {#exact}

`ca_check_is_zero` proves $x = 0$ by the following methods, in order of
increasing cost.

::: algorithm Zero test {#alg-zero}
**Input:** a `ca_t` value $x$.

**Output:** `T_TRUE`, `T_FALSE` or `T_UNKNOWN`.

1. *Representation.* If $x$ is a special value, return `T_FALSE`
   (`T_UNKNOWN` for *Unknown*). If $x \in \mathbb{Q}$ or $x$ is an element
   of a field with one algebraic generator, return whether its canonical
   representation is zero.
2. Otherwise $x = p/q$ in a field $\mathbb{Q}(a_1, \ldots, a_n)$; replace
   $x$ by the numerator $p$, which is zero if and only if $x$ is.
3. *Numerical evaluation.* For precisions $64, 128, \ldots$ up to
   `CA_OPT_PREC_LIMIT`: evaluate $p$ in ball arithmetic and return
   `T_FALSE` if the ball excludes zero. At the first precision, if every
   generator can be evaluated as an algebraic number, evaluate $p$ exactly
   with [](qqbar) (`ca_get_qqbar`) and return the answer, unless the
   degree limit `CA_OPT_QQBAR_DEG_LIMIT` is exceeded.
4. *Normal form.* Rewrite the generators recursively in a normal form
   (see below) and evaluate $p$ on the rewritten generators; if the
   result is represented as zero, return `T_TRUE`.
5. *Factoring.* Factor $p$ as a polynomial in the generators; if it has
   at least two nonconstant factors, apply steps 1 to 4 to each factor
   and return `T_TRUE` if one of them is proved zero.
6. Return `T_UNKNOWN`.
:::

::: theorem Correctness {#thm-zero}
If [](#alg-zero) returns `T_TRUE`, then $x = 0$; if it returns
`T_FALSE`, then $x \ne 0$.
:::

::: proof
`T_FALSE` is returned only for a special value or when a ball computed
in ball arithmetic, which is guaranteed to contain $p(a_1, \ldots,
a_n)$, excludes zero, or when the canonical representation in
$\mathbb{Q}$, in a number field, or of an exact algebraic evaluation is
nonzero. `T_TRUE` is returned when an exact representation is zero. The
representations of steps 1, 3 and 4 are elements of formal fields
reduced modulo ideals of proved relations, so by
[](algorithms/exact/field-towers#prop-sound) a zero representation
implies a zero value. The rewriting of step 4 replaces each generator by
an expression with the same value, so the rewritten $p$ has the value of
$p$. In step 5, $p$ is the product of its factors, so it vanishes if one
of them does.
:::

The normal form of step 4 (`ca_rewrite_complex_normal_form`, internal)
makes relations visible that the original generators hide:

- $\sin z$, $\cos z$, $\tan z$ become rational functions of $e^{iz}$, and
  $\operatorname{atan} z$, $\operatorname{asin} z$ and
  $\operatorname{acos} z$ become logarithms;
- $\operatorname{Re} z$, $\operatorname{Im} z$, $\overline{z}$, $|z|$,
  $\operatorname{sgn} z$ and $\arg z$ are expressed with $z$ and its deep
  conjugate $\overline{z}$, with $|z| = \sqrt{z \overline{z}}$ and $\arg z
  = -i \log(z/|z|)$;
- algebraic generators that lie in a cyclotomic field are expressed with
  a root of unity (`qqbar_try_as_cyclotomic`);
- the arguments of $\exp$, $\log$ and $\sqrt{\cdot}$ are rewritten
  recursively.

Creating the fields of the rewritten expression triggers the relation
searches of [](algorithms/exact/field-towers#ideal), which find
multiplicative relations between exponentials and linear relations
between logarithms. In this way, for example, $\cos^2 x + \sin^2 x - 1$
is proved zero even when the trigonometric functions are represented
directly (`CA_TRIG_DIRECT`), and $16 \operatorname{atan}(1/5) - 4
\operatorname{atan}(1/239) - \pi$ (Machin's formula) is proved zero.

## What is decided in practice {#in-practice}

With the methods above, the following questions are always decided,
within the limits set in the context (see [](ca/options)).

- **Inequality of distinct numbers**, given enough precision
  ([](#numerical)). The precision is bounded by the context option `CA_OPT_PREC_LIMIT` (4096 bits by default), so a
  nonzero number smaller than about $2^{-4096}$ in absolute value is not
  recognised as nonzero by default.
- **Equality of algebraic numbers.** A number all of whose generators
  are algebraic can be evaluated exactly as a [](qqbar) number, whose
  minimal polynomial decides equality. This is limited by
  `CA_OPT_QQBAR_DEG_LIMIT` (degree 120 by default). With [](qqbar)
  itself, all predicates are decided without limit, at a cost that grows
  quickly with the degree.
- **Elements of one field.** If $x$ and $y$ belong to the same number
  field, or more generally if their difference is reduced to zero by the
  known relations of their field, equality is recognised formally.

For numbers involving transcendental functions, equality is decided
only sometimes: `ca_check_is_zero` proves $x = 0$ when it finds the
relations that make $x$ vanish, that is, algebraic relations, linear
relations between logarithms, multiplicative relations between
exponentials, the functional equations of some special functions, and
the rewriting of trigonometric functions as exponentials. This succeeds
for many identities that occur in practice, such as
$\sqrt{5 + 2\sqrt{6}} = \sqrt{2} + \sqrt{3}$,
$\log(\sqrt 2 + 1) + \log(\sqrt 2 - 1) = 0$ or Machin's formula
$16 \operatorname{atan}(1/5) - 4 \operatorname{atan}(1/239) = \pi$, and
fails for identities that need deeper knowledge, such as
$\Gamma(1/3) \Gamma(2/3) = 2\pi/\sqrt{3}$. [](ca/check-is-zero) shows these cases in a program and explains what
to do when the answer is `T_UNKNOWN`.

## Consequences for the API {#api}

Since a zero test can fail, every function whose result depends on one
reports this: scalar operations return *Unknown* (for example $1/x$
when $x$ cannot be decided), predicates return `T_UNKNOWN`, and
functions on polynomials and matrices return a failure flag or
`T_UNKNOWN` (for example `ca_poly_gcd`, `ca_mat_rank`, `ca_mat_inv`).
Comparisons $x < y$ of real numbers use the same scheme: they compare
enclosures at increasing precision and, when $x$ and $y$ cannot be
separated, test $x = y$ exactly. Functions that need many zero tests,
such as Gaussian elimination, choose the operands to be tested in an
order that makes success likely (see
[](algorithms/exact/linear-algebra)).

## Use in FLINT {#use}

- `ca_check_is_zero` implements [](#alg-zero); `ca_is_zero_check_fast`
  is its first step, and `ca_check_equal`, `ca_check_is_one` and the
  other predicates are built on it.
- `ca_get_acb_raw` and `ca_ext_get_acb_raw` evaluate numbers and
  generators in ball arithmetic; `ca_get_acb` adds precision control and
  uses the zero test to return exact zeros.
- `ca_get_qqbar` performs the exact algebraic evaluation, and
  `ca_factor` the factorisation of step 5.
- [](qqbar) decides equality of algebraic numbers completely; see
  [](algorithms/exact/qqbar-representation#predicates).
