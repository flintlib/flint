---
title: Precision and representation of p-adic numbers
summary: How FLINT represents a p-adic number as p^v u modulo p^N, the difference between absolute and relative precision, how much precision each arithmetic operation preserves, and the two precision models used by padic/qadic and by padic_radix.
complexity: Representation only; reduction modulo $p^{N - v}$ costs one division of an integer of $O((N - v)\log p)$ bits.
order: 10
---

A $p$-adic number can only be stored approximately, as a residue class
modulo a power of $p$. FLINT uses two models: the fixed-precision model of `padic`, `qadic`,
`padic_poly` and `padic_mat`, where every variable has a precision and
inputs are treated as exact, and the ball model of `padic_radix`, where
every element carries its own error term. [@Car2017] compares these and
other models in detail.

## Valuation, absolute and relative precision {#notation}

Every nonzero $x \in \mathbf{Q}_p$ can be written uniquely as
$x = p^v u$ with $v = \operatorname{ord}_p(x) \in \Z$ and $u$ a $p$-adic
unit, that is, $u \in \mathbf{Z}_p^\times$. We write $x + O(p^N)$ for
the set $x + p^N \mathbf{Z}_p$, the ball of radius $p^{-N}$ around $x$.
An approximation $x + O(p^N)$ has *absolute precision* $N$; if
$\operatorname{ord}_p(x) = v < N$ it has *relative precision* $N - v$,
the number of known $p$-adic digits of the unit. If $v \ge N$ the
approximation is indistinguishable from zero.

In FLINT the unit is an integer $u$, so the approximation is
$p^v u + O(p^N)$ with $u \in \Z$, $p \nmid u$, and $u$ reduced modulo
$p^{N - v}$. Storing $u$ (an `fmpz`, or a `radix_integer_t` in
`padic_radix`) and $v$ separately has two advantages over storing the
integer $p^v u$: elements of negative valuation need no special
treatment, and the size of the stored integer is proportional to the
relative precision, not the absolute precision.

## Precision of arithmetic operations {#propagation}

The following statements describe what is determined about the result
when the inputs are known only modulo powers of $p$. They are the rules
that `padic_radix` applies to its error terms, and the rules a user of
`padic` must follow to know which digits of a result are meaningful.

::: proposition Precision of ring operations {#prop-ring}
Let $\operatorname{ord}_p(x) = v$, $\operatorname{ord}_p(y) = w$, and let
$\delta, \varepsilon$ have valuations at least $a$ and $b$, with
$a > v$ and $b > w$. Then

1. $(x + \delta) \pm (y + \varepsilon) \in x \pm y + O(p^{\min(a, b)})$;
2. $(x + \delta)(y + \varepsilon) \in xy + O(p^{\min(a + w, b + v)})$;
3. $1/(x + \delta) \in 1/x + O(p^{a - 2v})$;
4. $(x + \delta)^e \in x^e + O(p^{a + (e - 1)v})$ for $e \ge 1$.

In each case the precision on the right cannot be improved in general.
:::

::: proof
The first statement is the ultrametric inequality. For the second,
$(x + \delta)(y + \varepsilon) - xy = x\varepsilon + y\delta +
\delta\varepsilon$, whose terms have valuations at least $v + b$,
$w + a$ and $a + b > \max(v + b, w + a)$. For the third,
$1/(x + \delta) - 1/x = -\delta / (x(x + \delta))$ and
$\operatorname{ord}_p(x + \delta) = v$ since $a > v$. The fourth follows
from $(x + \delta)^e - x^e = \sum_{i \ge 1} \binom{e}{i} x^{e - i}
\delta^i$, whose terms have valuation at least $(e - i)v + ia \ge
(e - 1)v + a$. Taking $\delta = p^a$ and $\varepsilon = p^b$ shows that
the bounds are attained.
:::

In terms of relative precision the statements are simpler: a product or
quotient of elements known to relative precisions $r$ and $s$ is known
to relative precision $\min(r, s)$, a power to the relative precision of
its base, and a sum to an absolute precision that may correspond to a
much smaller relative precision if the leading digits cancel.

::: proposition Precision of roots and transcendental functions {#prop-functions}
Let $u$ be a unit of $\mathbf{Z}_p$ and $\operatorname{ord}_p(\delta)
\ge a \ge 1$.

1. If $p$ is odd and $u = s^2$, then $u + \delta$ has a square root in
   $s + O(p^a)$. If $p = 2$, $a \ge 3$ and $u = s^2$, then $u + \delta$
   has a square root in $s + O(2^{a - 1})$.
2. If $\exp(x)$ converges, then $\exp(x + \delta) \in \exp(x)(1 +
   O(p^a))$, provided $\exp(\delta)$ converges.
3. If $\log(u)$ converges, then $\log(u + \delta) \in \log(u) + O(p^a)$.
:::

::: proof
For the first statement, $(u + \delta)/u = 1 + \delta/u$ and
$\operatorname{ord}_p(\delta/u) = a$. For odd $p$ the binomial series
$\sqrt{1 + t} = \sum \binom{1/2}{k} t^k$ has coefficients in
$\mathbf{Z}_p$, so $\sqrt{1 + \delta/u} \in 1 + O(p^a)$; for $p = 2$ the
coefficients $\binom{1/2}{k}$ have 2-adic valuation at least $-2k + 1$
for $k \ge 1$, so for $a \ge 3$ the series converges and lies in
$1 + O(2^{a - 1})$. Multiplying by $s$ gives the claim. The second
statement is the functional equation $\exp(x + \delta) = \exp(x)
\exp(\delta)$ together with $\exp(\delta) \in 1 + O(p^a)$
([](algorithms/p-adic/exp-log#lem-exp-terms)). For the third,
$\log(u + \delta) = \log(u) + \log(1 + \delta/u)$, and the logarithm
of $1 + t$ with $\operatorname{ord}_p(t) \ge a$ has valuation at least
$a$, because $\operatorname{ord}_p(t^i / i) \ge ia - \log_p i \ge a$.
:::

## The fixed-precision model of padic and qadic {#fixed}

In `padic`, every variable has a precision $N$ that is set when it is
initialised. The value $p^v u$ is kept *reduced*: $u$ is a nonnegative
integer coprime to $p$ and smaller than $p^{N - v}$, or $u = v = 0$.
A function computes the exact result of the operation applied to the
stored values, which are elements of $\Z[1/p]$ (or, for a square root or
logarithm, of $\mathbf{Q}_p$), and reduces it modulo $p^N$ for the
precision $N$ of the output. This is the analogue of correct rounding in
MPFR, and it makes every function deterministic and every result as
precise as the output variable allows. It also means that the results
are exact functions of the stored representatives, so that, for
example, `padic_sub` of two variables that agree modulo $p^{20}$ returns
their exact difference with as many digits as the output can hold. The
precision rules above are not applied; the caller must apply them to
know which digits are meaningful.

The same model is used by `qadic` for elements $p^v g(X)$ of
$\mathbf{Q}_q$, where $g$ is a polynomial whose coefficients are
reduced modulo $p^{N - v}$ and not all divisible by $p$, and by
`padic_poly` and `padic_mat`, where one valuation and one precision are
shared by all coefficients. In these structures a coefficient of
valuation $w > v$ is stored to relative precision $N - w$, and
coefficients of valuation at least $N$ are zero.

::: summary algorithms/p-adic/overview :::

## The ball model of padic_radix {#balls}

In `padic_radix`, an element is a ball $p^v u + O(p^N)$ with its own
absolute precision $N$, or an exact element without error term. The
context fixes a default relative precision $r$ and absolute precision
$a$, either of which may be infinite. The functions propagate the error
terms by the rules of [](#prop-ring) and [](#prop-functions): a sum has
precision $\min(N_x, N_y)$, a product $\min(N_x + v_y, N_y + v_x)$, a
quotient the smaller relative precision of its operands, and a square
root for $p = 2$ one digit less than its operand. The result is then
truncated to precision $\min(a, v + r)$ for its valuation $v$, and the
truncation is recorded in the error term. An exact result that fits
within this precision is kept exact, so that, with both precisions
infinite, the arithmetic is exact arithmetic in $\Z[1/p]$.

Applying these rules operation by operation always gives valid error
terms, but not always the smallest ones: the guaranteed precision of a
result depends on the algorithm. For sums of products and block
products, see [](algorithms/p-adic/block-products); for linear algebra,
where full pivoting on valuations and the differential of the function
give the optimal precision, see [](algorithms/p-adic/linear-algebra).

Predicates follow the conventions of the generic rings: `padic_radix_equal`
returns `T_FALSE` if the balls are disjoint, `T_TRUE` only if both
elements are exact and equal, and `T_UNKNOWN` otherwise; division by an
element that may be zero returns `GR_UNABLE`.

## Use in FLINT {#use}

- `padic_reduce`, `_padic_canonicalise` and `_padic_reduce` implement
  canonicalisation and reduction for `padic_t`; `padic_poly_reduce`,
  `padic_poly_canonicalise`, `padic_mat_reduce` and `qadic_reduce` do
  the same for the other structures.
- `padic_pow_si` and `padic_poly_pow` state the precision loss of
  [](#prop-ring) for negative valuations but do not apply it.
- `padic_radix_add`, `padic_radix_mul`, `padic_radix_div`,
  `padic_radix_inv`, `padic_radix_sqrt`, `padic_radix_exp` and
  `padic_radix_log` apply the rules above to their error terms;
  `padic_radix_dot` applies them to whole sums of products
  ([](algorithms/p-adic/block-products#dot)).
