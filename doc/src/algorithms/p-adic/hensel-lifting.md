---
title: Hensel lifting of inverses, square roots and Teichmüller lifts
summary: Newton iteration in Z_p and Z_q with the precision sequence N, ⌈N/2⌉, ..., 1, applied to inverses, inverse square roots (with the final step producing the square root) and Teichmüller lifts, including the modifications needed for p = 2.
complexity: Inverse and square root to precision $N$: $O(M(N \log p))$ bit operations in $\mathbf{Z}_p$ and $O(M(d) M(N \log p))$ in $\mathbf{Z}_q$, plus the cost at precision 1; Teichmüller lift: $O(\log q)$ multiplications at each precision, $O(\log q \, M(N \log p))$ in total for $\mathbf{Z}_p$.
order: 20
---

Many $p$-adic quantities are roots of an equation whose reduction
modulo $p$ is easy to solve: the inverse of a unit, a square root, a
root of unity. Hensel's lemma [@Hen1908] lifts a simple root modulo $p$
to a root modulo $p^N$, and in the Newton form used by FLINT each step
doubles the number of correct digits, so the whole computation costs a
constant number of multiplications at the target precision. The same
pattern appears for power series in [](algorithms/polynomials/newton);
textbook treatments are [@GG2003, §9] and [@BZ2011, §4.2].

Throughout, $R$ is $\mathbf{Z}_p$ or the ring of integers
$\mathbf{Z}_q$ of an unramified extension of degree $d$
([](algorithms/p-adic/unramified-extensions)); in both rings $p$
generates the maximal ideal and $R/pR$ is a field. $M(n)$ is the cost of
multiplying $n$-bit integers.

## Newton's method in R {#newton}

::: lemma Quadratic convergence {#lem-newton}
Let $f \in R[t]$, let $k \ge 1$ and let $t, s \in R$ with
$f(t) \equiv 0 \pmod{p^k}$ and $s f'(t) \equiv 1 \pmod{p^k}$. Then
$t' = t - s f(t)$ satisfies $f(t') \equiv 0 \pmod{p^{2k}}$, and
$t' \equiv t \pmod{p^k}$.
:::

::: proof
By Taylor expansion, $f(t + h) = f(t) + f'(t) h + h^2 g(t, h)$ with
$g \in R[t, h]$. With $h = -s f(t) \in p^k R$ we get
$f(t') = f(t)(1 - s f'(t)) + h^2 g(t, h)$, and both terms are divisible
by $p^{2k}$.
:::

In particular, if $f'(t)$ is a unit at a root $t_1$ of $f$ modulo $p$,
there is a unique root of $f$ in $R$ congruent to $t_1$ modulo $p$
(Hensel's lemma), and Newton's iteration converges to it. The inverse
$s$ need only be known to the current precision, so it can be updated by
its own Newton step $s \gets s(2 - s f'(t))$ instead of being recomputed.

### The precision sequence {#sequence}

FLINT runs the iteration from the target precision downwards. For a
target $N \ge 1$, `_padic_lifts_exps` computes

$$
a_0 = N, \quad a_{i+1} = \lceil a_i / 2 \rceil, \quad \ldots, \quad
a_{n - 1} = 1, \qquad n = \lceil \log_2 N \rceil + 1,
$$

and `_padic_lifts_pows` computes the moduli $p^{a_i}$ with about one
multiplication each, since $p^{a_i}$ is the square of $p^{a_{i+1}}$,
possibly divided by $p$. The input is reduced modulo every $p^{a_i}$
once, from the top. Then the root is found modulo $p^{a_{n-1}} = p$ and
lifted through $a_{n-2}, \ldots, a_0$ by [](#lem-newton), which applies
because $a_i \le 2 a_{i+1}$. Lifting to $\lceil N/2 \rceil$ and then to
$N$, rather than doubling from 1, never computes more digits than
needed.

::: proposition {#prop-cost}
If one step at precision $a$ costs $c \, M(a \log p)$ and $M(n)/n$ is
nondecreasing, the whole iteration costs at most
$2c \, M(N \log p) + O(c \log N)$ plus the cost of the base case modulo
$p$.
:::

::: proof
Since $a_{i+1} \le (a_i + 1)/2$, the sum of $M(a_i \log p)$ over
$i \ge 1$ is bounded by a geometric series with ratio $1/2$ up to the
rounding terms.
:::

## Inverses {#inverse}

For a unit $u$ the inverse is the root of $f(z) = u - 1/z$; the Newton
step simplifies to

$$
z' = z + z(1 - uz) = 2z - uz^2,
\qquad 1 - uz' = (1 - uz)^2 ,
$$

so the number of correct digits doubles exactly, as the identity on the
right shows without appeal to [](#lem-newton). `_padic_inv_precomp`
computes $z \equiv u^{-1}$ modulo $p$ by the extended Euclidean
algorithm (`fmpz_invmod`) and applies the step modulo $p^{a_{n-2}},
\ldots, p^{a_0}$, using at each level the reduction of $u$ modulo
$p^{a_i}$. Each step costs two multiplications and a reduction.
`_padic_inv` uses `fmpz_invmod` alone for $N = 1$ and, for $p = 2$,
`fmpz_invmod_2exp`, which runs the same iteration on the limbs of $u$
without divisions. In $\mathbf{Z}_q$, `_qadic_inv` uses the same step
with products reduced modulo the defining polynomial, and computes the
inverse modulo $p$ by the extended Euclidean algorithm in
$\mathbf{F}_p[X]$.

## Square roots {#sqrt}

Let $u$ be a unit. FLINT computes the inverse square root $z = u^{-1/2}$
by Newton's method for $f(z) = u z^2 - 1$ and obtains the square root at
the last step, a variant of the Karp–Markstein trick. Write
$e = uz^2 - 1$.

::: lemma {#lem-rsqrt}
Let $p$ be odd. If $e \equiv 0 \pmod{p^k}$ and $z' = z - z e / 2$, then
$u z'^2 \equiv 1 \pmod{p^{2k}}$. If $p = 2$, $k \ge 3$ and
$e \equiv 0 \pmod{2^k}$, then $u z'^2 \equiv 1 \pmod{2^{2k - 2}}$.
:::

::: proof
We have $u z'^2 = (1 + e)(1 - e/2)^2 = 1 - \tfrac34 e^2 + \tfrac14 e^3$.
For odd $p$, $2$ and $4$ are units, and $e^2$ is divisible by $p^{2k}$.
For $p = 2$, $e^2/4$ and $e^3/4$ are divisible by $2^{2k - 2}$.
:::

::: lemma Final step {#lem-sqrt-final}
Let $p$ be odd, $s^2 = u$, $z \equiv s^{-1} \pmod{p^k}$, and
$b \equiv uz \pmod{p^k}$. Then $r = b + z(u - b^2)/2$ satisfies
$r \equiv s \pmod{p^{2k}}$ (for the appropriate sign of $s$).
:::

::: proof
Write $b = s + \beta$ and $z = s^{-1} + \eta$ with $\beta, \eta \in
p^k R$. Then $u - b^2 = -2s\beta - \beta^2$ and
$r = s + \beta - (s^{-1} + \eta)(2s\beta + \beta^2)/2 = s - \beta^2 /
(2s) - \eta s \beta - \eta\beta^2/2$, and all three correction terms are
divisible by $p^{2k}$.
:::

`_padic_sqrt` for odd $p$ computes a square root modulo $p$ with
`fmpz_sqrtmod` (which decides whether $u$ is a square: by Hensel's lemma
a unit is a square in $\mathbf{Z}_p$ if and only if it is a square
modulo $p$), inverts it, lifts the inverse square root by
[](#lem-rsqrt) through $a_{n-2}, \ldots, a_1$, and applies
[](#lem-sqrt-final) from $a_1$ to $a_0 = N$. The last step replaces an
inverse square root step and a multiplication by one step of about the
same cost.

For $p = 2$ the derivative $2uz$ of $f$ is not a unit, and the precision
grows from $k$ to $2k - 2$ only. A unit is a square in $\mathbf{Z}_2$
if and only if it is $\equiv 1 \pmod 8$; `_padic_sqrt` tests this and
then calls `fmpz_sqrtmod_2exp`. In $\mathbf{Z}_q$ with $p = 2$,
`qadic_sqrt` lifts the inverse square root with the sequence
$e_0 = N$, $e_{i+1} = \lfloor e_i/2 \rfloor + 1$ (stopping at 2), which
satisfies $2e_{i+1} - 2 \ge e_i$. Its starting value is an inverse
square root modulo 2, which always exists since $x \mapsto x^2$ is an
automorphism of $\mathbf{F}_q$, corrected by a root of an Artin–Schreier
equation $Y^2 + Y = c$ over $\mathbf{F}_q$. That equation is solvable if
and only if $u$ is a square, and is solved by linear algebra over
$\mathbf{F}_2$ ([](algorithms/p-adic/unramified-extensions#sqrt)).

## Teichmüller lifts {#teichmuller}

Let $q = |R/pR|$. Every unit $x \in R$ has a unique *Teichmüller lift*
$\omega(x)$: a $(q - 1)$-th root of unity congruent to $x$ modulo $p$.
It is the root of $f(t) = t^q - t$ congruent to $x$, which exists and is
unique by Hensel's lemma since $f'(t) = q t^{q - 1} - 1 \equiv -1
\pmod p$.

::: lemma {#lem-teich}
If $t \equiv \omega(x) \pmod{p^k}$ and $s \equiv (q - 1)^{-1}
\pmod{p^k}$, then $t' = t - s(t^q - t) \equiv \omega(x) \pmod{p^{2k}}$.
:::

::: proof
With $\omega = \omega(x)$ and $t = \omega + h$, $h \in p^k R$, we have
$f'(t) = q t^{q-1} - 1 \equiv q \omega^{q-1} - 1 = q - 1
\pmod{p^{k + 1}}$, because $q \equiv 0 \pmod p$ and
$t^{q - 1} \equiv \omega^{q - 1} \pmod{p^k}$. Hence $s f'(t) \equiv 1
\pmod{p^k}$ and [](#lem-newton) applies.
:::

`_padic_teichmuller` ($q = p$) and `_qadic_teichmuller` ($q = p^d$)
start from $t = x \bmod p$ and $s = p - 1 \equiv (q - 1)^{-1} \pmod p$,
and at each precision of the sequence apply [](#lem-teich) and then a
Newton step $s \gets 2s - (q - 1)s^2$ for the inverse. The power $t^q$
dominates the cost: $O(\log q)$ multiplications at each precision. For
$p = 2$ the only odd root of unity in $\mathbf{Z}_2$ is 1, and
`_padic_teichmuller` returns it directly.

## Use in FLINT {#use}

- Inverses: `_padic_inv`, `_padic_inv_precomp` (used repeatedly in
  `padic_log_balanced`), `padic_inv`, `padic_div`, `padic_poly_inv_series`
  (for the constant coefficient), `_qadic_inv`, `qadic_inv`.
- Square roots: `_padic_sqrt`, `padic_sqrt`, `qadic_sqrt`; the
  `padic_radix` functions `padic_radix_sqrt` and `padic_radix_rsqrt` use
  the analogous lifting in the radix representation
  (`radix_integer_sqrtmod_limbs`, `radix_integer_rsqrtmod_limbs`).
- Teichmüller lifts: `_padic_teichmuller`, `padic_teichmuller`,
  `_qadic_teichmuller`, `qadic_teichmuller`.
- The Frobenius lift in `_qadic_frobenius_a` is a Newton iteration of the
  same shape for the defining polynomial
  ([](algorithms/p-adic/frobenius#frobenius)).
