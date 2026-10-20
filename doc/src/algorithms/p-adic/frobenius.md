---
title: Frobenius automorphism and Teichmüller lifts in unramified extensions
summary: How qadic computes the Frobenius automorphism of Q_q by Newton lifting the image of the generator and evaluating by rectangular splitting, and how Teichmüller lifts are computed, with the relations between the two.
complexity: $\Sigma^e(X)$ modulo $p^N$: $O(e \log p)$ multiplications modulo $p$, then $O(\sqrt d)$ multiplications in $\mathbf{Z}_q$ and $O(d)$ scalar operations at each of $O(\log N)$ precisions; applying $\Sigma^e$ to an element: $O(\sqrt d)$ multiplications in $\mathbf{Z}_q$. Teichmüller lift: $O(d \log p)$ multiplications in $\mathbf{Z}_q$ at each precision.
order: 50
---

The extension $\mathbf{Q}_q/\mathbf{Q}_p$ of degree $d$, $q = p^d$, is
Galois with cyclic Galois group of order $d$, and reduction modulo $p$
identifies it with the Galois group of $\mathbf{F}_q/\mathbf{F}_p$
[@Rob2000, Ch. 3]. The generator $\Sigma$ corresponding to the Frobenius
$\sigma(x) = x^p$ of $\mathbf{F}_q$ is the *Frobenius automorphism* of
$\mathbf{Q}_q$; it is the basic tool of $p$-adic point counting
algorithms such as Kedlaya's [@Ked2001] and of the canonical lift
methods of Satoh, Skjernaa and Taguchi [@SST2003]. The representation
of $\mathbf{Q}_q = \mathbf{Q}_p[X]/(f)$ is described in
[](algorithms/p-adic/unramified-extensions).

## The Frobenius automorphism {#frobenius}

Since $\Sigma$ fixes $\mathbf{Q}_p$, it is determined by the image of
the generator: for $g = \sum_i c_i X^i$ with $c_i \in \mathbf{Q}_p$,
$\Sigma^e(g) = \sum_i c_i \, \Sigma^e(X)^i = g(\Sigma^e(X))$.

::: proposition {#prop-frob-x}
$\Sigma^e(X)$ is the unique root $\xi \in \mathbf{Z}_q$ of $f$ with
$\xi \equiv X^{p^e} \pmod p$.
:::

::: proof
$\Sigma^e(X)$ is a root of $f$ since $\Sigma^e$ is an automorphism, and
it reduces to $\sigma^e(\bar X) = \bar X^{p^e}$. The reduction $\bar f$
is irreducible over the finite field $\mathbf{F}_p$, hence separable,
so $f'(\xi) \not\equiv 0 \pmod p$ and Hensel's lemma
([](algorithms/p-adic/hensel-lifting#lem-newton)) gives uniqueness.
:::

`_qadic_frobenius_a` computes $\xi$ modulo $p^N$ by Newton's method.
It computes $X^{p^e}$ modulo $(f, p)$ by binary powering, inverts
$f'(\xi)$ modulo $(f, p)$, and then, along the precision sequence
$N, \lceil N/2 \rceil, \ldots$ of
[](algorithms/p-adic/hensel-lifting#sequence), applies

$$
\xi \gets \xi - f(\xi) \, s, \qquad
s \gets s\bigl(2 - f'(\xi) s\bigr),
$$

where $s$ approximates $f'(\xi)^{-1}$; by
[](algorithms/p-adic/hensel-lifting#lem-newton) the precision of $\xi$
and of $s$ doubles in each step. The polynomials $f(\xi)$ and
$f'(\xi)$ are compositions: polynomials of degree $d$ and $d - 1$ with
coefficients in $\mathbf{Z}_p$ evaluated at an element of $\mathbf{Z}_q$.

::: algorithm Composition modulo f by rectangular splitting {#alg-compose}
**Input:** $g = \sum_{i < \ell} c_i X^i$ with $c_i \in \Z$, an element
$\xi$ of $\mathbf{Z}_q$, a modulus $p^M$.

**Output:** $g(\xi) \bmod (f, p^M)$.

1. Let $B = \lfloor \sqrt{\ell} \rfloor$ and compute
   $\xi^2, \ldots, \xi^B$ modulo $(f, p^M)$.
2. For $h = \lceil \ell/B \rceil - 1, \ldots, 0$: set
   $r \gets r \xi^B + \sum_{k < B} c_{hB + k} \xi^k$, reducing modulo
   $(f, p^M)$.
:::

This is the Paterson–Stockmeyer scheme [@PS1973]: about $2\sqrt\ell$
multiplications in $\mathbf{Z}_q$ and $\ell$ scalar multiplications of
elements by integers. FLINT uses it for $\ell \ge 6$ and Horner's rule
for shorter $g$.

`_qadic_frobenius` then applies $\Sigma^e$ to an element $g$ by
evaluating $g(\Sigma^e(X))$ with the same algorithm. Two special cases
avoid the lifting: an element of $\mathbf{Q}_p$ (length 1) is fixed, and
for $N = 1$ the image is $g^{p^e}$ computed by binary powering modulo
$(f, p)$. `qadic_frobenius` reduces $e$ modulo $d$, so that negative
$e$ give the inverse automorphisms, and applies the map to the unit
part, the valuation being unchanged.

`qadic_frobenius` recomputes $\Sigma^e(X)$ on every call. A caller
that applies the same power of $\Sigma$ to many elements can call
`_qadic_frobenius_a` once and evaluate the elements at its result.

## Teichmüller lifts {#teichmuller}

The Teichmüller lift $\omega(x)$ of a unit $x \in \mathbf{Z}_q$ is the
unique $(q - 1)$-th root of unity congruent to $x$ modulo $p$. The map
$x \mapsto \omega(x)$ is multiplicative, the roots of unity of order
prime to $p$ in $\mathbf{Q}_q$ are exactly the Teichmüller lifts, and

$$
\Sigma(\omega(x)) = \omega(x)^p ,
$$

because both sides are $(q - 1)$-th roots of unity congruent to $x^p$
modulo $p$. `_qadic_teichmuller` computes $\omega(x)$ as the root of
$t^q - t$ congruent to $x$, by the iteration
$t \gets t - (q - 1)^{-1}(t^q - t)$ of
[](algorithms/p-adic/hensel-lifting#teichmuller), the inverse of
$q - 1$ being lifted alongside. Each step computes $t^q$ by binary
powering with exponent $q = p^d$, which takes about $d \log_2 p$
squarings in $\mathbf{Z}_q$, so the cost is $O(d \log p)$
multiplications at precision $N$, against $O(\sqrt d)$ for one
application of $\Sigma$. For an element of $\mathbf{Z}_p$ the lift is
computed in $\mathbf{Z}_p$ by `_padic_teichmuller`, where the exponent
is $p$.

## Use in FLINT {#use}

- `_qadic_frobenius_a` computes $\Sigma^e(X)$ modulo $p^N$ for
  $0 < e < d$.
- `_qadic_frobenius` and `qadic_frobenius` apply $\Sigma^e$ to an
  element.
- `_qadic_teichmuller` and `qadic_teichmuller` compute Teichmüller
  lifts; `_padic_teichmuller` and `padic_teichmuller` do the same in
  $\mathbf{Z}_p$.
- `qadic_trace` and `qadic_norm` do not use the Frobenius map; they use
  Newton's identities, resultants and logarithms
  ([](algorithms/p-adic/unramified-extensions#trace)).
