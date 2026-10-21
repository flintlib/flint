---
title: Rational reconstruction
functions: fmpq_mod_fmpz fmpq_reconstruct_fmpz_2 fmpq_reconstruct_fmpz
---

A multimodular or $p$-adic algorithm whose answer is a rational number
computes it as a residue $a$ modulo a large integer $m$, a product of
primes or a prime power. *Rational reconstruction* recovers the fraction
$n/d$ from $a = n d^{-1} \bmod m$, provided that $|n|$ and $d$ are small
compared with $m$. `fmpq_mod_fmpz` computes the residue of a fraction,
and the reconstruction functions invert it:

- `fmpq_reconstruct_fmpz_2` takes explicit bounds $N$ and $D$ with $2ND <
  m$ and finds the unique $n/d$ with $|n| \le N$ and $0 < d \le D$, if it
  exists;
- `fmpq_reconstruct_fmpz` uses the balanced bounds $N = D =
  \lfloor\sqrt{\lfloor (m-1)/2 \rfloor}\rfloor$.

The residue must be reduced, $0 \le a < m$, and $m > 2$. The result is
canonical: $\gcd(n, d) = 1$ and $d > 0$.

::: summary algorithms/number-theory/rational-reconstruction :::

## Success and failure {#failure}

The bound $2ND < m$ makes the answer unique: two fractions $n_1/d_1 \ne
n_2/d_2$ within the bounds have different residues, because $m$ would
otherwise divide $n_1 d_2 - n_2 d_1$, a nonzero number of absolute value
at most $2ND$. A return value of 0 therefore proves that no fraction
within the bounds has residue $a$.

A return value of 1 proves only that the fraction found is the unique one
*within the bounds*. If the true answer of a computation is larger than
the bounds, the function either fails or returns a different, small
fraction with the same residue. About $6/\pi^2 \approx 61\%$ of all
residues have a reconstruction with the balanced bounds, so a successful
reconstruction of a residue that is not known to come from a small
fraction is weak evidence. Multimodular algorithms therefore either use a
priori bounds on the answer (Hadamard's bound for determinants and
solutions of linear systems, for example) and choose $m$ larger than $2ND$,
or verify the reconstructed result, or add primes until the reconstruction
stabilises.

When the shape of the answer is known, unbalanced bounds help. For a
fraction whose denominator divides a known $d_0$, take $D = d_0$ and $N =
\lfloor (m-1)/(2 d_0) \rfloor$: almost all of $m$ is then available for
the numerator.

## Algorithm {#algorithm}

The fraction is found among the remainders and cofactors of the extended
Euclidean algorithm applied to $m$ and $a$, stopping at the first
remainder $r_j \le N$; then $n/d = \pm r_j / t_j$ where $t_j$ is the
cofactor of $a$, and the reconstruction succeeds if $|t_j| \le D$ and
$\gcd(r_j, t_j) = 1$ (Wang's algorithm [@WGD1982]). The implementation
works in one or two words or on fixed arrays of up to 12 limbs for small
$m$, with Lehmer's method for larger $m$, and with the half-GCD algorithm
when the modulus has at least `FMPQ_RECONSTRUCT_HGCD_CUTOFF` (500) more
limbs than $N$, for a total cost of $O(M(\log m) \log\log m)$. Since the
search stops at a remainder of the size of $N$, its cost depends on the
size of $m$ relative to $N$, not only on the size of $m$.
`_fmpq_reconstruct_fmpz_2_naive` is the plain Euclidean algorithm on
`fmpz`, used for testing.

## Example {#example}

The modulus is the prime $m = 10^9 + 7$, for which the balanced bound is
$N = D = 22360$. The fraction $-355/113$ is recovered. The fraction
$123456/789 = 41152/263$ has a numerator above the balanced bound; the
balanced reconstruction fails, and unbalanced bounds $N = 400000$, $D =
1000$ (with $2ND < m$) recover it. The residue of $1234567/89$ is
reconstructed with the balanced bounds, but to a different fraction of
small height, which shows that a successful reconstruction has to be
checked when the answer is not known to be within the bounds.

```c example name=fmpq_reconstruct_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>

static void
show(const fmpq_t x, const fmpz_t m, const fmpz_t N, const fmpz_t D)
{
    fmpz_t a;
    fmpq_t y;
    int ok;

    fmpz_init(a);
    fmpq_init(y);
    fmpq_mod_fmpz(a, x, m);
    if (N == NULL)
        ok = fmpq_reconstruct_fmpz(y, a, m);
    else
        ok = fmpq_reconstruct_fmpz_2(y, a, m, N, D);

    fmpq_print(x);
    flint_printf(" -> ");
    fmpz_print(a);
    if (ok)
    {
        flint_printf(" -> ");
        fmpq_print(y);
        flint_printf("\n");
    }
    else
        flint_printf(" -> no fraction within the bounds\n");

    fmpz_clear(a);
    fmpq_clear(y);
}

int main(void)
{
    fmpz_t m, N, D;
    fmpq_t x;

    fmpz_init_set_ui(m, 1000000007);
    fmpz_init_set_ui(N, 400000);
    fmpz_init_set_ui(D, 1000);
    fmpq_init(x);

    fmpq_set_si(x, -355, 113);
    show(x, m, NULL, NULL);

    fmpq_set_si(x, 123456, 789);
    show(x, m, NULL, NULL);
    show(x, m, N, D);

    fmpq_set_si(x, 1234567, 89);
    show(x, m, NULL, NULL);

    fmpz_clear(m);
    fmpz_clear(N);
    fmpz_clear(D);
    fmpq_clear(x);
    return 0;
}
```
```output
-355/113 -> 699115046 -> -355/113
41152/263 -> 726235903 -> no fraction within the bounds
41152/263 -> 726235903 -> 41152/263
1234567/89 -> 955070058 -> -14062/18629
```
