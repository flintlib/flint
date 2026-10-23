---
title: Finding roots that are expected to be distinct
functions: nmod_poly_find_distinct_nonzero_roots _nmod_poly_split_rabin
---

`nmod_poly_find_distinct_nonzero_roots` answers a specific question: does
$P$ of degree $d$ split over $\mathbb{F}_p$ into $d$ distinct linear
factors $x - r_i$ with $r_i \ne 0$, and if so, what are the $r_i$? This is
the situation in the decoding of Reed–Solomon and BCH codes (see
[](tutorials/algebra/reed-solomon)), in sparse interpolation, where the roots are
the values of the monomials, and in other algorithms that produce a
polynomial whose roots are known to be distinct powers of an element. For
general polynomials, `nmod_poly_roots` and `nmod_poly_factor` in
`nmod_poly_factor.h` find all roots with multiplicities
([](algorithms/polynomials/root-finding)).

## Algorithm {#algorithm}

The algorithm is the equal-degree splitting of Berlekamp and Rabin
[@Ber1970] [@Rab1980b] for factors of degree 1, combined with the test
that all factors are linear and distinct. For odd $p$, every $a \in
\mathbb{F}_p^\times$ satisfies $a^{(p-1)/2} = \pm 1$, so
$$
x^{p-1} - 1 = \bigl(x^{(p-1)/2} - 1\bigr)\bigl(x^{(p-1)/2} + 1\bigr) =
\prod_{a \ne 0} (x - a)
$$
separates the nonzero squares from the non-squares. The function computes
$t = x^{(p-1)/2} \bmod P$ with `nmod_poly_powmod_ui_binexp_preinv` and the
two GCDs $\gcd(t - 1, P)$ and $\gcd(t + 1, P)$. Their degrees add up to
$d$ exactly when $P$ divides $x^{p-1} - 1$, that is, when $P$ has $d$
distinct nonzero roots in $\mathbb{F}_p$; otherwise the function returns 0
at this point, at the cost of one exponentiation and two GCDs.

Each factor of degree at least 2 is then split by
`_nmod_poly_split_rabin`: the substitution $x \mapsto x + \delta$ for a
random $\delta$ permutes the roots, and $\gcd((x + \delta)^{(p-1)/2} - 1,
f)$ collects those $r$ for which $r + \delta$ is a nonzero square, which is
a proper factor with probability about $1 - 2^{1 - \deg f}$. The factors
are kept on a stack until all are linear. The expected cost is $O(M(d)
\log d \log p)$ operations. The random state is created inside the
function with a fixed seed, so the results are reproducible.

The function needs an odd prime modulus; for $p = 2$ it returns 0 for
every $P$ of degree at least 2, since $\mathbb{F}_2$ has only one nonzero
element. The roots are returned in no particular order.

## Example {#example}

A polynomial with 1000 prescribed distinct roots, the same polynomial
times $x^2 + 1$ (which has no roots modulo $p \equiv 3 \pmod 4$), and one
with a repeated root:

```c example name=nmod_poly_distinct_roots
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

static int cmp(const void * a, const void * b)
{
    ulong x = *(const ulong *) a, y = *(const ulong *) b;
    return (x > y) - (x < y);
}

int main(void)
{
    nmod_poly_t P, Q;
    ulong p = 1000003;                       /* p = 3 mod 4 */
    ulong xs[1000], roots[1002];
    slong i;

    nmod_poly_init(P, p);
    nmod_poly_init(Q, p);

    for (i = 0; i < 1000; i++)
        xs[i] = (UWORD(7919) * (i + 1)) % p;
    nmod_poly_product_roots_nmod_vec(P, xs, 1000);

    flint_printf("split: %d\n", nmod_poly_find_distinct_nonzero_roots(roots, P));
    qsort(roots, 1000, sizeof(ulong), cmp);
    qsort(xs, 1000, sizeof(ulong), cmp);
    flint_printf("roots found: %d\n", _nmod_vec_equal(roots, xs, 1000));

    nmod_poly_set_str(Q, "3 1000003  1 0 1");      /* x^2 + 1 */
    nmod_poly_mul(Q, Q, P);
    flint_printf("times x^2 + 1: %d\n", nmod_poly_find_distinct_nonzero_roots(roots, Q));

    nmod_poly_set_coeff_ui(Q, 0, p - xs[0]);       /* repeated root */
    nmod_poly_set_coeff_ui(Q, 1, 1);
    nmod_poly_truncate(Q, 2);
    nmod_poly_mul(Q, Q, P);
    flint_printf("repeated root: %d\n", nmod_poly_find_distinct_nonzero_roots(roots, Q));

    nmod_poly_clear(P);
    nmod_poly_clear(Q);
    return 0;
}
```
```output
split: 1
roots found: 1
times x^2 + 1: 0
repeated root: 0
```
