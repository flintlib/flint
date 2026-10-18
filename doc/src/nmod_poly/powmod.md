---
title: Exponentiation modulo a polynomial
functions: nmod_poly_powmod_ui_binexp nmod_poly_powmod_fmpz_binexp nmod_poly_powmod_ui_binexp_preinv nmod_poly_powmod_fmpz_binexp_preinv nmod_poly_powmod_x_ui_preinv nmod_poly_powmod_x_fmpz_preinv
---

Powers $a^e \bmod f$ with huge exponents are the basic operation of
computations in finite fields $\mathbb{F}_p[x]/(f)$: root finding and
equal-degree factorisation use $(x + \delta)^{(p-1)/2} \bmod f$,
distinct-degree factorisation and irreducibility tests use $x^{p^i} \bmod
f$, and primality tests for polynomials use $a^{p^n - 1}$. These functions
compute such powers by binary exponentiation, with $\lfloor \log_2 e
\rfloor$ squarings and at most as many multiplications, each followed by a
reduction modulo $f$, for $O(M(\deg f) \log e)$ operations in total.

## Choosing a function {#choice}

- The functions without `preinv` reduce with `_nmod_poly_divrem`, which
  for large $\deg f$ computes an inverse of the reversed $f$ in every
  reduction. For anything but small $f$, use the `preinv` versions with
  `finv` computed once, as explained in [](nmod_poly/modular-arithmetic).
- For the base $x$, the `powmod_x` functions are faster: they use a
  sliding window of about $\log_2 \deg f$ bits, and multiplying by $x^w$
  for a window $w$ is a shift followed by one reduction instead of a
  product. Requires $\deg f \ge 2$ for the underscore versions.
- The `fmpz` versions accept exponents of any size, such as $p^n$.
- The non-underscore functions reduce `poly` modulo $f$ first, accept $e =
  0$ (giving 1) and allow aliasing. The underscore functions require
  $\deg f \ge 1$, `poly` reduced and zero-padded to length $\deg f$, $e
  \ge 1$, and no aliasing.

For many powers of the same element, `nmod_poly_powers_mod_bsgs` computes
$a^0, \ldots, a^{k-1}$ in parallel, and for the iterates $x^{p^i}$ of the
Frobenius map, composition with $x^p$ is faster than repeated
exponentiation for large $p$; see [](nmod_poly/compose_mod).

## Example: counting roots {#example}

The distinct roots of $f$ in $\mathbb{F}_p$ are the roots of $\gcd(x^p -
x, f)$, since $x^p - x = \prod_{c \in \mathbb{F}_p} (x - c)$. One
exponentiation of $x$ modulo $f$ therefore counts them, and the factor
$\gcd(x^p - x, f)$ splits into linear factors, whose roots
`nmod_poly_find_distinct_nonzero_roots` finds:

```c example name=nmod_poly_count_roots
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t f, finv, t, g, x;
    ulong p = 1000003, roots[10];
    slong i, k;

    nmod_poly_init(f, p);
    nmod_poly_init(finv, p);
    nmod_poly_init(t, p);
    nmod_poly_init(g, p);
    nmod_poly_init(x, p);

    /* f = (x - 10)(x - 20)(x - 30)^2 (x^50 + x + 7) */
    nmod_poly_set_str(f, "51 1000003  7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1");
    for (i = 10; i <= 30; i += 10)
    {
        nmod_poly_set_coeff_ui(t, 1, 1);
        nmod_poly_set_coeff_ui(t, 0, p - i);
        nmod_poly_mul(f, f, t);
    }
    nmod_poly_mul(f, f, t);                       /* second factor x - 30 */

    nmod_poly_reverse(finv, f, f->length);
    nmod_poly_inv_series(finv, finv, f->length);

    nmod_poly_powmod_x_ui_preinv(t, p, f, finv);  /* x^p mod f */
    nmod_poly_set_coeff_ui(x, 1, 1);
    nmod_poly_sub(t, t, x);
    nmod_poly_gcd(g, t, f);
    flint_printf("deg f = %wd, distinct roots: %wd\n",
                 nmod_poly_degree(f), nmod_poly_degree(g));

    if (nmod_poly_find_distinct_nonzero_roots(roots, g))
    {
        k = nmod_poly_degree(g);
        flint_printf("roots:");
        for (i = 0; i < k; i++)
            flint_printf(" %wu", roots[i]);
        flint_printf("\n");
    }

    /* the same power by the generic function */
    nmod_poly_powmod_ui_binexp(g, x, p, f);
    nmod_poly_add(t, t, x);
    flint_printf("powmod agrees: %d\n", nmod_poly_equal(g, t));

    nmod_poly_clear(f);
    nmod_poly_clear(finv);
    nmod_poly_clear(t);
    nmod_poly_clear(g);
    nmod_poly_clear(x);
    return 0;
}
```
```output
deg f = 54, distinct roots: 3
roots: 10 20 30
powmod agrees: 1
```
