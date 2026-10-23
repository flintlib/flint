---
title: Getting started with nf_elem
summary: A first program with number field elements in Q(2^(1/3)): creating the field, the generator, sums, powers, inverses, norms and traces, and the rules for reduction, canonical form and delayed reduction.
order: 1
---

A number field is created from its defining polynomial with `nf_init`,
and its elements are `nf_elem_t` variables that must be initialised and
cleared with the field. Every function takes the field as its last
argument.

## A first program {#first}

The program works in $K = \mathbf{Q}(a)$ with $a^3 = 2$.

```c example name=nf_elem_first
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>

int main(void)
{
    fmpq_poly_t pol;
    nf_t nf;
    nf_elem_t a, b, c;
    fmpq_t t;

    /* K = Q(a) with a^3 = 2 */
    fmpq_poly_init(pol);
    fmpq_poly_set_str(pol, "4  -2 0 0 1");
    nf_init(nf, pol);

    nf_elem_init(a, nf);
    nf_elem_init(b, nf);
    nf_elem_init(c, nf);
    fmpq_init(t);

    nf_elem_gen(a, nf);
    nf_elem_add_si(b, a, 1, nf);                /* b = 1 + a */
    nf_elem_pow(c, b, 5, nf);                   /* reduced modulo a^3 - 2 */
    flint_printf("(1 + a)^5 = "); nf_elem_print_pretty(c, nf, "a");
    flint_printf("\n");

    nf_elem_inv(c, b, nf);
    flint_printf("1/(1 + a) = "); nf_elem_print_pretty(c, nf, "a");
    nf_elem_mul(c, c, b, nf);
    flint_printf("\ncheck: %d\n", nf_elem_is_one(c, nf));

    nf_elem_norm(t, b, nf);
    flint_printf("N(1 + a) = "); fmpq_print(t);
    nf_elem_trace(t, b, nf);
    flint_printf(", Tr(1 + a) = "); fmpq_print(t);
    flint_printf("\n");

    fmpq_set_si(t, 1, 2);
    nf_elem_scalar_mul_fmpq(c, b, t, nf);       /* (1 + a)/2 */
    nf_elem_norm(t, c, nf);
    flint_printf("N((1 + a)/2) = "); fmpq_print(t); flint_printf("\n");

    fmpq_clear(t);
    nf_elem_clear(a, nf); nf_elem_clear(b, nf); nf_elem_clear(c, nf);
    nf_clear(nf);
    fmpq_poly_clear(pol);
    return 0;
}
```
```output
(1 + a)^5 = 12*a^2+15*a+21
1/(1 + a) = 1/3*a^2-1/3*a+1/3
check: 1
N(1 + a) = 3, Tr(1 + a) = 3
N((1 + a)/2) = 3/8
```

The power is reduced with $a^3 = 2$: $(1 + a)^5 = 21 + 15a + 12a^2$.
The inverse follows from $(1 + a)(1 - a + a^2) = 1 + a^3 = 3$. The norm
of $1 + a$ is $-f(-1) = 3$ for $f = x^3 - 2$, the product of the three
conjugates $1 + a\zeta^k$, and the norm is multiplicative, so dividing
by 2 divides it by $2^3$.

## Representation and canonical form {#representation}

An element is stored as a numerator polynomial in $\Z[x]$ of degree
less than $d = \deg f$ and a positive denominator, with no common
factor; the zero element has denominator 1. Fields of degree 1 and 2
store the numerator in a fixed-size array instead of an `fmpq_poly_t`,
and $\mathbf{Q}(i)$, given by exactly $x^2 + 1$, uses a dedicated
multiplication. This is decided by `nf_init` and is invisible to the
user, except that the accessors `NF_ELEM`, `LNF_ELEM` and `QNF_ELEM`
must match the field.

All functions without a leading underscore return canonical, reduced
elements. The functions with an underscore skip the canonicalisation,
and `nf_elem_mul_red` with `red` $= 0$ also skips the reduction modulo
$f$, so that a sum of products can be accumulated in unreduced form and
reduced once with `nf_elem_reduce`; such unreduced elements must not be
passed to other functions, which assume reduced inputs.

## Fields given by reducible polynomials {#reducible}

`nf_init` does not check that $f$ is irreducible. With a reducible $f$
the functions compute in the ring $\mathbf{Q}[x]/(f)$, which is fine for
addition and multiplication, but inverses of zero divisors are
meaningless and norms of zero divisors are zero.
`_nf_elem_invertible_check` tests whether an element is invertible.

## The generic ring interface {#generic}

`gr_ctx_init_nf` creates a generic ring for a number field, with
elements of type `nf_elem_t`, so that `gr_poly`, `gr_mat` and the other
generic modules work over number fields.
