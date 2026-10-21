---
title: Getting started with finite fields
summary: How to create a finite field with fq_default, compute with its elements, choose or inspect the representation, give your own defining polynomial, convert to and from polynomials over the prime field, and use polynomials, matrices and the generic ring interface over the same field.
order: 1
---

FLINT represents a finite field $\mathbb{F}_q$, $q = p^d$, by a context
object, and its elements by values that are only meaningful together with
the context. The module `fq_default` hides the choice between five
representations behind one interface: it is the module to start with,
and the only one most programs need.

## A first program {#first}

The program creates $\mathbb{F}_{5^3}$, takes its generator $a$ (a root of
the defining polynomial), and computes with $a + 3$.

```c example name=fq_default_first
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq_default.h>

int main(void)
{
    fmpz_t p, t;
    fq_default_ctx_t K;
    fq_default_t a, b, c;

    fmpz_init_set_ui(p, 5);
    fmpz_init(t);
    fq_default_ctx_init(K, p, 3, "a");         /* GF(5^3) */
    fq_default_ctx_print(K);

    fq_default_init(a, K);
    fq_default_init(b, K);
    fq_default_init(c, K);

    fq_default_gen(a, K);                      /* a */
    fq_default_set_ui(b, 3, K);
    fq_default_add(b, b, a, K);                /* b = a + 3 */
    fq_default_mul(c, a, b, K);
    flint_printf("a*(a+3) = "); fq_default_print_pretty(c, K); flint_printf("\n");
    fq_default_inv(c, b, K);
    fq_default_mul(c, c, b, K);
    flint_printf("(a+3)^-1 * (a+3) is one: %d\n", fq_default_is_one(c, K));
    fq_default_pow_ui(c, a, 124, K);
    flint_printf("a^124 = "); fq_default_print_pretty(c, K); flint_printf("\n");
    fq_default_norm(t, b, K);
    flint_printf("N(a+3) = %{fmpz}\n", t);
    fq_default_trace(t, b, K);
    flint_printf("Tr(a+3) = %{fmpz}\n", t);

    fq_default_clear(a, K);
    fq_default_clear(b, K);
    fq_default_clear(c, K);
    fq_default_ctx_clear(K);
    fmpz_clear(p);
    fmpz_clear(t);
    return 0;
}
```
```output
Zech Representation:
p = 5
d = 3
f(X) = 3 + 3*X + X^3
a*(a+3) = a^120
(a+3)^-1 * (a+3) is one: 1
a^124 = a^0
N(a+3) = 3
Tr(a+3) = 4
```

The field is defined by the Conway polynomial $f = X^3 + 3X + 3$, and
since it has only 125 elements, `fq_default` chose the Zech logarithm
representation, in which every nonzero element is stored as its discrete
logarithm with respect to $a$. This is why the pretty printer writes
$a(a + 3)$ as $a^{120}$, and $a^{124} = a^{q-1}$ as $a^0 = 1$. In the
other representations elements print as polynomials in $a$. The norm and
the trace are elements of $\mathbb{F}_5$, returned as integers in
$[0, 5)$: indeed $N(a + 3) = -f(-3) = 33 \equiv 3$ and
$\operatorname{Tr}(a + 3) = \operatorname{Tr}(a) + 3 \cdot 3 = 0 + 9
\equiv 4$.

Every element must be initialised with `fq_default_init` for its
context and cleared with `fq_default_clear`; the context is cleared last.
Outputs may coincide with inputs, as in `fq_default_add(b, b, a, K)`.

## Which representation? {#representation}

`fq_default_ctx_init` picks the representation from $p$ and $d$:

| Field | Representation | Type code |
|---|---|---|
| $d > 1$, $q < 2^{16}$ (more precisely $\operatorname{bits}(p) \cdot d \le 16$) | Zech logarithms, `fq_zech` | `FQ_DEFAULT_FQ_ZECH` |
| $d > 1$, $p < 2^{\bits}$ | polynomials over `nmod`, `fq_nmod` | `FQ_DEFAULT_FQ_NMOD` |
| $d > 1$, larger $p$ | polynomials over `fmpz_mod`, `fq` | `FQ_DEFAULT_FQ` |
| $d = 1$, $p < 2^{\bits}$ | residues, `nmod` | `FQ_DEFAULT_NMOD` |
| $d = 1$, larger $p$ | residues, `fmpz_mod` | `FQ_DEFAULT_FMPZ_MOD` |

Zech logarithms make multiplication and division a single addition or
subtraction modulo $q - 1$, and addition one table lookup, but they need
tables with $2q$ words. The polynomial representations reduce a product of
two polynomials of degree less than $d$ modulo the defining polynomial;
a sparse defining polynomial, such as most Conway polynomials and the
trinomials chosen when no Conway polynomial is known, makes this
reduction cheap. [](algorithms/finite-fields/overview) compares the
representations in detail.

The program below runs the same code in all five representations and
checks two identities, $x^q = x$ and $\sigma^d = \mathrm{id}$ for the
Frobenius automorphism $\sigma(x) = x^p$.

```c example name=fq_default_all_types
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq_default.h>

void check_field(const fmpz_t p, slong d)
{
    fq_default_ctx_t K;
    fq_default_t a, b;
    fmpz_t q;
    slong i;
    int ok;

    fq_default_ctx_init(K, p, d, "a");
    fq_default_init(a, K);
    fq_default_init(b, K);
    fmpz_init(q);
    fq_default_ctx_order(q, K);

    fq_default_gen(a, K);
    fq_default_set_ui(b, 1, K);
    fq_default_add(a, a, b, K);              /* a = gen + 1 */

    fq_default_pow(b, a, q, K);              /* a^q = a */
    ok = fq_default_equal(a, b, K);

    fq_default_set(b, a, K);
    for (i = 0; i < d; i++)                  /* sigma^d = identity */
        fq_default_frobenius(b, b, 1, K);
    ok = ok && fq_default_equal(a, b, K);

    flint_printf("p = %{fmpz}, d = %wd, type %d: %s\n",
                 p, d, fq_default_ctx_type(K), ok ? "ok" : "FAIL");

    fq_default_clear(a, K);
    fq_default_clear(b, K);
    fq_default_ctx_clear(K);
    fmpz_clear(q);
}

int main(void)
{
    fmpz_t p;
    fmpz_init(p);

    fmpz_set_ui(p, 3);
    check_field(p, 5);
    check_field(p, 20);
    fmpz_set_ui(p, 1000003);
    check_field(p, 1);
    fmpz_set_str(p, "170141183460469231731687303715884105727", 10);
    check_field(p, 1);                       /* 2^127 - 1 */
    check_field(p, 3);

    fmpz_clear(p);
    return 0;
}
```
```output
p = 3, d = 5, type 1: ok
p = 3, d = 20, type 2: ok
p = 1000003, d = 1, type 4: ok
p = 170141183460469231731687303715884105727, d = 1, type 5: ok
p = 170141183460469231731687303715884105727, d = 3, type 3: ok
```

`fq_default_ctx_init_type` forces a representation, for example to
compare timings. The type must suit the field: a Zech context for
$\mathbb{F}_{3^{20}}$ would need tables of $2 \cdot 3^{20}$ words (55
gigabytes), and the prime field types ignore $d$.

## Your own defining polynomial {#modulus}

`fq_default_ctx_init_modulus` defines the field as
$\mathbb{F}_p[x]/(f)$ for a given irreducible $f$. The program uses the
field of the AES block cipher, $\mathbb{F}_2[x]/(x^8 + x^4 + x^3 + x +
1)$, and inverts the byte `0x53`, that is $x^6 + x^4 + x + 1$. Elements
are created from polynomials over $\mathbb{F}_p$ with
`fq_default_set_fmpz_mod_poly`, which reduces modulo $f$.

```c example name=fq_default_aes
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fq_default.h>

int main(void)
{
    fmpz_t p;
    fmpz_mod_ctx_t F2;
    fmpz_mod_poly_t f;
    fq_default_ctx_t K;
    fq_default_t a, b;

    fmpz_init_set_ui(p, 2);
    fmpz_mod_ctx_init(F2, p);
    fmpz_mod_poly_init(f, F2);

    /* x^8 + x^4 + x^3 + x + 1 */
    fmpz_mod_poly_set_coeff_ui(f, 8, 1, F2);
    fmpz_mod_poly_set_coeff_ui(f, 4, 1, F2);
    fmpz_mod_poly_set_coeff_ui(f, 3, 1, F2);
    fmpz_mod_poly_set_coeff_ui(f, 1, 1, F2);
    fmpz_mod_poly_set_coeff_ui(f, 0, 1, F2);

    fq_default_ctx_init_modulus(K, f, F2, "x");
    flint_printf("type %d\n", fq_default_ctx_type(K));

    fq_default_init(a, K);
    fq_default_init(b, K);

    /* 0x53 = x^6 + x^4 + x + 1 */
    fmpz_mod_poly_zero(f, F2);
    fmpz_mod_poly_set_coeff_ui(f, 6, 1, F2);
    fmpz_mod_poly_set_coeff_ui(f, 4, 1, F2);
    fmpz_mod_poly_set_coeff_ui(f, 1, 1, F2);
    fmpz_mod_poly_set_coeff_ui(f, 0, 1, F2);
    fq_default_set_fmpz_mod_poly(a, f, K);
    fq_default_inv(b, a, K);
    flint_printf("1/("); fq_default_print_pretty(a, K);
    flint_printf(") = "); fq_default_print_pretty(b, K); flint_printf("\n");

    fq_default_clear(a, K);
    fq_default_clear(b, K);
    fq_default_ctx_clear(K);
    fmpz_mod_poly_clear(f, F2);
    fmpz_mod_ctx_clear(F2);
    fmpz_clear(p);
    return 0;
}
```
```output
type 2
1/(x^6+x^4+x+1) = x^7+x^6+x^3+x
```

The inverse is `0xCA`, as in the AES specification. The field has only
256 elements, but the context uses `fq_nmod` (type 2), not Zech
logarithms: the AES polynomial is irreducible but not primitive (the
class of $x$ has order 51, not 255), and a Zech context needs the
generator to be a primitive element, so the initialisation falls back to
the next representation. Irreducibility is not checked; a reducible $f$
gives a ring in which some inversions fail or return wrong results.

## Conversions {#conversions}

An element of $\mathbb{F}_q = \mathbb{F}_p[x]/(f)$ has a unique
representative of degree less than $d$ in $\mathbb{F}_p[x]$:

- `fq_default_get_fmpz_mod_poly`, `fq_default_get_nmod_poly` and
  `fq_default_get_fmpz_poly` return it, and the corresponding `set`
  functions convert back, reducing modulo $f$;
- `fq_default_get_coeff_fmpz` extracts a single coefficient;
- `fq_default_set_ui`, `fq_default_set_si` and `fq_default_set_fmpz` map
  integers into the prime field, and `fq_default_get_fmpz` maps back,
  returning 0 if the element is not in $\mathbb{F}_p$;
- `fq_default_get_str_pretty` and `fq_default_print_pretty` give a
  readable form, and `fq_default_get_str` and `fq_default_print` the
  plain format of the representation.

In the Zech representation the conversions to and from polynomials go
through tables, which is fast; the readable form is a power of the
generator, as in the first example.

## Polynomials, matrices and factorisation {#poly-mat}

The modules `fq_default_poly`, `fq_default_mat` and
`fq_default_poly_factor` provide polynomials and matrices over an
`fq_default` field, with the same context. Factoring $2x^4 + 2$ over
$\mathbb{F}_9$:

```c example name=fq_default_factor
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq_default_poly_factor.h>

int main(void)
{
    fmpz_t p;
    fq_default_ctx_t K;
    fq_default_poly_t f;
    fq_default_poly_factor_t fac;
    fq_default_t c;

    fmpz_init_set_ui(p, 3);
    fq_default_ctx_init(K, p, 2, "a");         /* GF(9) */
    fq_default_poly_init(f, K);
    fq_default_poly_factor_init(fac, K);
    fq_default_init(c, K);

    fq_default_set_ui(c, 2, K);
    fq_default_poly_set_coeff(f, 4, c, K);
    fq_default_poly_set_coeff(f, 0, c, K);     /* f = 2x^4 + 2 */
    fq_default_poly_print_pretty(f, "x", K);
    flint_printf("\n");

    fq_default_poly_factor(fac, c, f, K);
    flint_printf("leading coefficient "); fq_default_print_pretty(c, K);
    flint_printf("\n");
    fq_default_poly_factor_print_pretty(fac, "x", K);

    fq_default_clear(c, K);
    fq_default_poly_factor_clear(fac, K);
    fq_default_poly_clear(f, K);
    fq_default_ctx_clear(K);
    fmpz_clear(p);
    return 0;
}
```
```output
(a^4)*x^4+(a^4)
leading coefficient a^4
x+(a^7) ^ 1
x+(a^3) ^ 1
x+(a^5) ^ 1
x+(a^1) ^ 1
```

In $\mathbb{F}_9$ the element $2 = -1$ is $a^4$, and $x^4 + 1$ splits into
the four linear factors $x - \zeta$ for the primitive 8th roots of unity
$\zeta = a^1, a^3, a^5, a^7$ (note $-a^7 = a^4 a^7 = a^3$).

## The generic ring interface {#gr}

An `fq_default_ctx_t` is a generic ring context, and an `fq_default_t`
an element of that ring. `FQ_DEFAULT_GR_CTX` gives access to all
functions of [](gr) and the modules built on it (`gr_poly`, `gr_mat`,
`gr_mpoly`, ...), which return a status code instead of throwing. For
speed, `fq_default_ctx_inner` and the members of the element union give
direct access to the underlying module.

```c example name=fq_default_gr
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq_default.h>
#include <flint/fq_nmod.h>
#include <flint/gr.h>

int main(void)
{
    fmpz_t p;
    fq_default_ctx_t K;
    fq_default_t a;
    gr_ctx_struct * R;
    int status = GR_SUCCESS;

    fmpz_init_set_ui(p, 101);
    fq_default_ctx_init(K, p, 4, "z");
    fq_default_init(a, K);

    /* the context as a generic ring */
    R = FQ_DEFAULT_GR_CTX(K);
    gr_ctx_println(R);
    status |= gr_gen(a, R);
    status |= gr_add_ui(a, a, 1, R);
    status |= gr_pow_ui(a, a, 5, R);
    status |= gr_println(a, R);

    /* direct access to the fq_nmod representation */
    if (fq_default_ctx_type(K) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_ctx_struct * L = fq_default_ctx_inner(K);
        fq_nmod_sqr(a->fq_nmod, a->fq_nmod, L);
        fq_nmod_print_pretty(a->fq_nmod, L);
        flint_printf("\n");
    }

    fq_default_clear(a, K);
    fq_default_ctx_clear(K);
    fmpz_clear(p);
    return status;
}
```
```output
GF(101^4) (fq_nmod)
9*z^3+28*z^2+17*z+92
28*z^3+5*z^2+77*z+83
```

Code that works with one known representation can use `fq_nmod`,
`fq_zech` or `fq` directly; their interfaces are almost identical to that
of `fq_default`, without the dispatch, and add low-level functions on
coefficient arrays.
