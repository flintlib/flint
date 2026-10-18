---
title: Computing in number fields
summary: Exact arithmetic with algebraic numbers in one or several quadratic fields: Fibonacci numbers from Binet's formula, the precision limit that keeps large powers unexpanded, and the Swinnerton-Dyer polynomials; the example programs binet and swinnerton_dyer_poly.
order: 2
---

When all generators of a field are algebraic, [](ca) computes with
exact number field elements, and results that are rational come out as
rational numbers without a zero test. This tutorial shows this with two
classical constructions, following the example programs `binet` and
`swinnerton_dyer_poly` in the `examples` directory.

## Binet's formula {#binet}

The Fibonacci numbers are given by

$$
F_n = \frac{\varphi^n - (1 - \varphi)^n}{\sqrt 5}, \qquad
\varphi = \frac{1 + \sqrt 5}{2}.
$$

The numbers $\varphi$, $1 - \varphi$ and $\sqrt 5$ all lie in the
quadratic field $\mathbb{Q}(\sqrt 5)$, which [](ca) represents as a
number field with [](nf_elem) arithmetic. The powers are computed by
repeated squaring in the number field, and the final division leaves an
element with zero irrational part, which is printed as an integer.

Large powers are not always expanded. `ca_pow_ui` computes $x^n$ in a
number field only if the bit size of $x$ times $n$ is below the context
option `CA_OPT_PREC_LIMIT` (4096 by default); otherwise $x^n$ becomes
a new generator, a symbolic power. This protects against accidental
computations with gigantic numbers, but it also means that $F_{10^6}$,
an integer with 208988 digits, is not recognised as an integer by
default. In a context with a larger limit, the powers are expanded and
`ca_get_fmpz` converts the result to an `fmpz_t`.

## Swinnerton-Dyer polynomials {#swinnerton-dyer}

The Swinnerton-Dyer polynomial

$$
S_n = \prod (x \pm \sqrt{2} \pm \sqrt{3} \pm \sqrt{5} \pm \cdots \pm \sqrt{p_n}),
$$

with all $2^n$ combinations of signs, has integer coefficients, degree
$2^n$, and is irreducible over $\mathbb{Q}$ although it factors into
polynomials of degree at most 2 modulo every prime; it is a standard
worst case for polynomial factorisation algorithms. The coefficients
can be computed by multiplying out the linear factors over the field
$\mathbb{Q}(\sqrt 2, \sqrt 3, \ldots, \sqrt{p_n})$, which [](ca)
represents with $n$ generators and the relations $a_k^2 = p_k$. The
program below does this for $n = 3$ with `ca_poly_set_roots`, which
multiplies the factors in a balanced product tree.

```c example name=tutorial_exact_number_fields
#include <flint/ca.h>
#include <flint/ca_vec.h>
#include <flint/ca_poly.h>

/* F_n = (phi^n - (1 - phi)^n) / sqrt(5) */
static void
binet(ca_t res, ulong n, ca_ctx_t ctx)
{
    ca_t sqrt5, phi, psi;

    ca_init(sqrt5, ctx);
    ca_init(phi, ctx);
    ca_init(psi, ctx);

    ca_sqrt_ui(sqrt5, 5, ctx);
    ca_add_ui(phi, sqrt5, 1, ctx);
    ca_div_ui(phi, phi, 2, ctx);
    ca_ui_sub(psi, 1, phi, ctx);
    ca_pow_ui(phi, phi, n, ctx);
    ca_pow_ui(psi, psi, n, ctx);
    ca_sub(res, phi, psi, ctx);
    ca_div(res, res, sqrt5, ctx);

    ca_clear(sqrt5, ctx);
    ca_clear(phi, ctx);
    ca_clear(psi, ctx);
}

int main(void)
{
    ca_ctx_t ctx, ctx2;
    ca_t f;
    fmpz_t z;
    ca_vec_t roots;
    ca_poly_t S;
    ulong primes[3] = {2, 3, 5};
    ulong mult[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    slong i, j;

    ca_ctx_init(ctx);
    ca_init(f, ctx);
    fmpz_init(z);

    binet(f, 250, ctx);
    flint_printf("F_250 = "); ca_print(f, ctx); flint_printf("\n");

    binet(f, 1000000, ctx);
    flint_printf("F_1000000 = "); ca_print(f, ctx); flint_printf("\n");
    ca_clear(f, ctx);

    /* a context with a larger precision limit expands the powers */
    ca_ctx_init(ctx2);
    ca_ctx_set_option(ctx2, CA_OPT_PREC_LIMIT, 10000000);
    ca_init(f, ctx2);
    binet(f, 1000000, ctx2);
    if (ca_get_fmpz(z, f, ctx2))
        flint_printf("F_1000000 is an integer with %wd digits\n",
            (slong) fmpz_sizeinbase(z, 10));
    ca_clear(f, ctx2);
    ca_ctx_clear(ctx2);

    /* Swinnerton-Dyer polynomial: roots +/- sqrt(2) +/- sqrt(3) +/- sqrt(5) */
    ca_vec_init(roots, 8, ctx);
    ca_init(f, ctx);
    for (i = 0; i < 8; i++)
    {
        ca_zero(ca_vec_entry_ptr(roots, i), ctx);
        for (j = 0; j < 3; j++)
        {
            ca_sqrt_ui(f, primes[j], ctx);
            if ((i >> j) & 1)
                ca_neg(f, f, ctx);
            ca_add(ca_vec_entry_ptr(roots, i), ca_vec_entry_ptr(roots, i), f, ctx);
        }
    }
    ca_poly_init(S, ctx);
    ca_poly_set_roots(S, roots, mult, ctx);
    flint_printf("coefficients of S_3:");
    for (i = 0; i < S->length; i++)
    {
        flint_printf(" ");
        ca_print(ca_poly_coeff_ptr(S, i), ctx);
    }
    flint_printf("\n");

    ca_poly_clear(S, ctx);
    ca_vec_clear(roots, ctx);
    ca_clear(f, ctx);
    fmpz_clear(z);
    ca_ctx_clear(ctx);
    return 0;
}
```

```output
F_250 = 7.89633e+51 {7896325826131730509282738943634332893686268675876375}
F_1000000 = 1.95328e+208987 {(a*c-b*c)/5 where a = 4.36767e+208987 [Pow(1.61803 {(c+1)/2}, 1.00000e+6 {1000000})], b = 2.28955e-208988 [Pow(-0.618034 {(-c+1)/2}, 1.00000e+6 {1000000})], c = 2.23607 [c^2-5=0]}
F_1000000 is an integer with 208988 digits
coefficients of S_3: 576 0 -960 0 352 0 -40 0 1
```

The example program `binet n` prints $F_n$, and `binet -limit B n`
sets the precision limit to $B$ first. The program
`swinnerton_dyer_poly n` computes all coefficients of $S_n$ with a
more economical product: since $S_n(x)$ is a polynomial in $x^2$, it
multiplies the $2^{n-1}$ quadratic factors $x^2 - (\pm\sqrt 2 \pm \cdots
\pm \sqrt{p_n})^2$ as polynomials in $x^2$ and inflates the result.
For $n = 10$, with integer coefficients of up to about 800 digits,
this takes a few seconds. Much faster ways to compute $S_n$ exist
(through numerical evaluation of the roots in [](arb) and rounding the
coefficients, see `fmpz_poly_swinnerton_dyer`); the point of the
example is that arithmetic in fields with many algebraic generators
works without any special setup.
