---
title: Computing in the class group
summary: How to list the reduced forms of a discriminant, compose and power forms with NUCOMP followed by reduction, and compute the order of a form and the exponent of the class group, with an example.
order: 1
---

The primitive positive definite binary quadratic forms of a discriminant
$D < 0$, up to $\mathrm{SL}_2(\Z)$-equivalence, form the *form class
group* $\mathrm{Cl}(D)$, a finite abelian group whose order is the class
number $h(D)$. For a fundamental discriminant it is isomorphic to the
ideal class group of the imaginary quadratic field $\Q(\sqrt D)$
[@Coh1996, §5.2], [@Buell1989].

## Representing classes {#classes}

Each class contains exactly one reduced form $(a, b, c)$, with $|b| \le a
\le c$ and $b \ge 0$ if $|b| = a$ or $a = c$; `qfb_reduce` finds it, and
the reduced form serves as the canonical representative, so that two
classes are equal if and only if their reduced forms are equal
(`qfb_equal`). The identity is the principal form (`qfb_principal_form`),
and the inverse of $(a, b, c)$ is $(a, -b, c)$ (`qfb_inverse`). Since $a
\le \sqrt{|D|/3}$ for a reduced form, all reduced forms can be listed for
moderate $|D|$ (`qfb_reduced_forms`), which gives $h(D)$ in
$O(|D|^{1/2 + o(1)})$ time.

## Composition and powers {#composition}

The group operation is Gauss composition. `qfb_nucomp` computes it with
Shanks' NUCOMP algorithm, which interleaves the composition with a
partial reduction so that the intermediate numbers stay of the size of
$|D|^{1/2}$ rather than $|D|$; its output is nearly reduced and must be
passed to `qfb_reduce`:

```
fmpz_abs(L, D);
fmpz_root(L, L, 4);            /* L = floor(|D|^(1/4)) */
qfb_nucomp(r, f, g, D, L);
qfb_reduce(r, r, D);
```

`qfb_nudupl` squares a form, and `qfb_pow_ui` and `qfb_pow` compute
reduced powers by binary powering; their exponents must be nonnegative.
The inputs of the composition functions must be primitive.

## Orders and the exponent {#orders}

For large $|D|$ the reduced forms cannot be listed, but the structure of
the group can still be explored. `qfb_exponent_element` finds the order
of a form whose order has only small prime factors, except possibly one
larger prime (a two-stage method like Pollard's $p - 1$), and
`qfb_exponent` and `qfb_exponent_grh` compute the exponent of the class
group from the orders of prime forms (`qfb_prime_form`) for small primes;
the GRH version uses enough primes to generate the group if the
generalised Riemann hypothesis holds. The exponent is a multiple of the
order of every element and divides $h(D)$; for most discriminants the
class group is cyclic or nearly so, and the exponent is close to $h(D)$.
The stage 1 bound `B1` should be at least about 8200 (see the warning at
`qfb_exponent_element`).

## Example {#example}

The program lists the four reduced forms of discriminant $-56$, shows
that $(3, 2, 5)$ generates the cyclic group $\mathrm{Cl}(-56)$ of order
4, compares the class number of $-1000003$ with the exponent found from
prime forms, and computes the exponent of the class group of
discriminant $-(10^{18} + 3)$, whose class number is near $10^8$.

```c example name=qfb_class_group
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/qfb.h>

int main(void)
{
    qfb * forms;
    qfb_t f, g;
    fmpz_t D, e;
    slong i, num;

    fmpz_init(D);
    fmpz_init(e);
    qfb_init(f);
    qfb_init(g);

    num = qfb_reduced_forms(&forms, -56);
    flint_printf("h(-56) = %wd:", num);
    for (i = 0; i < num; i++)
    {
        flint_printf(" ");
        qfb_print(forms + i);
    }
    flint_printf("\n");

    fmpz_set_si(D, -56);
    qfb_set(f, forms + 2);
    flint_printf("powers of ");
    qfb_print(f);
    flint_printf(":");
    for (i = 1; i <= 4; i++)
    {
        qfb_pow_ui(g, f, D, i);
        flint_printf(" ");
        qfb_print(g);
    }
    flint_printf("\n");
    qfb_array_clear(&forms, num);

    num = qfb_reduced_forms(&forms, -1000003);
    qfb_array_clear(&forms, num);
    fmpz_set_si(D, -1000003);
    qfb_exponent_grh(e, D, 10000, 100);
    flint_printf("h(-1000003) = %wd, exponent ", num);
    fmpz_print(e);
    flint_printf("\n");

    fmpz_set_str(D, "-1000000000000000003", 10);
    if (qfb_exponent_grh(e, D, 100000, 10000))
    {
        flint_printf("exponent of Cl(-(10^18 + 3)) = ");
        fmpz_print(e);
        flint_printf("\n");
    }

    qfb_clear(f);
    qfb_clear(g);
    fmpz_clear(D);
    fmpz_clear(e);
    return 0;
}
```
```output
h(-56) = 4: (1, 0, 14) (2, 0, 7) (3, 2, 5) (3, -2, 5)
powers of (3, 2, 5): (3, 2, 5) (2, 0, 7) (3, -2, 5) (1, 0, 14)
h(-1000003) = 105, exponent 105
exponent of Cl(-(10^18 + 3)) = 102764373
```
