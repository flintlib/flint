---
title: Getting started with fmpz
summary: How to create, compute with, print and free integers, the rules for aliasing, how to work with arrays of integers, and how to convert to and from machine types, strings and GMP.
order: 1
---

An `fmpz_t` holds an integer of any size. Like GMP's `mpz_t`, it must be
initialised before use and cleared afterwards, every function writes its
result to its first argument, and outputs may be the same objects as
inputs. Unlike an `mpz_t`, an `fmpz_t` holding a small value is a single
machine word with no allocated memory, which makes FLINT's integers cheap
to create, copy and store in large arrays.

## A first program {#first}

The program computes with $a = 3^{50}$ and $b = -10^{21}$: it prints $a$
and its size, the product $ab$, the quotient and remainder of floor
division of $a$ by $b$, and $2a$ in hexadecimal.

```c example name=fmpz_first
#include <flint/flint.h>
#include <flint/fmpz.h>

int main(void)
{
    fmpz_t a, b, q, r;
    char * s;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(q);
    fmpz_init(r);

    fmpz_set_ui(a, 3);
    fmpz_pow_ui(a, a, 50);                  /* a = 3^50 */
    fmpz_set_str(b, "-1000000000000000000000", 10);

    flint_printf("a = "); fmpz_print(a); flint_printf("\n");
    flint_printf("bits(a) = %wu, digits(a) = %zu\n",
                 fmpz_bits(a), fmpz_sizeinbase(a, 10));

    fmpz_mul(r, a, b);
    flint_printf("a * b = "); fmpz_print(r); flint_printf("\n");

    fmpz_fdiv_qr(q, r, a, b);               /* floor division */
    flint_printf("q = "); fmpz_print(q);
    flint_printf(", r = "); fmpz_print(r); flint_printf("\n");

    fmpz_add(a, a, a);                      /* outputs may alias inputs */
    s = fmpz_get_str(NULL, 16, a);
    flint_printf("2a in hex = %s\n", s);
    flint_free(s);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(q);
    fmpz_clear(r);
    return 0;
}
```
```output
a = 717897987691852588770249
bits(a) = 80, digits(a) = 24
a * b = -717897987691852588770249000000000000000000000
q = -718, r = -102012308147411229751
2a in hex = 1300aa7e1b65fa13bc792
```

Floor division rounds the quotient $a/b \approx -717.9$ down to $-718$,
so the remainder $a - qb$ has the sign of the divisor. The other rounding
conventions are explained in [](fmpz/division).

## Conventions {#conventions}

**Initialisation.** `fmpz_init` sets a variable to zero; it only writes a
word and cannot fail. `fmpz_clear` releases the memory of a large value.
A variable that is never cleared leaks memory only if it ever held a
large value, but every variable should be cleared. `fmpz_init_set`,
`fmpz_init_set_ui` and `fmpz_init_set_si` combine initialisation with an
assignment.

**Outputs first.** Functions write their results to the first
arguments: `fmpz_sub(r, a, b)` sets $r = a - b$. Any output may be the
same variable as any input, so `fmpz_mul(a, a, a)` squares $a$ in place.
A function with two outputs, such as `fmpz_tdiv_qr`, needs two distinct
variables for them.

**Variants.** Most operations have versions taking an `ulong` (`_ui`) or
an `slong` (`_si`) operand, which avoid creating an `fmpz_t` for a
machine-size value: `fmpz_add_ui`, `fmpz_mul_si`, `fmpz_divexact_ui`,
`fmpz_cmp_si`, `fmpz_pow_ui`. Multiplication and division by powers of
two have `_2exp` versions.

**Exact operations.** When a division is known to be exact, use
`fmpz_divexact`, which is considerably faster than a general division; it
does not check exactness. `fmpz_div` checks it and returns a flag.
`fmpz_addmul` and `fmpz_submul` update an accumulator without a
temporary for the product.

**Errors.** Division by zero, roots of negative numbers and similar
errors raise an exception through `flint_throw`, which by default prints
a message and aborts. A few functions document unchecked preconditions
instead, such as `fmpz_mod` (nonzero modulus) and `fmpz_divexact`
(exactness).

## Arrays of integers {#arrays}

An array of integers is an array of `fmpz`, and `v + i` (a pointer to an
entry) can be passed wherever an `fmpz_t` is expected. Since an `fmpz`
with value zero is a zero word, `_fmpz_vec_init` allocates an array with
`flint_calloc`, and `_fmpz_vec_clear` clears every entry and frees it.
The vector functions of [](fmpz_vec) work on such arrays.

```c example name=fmpz_arrays
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>

int main(void)
{
    slong i, n = 10;
    fmpz * v = _fmpz_vec_init(n);       /* n integers, all zero */
    fmpz_t s;

    fmpz_init(s);

    for (i = 0; i < n; i++)
        fmpz_fib_ui(v + i, 10 * i);     /* v + i is an fmpz_t */

    for (i = 0; i < n; i++)
        fmpz_addmul(s, v + i, v + i);   /* s += v[i]^2 */

    flint_printf("F_90 = "); fmpz_print(v + 9); flint_printf("\n");
    flint_printf("sum of squares = "); fmpz_print(s); flint_printf("\n");
    flint_printf("v[9] is %s\n", COEFF_IS_MPZ(v[9]) ? "an mpz" : "inline");

    fmpz_clear(s);
    _fmpz_vec_clear(v, n);              /* clears every entry, then frees */
    return 0;
}
```
```output
F_90 = 2880067194370816120
sum of squares = 8295335423509227258649641464938952750
v[9] is inline
```

$F_{90} \approx 2.9 \cdot 10^{18}$ is below `COEFF_MAX` $= 2^{62} - 1
\approx 4.6 \cdot 10^{18}$, so it occupies just the eight bytes of its
array entry; the sum of squares needs an `mpz`. An array declared as
`fmpz v[10]` or obtained from `flint_malloc` works the same way after each
entry is initialised, and each entry must be cleared before the memory is
released. Copying an `fmpz` word with `=` copies a pointer for a large
value, so two entries would share one `mpz`; use `fmpz_set` or
`fmpz_swap` instead.

## Conversions {#conversions}

`fmpz_set_si`, `fmpz_set_ui`, `fmpz_set_d` and `fmpz_set_str` create
integers from machine types and strings. In the other direction,
`fmpz_get_si`, `fmpz_get_ui` and `fmpz_get_d` give unspecified results
when the value does not fit, so test first with `fmpz_fits_si` or
`fmpz_abs_fits_ui`, or use the functions for two-word values
(`fmpz_get_uiui`, `fmpz_get_signed_uiui`), limb arrays
(`fmpz_get_ui_array`) or a double with a separate exponent
(`fmpz_get_d_2exp`).

```c example name=fmpz_conversions
#include <flint/flint.h>
#include <flint/fmpz.h>

static void show(const fmpz_t x)
{
    fmpz_print(x);
    if (fmpz_fits_si(x))
        flint_printf(": slong %wd", fmpz_get_si(x));
    else if (fmpz_abs_fits_ui(x) && fmpz_sgn(x) > 0)
        flint_printf(": ulong %wu", fmpz_get_ui(x));
    else
        flint_printf(": too large, double %g", fmpz_get_d(x));
    flint_printf("\n");
}

int main(void)
{
    fmpz_t x;
    fmpz_init(x);

    fmpz_set_si(x, -42);
    show(x);
    fmpz_set_ui(x, UWORD_MAX);
    show(x);
    fmpz_add_ui(x, x, 1);
    show(x);
    fmpz_set_d(x, -2.75);                   /* rounds towards zero */
    show(x);

    fmpz_clear(x);
    return 0;
}
```
```output
-42: slong -42
18446744073709551615: ulong 18446744073709551615
18446744073709551616: too large, double 1.84467e+19
-2: slong -2
```

Strings in bases 2 to 62 are handled by `fmpz_set_str` and
`fmpz_get_str`; [](fmpz/strings) describes them and the conversion
algorithms for large numbers.

## Working with GMP {#gmp}

FLINT is built on GMP, and its headers declare the conversion functions
`fmpz_get_mpz` and `fmpz_set_mpz` (and `fmpz_get_mpfr` for MPFR) only if
`gmp.h` (respectively `mpfr.h`) is included before them. These functions
copy the value. To pass an `fmpz_t` to a GMP function without copying,
use `flint_mpz_init_set_readonly`, described with the representation in
[](fmpz/representation#gmp).

For a computation that is mostly integer arithmetic, `fmpz_t` is usually
faster than `mpz_t`: small values need no allocation and many operations
on them are inline, and for large values FLINT uses its own
multiplication, division and square root code where it is faster than GMP's
([](algorithms/integers/multiplication)).
