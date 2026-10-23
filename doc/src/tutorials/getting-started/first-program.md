---
title: A first program
summary: Integers and polynomials with fmpz and fmpz_poly: declaring, initialising and clearing variables, arithmetic, printing and factoring, and how to compile, link and run the program.
order: 1
---

This tutorial writes two short programs, one with integers and one with
polynomials, and shows how to compile and run them. It assumes that FLINT
is installed as described in [](guide/getting-started/building).

## Integers {#integers}

FLINT's integers have the type `fmpz_t` and are provided by the module
[](fmpz). Unlike C integers, they have no size limit. The program below
computes $30!$, a greatest common divisor, and two factorisations.

```c example name=tutorial_first_fmpz
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>

int main(void)
{
    fmpz_t a, b, g;
    fmpz_factor_t fac;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(g);
    fmpz_factor_init(fac);

    fmpz_fac_ui(a, 30);                 /* a = 30! */
    flint_printf("30! = %{fmpz}\n", a);
    flint_printf("30! has %wu bits\n", fmpz_bits(a));

    fmpz_one(b);
    fmpz_mul_2exp(b, b, 60);
    fmpz_sub_ui(b, b, 1);               /* b = 2^60 - 1 */
    fmpz_gcd(g, a, b);
    flint_printf("gcd(30!, 2^60 - 1) = %{fmpz}\n", g);

    fmpz_factor(fac, b);
    flint_printf("2^60 - 1 = ");
    fmpz_factor_print(fac);
    flint_printf("\n");

    fmpz_add_ui(b, b, 1);
    fmpz_mul_2exp(b, b, 4);
    fmpz_add_ui(b, b, 1);               /* b = 2^64 + 1 */
    fmpz_factor(fac, b);
    flint_printf("2^64 + 1 = ");
    fmpz_factor_print(fac);
    flint_printf("\n");

    fmpz_factor_clear(fac);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(g);
    return 0;
}
```

```output
30! = 265252859812191058636308480000000
30! has 108 bits
gcd(30!, 2^60 - 1) = 225225
2^60 - 1 = 3^2 * 5^2 * 7 * 11 * 13 * 31 * 41 * 61 * 151 * 331 * 1321
2^64 + 1 = 274177 * 67280421310721
```

A few things to notice:

- Each module has its own header, included as `<flint/module.h>`;
  `flint.h` is included by all of them but is often listed for clarity.
- Every variable is initialised with an `init` function before it is used
  and cleared with the matching `clear` function at the end. The `fmpz_t`
  type is an array of length one, so variables are passed to functions
  without `&`.
- Functions write their result to the first argument, as in
  `fmpz_gcd(g, a, b)`, which sets $g = \gcd(a, b)$. An output may also be
  an input: `fmpz_sub_ui(b, b, 1)` subtracts 1 from $b$ in place.
- Function names tell the types of the arguments: `fmpz_sub_ui`
  subtracts a `ulong`, and `fmpz_fac_ui` takes a `ulong` argument.
- `flint_printf` is `printf` with two additions: `%wu` and `%wd` print
  the word types `ulong` and `slong` (which `fmpz_bits` returns), and
  `%{fmpz}` prints an integer; see [](flint/printf).
- An `fmpz_factor_t` holds a factorisation. `fmpz_factor` overwrites
  it, so the same variable can be reused.

The conventions behind these rules are explained in
[](guide/conventions/types-and-naming) and
[](guide/conventions/memory).

## Polynomials {#polynomials}

Polynomials with integer coefficients have the type `fmpz_poly_t`
([](fmpz_poly)). This program builds $f = x^{12} - 1$ and
$g = (x + 1)^3$, divides one by the other, computes their GCD, evaluates
$f$ at 10 and factors $f$ into irreducible polynomials with
[](fmpz_poly_factor).

```c example name=tutorial_first_fmpz_poly
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

int main(void)
{
    fmpz_poly_t f, g, q, r;
    fmpz_poly_factor_t fac;
    fmpz_t y;
    slong i;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(q);
    fmpz_poly_init(r);
    fmpz_poly_factor_init(fac);
    fmpz_init(y);

    /* f = x^12 - 1 */
    fmpz_poly_set_coeff_si(f, 12, 1);
    fmpz_poly_set_coeff_si(f, 0, -1);

    /* g = (x + 1)^3 */
    fmpz_poly_set_str(g, "2  1 1");
    fmpz_poly_pow(g, g, 3);

    flint_printf("f = "); fmpz_poly_print_pretty(f, "x"); flint_printf("\n");
    flint_printf("g = "); fmpz_poly_print_pretty(g, "x"); flint_printf("\n");

    fmpz_poly_divrem(q, r, f, g);
    flint_printf("f mod g = "); fmpz_poly_print_pretty(r, "x"); flint_printf("\n");

    fmpz_poly_gcd(r, f, g);
    flint_printf("gcd(f, g) = "); fmpz_poly_print_pretty(r, "x"); flint_printf("\n");

    fmpz_set_ui(y, 10);
    fmpz_poly_evaluate_fmpz(y, f, y);
    flint_printf("f(10) = %{fmpz}\n", y);

    fmpz_poly_factor(fac, f);
    flint_printf("f has %wd irreducible factors:\n", fac->num);
    for (i = 0; i < fac->num; i++)
    {
        flint_printf("  ");
        fmpz_poly_print_pretty(fac->p + i, "x");
        flint_printf("\n");
    }

    fmpz_clear(y);
    fmpz_poly_factor_clear(fac);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(q);
    fmpz_poly_clear(r);
    return 0;
}
```

```output
f = x^12-1
g = x^3+3*x^2+3*x+1
f mod g = 66*x^2+120*x+54
gcd(f, g) = x+1
f(10) = 999999999999
f has 6 irreducible factors:
  x-1
  x^2+x+1
  x+1
  x^2-x+1
  x^2+1
  x^4-x^2+1
```

Coefficients are set one at a time with `fmpz_poly_set_coeff_si`, or the
whole polynomial is read from a string with `fmpz_poly_set_str`, whose
format is the length, two spaces, and the coefficients from the constant
term up: `"2  1 1"` is $1 + x$. Polynomials grow automatically as
coefficients are set. The factorisation found by `fmpz_poly_factor`
consists of the content and the irreducible factors with
multiplicities `fac->exp[i]`; here they are the cyclotomic polynomials
$\Phi_d$ for the divisors $d$ of 12, each with multiplicity 1. The
factors are entries of an array of `fmpz_poly_struct`, so the $i$-th one
is passed as `fac->p + i`.

## Compiling and running {#compiling}

Save a program as `first.c` and compile it with

```
gcc first.c -o first -lflint
```

then run it with `./first`. If FLINT was installed under a prefix that
the compiler does not search by default, such as `$HOME/local`, give the
directories of the headers and the library, and tell the dynamic linker
where to find the library at run time:

```
gcc first.c -o first -I$HOME/local/include -L$HOME/local/lib -lflint \
    -Wl,-rpath,$HOME/local/lib
```

Instead of `-Wl,-rpath`, the directory can be added to the environment
variable `LD_LIBRARY_PATH` (`DYLD_LIBRARY_PATH` on macOS) before running
the program. With pkg-config, `gcc first.c -o first $(pkg-config
--cflags --libs flint)` supplies the flags. If FLINT was built only as a
static library, add the libraries it depends on:
`-lflint -lmpfr -lgmp -lpthread -lm`. Common errors are:

- `fatal error: flint/fmpz.h: No such file or directory`: the include
  directory is missing (`-I`);
- `undefined reference to fmpz_init`: `-lflint` is missing or comes
  before `first.c` on the command line, or the library directory is
  missing (`-L`);
- `error while loading shared libraries: libflint.so`: the program was
  linked but the dynamic linker cannot find the library at run time
  (`-Wl,-rpath` or `LD_LIBRARY_PATH`).

While developing, it is worth checking programs with a memory checker:
`valgrind ./first` reports missing `clear` calls as leaks. Calling
`flint_cleanup_master()` at the end of `main` frees FLINT's internal
caches, so that only genuine leaks are reported.

The next tutorial, [](tutorials/getting-started/multimodular), combines
integer matrices, modular arithmetic and the Chinese remainder theorem.
