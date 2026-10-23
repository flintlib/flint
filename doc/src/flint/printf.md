---
title: Formatted output
functions: flint_printf
---

`flint_printf` and its variants work like the C library's `printf`
family, with two extensions: the length modifier `w` for the word types
`ulong` and `slong`, and the conversion `%{type}`, which prints a FLINT
or GMP object. Ordinary conversions such as `%d`, `%s`, `%f` and `%p`,
with flags, field widths and precisions, are passed to the C library and
behave as usual.

## Words {#words}

`%wd`, `%wi`, `%wu`, `%wx`, `%wX` and `%wo` print an `slong` or a `ulong`
with the conversion `d`, `i`, `u`, `x`, `X` or `o`. Since `ulong` is
`unsigned long` on most platforms but `unsigned long long` on 64-bit
Windows, `%lu` is not portable, while `%wu` is. The macros `WORD_FMT` and
`UWORD` serve the same purpose for the standard `printf` and for
constants.

## FLINT types {#types}

A conversion `%{type}` consumes one argument, the object to print, passed
as for any function taking that type (for example an `fmpz_t`, which
decays to a pointer). The supported types are:

| Conversion | Argument | Prints |
|---|---|---|
| `%{ulong}`, `%{slong}` | `ulong`, `slong` | the integer |
| `%{fmpz}`, `%{fmpq}` | `fmpz_t`, `fmpq_t` | the number, a fraction as `p / q` |
| `%{mpz}`, `%{mpq}` | `mpz_t`, `mpq_t` | the number |
| `%{mag}`, `%{arf}`, `%{arb}`, `%{acb}` | `mag_t`, `arf_t`, `arb_t`, `acb_t` | the number, as by the module's `_print` function |
| `%{nmod}` | `nmod_t` | the modulus, as `mod n` |
| `%{fmpz_mod_ctx}` | `fmpz_mod_ctx_t` | the modulus |
| `%{truth}` | `truth_t` | `T_TRUE`, `T_FALSE` or `T_UNKNOWN` |
| `%{nmod_poly}`, `%{fmpz_poly}`, `%{fmpz_mod_poly}`, `%{fmpq_poly}`, `%{arb_poly}`, `%{acb_poly}` | the polynomial | the polynomial in `x` |
| `%{nmod_mat}`, `%{fmpz_mat}`, `%{fmpz_mod_mat}`, `%{fmpq_mat}`, `%{arb_mat}`, `%{acb_mat}` | the matrix | the matrix, row by row |

Vectors are printed with `%{type*}`, which consumes two arguments, the
pointer to the first entry and the length as an `slong`, and prints
`[a, b, ...]`. It is available for `ulong` (also for vectors of residues
of `nmod` type), `slong`, `fmpz` (also for `fmpz_mod` vectors), `fmpq`,
`mag`, `arf`, `arb` and `acb`. The length must be an `slong`: passing an
`int` such as the literal `3` reads garbage on platforms where `int` is
narrower than `slong`, so write `(slong) 3`.

Elements of [generic rings](gr) are printed with `%{gr}` (element and
context), `%{gr*}` (pointer, length and context), `%{gr_poly}`,
`%{gr_ore_poly}` and `%{gr_mat}` (object and context), and the context
itself with `%{gr_ctx}`. The context argument always comes after the
object, so printing an element and its ring takes the context twice:
`flint_printf("%{gr} in %{gr_ctx}", x, ctx, ctx)`.

## Example {#example}

```c example name=flint_printf_types
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/nmod.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_poly.h>
#include <flint/arb.h>
#include <flint/gr.h>

int main(void)
{
    ulong u = UWORD_MAX;
    slong s = WORD_MIN;
    fmpz_t a;
    fmpq_t q;
    nmod_t mod;
    fmpz v[3];
    slong w[4] = {3, 1, 4, 1};
    fmpz_poly_t f;
    nmod_poly_t g;
    arb_t x;
    gr_ctx_t ZZ;
    char buf[32];
    int n;

    fmpz_init(a);
    fmpz_ui_pow_ui(a, 3, 50);
    fmpq_init(q);
    fmpq_set_si(q, -22, 7);
    nmod_init(&mod, 17);
    fmpz_init_set_si(v + 0, 1);
    fmpz_init_set_si(v + 1, -2);
    fmpz_init_set_si(v + 2, 3);
    fmpz_poly_init(f);
    fmpz_poly_set_str(f, "3  1 0 -2");
    nmod_poly_init(g, 17);
    nmod_poly_set_coeff_ui(g, 2, 5);
    nmod_poly_set_coeff_ui(g, 0, 1);
    arb_init(x);
    arb_set_ui(x, 2);
    arb_sqrt(x, x, 30);
    gr_ctx_init_fmpz(ZZ);

    flint_printf("%wu %wd %wx\n", u, s, u);
    flint_printf("%{ulong} %{slong}\n", u, s);
    flint_printf("a = %{fmpz}, q = %{fmpq}, %{nmod}\n", a, q, mod);
    flint_printf("v = %{fmpz*}, w = %{slong*}\n", v, (slong) 3, w, (slong) 4);
    flint_printf("f = %{fmpz_poly}, g = %{nmod_poly}\n", f, g);
    flint_printf("x = %{arb}\n", x);
    flint_printf("%{gr} in %{gr_ctx}\n", a, ZZ, ZZ);
    flint_printf("%-8s|%5.2f|%{fmpz}|\n", "text", 3.14159, v + 1);
    n = flint_snprintf(buf, sizeof(buf), "%{fmpz}", a);
    flint_printf("%d characters: %s\n", n, buf);

    gr_ctx_clear(ZZ);
    arb_clear(x);
    nmod_poly_clear(g);
    fmpz_poly_clear(f);
    fmpz_clear(v + 0);
    fmpz_clear(v + 1);
    fmpz_clear(v + 2);
    fmpq_clear(q);
    fmpz_clear(a);
    return 0;
}
```

```output
18446744073709551615 -9223372036854775808 ffffffffffffffff
18446744073709551615 -9223372036854775808
a = 717897987691852588770249, q = -22 / 7, mod 17
v = [1, -2, 3], w = [3, 1, 4, 1]
f = -2 * x^2 + 1, g = 5 * x^2 + 1
x = [1.41421 +/- 3.57e-6]
717897987691852588770249 in Integer ring (fmpz)
text    | 3.14|-2|
24 characters: 717897987691852588770249
```

## Limitations {#limitations}

- Flags, field widths and precisions are not supported for `%{type}`
  conversions: `%5{slong}` is printed literally, like an unknown type
  name such as `%{foo}`, and no argument is consumed.
- Formats with a `*` width or precision mixed with unusual combinations of
  flags (such as `"%* p"`) may be parsed differently from the C library.
- `%n` is not supported and is not reported as an error.
- The compiler cannot check the arguments against the format string, so
  mismatched types are not diagnosed.
