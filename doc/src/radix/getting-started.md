---
title: Arithmetic in radix b^e
summary: Choosing the radix, working with radix_integer_t, and cheap operations on digits, with decimal and 7-adic examples.
order: 1
---

The `radix` module does exact integer arithmetic in a radix other than a
power of two. A radix is described by
a `radix_t`, created with `radix_init(radix, b, e)` for the digit radix
$b$ and $e$ digits per limb. Passing $e = 0$ chooses the largest $e$ with
$b^e < 2^{\bits}$, which is almost always the best choice: $10^{19}$ for
decimal arithmetic and $7^{22}$ for 7-adic arithmetic on 64-bit
machines.

Integers are of type `radix_integer_t`, which works like `fmpz_t` but
needs the radix as an extra argument to every function (the radix is not
stored in the integer). They can also be used through the generic ring
interface with `gr_ctx_init_radix_integer`.

## What is cheap {#cheap}

Multiplication and division cost a little more than in binary, because
the carries are propagated in radix $B$; `radix_mul` uses the
small-prime FFT for large operands, so the asymptotic cost is the same.
What becomes cheap is anything expressed in digits:

- the number of digits (`radix_integer_size_digits`), a given digit
  (`radix_integer_get_digit`), and the number of trailing zero digits
  (`radix_integer_valuation_digits`);
- truncation to the last $n$ digits (`radix_integer_mod_digits`,
  `radix_integer_trunc_digits`) and shifts by digits
  (`radix_integer_lshift_digits`, `radix_integer_rshift_digits`), which
  in binary would need a division by $10^n$;
- conversion to a decimal string, in linear time for $b = 10$
  (`radix_get_str_decimal`);
- arithmetic modulo $b^n$: low products (`radix_integer_mullow_limbs`)
  and inverses, quotients and square roots modulo $B^n$ by Hensel lifting
  (`radix_integer_invmod_limbs` and the following functions), which for a
  prime $b = p$ is arithmetic in the $p$-adic integers.

Conversion to and from `fmpz_t` takes quasi-linear time
([](radix/conversion)).

## Example {#example}

Compute $2^{200}$ in decimal, inspect its digits, and compute a square
root of 2 in the 7-adic integers modulo $7^{22}$:

```c example name=radix_getting_started
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/radix.h>

int main(void)
{
    radix_t dec, sev;
    radix_integer_t x, y, t;
    fmpz_t z;
    char * s;
    slong i;

    radix_init(dec, 10, 0);
    radix_integer_init(x, dec);
    radix_integer_init(y, dec);

    radix_integer_set_ui(x, 1, dec);
    radix_integer_set_ui(y, 2, dec);
    for (i = 0; i < 200; i++)
        radix_integer_mul(x, x, y, dec);

    s = radix_get_str_decimal(NULL, x->d, x->size, 0, dec);
    flint_printf("2^200 = %s\n", s);
    flint_free(s);
    flint_printf("%wd digits in %wd limbs\n",
        radix_integer_size_digits(x, dec), radix_integer_size_limbs(x, dec));

    radix_integer_mod_digits(y, x, 10, dec);
    s = radix_get_str_decimal(NULL, y->d, y->size, 0, dec);
    flint_printf("last 10 digits: %s\n", s);
    flint_free(s);

    radix_integer_rshift_digits(y, x, 50, dec);
    s = radix_get_str_decimal(NULL, y->d, y->size, 0, dec);
    flint_printf("first 11 digits: %s\n", s);
    flint_free(s);

    /* 7-adic square root of 2 to 22 digits */
    radix_init(sev, 7, 0);
    radix_integer_init(t, sev);
    fmpz_init(z);
    radix_integer_set_ui(t, 2, sev);
    radix_integer_sqrtmod_limbs(t, t, 1, sev);
    radix_integer_get_fmpz(z, t, sev);
    flint_printf("sqrt(2) mod 7^%wu = %{fmpz}\n", radix_limb_exponent(sev), z);
    flint_printf("lowest 7-adic digit: %wu\n", radix_integer_get_digit(t, 0, sev));

    radix_integer_mullow_limbs(t, t, t, 1, sev);
    radix_integer_get_fmpz(z, t, sev);
    flint_printf("its square: %{fmpz}\n", z);

    fmpz_clear(z);
    radix_integer_clear(x, dec);
    radix_integer_clear(y, dec);
    radix_integer_clear(t, sev);
    radix_clear(dec);
    radix_clear(sev);
    return 0;
}
```
```output
2^200 = 1606938044258990275541962092341162602522202993782792835301376
61 digits in 4 limbs
last 10 digits: 2835301376
first 11 digits: 16069380442
sqrt(2) mod 7^22 = 1830612359265707720
lowest 7-adic digit: 3
its square: 2
```

Of the two square roots of 2 modulo $7^{22}$, the function returns the
one congruent to 3 (rather than 4) modulo 7, the smaller residue.
