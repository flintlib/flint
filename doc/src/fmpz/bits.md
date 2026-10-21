---
title: Bit operations
functions: fmpz_tstbit fmpz_setbit fmpz_complement fmpz_and fmpz_popcnt
---

The bit operations treat an integer as an infinite sequence of bits in
two's complement, as GMP does: a nonnegative integer has finitely many one
bits followed by zeros, and a negative integer $x$ is the bit sequence of
$2^N + x$ for any large enough $N$, that is, it has infinitely many
leading ones. With this convention $-x = \mathrm{not}(x) + 1$, the
bitwise operations satisfy the usual identities, and `fmpz_fdiv_q_2exp`
and `fmpz_fdiv_r_2exp` are the right shift and the extraction of the low
bits.

## Semantics {#semantics}

- `fmpz_tstbit(f, i)` returns bit $i$, which for $i$ beyond the length of
  $|f|$ is 0 for $f \ge 0$ and 1 for $f < 0$.
- `fmpz_setbit`, `fmpz_clrbit` and `fmpz_combit` change one bit, so they
  add or subtract $2^i$. Setting a bit of a negative number beyond its
  length has no effect, and clearing one subtracts $2^i$.
- `fmpz_complement(r, f)` is the bitwise not, $r = -f - 1$.
- `fmpz_and`, `fmpz_or` and `fmpz_xor` combine the infinite sequences;
  the result is negative if the corresponding operation on the sign bits
  gives 1.
- `fmpz_popcnt(c)` returns the number of one bits of $c \ge 0$; it
  returns 0 for negative $c$, which has infinitely many (GMP's
  `mpz_popcount` returns the largest `mp_bitcnt_t` instead).

## Example {#example}

The program shows the low eight bits of the results; the dots stand for
the infinitely many copies of the leftmost bit shown.

```c example name=fmpz_bits_ops
#include <flint/flint.h>
#include <flint/fmpz.h>

static void show(const char * name, const fmpz_t x)
{
    slong i;
    flint_printf("%-12s = %5wd = ...", name, fmpz_get_si(x));
    for (i = 7; i >= 0; i--)
        flint_printf("%d", fmpz_tstbit(x, i));
    flint_printf("\n");
}

int main(void)
{
    fmpz_t a, b, r;

    fmpz_init(a); fmpz_init(b); fmpz_init(r);

    fmpz_set_si(a, -6);
    fmpz_set_si(b, 13);
    show("a", a);
    show("b", b);
    fmpz_and(r, a, b);    show("a and b", r);
    fmpz_or(r, a, b);     show("a or b", r);
    fmpz_xor(r, a, b);    show("a xor b", r);
    fmpz_complement(r, a); show("not a", r);

    fmpz_set(r, a);
    fmpz_clrbit(r, 1);    show("clrbit(a, 1)", r);
    fmpz_set(r, b);
    fmpz_setbit(r, 7);    show("setbit(b, 7)", r);

    flint_printf("bit 1000 of a = %d, of b = %d\n",
                 fmpz_tstbit(a, 1000), fmpz_tstbit(b, 1000));
    flint_printf("popcnt(b) = %wu, popcnt(a) = %wu\n",
                 fmpz_popcnt(b), fmpz_popcnt(a));

    fmpz_set(r, a);
    fmpz_combit(r, 64);   /* -6 - 2^64 */
    flint_printf("combit(a, 64) = "); fmpz_print(r); flint_printf("\n");

    fmpz_clear(a); fmpz_clear(b); fmpz_clear(r);
    return 0;
}
```
```output
a            =    -6 = ...11111010
b            =    13 = ...00001101
a and b      =     8 = ...00001000
a or b       =    -1 = ...11111111
a xor b      =    -9 = ...11110111
not a        =     5 = ...00000101
clrbit(a, 1) =    -8 = ...11111000
setbit(b, 7) =   141 = ...10001101
bit 1000 of a = 1, of b = 0
popcnt(b) = 3, popcnt(a) = 0
combit(a, 64) = -18446744073709551622
```

Complementing bit 64 of $-6$ turns a leading one into a zero, which
subtracts $2^{64}$.
