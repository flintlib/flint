---
title: Decimal input and output
functions: arb_set_str arb_get_str arb_printn arb_printd arb_print
---

Balls are stored in binary, but are usually read and written in decimal.
Converting between the two is inexact in general, and FLINT provides
conversions that take the conversion error into account, so that a
decimal string obtained from a ball, or a ball obtained from a decimal
string, is still a rigorous enclosure.

## The decimal ball format {#format}

`arb_get_str` and `arb_printn` produce a *decimal ball* `[m +/- r]`
whose midpoint $m$ has at most `n` significant digits and whose radius
$r$ has three. The guarantee is that the decimal interval $[m - r, m +
r]$ contains the binary ball, so it can be read back with `arb_set_str`
without loss of rigour. By default only digits that are correct up to
one unit in the last place are printed: the number of digits is the
largest for which the printed radius, including the error of rounding
the midpoint to that many digits, is at most one unit in the last digit.
Consequently

- an exact ball whose value has at most `n` significant digits is
  printed without radius, for example `0` or `1.5` or
  `1.000000000000000000000000000000e+100`;
- a ball with a tiny radius is printed with the radius it has, not with
  the rounding error of the last printed digit: $[1 \pm 2^{-1000}]$ at 30
  digits gives `[1.00000000000000000000000000000 +/- 9.34e-302]`;
- a ball without a single correct digit is printed as `[+/- r]`;
- an indeterminate ball is printed as `nan`, and any ball with an
  infinite midpoint or radius as `[+/- inf]`, which drops the sign of an
  infinite midpoint (the string is still a valid enclosure).

The flags change this behaviour. `ARB_STR_MORE` prints up to `n` digits
even if they are not all correct (the radius still makes the output
rigorous). `ARB_STR_NO_RADIUS` omits the radius; the digits printed are
still correct to one unit in the last place, and a ball with no correct
digit is printed as `0e+n`, meaning a number of absolute value below
$10^n$. A multiple $k$ of `ARB_STR_CONDENSE` abbreviates runs of more
than $3k$ digits.

```c example name=arb_printing
#include <flint/flint.h>
#include <flint/arb.h>

int main(void)
{
    arb_t x, y;
    char * s;

    arb_init(x); arb_init(y);

    arb_const_pi(x, 100);
    arb_div_ui(x, x, 1000, 100);

    arb_printn(x, 10, 0);                  flint_printf("\n");
    arb_printn(x, 50, 0);                  flint_printf("\n");
    arb_printn(x, 50, ARB_STR_MORE);       flint_printf("\n");
    arb_printn(x, 10, ARB_STR_NO_RADIUS);  flint_printf("\n");
    arb_printd(x, 10);                     flint_printf("\n");
    arb_print(x);                          flint_printf("\n");

    /* a wide ball: only the first digit is certain */
    arb_set_str(y, "[3.14 +/- 0.2]", 53);
    arb_printn(y, 10, 0);                  flint_printf("\n");
    arb_printn(y, 10, ARB_STR_NO_RADIUS);  flint_printf("\n");

    /* round trip through a decimal string */
    s = arb_get_str(x, 20, 0);
    arb_set_str(y, s, 100);
    flint_printf("%s -> contains x: %d\n", s, arb_contains(y, x));
    flint_free(s);

    /* condensed output of a long number */
    arb_const_pi(x, 3000);
    arb_printn(x, 900, ARB_STR_CONDENSE * 10); flint_printf("\n");

    arb_clear(x); arb_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
[0.003141592654 +/- 4.11e-13]
[0.00314159265358979323846264338328 +/- 7.52e-33]
[0.0031415926535897932384626433832771392962513160822984 +/- 4.66e-33]
0.003141592654
0.003141592654 +/- 4.6592e-33
(1019505104126898525104217885171 * 2^-108) +/- (811748820 * 2^-137)
[3e+0 +/- 0.341]
3e+0
[0.0031415926535897932385 +/- 3.74e-23] -> contains x: 1
[3.1415926535{...879 digits...}7766914730 +/- 3.61e-900]
```

The first line shows how the radius absorbs the rounding of the midpoint
to 10 digits: the binary ball has radius $4.7 \cdot 10^{-33}$, but the
decimal ball needs $4.11 \cdot 10^{-13}$. At 50 digits only 33 are
certain, and `ARB_STR_MORE` shows the remaining digits of the binary
midpoint. The wide ball $[3.14 \pm 0.2]$ has a single certain digit; the
radius $0.341$ covers both the original radius and the distance from
3.14 to 3.

## Other output functions {#other}

`arb_printd` prints decimal approximations of the midpoint and the radius
(with 5 digits). It is cheap and shows the binary radius directly, but
the printed interval is not guaranteed to contain the ball. `arb_print`
prints the exact binary representation `(m * 2^e) +/- (r * 2^f)` with
integer mantissas, which is what to look at when debugging. `arb_dump_str`
and `arb_load_str` serialise a ball exactly in a compact hexadecimal
format, for storing values and reading them back without any rounding.

## Reading decimal input {#input}

`arb_set_str` accepts plain decimal numbers (`"0.1"`, `"-7e+141"`),
balls in the output format above (`"[0.1 +/- 1e-20]"`, also without the
brackets, or `"[+/- 1e-20]"`), and `"inf"`, `"-inf"` and `"nan"`. The
midpoint and the radius are converted rigorously: the result is a ball
with a `prec`-bit midpoint that contains the decimal ball. It returns 0
on success; on a syntax error it returns nonzero and sets the ball to an
indeterminate value.

To enter a decimal constant, `arb_set_str` is the right tool, while
`arb_set_d` with a C literal such as `0.1` gives the binary `double`
nearest to it, a different number:

```c
arb_set_str(x, "0.1", prec);   /* a ball containing 1/10 */
arb_set_d(x, 0.1);             /* exactly 0.1000000000000000055511151231257827... */
```

Rational numbers can also be entered exactly with `arb_set_fmpq` or by
dividing integers, for example `arb_set_ui(x, 1); arb_div_ui(x, x, 10,
prec)`.
