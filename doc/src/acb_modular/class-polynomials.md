---
title: Hilbert class polynomials
functions: acb_modular_hilbert_class_poly
---

For a negative discriminant $D$ (that is, $D < 0$ and $D \equiv 0, 1 \pmod 4$),
the Hilbert class polynomial is

$$
H_D(x) = \prod_{(a,b,c)} \left(x - j\left(\frac{-b+\sqrt{D}}{2a}\right)\right)
$$

where $(a,b,c)$ ranges over the primitive reduced positive definite
binary quadratic forms $ax^2 + bxy + cy^2$ of discriminant
$b^2 - 4ac = D$. It has integer coefficients and degree equal to the
class number $h(D)$, and its roots generate the ring class field of the
order of discriminant $D$ [@Coh1996]. It is used to construct
elliptic curves with complex multiplication, for instance in elliptic
curve primality proving.

`acb_modular_hilbert_class_poly` computes $H_D$ numerically
[@Enge2009]:

1. It enumerates the reduced forms with $|b| \le a \le c$
   [@Coh1996, Algorithm 5.3.5], pairing $(a, b, c)$ with $(a, -b, c)$,
   whose roots are complex conjugate.
2. It estimates the bit size of the coefficients from
   $\log |H_D(0)| \lesssim \pi \sqrt{|D|} \sum 1/a$ and chooses a working
   precision 0.5% above it, plus 20 bits.
3. It evaluates the roots with `acb_modular_j`, forming a real quadratic
   factor $x^2 - 2\operatorname{Re}(j) x + |j|^2$ for each conjugate pair,
   and multiplies the factors by binary splitting with `arb_poly_mul`,
   in parallel with `flint_parallel_binary_splitting` when several
   threads are available.
4. It converts the result with `arb_poly_get_unique_fmpz_poly`. Since
   the real polynomial is an enclosure of $H_D$, this succeeds only when
   every coefficient ball contains a unique integer, and then the result
   is proven correct. If the enclosure is too wide, it prints a message
   to standard output and retries with 1.2 times the precision.

Each $j$-value costs $O(\sqrt{p})$ multiplications at precision $p$, and
$p$ grows like $\sqrt{|D|} \log |D|$; this simple method is practical
for discriminants up to several million. Faster methods for large $D$ use
class invariants smaller than $j$ (Weber functions and other modular
functions) and the Chinese remainder theorem [@EngeCM].

If $D \ge 0$ or $D \equiv 2, 3 \pmod 4$, the output is the zero
polynomial.

## Example {#example}

The class polynomial of discriminant $-23$ (class number 3), the size of
the class polynomial of discriminant $-1000003$, and the near-integer
$j((1 + \sqrt{-163})/2) = -640320^3$ for the class number one
discriminant $-163$.

```c example name=acb_modular_hilbert
#include <flint/fmpz_poly.h>
#include <flint/acb.h>
#include <flint/acb_modular.h>

int main(void)
{
    fmpz_poly_t H;
    acb_t tau, j;

    fmpz_poly_init(H);
    acb_init(tau);
    acb_init(j);

    acb_modular_hilbert_class_poly(H, -23);
    fmpz_poly_print_pretty(H, "x");
    flint_printf("\n");

    acb_modular_hilbert_class_poly(H, -1000003);
    flint_printf("D = -1000003: degree %wd, %wd bits\n",
        fmpz_poly_degree(H), FLINT_ABS(fmpz_poly_max_bits(H)));

    arb_set_ui(acb_imagref(tau), 163);
    arb_sqrt(acb_imagref(tau), acb_imagref(tau), 128);
    arb_one(acb_realref(tau));
    acb_mul_2exp_si(tau, tau, -1);
    acb_modular_j(j, tau, 128);
    acb_printn(j, 30, 0);
    flint_printf("\n");

    fmpz_poly_clear(H);
    acb_clear(tau);
    acb_clear(j);
    return 0;
}
```
```output
x^3+3491750*x^2-5151296875*x+12771880859375
D = -1000003: degree 105, 8527 bits
[-262537412640768000.000000000000 +/- 2.85e-19]
```
