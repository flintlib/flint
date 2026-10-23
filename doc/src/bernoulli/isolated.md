---
title: Isolated Bernoulli numbers
functions: bernoulli_fmpq_ui _bernoulli_fmpq_ui_zeta bernoulli_mod_p_harvey
---

`bernoulli_fmpq_ui` is the function to use for a single Bernoulli number
$B_n$. For a range of consecutive values, `bernoulli_fmpq_vec_no_cache`
or the cache (`bernoulli_cache_compute`) is much faster than repeated
calls, and for the residue of $B_n$ modulo a prime,
`bernoulli_mod_p_harvey` avoids computing $B_n$ at all.

::: summary algorithms/number-theory/bernoulli :::

## How $B_n$ is computed {#method}

The denominator of $B_n$ is known in advance: by the theorem of von
Staudt and Clausen it is the product of the primes $p$ with $(p - 1) \mid
n$ (`arith_bernoulli_number_denom`). It has few bits, about $O(\log n
\cdot n^{\epsilon})$ in the worst case, so the work is in the numerator,
which has about $n \log_2 n$ bits: $\log_2 |B_n| \approx n \log_2(n /
(2 \pi e))$. Two methods find it.

- **The zeta function.** For even $n \ge 2$,
  $$
  B_n = (-1)^{n/2+1} \frac{2\, n!}{(2\pi)^n} \zeta(n),
  $$
  and $\zeta(n) = 1 + 2^{-n} + 3^{-n} + \cdots$ converges so fast that
  about $n^{1 + o(1)}$ bits of precision need only $O(n / \log n)$ terms
  of the series, or an Euler product. Multiplying the approximation by the
  denominator gives a ball containing one integer, the numerator.
  `_bernoulli_fmpq_ui_zeta` does this with `arb_bernoulli_ui_zeta`; the
  cost is $O(n^2 \log^{O(1)} n)$ bit operations.
- **Multimodular.** Harvey's algorithm computes $B_n \bmod p$ in time
  $O(p \log^{O(1)} p)$, independent of $n$. With primes up to about $n
  \log n$ this gives the numerator modulo a number larger than itself in
  total time $O(n^2 \log^{O(1)} n)$ [@Har2010], with a better constant
  and trivial parallelism. `_bernoulli_fmpq_ui_multi_mod` computes only
  the low bits this way, a fraction `alpha` of the total, and the high
  bits from the zeta function, which is the cheapest combination in
  practice.

`bernoulli_fmpq_ui` first looks in the cache, then uses the zeta method
for $n < 18000$ and the combined method above, where most of the time
is spent on the residues, which run in parallel with
`flint_set_num_threads`.

## Residues modulo primes {#mod-p}

`bernoulli_mod_p_harvey` returns $B_n \bmod p$ for a prime $p <
2^{53}$. It reduces $n$ modulo $p - 1$ with Kummer's congruence $B_n/n
\equiv B_m/m \pmod p$ for $n \equiv m \not\equiv 0 \pmod{p-1}$, then
evaluates a sum of $O(p)$ terms whose inner loop is a table lookup on the
binary expansion of $1/p$. When $p - 1$ divides $n$, $p$ divides the
denominator of $B_n$ and the function returns `UWORD_MAX`. A result of 0
means that $p$ divides the numerator, as for $p = 691$ and $B_{12} =
-691/2730$; primes with this property for some $n < p - 1$ are the
*irregular primes*.

## Example {#example}

```c example name=bernoulli_isolated
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/bernoulli.h>

int main(void)
{
    fmpq_t b;
    fmpz_t p, r;

    fmpq_init(b);
    fmpz_init(p);
    fmpz_init(r);

    bernoulli_fmpq_ui(b, 20);
    flint_printf("B_20 = ");
    fmpq_print(b);
    flint_printf("\n");

    bernoulli_fmpq_ui(b, 10000);
    flint_printf("B_10000: numerator of %wu bits, denominator ",
        fmpz_bits(fmpq_numref(b)));
    fmpz_print(fmpq_denref(b));
    flint_printf("\n");

    /* compare with Harvey's algorithm modulo p = 1000003 */
    fmpz_set_ui(p, 1000003);
    fmpq_mod_fmpz(r, b, p);
    flint_printf("B_10000 mod p: %wu (Harvey), ",
        bernoulli_mod_p_harvey(10000, 1000003));
    fmpz_print(r);
    flint_printf(" (exact)\n");

    flint_printf("B_12 mod 691 = %wu\n", bernoulli_mod_p_harvey(12, 691));
    flint_printf("B_12 mod 13 is undefined: %d\n",
        bernoulli_mod_p_harvey(12, 13) == UWORD_MAX);

    fmpq_clear(b);
    fmpz_clear(p);
    fmpz_clear(r);
    return 0;
}
```
```output
B_20 = -174611/330
B_10000: numerator of 91986 bits, denominator 2338224387510
B_10000 mod p: 932435 (Harvey), 932435 (exact)
B_12 mod 691 = 0
B_12 mod 13 is undefined: 1
```
