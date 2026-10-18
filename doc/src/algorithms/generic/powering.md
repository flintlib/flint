---
title: Generic powering
summary: Binary exponentiation and the sliding window method as implemented for any ring, with operation counts, and when each is preferable.
complexity: $\lfloor \log_2 e \rfloor$ squarings plus $\nu(e) - 1$ multiplications by $x$ (binary); about $\log_2 e$ squarings plus $2^{k-1} + \log_2(e)/(k+1)$ multiplications (sliding window of width $k$).
order: 50
---

The default integer powers of the generic interface, `gr_pow_ui`,
`gr_pow_si` and `gr_pow_fmpz`, are computed by binary exponentiation in
any ring, using only the ring's `mul` and `sqr` methods
(`gr_generic_pow_ui`, `gr_generic_pow_fmpz`). A sliding window variant
is also available for exponents and rings in which it is faster. Both
are classical; see Knuth [@Knu1997, §4.6.3] and Cohen [@Coh1996, §1.2]. Throughout, $x$ is an
element of a ring, $e \ge 1$ an integer with binary expansion
$e = \sum_{i=0}^{\ell} e_i 2^i$, $e_\ell = 1$, and $\nu(e)$ the number of
one bits.

## Binary exponentiation {#binary}

::: algorithm Left-to-right binary exponentiation {#alg-binary}
1. Set $r \gets x$.
2. For $i = \ell - 1, \ell - 2, \ldots, 0$: set $r \gets r^2$, and if
   $e_i = 1$ set $r \gets r x$.
3. Return $r$.
:::

::: proposition {#prop-binary}
[](#alg-binary) returns $x^e$ using $\ell = \lfloor \log_2 e \rfloor$
squarings and $\nu(e) - 1$ multiplications by $x$.
:::

::: proof
After the step for index $i$, $r = x^{\lfloor e / 2^i \rfloor}$, since
$\lfloor e / 2^i \rfloor = 2 \lfloor e / 2^{i+1} \rfloor + e_i$; at
$i = 0$ this is $x^e$. Each index costs one squaring and, for each of
the $\nu(e) - 1$ one bits below the leading one, a multiplication.
:::

`gr_generic_pow_ui_binexp` and `gr_generic_pow_fmpz_binexp` implement
this, with direct formulas for $e \le 4$ (for example $x^4 = (x^2)^2$
and $x^3 = x^2 \cdot x$). To avoid copying, the `ulong` version
alternates between the output and one temporary, choosing the initial
buffer from the parity of the number of zero bits so that the final
result lands in the output; the output may alias $x$, in which case $x$
is copied first. All multiplications are by the base $x$, which matters
when $x$ is small compared with the intermediate powers: multiplying a
large element by a small one costs much less than a general product.

For negative exponents, `gr_generic_pow_si` and `gr_generic_pow_fmpz`
invert $x$ with `gr_inv` and raise the inverse to $-e$, so the power
fails with the inverse (`GR_DOMAIN` for a non-unit). `gr_pow_fmpq` with
exponent $p/2$ computes a square root (or reciprocal square root) and
then an integer power.

## Sliding windows {#sliding}

The window method processes several bits of $e$ at a time using
precomputed odd powers of $x$.

::: algorithm Sliding window exponentiation {#alg-sliding}
1. Choose a window width $k \ge 1$. Compute $x^2$ and the odd powers
   $x, x^3, x^5, \ldots, x^{2^k - 1}$, each from the previous one by a
   multiplication by $x^2$.
2. Set $r \gets 1$ and $i \gets \ell$. While $i \ge 0$:
    - if $e_i = 0$, set $r \gets r^2$ and $i \gets i - 1$;
    - otherwise let $j \ge \max(i - k + 1, 0)$ be the smallest index
      with $e_j = 1$, let $w = (e_i \ldots e_j)_2$ (an odd number below
      $2^k$), set $r \gets r^{2^{i - j + 1}} x^w$ and $i \gets j - 1$.
3. Return $r$.
:::

::: proposition {#prop-sliding}
[](#alg-sliding) returns $x^e$. It uses $\ell + 2$ squarings at most
(one in the precomputation), $2^{k-1} - 1$ multiplications in the
precomputation, and one multiplication per window; the number of
windows is at most $\nu(e)$ and, for random exponents, about
$(\ell + 1) / (k + 1)$ on average.
:::

::: proof
Each step consumes the bits $e_i, \ldots, e_j$ and maintains
$r = x^{\lfloor e / 2^{j} \rfloor}$ after a window (and
$r = x^{\lfloor e / 2^{i} \rfloor}$ after a zero bit), as in
[](#prop-binary), since appending the bits $e_i \ldots e_j$ to the
exponent multiplies it by $2^{i - j + 1}$ and adds $w$. Every bit
position causes one squaring of $r$ (the first squarings of $r = 1$ are
trivial). A window starts at a one bit and covers at least that bit, so
there are at most $\nu(e)$ windows; for random bits a window covers
$k$ bits and is followed on average by one zero bit, which gives the
average count.
:::

The multiplication count $2^{k-1} + (\ell + 1)/(k + 1)$ is minimised by
$k \approx \log_2 \ell - 2 \log_2 \log_2 \ell$;
`gr_generic_pow_ui_sliding` and `gr_generic_pow_fmpz_sliding` choose $k$
from a table, from $k = 1$ for exponents of at most 8 bits to $k = 11$
above 56651 bits. Unlike binary exponentiation, the window method
multiplies by precomputed powers $x^w$ that are as large as $r$, so it
saves work only when every multiplication by the base is a full-size
product anyway.

## Use in FLINT {#use}

- `gr_generic_pow_ui`, `gr_generic_pow_si` and `gr_generic_pow_fmpz`,
  the defaults of `gr_pow_ui`, `gr_pow_si` and `gr_pow_fmpz`, use
  binary exponentiation ([](#alg-binary)). Rings with a better method
  override them; for example `nmod` calls `nmod_pow_ui` and `arb`
  calls `arb_pow_ui`.
- `fmpz_powm` for moduli of at least 70000 bits (in builds with
  `fft_small`), and `fmpz_mod_pow_ui` and `fmpz_mod_pow_fmpz` for such
  moduli, run the generic
  algorithms over an `fmpz_mod` context: binary exponentiation when the
  base is a single-word integer, and the sliding window method when the
  base is large (and, in `fmpz_mod_pow_fmpz`, the exponent has at least
  20 bits). Smaller moduli use GMP's `mpz_powm`.
- `gr_generic_mul_2exp_si`, `gr_generic_set_fmpz_2exp_fmpz` and
  `gr_generic_set_fmpz_10exp_fmpz` compute powers of 2 and 10 in the
  ring with these functions.
