---
title: Radix conversion
summary: Conversion between binary and decimal (or any radix $b^e$) in quasi-linear time by divide and conquer with multiplications only, carrying out the arithmetic in the target radix.
complexity: $O(M(n) \log n)$ for $n$-limb numbers, where $M(n)$ is the cost of multiplication.
order: 40
---

Converting an integer from binary to decimal and back is needed for
input and output, and for arithmetic in a decimal (or $p$-adic) radix.
The schoolbook methods (repeated division by the radix, or Horner's rule)
take quadratic time, which is far too slow for numbers of millions of
digits. The [`radix`](radix) module uses a quasi-linear variant that
avoids divisions entirely by doing the arithmetic in the target radix.

Throughout, $\beta = 2^{\bits}$ is the binary limb radix and $B = b^e < \beta$
the target limb radix, for example $B = 10^{19}$ for decimal arithmetic on
a 64-bit machine. Converting an integer of $n$ limbs changes the number of
limbs by the constant factor $\log \beta / \log B$.

## Quadratic methods {#basecase}

From radix $B$ to binary, Horner's rule computes
$x = (\cdots(a_{n-1} B + a_{n-2}) B + \cdots) B + a_0$ with $n - 1$
multiplications of a growing binary number by the word $B$
(`radix_get_mpn_basecase`). From binary to radix $B$, repeated division
by $B$ produces the limbs $a_0, a_1, \ldots$ as remainders
(`radix_set_mpn_basecase`, which uses GMP's assembly `mpn_divrem_1`
while the number is long and a precomputed inverse of $B$ for the last
dozen limbs). Both cost $O(n^2)$ and are used below 40 to 50 limbs.

## Divide and conquer {#divide-and-conquer}

The classical subquadratic algorithm [@Knu1997, §4.4; @BZ2011, §1.7]
converts from binary to radix $B$ by splitting
$x = q B^k + r$ with $0 \le r < B^k$, $k \approx n/2$, by one division by
the precomputed power $B^k$, and converting $q$ and $r$ recursively; the
conversion from radix $B$ to binary splits the digits instead and
computes $x = x_{\mathrm{hi}} B^k + x_{\mathrm{lo}}$ with one
multiplication. With $T(n) = 2T(n/2) + O(M(n))$, both cost
$O(M(n) \log n)$, but the division-based direction is several times
slower in practice because a division costs a few multiplications.

The key observation of the method used by FLINT is that splitting is
free in the *source* radix and combining needs only multiplication and
addition, so it is enough to always split in the source radix and to
combine in the *target* radix.

::: algorithm Conversion by multiplication in the target radix {#alg-convert}
**Input:** the digits $a_0, \ldots, a_{m-1}$ of $x = \sum a_i c^i$ in a
source radix $c$; a target radix $B$ with arithmetic on numbers written in
radix $B$; the powers $c^{k}$, written in radix $B$, for the exponents $k$
used below.

**Output:** $x$ written in radix $B$.

1. If $m$ is small, evaluate by Horner's rule in radix-$B$ arithmetic.
2. Choose $k \approx m/2$. Recursively convert $x_{\mathrm{lo}} =
   \sum_{i<k} a_i c^i$ and $x_{\mathrm{hi}} = \sum_{i \ge k} a_i c^{i-k}$.
3. Return $x_{\mathrm{lo}} + x_{\mathrm{hi}} \cdot c^k$, computed in radix
   $B$.
:::

::: theorem {#thm-convert}
[](#alg-convert) is correct. If multiplication of $n$-limb numbers in
radix $B$ costs $M_B(n)$ with $M_B(n)/n$ nondecreasing, the conversion
costs $O(M_B(n) \log n)$ for an $n$-limb result, plus the cost of the
powers.
:::

::: proof
Correctness is the identity $x = x_{\mathrm{lo}} + x_{\mathrm{hi}}
c^k$. The two halves have about half the size of $x$ in either radix,
so the cost satisfies $T(n) = 2T(n/2) + O(M_B(n))$, whence
$T(n) = O(M_B(n) \log n)$.
:::

To convert from radix $B$ to binary, take $c = B$ and do the arithmetic
in binary (`radix_get_mpn`). To convert from binary to radix $B$, take
$c = 2^s$ with $s$ the largest integer such that $2^s < B$
(so $s = 63$ for $B = 10^{19}$): reading the binary input in chunks of
$s$ bits gives its digits in radix $c$ at no cost, and the arithmetic is
done in radix $B$ with `radix_mul`, which uses the small-prime FFT and
propagates the carries in radix $B$ (`radix_set_mpn`). In either
direction no division by a large number is needed.

## The table of powers {#powers}

The recursion uses the powers $c^{k_0}, c^{k_1}, \ldots$ with
$k_0 = \lceil m/2 \rceil$ and $k_{i+1} = \lceil k_i / 2 \rceil$ down to
$k = 1$. At every node the split is made at the largest $k_i$ that is at
most half the length of the current piece, so that the recursion stays
balanced and only these $O(\log m)$ powers are needed.

::: lemma {#lem-powers}
$k_i = 2 k_{i+1}$ or $k_i = 2 k_{i+1} - 1$. Hence $c^{k_i}$ is either the
square of $c^{k_{i+1}}$ or that square divided exactly by $c$.
:::

::: proof
$k_{i+1} = \lceil k_i/2 \rceil$ is $k_i/2$ for even $k_i$ and
$(k_i + 1)/2$ for odd $k_i$.
:::

The table therefore costs one squaring and possibly one exact division by
a single word per level, $O(M_B(n))$ in total. In binary, a power of
$B = 10^{19} = 2^{19} \cdot 5^{19}$ is divisible by $2^{19k}$, so the
binary powers $B^k$ end in many zero limbs; these are recorded and
skipped in the multiplications.
Since the table depends only on the size of the input, `radix_set_mpn`
keeps the most recent one in a per-thread cache, which speeds up the
conversion of many numbers of similar size, as when printing a vector.

## Decimal strings {#strings}

Once a number is in radix $10^{19}$, writing its decimal digits takes
linear time: each limb is split into groups of digits with a few
divisions by constants, and pairs of digits are copied from a 200-byte
table (`radix_get_str_decimal`). Accordingly:

- `flint_mpn_get_str` and `fmpz_get_str` in base 10 convert with
  `radix_set_mpn` and `radix_get_str_decimal` when `fft_small` is
  available, except for numbers of 40 to 1400 limbs, where GMP's
  `mpz_get_str` (a division-based subquadratic algorithm with fast
  basecases) is currently faster; other bases use GMP.
- `fmpz_set_str` in base 10 parses short strings with GMP's
  `mpn_set_str` and longer ones (from 24000 digits) by divide and conquer
  over the string, combining the halves as $x_{\mathrm{hi}} \cdot 10^k +
  x_{\mathrm{lo}} = (x_{\mathrm{hi}} \cdot 5^k) 2^k + x_{\mathrm{lo}}$
  with precomputed powers of 5, in parallel when threads are available.

Above 500 limbs, the two halves of the recursion of `radix_get_mpn` and
`radix_set_mpn` are converted in parallel when threads are available,
leaving the remaining threads to the multithreaded multiplication near
the root.

## Use in FLINT {#use}

`radix_get_mpn` and `radix_set_mpn` (with the divide and conquer
versions `radix_get_mpn_divconquer` and `radix_set_mpn_divconquer`),
`radix_get_str_decimal`, `flint_mpn_get_str`, and through them
`fmpz_get_str` and the conversions between `fmpz_t` and
`radix_integer_t` (`radix_integer_set_fmpz`, `radix_integer_get_fmpz`);
the [`decimal`](decimal) module builds on the `radix` arithmetic. See also
[](radix/conversion).

::: functions radix_get_mpn radix_set_mpn radix_get_str_decimal flint_mpn_get_str :::
