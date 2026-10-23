---
title: Radix-2 and truncated transforms
functions: fft_radix2 ifft_radix2 fft_truncate fft_truncate1 ifft_truncate ifft_truncate1 fft_truncate_sqrt2 ifft_truncate_sqrt2 fft_negacyclic ifft_negacyclic
---

These are the transforms from which the convolutions of the module are
built. They work in place on an array of coefficient pointers, in the
ring $\Z/p\Z$, $p = 2^{nw} + 1$, where multiplication by the root of
unity $2^w$ is a shift. Most applications should use `fft_convolution`
instead, which combines them correctly; the functions here are for
building other transforms, for example a convolution with a different
pointwise operation.

## What is computed {#definition}

Let $L$ be the transform length, $L = 2n$ for `fft_radix2` and the
truncated functions and $L = 4n$ for the sqrt2 functions, let
$\omega = 2^w$, respectively $\omega = \sqrt 2^{\,w}$, a root of unity of
order exactly $L$ in $\Z/p\Z$, and let $\mathrm{rev}(k)$ be the integer
whose $\log_2 L$ binary digits are those of $k$ in reverse order. For
$a(x) = \sum_{j < L} a_j x^j$, the forward transforms compute

$$
\hat a_k = a\bigl(\omega^{\mathrm{rev}(k)}\bigr), \qquad 0 \le k < L,
$$

so that $\hat a_0, \ldots, \hat a_{L/2 - 1}$ are the values at the even
powers of $\omega$ and the second half the values at the odd powers. The
inverse transforms take the values in this order and return
$L \cdot a_j$; dividing by $L$ is a shift by $\log_2 L$ bits
(`mpn_div_2expmod_2expp1` with $d = \log_2 L$). Since
$\omega^{\mathrm{rev}(k)}$ runs over all powers, pointwise products of
two transforms are the transform of the cyclic product modulo
$x^L - 1$, whatever the order.

The recursion is that of Gentleman and Sande [@GS1966] (decimation in
frequency): a layer of butterflies
$(a_i, a_{i + L/2}) \mapsto (a_i + a_{i + L/2}, \omega^i (a_i - a_{i + L/2}))$
followed by two transforms of half the length with root $\omega^2$. The
inverse transforms reverse these steps with the butterflies
$(x, y) \mapsto (x + \omega^{-i} y, x - \omega^{-i} y)$, which produce
twice the inputs of the forward butterfly; this is where the factor $L$
comes from.

## Truncation {#truncation}

`fft_truncate(ii, n, w, t1, t2, trunc)` computes only
$\hat a_0, \ldots, \hat a_{\mathit{trunc} - 1}$ of a polynomial with
$a_j = 0$ for $j \ge \mathit{trunc}$. `ifft_truncate` inverts it: from
these `trunc` values it recovers $L a_0, \ldots, L a_{\mathit{trunc} - 1}$.
Together they multiply polynomials whose product has at most `trunc`
coefficients with about $\mathit{trunc} \log_2 L$ butterflies instead of
$L \log_2 L$. `fft_truncate1` computes the same outputs from a full input,
and `ifft_truncate1` is the general inverse step of
[van der Hoeven's algorithm](algorithms/fft/truncated-fft): given
$\hat a_k$ for $k < \mathit{trunc}$ and $L a_j$ for $j \ge \mathit{trunc}$
in the same array, it computes $L a_j$ for $j < \mathit{trunc}$. In all
four functions `trunc` must be even and between 2 and $L$.

The entries at and beyond `trunc` are scratch space for `fft_truncate`
and `ifft_truncate`, and their initial values are ignored, but they must
exist. `fft_truncate_sqrt2` and `ifft_truncate_sqrt2` require
$L/2 < \mathit{trunc} \le L$; for odd $w$ the first layer uses the
twiddles $\sqrt 2^{\,iw}$ for odd $i$ (`fft_butterfly_sqrt2`), and for
even $w$ the functions simply call `fft_truncate` and `ifft_truncate`
with $2n$ and $w/2$. For odd $w$, a value $\mathit{trunc} \le L/2$
makes them access entries before the start of the array.

## Negacyclic transforms {#negacyclic}

`fft_negacyclic` evaluates $a(x)$ at the roots of $x^{2n} + 1$: it
multiplies $a_j$ by $\varphi^j$, $\varphi = \sqrt 2^{\,w}$, which has
order $4n$, and then applies the transform of length $2n$ with root
$\varphi^2 = 2^w$. Pointwise products of two such transforms are the
transform of the product modulo $x^{2n} + 1$, which
`ifft_negacyclic` recovers, multiplied by $2n$. With pieces of an integer
as coefficients this is multiplication modulo $2^N + 1$; see
`_fft_mulmod_2expp1`. These are the weighted transforms of Crandall and
Fagin [@CF1994] in a ring where the weights are shifts.

## Data movement {#memory}

The functions never copy coefficients: a butterfly writes its results to
the two scratch blocks `*t1` and `*t2` and then swaps the pointers with
the two entries of `ii`. After a transform, the entries of `ii` and the
scratch pointers therefore point into each other's original blocks. All
blocks should be part of one allocation, and code must always access the
coefficients through `ii` after a call, never through pointers saved
before it. The sqrt2 and negacyclic functions need a third scratch block
`*temp`.

The coefficients grow during a transform, since sums are not reduced:
after $k$ layers of butterflies the top limb holds values up to about
$2^k$ in absolute value, far within a limb. Reduce with
`mpn_normmod_2expp1` before using the values as reduced residues, for
example before `fft_mulmod_2expp1`.

## Example {#example}

The transform of the polynomial $x$ of length 8 modulo $p = 2^{64} + 1$
with root $2^{16}$. Position $k$ receives $2^{16 \, \mathrm{rev}(k)}$
with $\mathrm{rev}(0), \ldots, \mathrm{rev}(7) = 0, 4, 2, 6, 1, 5, 3, 7$,
so for example $2^{64} \equiv -1$ (printed as $2^{64}$, top limb 1),
$2^{96} \equiv -2^{32}$ and $2^{80} \equiv -2^{16}$. The inverse
transform returns $8x$.

```c example name=fft_radix2_order
#include <flint/flint.h>
#include <flint/fft.h>
#include <flint/mpn_extras.h>

int main(void)
{
    /* length 2n = 8 transform modulo p = 2^64 + 1 with root 2^w = 2^16 */
    mp_size_t n = 4, limbs = 1, size = limbs + 1, i;
    flint_bitcnt_t w = 16;
    mp_limb_t mem[10 * 2], * ii[8], * t1, * t2;

    for (i = 0; i < 8; i++)
        ii[i] = mem + i * size;
    t1 = mem + 8 * size;
    t2 = mem + 9 * size;

    /* the polynomial x */
    for (i = 0; i < 8; i++)
        flint_mpn_zero(ii[i], size);
    ii[1][0] = 1;

    fft_radix2(ii, n, w, &t1, &t2);

    for (i = 0; i < 8; i++)
    {
        mpn_normmod_2expp1(ii[i], limbs);
        flint_printf("position %wd: %016wx %wu\n", i, ii[i][0], ii[i][1]);
    }

    /* the inverse transform returns 2n times the input */
    ifft_radix2(ii, n, w, &t1, &t2);
    for (i = 0; i < 8; i++)
    {
        mpn_normmod_2expp1(ii[i], limbs);
        flint_printf("%wu%s", ii[i][0], i < 7 ? " " : "\n");
    }
    return 0;
}
```
```output
position 0: 0000000000000001 0
position 1: 0000000000000000 1
position 2: 0000000100000000 0
position 3: ffffffff00000001 0
position 4: 0000000000010000 0
position 5: ffffffffffff0001 0
position 6: 0001000000000000 0
position 7: ffff000000000001 0
0 8 0 0 0 0 0 0
```
