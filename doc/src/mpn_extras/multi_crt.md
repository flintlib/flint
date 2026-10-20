---
title: Multimodular reduction and Chinese remaindering
functions: flint_mpn_crt_init flint_mpn_multi_mod flint_mpn_multi_crt flint_mpn_multi_mod_vec flint_mpn_multi_crt_once
---

Multimodular algorithms map a large integer to its residues modulo many
word-size primes, do their work modulo each prime, and recover the integer
result by Chinese remaindering. A `flint_mpn_crt_t` holds the precomputed data for a fixed list of
moduli $m_0, \ldots, m_{k-1}$ (pairwise coprime for Chinese remaindering,
not necessarily prime), with product $P$; `flint_mpn_multi_mod` and
`flint_mpn_multi_crt` then convert in quasi-linear time
$O(M(N) \log N)$, where $N$ is the bit size of $P$ and $M(N)$ the cost of
multiplying $N$-bit integers.

These functions are the engine behind `fmpz_multi_mod_ui`,
`fmpz_multi_CRT_ui` and `fmpz_comb_t`, and behind the multimodular
algorithms of `fmpz_mat`, `fmpz_poly` and `arb_mat`.

::: summary algorithms/integers/multimodular-reduction :::

## Usage {#usage}

- Initialise once with `flint_mpn_crt_init`, or with `flint_mpn_crt_init2`
  and the flags `FLINT_MPN_CRT_MOD` and `FLINT_MPN_CRT_CRT` when only
  one direction is needed (the data for the other is then neither
  computed nor stored). The initialisation fails with an exception if the
  moduli are not pairwise coprime and Chinese remaindering was requested.
- A conversion needs `C->tmp_limbs` limbs of workspace; pass it as `tmp`
  when converting many values, or pass `NULL`.
- The reconstruction writes `C->prod_len` limbs; with `sign` nonzero it
  returns the symmetric residue, as a magnitude and a sign.
- For a single conversion, `flint_mpn_multi_mod_once` and
  `flint_mpn_multi_crt_once` are faster than initialising a structure,
  and the latter uses much less memory.
- For many small values (such as the entries of a matrix), the vector
  versions `flint_mpn_multi_mod_vec` and `flint_mpn_multi_crt_vec`
  traverse the precomputed tables once per block of values.

The reduction accepts inputs of any size, including inputs larger than
$P$, so the same structure serves for reducing the operands of a product
and for reconstructing the (smaller than $P$) result.

## Performance notes {#performance}

Moduli of at most 32 bits are batched in pairs (and 20-bit moduli in
triples) so that the tree works with full limbs, and when $P$ is small (a
single chunk of up to a few limbs), the reconstruction uses fixed-length
code from [`crt_helpers.h`](crt_helpers). The reconstruction also checks
first whether the result is much smaller than $P$ by reconstructing it
modulo the product $Q$ of the first few moduli: if the value obtained has
at least a limb of slack below $Q$ and has the right residues modulo all
the other moduli, it is the answer, and the tree is not traversed. This
makes the reconstruction of small entries of a matrix whose bound is
dominated by a few large entries almost free.

## Example {#example}

Reduce a 1000-bit number modulo the first 20 primes above $2^{60}$ and
reconstruct it, and also the negative number $-12345$:

```c example name=mpn_multi_crt_basic
#include <gmp.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/mpn_extras.h>

#define K 20

int main(void)
{
    flint_rand_t state;
    flint_mpn_crt_t C;
    ulong primes[K], res[K], x[16], y[K];
    slong i, xn = 1000 / FLINT_BITS;
    int neg;

    primes[0] = n_nextprime(UWORD(1) << 60, 1);
    for (i = 1; i < K; i++)
        primes[i] = n_nextprime(primes[i - 1], 1);

    flint_mpn_crt_init(C, primes, K);
    flint_printf("product of the moduli: %wd limbs\n", C->prod_len);

    flint_rand_init(state);
    flint_mpn_urandomb(x, state, xn * FLINT_BITS);

    flint_mpn_multi_mod(res, x, xn, C, NULL);
    flint_printf("x mod p_0 correct: %d\n",
        res[0] == mpn_mod_1(x, xn, primes[0]));

    neg = flint_mpn_multi_crt(y, res, C, 0, NULL);
    flint_printf("reconstructed: %d (sign %d)\n",
        mpn_cmp(y, x, xn) == 0 && flint_mpn_zero_p(y + xn, C->prod_len - xn), neg);

    /* the residues of -12345, reconstructed symmetrically */
    for (i = 0; i < K; i++)
        res[i] = primes[i] - 12345;
    neg = flint_mpn_multi_crt(y, res, C, 1, NULL);
    flint_printf("symmetric: %s%wu\n", neg ? "-" : "", y[0]);

    flint_mpn_crt_clear(C);
    flint_rand_clear(state);
    return 0;
}
```
```output
product of the moduli: 19 limbs
x mod p_0 correct: 1
reconstructed: 1 (sign 0)
symmetric: -12345
```
