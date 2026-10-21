---
title: Chinese remaindering
functions: fmpz_CRT_ui fmpz_comb_init fmpz_comb_temp_init fmpz_multi_mod_ui fmpz_multi_CRT_ui fmpz_multi_mod_ui_once fmpz_multi_CRT fmpz_multi_CRT_precompute fmpz_multi_CRT_precomp fmpz_multi_mod_precompute fmpz_multi_mod_precomp
---

Multimodular algorithms compute an integer result modulo many primes and
reconstruct it by the Chinese remainder theorem: if $m_1, \ldots, m_k$
are pairwise coprime with product $M$, an integer $x$ with $|x| < M/2$ is
determined by its residues $x \bmod m_i$. The module offers three
interfaces, which differ in what is known in advance:

| Interface | Moduli | Use |
|---|---|---|
| `fmpz_CRT`, `fmpz_CRT_ui` | two at a time | adding one modulus at a time, when the number of primes needed is not known beforehand |
| `fmpz_comb_t` with `fmpz_multi_mod_ui`, `fmpz_multi_CRT_ui` | a fixed list of word-size moduli | many reductions or reconstructions with the same primes, as in multimodular polynomial and matrix arithmetic |
| `fmpz_multi_CRT_t`, `fmpz_multi_mod_t` | a fixed list of moduli of any size | moduli that are themselves large, such as prime powers |

In all of them, the flag `sign` selects the range of the result: $[0, M)$
if it is 0, and the symmetric range $(-M/2, M/2]$ otherwise, which is the
right choice when the result may be negative.

## Two moduli at a time {#incremental}

`fmpz_CRT(out, r1, m1, r2, m2, sign)` returns the $x$ with $x \equiv r_1
\pmod{m_1}$ and $x \equiv r_2 \pmod{m_2}$ by Garner's formula
$$
x = r_1 + m_1 \bigl((r_2 - r_1) m_1^{-1} \bmod m_2\bigr),
$$
which needs one inverse modulo $m_2$ and a few products. `fmpz_CRT_ui`
does the same for a word-size $m_2$, with word arithmetic for the
inverse and the reduction, so that adding a prime to a product $m_1$ of
$n$ limbs costs $O(n)$. A loop over $k$ primes thus costs $O(k^2)$ word
operations; `_fmpz_CRT_ui_precomp` saves the inverse when the same pair of
moduli recurs. The residue $r_1$ may lie in $[-m_1, m_1)$, so the signed
output of one step can be the input of the next. The moduli must be
coprime; `fmpz_CRT` checks this (it raises an exception) and
`fmpz_CRT_ui` does not.

## A fixed list of word-size moduli {#comb}

`fmpz_comb_init` precomputes a product tree of the moduli together with
the data for both directions; `fmpz_comb_init2` restricts it to one
direction. A `fmpz_comb_temp_t` holds scratch space and must not be
shared between threads. `fmpz_multi_mod_ui` then reduces an integer
modulo all primes at once, and `fmpz_multi_CRT_ui` reconstructs it. Both
work in $O(M(n) \log n)$ word operations for $n$ limbs of input or
output, compared with $O(n^2)$ for one prime at a time, using remainder
trees, a linear basecase on chunks of the tree, and batching of small
primes into one word ([](algorithms/integers/binary-splitting#remainder-trees),
[](algorithms/integers/multimodular-reduction)).
The work is done by `flint_mpn_multi_mod` and `flint_mpn_multi_crt` of
[](mpn_extras), of which `fmpz_comb_t` is a thin wrapper. The functions
`fmpz_multi_mod_ui_once` and `fmpz_multi_CRT_ui_once` build and discard
the precomputation for a single call, and `_fmpz_vec_multi_CRT_ui`
reconstructs a whole vector of integers.

The moduli need not be prime, but they must be pairwise coprime and
nonzero; this is not checked. The residues passed to `fmpz_multi_CRT_ui`
must be reduced.

## Moduli of any size {#multi}

`fmpz_multi_CRT_precompute` takes moduli of any size, checks that they
are pairwise coprime, and compiles a straight-line program that combines
them in a binary tree; as in Huffman coding, the two smallest current
moduli are always combined first, so that unbalanced lists of moduli
lead to cheap trees. `fmpz_multi_CRT_precomp` runs the program on a list
of residues; the residues need not be reduced. `fmpz_multi_CRT` does both
steps for a single reconstruction. `fmpz_multi_mod_precompute` and
`fmpz_multi_mod_precomp` are the analogous pair for reducing an integer
modulo moduli of any size.

## Example {#example}

The program reduces $-100!$, a number of 525 bits, modulo ten 60-bit
primes and reconstructs it with a comb and incrementally, and solves a
system of three congruences with moduli 9, $10^{20}$ and $7^{30}$.

```c example name=fmpz_crt
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>

#define NUM 10

int main(void)
{
    ulong primes[NUM], res[NUM];
    fmpz_comb_t comb;
    fmpz_comb_temp_t temp;
    fmpz_t x, y, m, r;
    fmpz moduli[3], values[3];
    slong i;

    fmpz_init(x); fmpz_init(y); fmpz_init(m); fmpz_init(r);

    /* ten 60-bit primes, product of 600 bits */
    primes[0] = n_nextprime(UWORD(1) << 59, 1);
    for (i = 1; i < NUM; i++)
        primes[i] = n_nextprime(primes[i - 1], 1);

    /* x = -100!, reduced modulo each prime and reconstructed */
    fmpz_fac_ui(x, 100);
    fmpz_neg(x, x);
    fmpz_comb_init(comb, primes, NUM);
    fmpz_comb_temp_init(temp, comb);
    fmpz_multi_mod_ui(res, x, comb, temp);
    fmpz_multi_CRT_ui(y, res, comb, temp, 1);       /* signed range */
    flint_printf("comb: %s\n", fmpz_equal(x, y) ? "recovered -100!" : "wrong");
    fmpz_multi_CRT_ui(y, res, comb, temp, 0);       /* range [0, M) */
    flint_printf("unsigned result has %wu bits\n", fmpz_bits(y));

    /* the same, one prime at a time */
    fmpz_set_ui(y, res[0]);
    fmpz_set_ui(m, primes[0]);
    for (i = 1; i < NUM; i++)
    {
        fmpz_CRT_ui(y, y, m, res[i], primes[i], 1);
        fmpz_mul_ui(m, m, primes[i]);
    }
    flint_printf("incremental: %s\n", fmpz_equal(x, y) ? "recovered -100!" : "wrong");

    /* arbitrary moduli: x = 2 mod 9, x = 3 mod 10^20, x = -1 mod 7^30 */
    for (i = 0; i < 3; i++)
    {
        fmpz_init(moduli + i);
        fmpz_init(values + i);
    }
    fmpz_set_ui(moduli + 0, 9);
    fmpz_ui_pow_ui(moduli + 1, 10, 20);
    fmpz_ui_pow_ui(moduli + 2, 7, 30);
    fmpz_set_ui(values + 0, 2);
    fmpz_set_ui(values + 1, 3);
    fmpz_set_si(values + 2, -1);
    if (fmpz_multi_CRT(r, moduli, values, 3, 0))
    {
        flint_printf("multi_CRT: "); fmpz_print(r); flint_printf("\n");
        fmpz_mod(y, r, moduli + 1);
        fmpz_add_ui(x, r, 1);
        flint_printf("r mod 9 = %wu, r mod 10^20 = %wd, 7^30 | r + 1: %d\n",
            fmpz_fdiv_ui(r, 9), fmpz_get_si(y), fmpz_divisible(x, moduli + 2));
    }

    for (i = 0; i < 3; i++)
    {
        fmpz_clear(moduli + i);
        fmpz_clear(values + i);
    }
    fmpz_comb_temp_clear(temp);
    fmpz_comb_clear(comb);
    fmpz_clear(x); fmpz_clear(y); fmpz_clear(m); fmpz_clear(r);
    return 0;
}
```
```output
comb: recovered -100!
unsigned result has 591 bits
incremental: recovered -100!
multi_CRT: 3093485019800846923426157000000000000000000003
r mod 9 = 2, r mod 10^20 = 3, 7^30 | r + 1: 1
```

In the unsigned range, $-100!$ is represented by $M - 100!$, which has
591 bits like the product $M$ of the primes.
