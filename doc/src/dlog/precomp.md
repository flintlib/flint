---
title: Precomputed discrete logarithms
functions: dlog_precomp_n_init dlog_precomp_p_init dlog_precomp_pe_init dlog_precomp_modpe_init dlog_precomp_small_init dlog_precomp dlog_precomp_clear
---

A `dlog_precomp_t` holds the data for computing many discrete logarithms
to a fixed base $a$ modulo a word-size $m$. It is initialised with the
base, the modulus, what is known about the order of $a$, and the
expected number `num` of logarithms; then each call to `dlog_precomp`
returns one logarithm, and `dlog_precomp_clear` frees the data. The
choice of the initialisation function tells the module what it may
assume:

| Function | Assumption |
|---|---|
| `dlog_precomp_n_init(pre, a, mod, n, num)` | $a$ has order $n$ modulo `mod` |
| `dlog_precomp_p_init(pre, a, mod, p, num)` | $a$ has prime order $p$ |
| `dlog_precomp_pe_init(pre, a, mod, p, e, pe, num)` | $a$ has order $\mathit{pe} = p^e$ |
| `dlog_precomp_modpe_init(pre, a, p, e, pe, num)` | $a$ is a primitive root modulo $\mathit{pe} = p^e$ |
| `dlog_precomp_small_init(pre, a, mod, n, num)` | the order $n$ is small: tabulate everything |

The order must be exact and the modulus a word; neither is checked.
`dlog_precomp_n_init` factors $n$ when it is large, so passing a known
structure saves that work and lets the function pick the best method.

::: summary algorithms/number-theory/discrete-log-tables :::

## Cost and memory {#cost}

The precomputation is sized for `num` queries. In a group of prime order
$p \ge 50$ it stores a baby-step giant-step table of about $\sqrt{p \cdot
\mathit{num}}$ entries (16 bytes each), and each query costs about
$\sqrt{p/\mathit{num}}$ multiplications and binary searches. For
$\mathit{num} \ge p/2$ the whole group is tabulated and a query is one
search. For composite orders the work is split over the prime factors by
the Pohlig–Hellman decomposition, so a logarithm in $(\Z/p\Z)^\times$
with smooth $p - 1$ is cheap, and for prime powers $p^e$ the logarithm
reduces to one modulo $p$ plus $O(e)$ modular powers. The field
`pre->cost` estimates the number of multiplications of one query.

For a single logarithm use `dlog_once`, and for the logarithms of all
integers below a bound `dlog_vec`, which sieves with the logarithms of
the primes.

## Example {#example}

The program computes 1000 logarithms modulo the prime $10^6 + 3$, whose
group order $p - 1 = 2 \cdot 3 \cdot 166667$ is decomposed (type
`DLOG_CRT`), one logarithm modulo $7^5$ with the $p$-adic method, and the
logarithms to the base 2 of the integers below 12 modulo 13 with
`dlog_vec`.

```c example name=dlog_precomp_basic
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/dlog.h>

int main(void)
{
    dlog_precomp_t pre;
    nmod_t mod, order;
    ulong p = 1000003, g, b, x, k;
    ulong v[12];
    int ok = 1;

    g = n_primitive_root_prime(p);
    nmod_init(&mod, p);

    /* tables for about 1000 logarithms in a group of order p - 1 */
    dlog_precomp_n_init(pre, g, p, p - 1, 1000);
    for (k = 1; k <= 1000; k++)
    {
        b = nmod_pow_ui(k, 3, mod);
        x = dlog_precomp(pre, b);
        ok = ok && (nmod_pow_ui(g, x, mod) == b);
    }
    flint_printf("g = %wu, CRT: %d, cost %wu, 1000 logarithms ok: %d\n",
        g, pre->type == DLOG_CRT, pre->cost, ok);
    dlog_precomp_clear(pre);

    /* modulus 7^5 with the primitive root 3 */
    dlog_precomp_modpe_init(pre, 3, 7, 5, 16807, 10);
    x = dlog_precomp(pre, 2);
    nmod_init(&mod, 16807);
    flint_printf("log_3(2) mod 7^5 = %wu, 3^%wu = %wu\n",
        x, x, nmod_pow_ui(3, x, mod));
    dlog_precomp_clear(pre);

    /* logarithms of 0, ..., 11 modulo 13 to the base 2, of order 12 */
    nmod_init(&mod, 13);
    nmod_init(&order, 12);
    dlog_vec(v, 12, 2, 1, mod, 12, order);
    flint_printf("log_2 mod 13:");
    for (k = 0; k < 12; k++)
    {
        if (v[k] == DLOG_NONE)
            flint_printf(" -");
        else
            flint_printf(" %wu", v[k]);
    }
    flint_printf("\n");
    return 0;
}
```
```output
g = 2, CRT: 1, cost 13, 1000 logarithms ok: 1
log_3(2) mod 7^5 = 9308, 3^9308 = 2
log_2 mod 13: - 0 1 4 2 9 5 11 3 8 10 7
```
