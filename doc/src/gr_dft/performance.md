---
title: Algorithms, complex mode and threads
summary: How a plan chooses its algorithm, what each algorithm costs, which constants control the implementation, and how transforms use threads.
order: 2
---

A plan made with `GR_DFT_ALG_AUTO` chooses an algorithm from the
factorisation of the length $n$, and composite algorithms recursively
choose the algorithms of their sub-plans. The mathematics is in
[](algorithms/fft/dft-algorithms).

## The automatic choice {#auto}

| Length | Algorithm | Sub-plans |
|---|---|---|
| $n = 1$ | radix-2 (nothing to do) | |
| $n = 2^k$, complex mode, no `GR_DFT_SCRAMBLED` | split-radix | none |
| $n = 2^k$ otherwise | radix-2 Cooley–Tukey | none |
| $n = p^e$, $p$ odd prime | mixed-radix, radices $p, \ldots, p$ | for $p > 19$, Bluestein of length $p$ if the ring has default roots |
| $n$ with $r \ge 2$ distinct prime factors | prime factor algorithm, $n = n_1 n_2$ with $n_1$ the power of the smallest prime | automatic, for $n_1$ and $n_2$ |

Bluestein's algorithm embeds the transform of length $p$ in a cyclic convolution of
power-of-two length $N$, the least power of two with $N \ge 2p - 1$,
computed with a scrambled radix-2 sub-plan. Over rings without default
roots, such as $\Z/m\Z$, this sub-plan cannot be built, and primes above
19 use the direct kernel, at a cost of $(p-1)(p-2)$ multiplications per
transform of length $p$.

The four-step algorithm `GR_DFT_ALG_BAILEY` is never chosen
automatically: in measurements up to $n = 2^{20}$, over balls, fixed-point
numbers and $\Z/m\Z$, its blocking did not beat the radix-2 and
split-radix transforms, serially or with threads. Algorithms can be
requested explicitly, subject to the length restrictions of
`GR_DFT_ALG_AUTO`; a request for four-step or split-radix with $n < 4$
gives radix-2, and mixed-radix with $n = 1$ gives the naive algorithm.

The function `gr_dft_precomp_print` shows the result. For three lengths
near 1000 over complex balls in complex mode:

```c example name=gr_dft_plan_print
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_dft.h>

int main(void)
{
    gr_ctx_t R, C;
    gr_dft_pre_t P;
    ulong lengths[3] = { 1024, 1000, 1013 };
    int i;

    gr_ctx_init_real_arb(R, 128);
    gr_ctx_init_complex_acb(C, 128);

    for (i = 0; i < 3; i++)
    {
        GR_MUST_SUCCEED(gr_dft_precomp_init_karatsuba(P, lengths[i],
                        GR_DFT_ALG_AUTO, 0, R, C));
        gr_dft_precomp_print(P);
        gr_dft_precomp_clear(P);
    }

    gr_ctx_clear(C);
    gr_ctx_clear(R);
    return 0;
}
```
```output
plan: n = 1024, alg = split, depth = 10
  tables: roots packed stages, root err 1.07e+04 ulp (est.)
plan: n = 1000, alg = pfa
  tables: roots no, root err 1.05e+04 ulp (est.)
  n1 = 8, n2 = 125, crt coefficients 625, 376
  sub-plan 1: n = 8, alg = split, depth = 3
    tables: roots packed stages, root err 66.9 ulp (est.)
  sub-plan 2: n = 125, alg = mixed
    tables: roots yes, complex mode, root err 1.29e+03 ulp
    radices: 5 5 5
plan: n = 1013, alg = mixed
  tables: roots yes, complex mode, root err 1.06e+04 ulp
  radices: 1013
  child: n = 1013, alg = bluestein
    tables: roots yes, complex mode, root err 1.06e+04 ulp, kernel (conv_len 2048)
    child: n = 2048, alg = ct, depth = 11, scrambled
      tables: roots packed stages, root err 2.15e+04 ulp (est.)
```

The "root err" figures are the fixed-point error estimates of
[](algorithms/fft/fixed-point-dft) and are irrelevant over `acb`. The
prime factor plan stores no root table, since it has no twiddle
factors; its sub-plans build their own.

## Costs {#costs}

The counts below are multiplications by roots of unity, the dominant
cost over most rings. The number of additions is of the same order.

| Algorithm | Multiplications by roots |
|---|---|
| naive | $(n - 1)^2$ |
| radix-2 | $\frac{n}{2} \log_2 n$ (fewer nontrivial ones) |
| split-radix | about $\frac{n}{3} \log_2 n$ nontrivial; in complex mode $n \log_2 n - 3n + 4$ real multiplications with the three-multiplication product |
| four-step, $n = n_1 n_2$ | those of the sub-transforms plus $n$ twiddles |
| mixed-radix, $n = p^e$ | $e \, n \, \frac{(p-1)^2}{p}$ (twiddles and direct kernels) |
| prime factor, $n = n_1 n_2$ | $n_2 C(n_1) + n_1 C(n_2)$, no twiddles |
| Bluestein, length $p$ | $2p$ chirp multiplications, $N$ pointwise products and two transforms of length $N < 4p$ |

Here $C(m)$ is the cost of a transform of length $m$. Plan construction
computes a table of $n$ roots (by `gr_exp_pi_i` for each entry with the
default root, by a running product over `acb`, and by successive
multiplications with a supplied root), and for Bluestein's algorithm one
transform of length $N$ for the kernel.

## Complex mode and arithmetic {#complex}

Over complex rings stored as two real parts, a plan made by
`gr_dft_precomp_init_karatsuba` with the real ring is cheaper: the powers
$\pm 1, \pm i$ cost nothing and the eighth roots of unity two real
multiplications. Generic roots take three real multiplications instead of
four over exact rings and at precisions of at least 1024 bits; below
that, one native complex multiplication is faster than three real
multiplications and three additions in `arb`. The one-shot functions
`gr_dft` and `gr_dft_inverse` do not use complex mode.

Radix-2 plans, and split-radix plans in complex mode, store the roots
used by each pass contiguously (`n - 1` or `n - 2` elements in total)
instead of the table $w^0, \ldots, w^{n-1}$, so that every pass reads its
roots sequentially; the inverse transform reads the same tables
backwards. Plans with three-multiplication tables keep the plain table.

## Tuning constants {#constants}

| Constant | Value | Effect |
|---|---|---|
| `GR_DFT_PRIME_CUTOFF` (internal) | 19 | largest prime radix that always uses the direct kernel |
| `GR_DFT_KARATSUBA_CUTOFF_PREC` (internal) | 1024 | precision from which complex mode uses three real multiplications |
| `GR_DFT_NFIXED_KARATSUBA_CUTOFF` (internal) | 22 | limb count from which fixed-point complex products use three real multiplications |
| `GR_DFT_CT_BLOCK_BYTES` | 128 KiB | block size for the depth-first order of the radix-2 transform |
| `GR_DFT_SERIAL_BLOCK_DEFAULT` | 256 | default threading granularity |

The radix-2 transform processes blocks larger than `GR_DFT_CT_BLOCK_BYTES`
one level at a time and then recurses into their halves, so that each
block that fits in the cache completes all its remaining levels while it
is resident. The butterflies and their order of dependence are those of
the level-by-level transform, so the results are the same.

## Threads {#threads}

A transform uses threads when the ring is thread-safe
(`gr_ctx_is_threadsafe` returns `T_TRUE`, for the real ring too in
complex mode) and either threads are attached to the plan with
`gr_dft_precomp_set_threads`, or `flint_get_num_threads` is greater than
1, in which case threads are requested from the global pool for the
duration of the transform. Work is divided into at most
$n / \mathit{serial\_block}$ chunks, so lengths below
$2 \cdot \mathit{serial\_block}$ run serially. The default granularity of
256 elements is small; for cheap rings such as $\Z/m\Z$ a larger value,
set with `gr_dft_precomp_set_serial_block`, can be faster. Each threaded
algorithm uses the same operations as the serial one, so the results do
not depend on the number of threads.

- *Radix-2.* The first $t$ levels of the transform, where $2^t$ is at
  least the number of chunks, are done pass by pass with the butterflies
  of each pass split evenly between the threads. After them, the array
  consists of $2^t$ independent blocks, and each thread completes all
  remaining levels of its own blocks without further synchronisation.
- *Split-radix.* The two transforms of length $n/4$ run on one thread
  while the transform of length $n/2$ runs on another, recursively
  dividing the available threads, and the $n/4$ independent steps of
  each combining pass are split between all threads.
- *Four-step.* The column transforms, the twiddle multiplications and
  the row transforms are each split between the threads; the final
  transposition is serial.
- *Product transforms.* The lines of each axis are split between the
  threads; an axis with a single line leaves the threads to its
  component plan.

The mixed-radix, prime factor, Bluestein and naive algorithms are not
threaded themselves, but their sub-plans of radix-2 and split-radix type
are, when they are long enough. Threads attached to a plan are borrowed:
`gr_dft_precomp_clear` does not return them to the pool.

The profile programs `p-gr_dft_threads` (speedup per algorithm and thread
count), `p-gr_dft_nfixed` (fixed-point arithmetic kernels), `p-gr_dft_acb`
and `p-gr_dft_acb_accuracy` in `src/gr_dft/profile` measure these
choices.
