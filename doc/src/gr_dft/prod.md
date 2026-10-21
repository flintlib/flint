---
title: Multidimensional transforms
functions: gr_dft_prod_pre_t gr_dft_prod_precomp_init gr_dft_prod_precomp_clear gr_dft_prod_precomp gr_dft_prod_precomp_set_serial_block
---

The DFT on a finite abelian group $G = \Z/c_0\Z \times \cdots \times
\Z/c_{r-1}\Z$ is the $r$-dimensional DFT of an array with $c_0 \cdots
c_{r-1}$ entries. These functions compute it with one cyclic plan per
distinct component length. They are the basis of the transforms on
groups of Dirichlet characters (`gr_dft_dirichlet`), whose components are
given by the Conrey decomposition and often repeat.

## Layout and algorithm {#layout}

The array is stored in row-major order: the entry with index
$(j_0, \ldots, j_{r-1})$ is at position
$j_{r-1} + c_{r-1}(j_{r-2} + c_{r-2}(\cdots + c_1 j_0))$, so the last index
varies fastest. The transform is computed axis by axis: for each axis
$a$, every line of $c_a$ entries along that axis is transformed with the
cyclic plan of length $c_a$, in place for the last axis and through a
buffer for the others, which are strided. Since

$$
\prod_{a} w_a^{j_a k_a}
$$

factors over the axes, there are no twiddle factors between the axes, in
contrast to the Cooley–Tukey decomposition of a cyclic transform of the
same total length. The inverse transform applies the inverse cyclic
transforms without their divisions, and divides by $n = \prod c_a$ once at
the end. Axes of length 1 are skipped.

The order of the axes does not matter mathematically, and the product
transform is not the cyclic transform of length $n$, even when the $c_a$
are pairwise coprime: the Good–Thomas algorithm reduces the cyclic
transform to a product transform only after permuting the input and the
output ([](algorithms/fft/dft-algorithms#pfa)).

## Roots {#roots}

`gr_dft_prod_precomp_init` and the one-shot functions use the default
roots $e^{-2\pi i/c_a}$, so they need a ring such as `acb`. Over other
rings, `gr_dft_prod_precomp_init_root` takes one principal root $w$ of
some order $N$ divisible by every $c_a$, for example
$N = \mathrm{lcm}(c_0, \ldots, c_{r-1})$, and uses $w^{N/c_a}$ for the
axis of length $c_a$. Each of these is then a principal $c_a$-th root of
unity, and the component plans check them as described in
[](gr_dft/plans#roots).

## Threads {#threads}

The lines of one axis are independent, and when the ring is thread-safe
and `flint_get_num_threads` is greater than 1, they are distributed over
at most $n / \mathit{serial\_block}$ chunks. Threads are requested from
the global pool for each axis separately, and only for an axis with at
least two lines; an axis with a single line, as in a product with one
long axis, is transformed by its component plan, which can then use the
threads of the pool itself. The granularity is set with
`gr_dft_prod_precomp_set_serial_block` (default
`GR_DFT_SERIAL_BLOCK_DEFAULT`).

## Example {#example}

The transform on $\Z/2\Z \times \Z/3\Z$ over $\Z/7\Z$, with the root 3 of
order 6: the axes use $3^3 = -1$ and $3^2 = 2$. The input is the
$2 \times 3$ array with rows $(1, 2, 3)$ and $(4, 5, 6)$.

```c example name=gr_dft_prod_nmod
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>
#include <flint/gr_dft.h>

/* DFT on Z/2 x Z/3 over Z/7Z, using w = 3 of order 6 */
int main(void)
{
    gr_ctx_t ctx;
    gr_dft_prod_pre_t P;
    gr_ptr x, y, w;
    ulong cyc[2] = { 2, 3 };
    slong j, n = 6;

    gr_ctx_init_nmod(ctx, 7);
    x = gr_heap_init_vec(n, ctx);
    y = gr_heap_init_vec(n, ctx);
    w = gr_heap_init(ctx);

    GR_MUST_SUCCEED(gr_set_ui(w, 3, ctx));
    GR_MUST_SUCCEED(gr_dft_prod_precomp_init_root(P, w, 6, cyc, 2, 0, ctx));

    /* entry 3 j0 + j1 holds x[j0][j1] */
    for (j = 0; j < n; j++)
        GR_MUST_SUCCEED(gr_set_si(GR_ENTRY(x, j, ctx->sizeof_elem), j + 1, ctx));

    GR_MUST_SUCCEED(gr_dft_prod_precomp(y, x, P, ctx));
    _gr_vec_print(y, n, ctx);
    flint_printf("\n");

    GR_MUST_SUCCEED(gr_dft_prod_inverse_precomp(y, y, P, ctx));
    _gr_vec_print(y, n, ctx);
    flint_printf("\n");

    gr_dft_prod_precomp_clear(P);
    gr_heap_clear(w, ctx);
    gr_heap_clear_vec(x, n, ctx);
    gr_heap_clear_vec(y, n, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```
```output
[0, 6, 2, 5, 0, 0]
[1, 2, 3, 4, 5, 6]
```

For instance, the entry $(1, 0)$ is $(1 + 2 + 3) - (4 + 5 + 6) = -9 \equiv 5$.
