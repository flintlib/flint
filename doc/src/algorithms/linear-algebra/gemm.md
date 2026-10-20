---
title: Floating-point matrix multiplication (gemm)
summary: How FLINT's own sgemm and dgemm kernels block, pack and thread a matrix product around a register-tile micro-kernel, and when floating-point matrix products of integers are exact.
complexity: $2mkn$ floating-point operations; packing copies $O(mk \lceil n/n_c \rceil + kn)$ entries.
order: 30
---

One of the ways FLINT multiplies matrices over $\Z/n\Z$ and $\Z$ is to
convert the entries to floating point, call a `gemm` (general matrix
multiplication) routine and convert back: a well-tuned floating-point
matrix product runs close to the peak arithmetic rate of the processor,
and it is exact when the entries are small enough ([](#exactness)).
Larger entries are first reduced modulo several small numbers, each
product is computed by `gemm`, and the result is reconstructed by
Chinese remaindering. This route is used only in some ranges of
moduli, entry sizes and dimensions; the topic pages [](algorithms/linear-algebra/matrix-multiplication) and
[](algorithms/linear-algebra/integer-matrix-multiplication) give the
complete dispatch. In outline:

- **Over $\Z/n\Z$**, `nmod_mat_mul` uses floating-point `gemm` in
  `nmod_mat_mul_blas` and in the byte kernels of `nmod_mat_mul_u8`,
  which call `flint_sgemm` on blocks for $17 \le n \le 255$.
  `nmod_mat_mul_blas` is chosen for $n < 2^{32}$ only for large matrices
  for which one `dgemm` is exact (from a tuned dimension,
  `FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF`), and otherwise mainly for
  moduli above $2^{52}$ (and on targets without the SIMD kernels) from a
  dimension between 100 and several hundred, depending on the modulus
  and the number of threads; with an external BLAS, also for large
  matrices with moduli of 33 to 52 bits. In
  the other cases `nmod_mat_mul` uses classical multiplication with
  delayed reduction ([](algorithms/linear-algebra/dot-products)), the
  blocked SIMD kernels `nmod_mat_mul_u32`, `nmod_mat_mul_u52`,
  `nmod_mat_mul_k52` and `nmod_mat_mul_fp50`
  ([](algorithms/linear-algebra/modular-matrix-kernels)), table lookups
  and bit-packing for tiny moduli
  ([](algorithms/linear-algebra/small-modulus-matrix-multiplication)),
  and Strassen–Winograd above them
  ([](algorithms/linear-algebra/strassen-winograd)).
- **Over $\Z$**, `fmpz_mat_mul` uses `gemm` (`fmpz_mat_mul_blas`) on
  64-bit machines only for dimensions above 50 and above a limit that
  depends on the bit size $c$ of the output entries: none for
  $c \le 53$, where one `dgemm` is exact, 1000 for $c \le 64$, 1300 for
  $c \le 128$, 250 for $c \le 256$ and about $200 + 8 \log_2 c$ above.
  Otherwise it uses dot products in one to five machine words, Waksman's
  and Strassen–Winograd's algorithms, multimodular multiplication through
  `nmod_mat_mul`, and, for entries of many thousands of bits,
  `fmpz_mat_mul_fft_small`, which transforms every entry once.
- **Other rings.** `fmpz_mod_mat_mul`, `fmpq_mat_mul` and the block
  products of `arb_mat` reduce to `fmpz_mat_mul`; `mpn_mod` matrices
  (moduli of several words) use classical multiplication with delayed
  reduction, Waksman's algorithm or multimodular multiplication through
  `nmod_mat_mul` ([](mpn_mod/mat-mul)); and `gr_mat_mul` calls the method
  of each ring, by default the classical algorithm
  ([](gr_mat/mul)). They reach `gemm` only through `nmod_mat_mul` and
  `fmpz_mat_mul`.

When FLINT is built with BLAS, `flint_sgemm` and `flint_dgemm` call it
(unless `flint_gemm_use_blas` is set to 0). Otherwise they use FLINT's
own kernels in `src/machine_vectors/gemm.c`, which are described
below. Their organisation follows Goto and van de Geijn
[@GvdG2008] in the loop ordering of BLIS [@VZvdG2015].

## Exactness {#exactness}

The products FLINT computes are exact: the entries of $A$ and $B$ are
integers, and they are bounded so that no rounding ever occurs.

::: theorem {#thm-exact}
Let $A$ and $B$ have integer entries and suppose that for all $i, j$

$$
\sum_{l} |A_{il}| \, |B_{lj}| \le 2^{53}.
$$

Then any algorithm that computes each $C_{ij}$ by additions,
multiplications and fused multiply-adds of doubles, starting from the
entries of $A$, $B$ and zeros, and combining products $A_{il} B_{lj}$ for
distinct $l$, returns the exact integer $C_{ij} = \sum_l A_{il} B_{lj}$.
The same holds for floats with $2^{24}$ in place of $2^{53}$.
:::

::: proof
Every intermediate value is a sum of a subset of the products
$A_{il} B_{lj}$, which is an integer of absolute value at most
$\sum_l |A_{il}||B_{lj}| \le 2^{53}$. Such integers are doubles, so by
induction every operation has an exactly representable result and
returns it without rounding.
:::

The theorem does not depend on the order of summation, so it applies to
FLINT's kernels including the zero padding in the packed buffers below,
and to any implementation that evaluates the sums directly, which is what
FLINT assumes of an external BLAS. A sufficient condition is
$k \max |A_{il}| \max |B_{lj}| \le 2^{53}$. For entries reduced modulo
$p$ and lifted to the symmetric range $|x| \le \lfloor p/2 \rfloor$,
this becomes $k \lfloor p/2 \rfloor^2 \le 2^{53}$, which is four times
weaker than $k (p - 1)^2 \le 2^{53}$ for entries in $[0, p)$.

## Structure {#structure}

The kernels compute $C = AB$ for row-major $A$ ($m \times k$), $B$
($k \times n$) and $C$ ($m \times n$) with arbitrary leading dimensions.
One source file, `gemm_templ.h`, is instantiated for `float` and
`double`, and all arithmetic goes through the vector types of
`machine_vectors.h`. The code has no target-specific branches beyond the
choice of vector type and tile height.

### The micro-kernel {#micro}

The innermost routine updates an $m_r \times n_r$ tile of $C$, held in
$2 m_r$ vector registers, with $n_r = 2v$ where $v$ is the number of
lanes per vector:

::: algorithm Micro-kernel {#alg-micro}
**Input:** a packed $m_r \times k_c$ block $\tilde A$ of $A$, a packed
$k_c \times n_r$ block $\tilde B$ of $B$, and an $m_r \times n_r$ tile of
$C$ (or zero).

**Output:** the tile plus $\tilde A \tilde B$.

1. Load the tile into accumulators $c_{r,0}, c_{r,1}$ for
   $0 \le r < m_r$, or set them to zero for the first block of $k$.
2. For $l = 0, \ldots, k_c - 1$:
   1. load the two vectors $b_0, b_1$ of row $l$ of $\tilde B$;
   2. for $r = 0, \ldots, m_r - 1$: broadcast $\tilde A_{rl}$ to a
      vector $a$, then set $c_{r,0} \gets a b_0 + c_{r,0}$ and
      $c_{r,1} \gets a b_1 + c_{r,1}$ with fused multiply-adds.
3. Store the accumulators.
:::

Each step of the loop performs $2 m_r$ independent fused multiply-adds
on two loaded vectors and $m_r$ broadcasts, enough to hide the latency
of the multiply-add units. The tile sizes chosen in `gemm.c` are:

| Backend | Double: vector, $m_r \times n_r$ | Float: vector, $m_r \times n_r$ | Registers used |
|---|---|---|---|
| AVX-512 | `vec8dz`, $8 \times 16$ | `vec16fz`, $8 \times 32$ | 18 of 32 |
| AVX2 | `vec4d`, $6 \times 8$ | `vec8f`, $6 \times 16$ | 14 of 16 |
| NEON | `vec2d`, $8 \times 4$ | `vec4f`, $8 \times 8$ | 18 of 32 |
| 64-bit, no AVX2 or NEON | `vec2d`, $6 \times 4$ | `vec4f`, $6 \times 8$ | 14 of 16 |
| Forced generic, GNU vectors with AVX | `vec4d`, $4 \times 8$ | `vec8f`, $4 \times 16$ | 10 |
| Forced generic, otherwise | `vec2d`, $4 \times 4$ | `vec4f`, $4 \times 8$ | 10 |
| Strict ISO C | `vec4d`, $4 \times 8$ | `vec4f`, $4 \times 8$ | |

The forced generic rows also apply to 32-bit builds, which always use a
generic backend.

### Blocking and packing {#blocking}

Around the micro-kernel, five loops partition the matrices so that the
data each level reuses stays in a level of the memory hierarchy:

::: algorithm Blocked product {#alg-blocked}
**Input:** $A$, $B$ and the block sizes $n_c = 2048$, $k_c = 256$,
$m_c = 96$ (in elements).

**Output:** $C = AB$.

1. For each block of $n_c$ columns of $B$ and $C$:
2. For each block of $k_c$ rows of this column block of $B$:
   1. pack the $k_c \times n_c$ block of $B$ into column panels of width
      $n_r$;
   2. for each block of $m_c$ rows of $A$ and $C$:
      1. pack the $m_c \times k_c$ block of $A$ into row panels of
         height $m_r$;
      2. for each panel of $B$ (step $n_r$), and for each panel of $A$
         (step $m_r$): run [](#alg-micro) on the corresponding tile of
         $C$, accumulating unless this is the first block of $k$.
:::

Packing copies a panel into a contiguous buffer in the order in which
the micro-kernel reads it: for $B$, row $l$ of the panel is the $n_r$
entries the kernel loads at step $l$, and for $A$, the $m_r$ entries of
column $l$ are adjacent. A panel that extends past the matrix is padded
with zeros, so the micro-kernel always works on full tiles. An edge tile
of $C$ is copied into a local $m_r \times n_r$ buffer, updated and copied
back, so no memory outside $C$ is touched. The packed buffers are
aligned to 64 bytes, and the scratch space is a thread-local buffer that
grows as needed and is kept between calls, or allocated per call without
thread-local storage.

In double precision the packed $A$ block takes
$96 \cdot 256 \cdot 8$ bytes $= 192$ KiB and a $k_c \times n_r$ panel
of $B$ takes 16 KiB with AVX2. In the scheme of [@GvdG2008] the former
stays in the level 2 cache and the latter in the level 1 cache while it
is multiplied by every panel of $A$. The block sizes are the same on all
targets.

### Special shapes {#shapes}

Packing costs $O(mk + kn)$ copies per call. For small or thin matrices
this is not repaid, and the serial driver dispatches as follows when
$k > 0$:

- $m \le m_r$: $B$ is read once, in strips of $4v$ columns, while the
  $m$ rows of a strip of $C$ stay in registers. Remaining columns are
  computed by scalar loops.
- all of $m, k, n \le 96$: the whole problem fits in the level 2 cache,
  and a direct kernel loads $A$ and $B$ in place, sweeping strips of $v$
  columns with $m_r$ accumulators; the last $n \bmod v$ columns of $B$
  are copied once into a zero-padded buffer. This path is disabled for
  the strict ISO C backend, where it was measured to be slower.
- $n \le n_r$: $B$ is packed once into a single $k \times n_r$ panel and
  $A$ is streamed through it.
- otherwise: [](#alg-blocked).

When $m \le 0$ or $n \le 0$ nothing is done, and when $k \le 0$ the
result is set to zero.

## Threading {#threading}

`flint_sgemm_fallback` and `flint_dgemm_fallback` use FLINT's thread
pool. Thin products ($m \le m_r$ or $n \le n_r$) are memory-bound and
always run serially. Otherwise the number of threads is at most
`flint_get_num_threads()` and at most $1 + 2mkn / 10^6$, so that each
thread gets about a million floating-point operations, and is further
limited by the number of free threads in the pool. Each worker is woken
once per call. There are two drivers:

- **Cooperative** (x86-64): all threads run the loops of
  [](#alg-blocked) together. For each block of $B$, the threads pack its
  panels in parallel, taking panel indices from an atomic counter, wait at
  a barrier, then take blocks of rows of $A$ from a second counter, each
  packing its own $A$ block and computing the corresponding rows of $C$.
  The row block size is $m/(4T)$ for $T$ threads, rounded up to a
  multiple of $m_r$ and capped at $m_c$. The barriers are
  sense-reversing spin barriers; except on Windows, a waiting thread
  yields the processor after every 2048 spins.
- **Split** (other targets): the rows of $C$ are divided into $T$ nearly
  equal blocks and each thread runs the serial algorithm on its rows of
  $A$ and all of $B$. Nothing is shared except read-only input, but each
  thread packs $B$ itself. The source notes that this driver is slightly
  slower at high thread counts and that the cooperative one, which relies
  on C11 atomics with acquire and release ordering, has only been stress
  tested on x86-64.

Defining `GEMM_MT_FORCE_ATOMIC` or `GEMM_MT_FORCE_SPLIT` when compiling
`gemm.c` overrides the choice.

## Use in FLINT {#use}

All conversions between FLINT's integer types and floating point happen
outside the kernels, which compute plain floating-point products.

- `nmod_mat_mul_blas` lifts entries modulo $p$ to the symmetric range.
  If $k \lfloor p/2 \rfloor^2 < 2^{24}$ it uses a single `flint_sgemm`,
  and if $k \lfloor p/2 \rfloor^2 < 2^{53}$ a single `flint_dgemm`.
  Otherwise it chooses up to 12 pairwise coprime moduli $t$ below
  $\sqrt{2^{55}/k}$, so that each product modulo $t$ is exact by
  [](#thm-exact), computes one `flint_dgemm` per modulus and recombines
  the results by Chinese remaindering
  ([](algorithms/linear-algebra/modular-matrix-kernels#blas)). Before each
  call it gives its worker threads back to the pool so that the kernels
  can use them.
- `fmpz_mat_mul_blas` uses a single `flint_dgemm` when the entries of $A$
  and $B$ have at most $a$ and $b$ bits with $a + b + \lceil \log_2(k+1)
  \rceil \le 53$. Otherwise it reduces the entries modulo primes $p$ with
  $k ((p - 1)/2)^2 < 2^{53}$, multiplies modulo each prime with
  `flint_dgemm` and reconstructs the result with a multi-modular
  algorithm.
- `nmod_mat_mul_u8` handles moduli $p \le 255$ with entries stored in
  bytes. For moduli above 16, the range of its byte-shuffle kernels, it
  converts entries in $[0, p)$ to floats and calls
  `flint_sgemm` on blocks of $k_{\max} = \lfloor (2^{24} - 1)/(p - 1)^2
  \rfloor$ columns, accumulating the exact block results in 32-bit
  integers ([](algorithms/linear-algebra/small-modulus-matrix-multiplication#sgemm)).

Exact linear algebra by floating-point matrix multiplication, with the
modular reduction delayed to the end, is the approach of the FFLAS and
FFPACK libraries [@DGP2008].
