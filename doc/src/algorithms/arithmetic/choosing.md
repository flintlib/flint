---
title: Choosing a method for modular arithmetic
kind: topic
summary: FLINT has a dozen ways to compute modulo a word-size $n$, for single values, vectors, SIMD lanes and matrices. This page compares their requirements and costs and says which to use.
order: 0
---

A hardware division of a two-word number by a word has a latency of tens
of cycles on current x86-64 processors, against three or four for a
multiplication. When the same modulus is used many times, FLINT replaces
the division by multiplications with precomputed data. The methods differ
in what they precompute, how they represent residues, which moduli and
inputs they accept, and whether they work on one value at a time or on
many. Moduli of more than one word are treated in
[](algorithms/arithmetic/multiword-overview).

Throughout, $B = 2^{\bits}$. A *full* product returns both words of a
$1 \times 1$ word product and a *high* product only the upper word; a
*low* product, the lower word only, is cheaper than either on most
processors.

## One value at a time {#scalar}

| Method | Entry points | Precomputed data | Requirements | Cost of $ab \bmod n$ |
|---|---|---|---|---|
| [Möller–Granlund](algorithms/arithmetic/moller-granlund) | `nmod_mul`, `NMOD_RED2`, `n_mulmod2_preinv` | one word per modulus, stored in `nmod_t` | none | 2 full + 1 low products |
| [Shoup](algorithms/arithmetic/shoup) | `n_mulmod_shoup` | one word per fixed *multiplier* | $n < 2^{\bits - 1}$, one operand fixed | 1 high + 2 low products |
| [Montgomery](algorithms/arithmetic/montgomery) | `nmod_redc_mul`, `nmod_redc_fast_mul` | $-n^{-1} \bmod B$, in `nmod_redc_ctx_t` | odd $n$; values in Montgomery form | 1 full + 1 low + 1 high product |
| Montgomery, half word | `nmod_redc_half_mul` | $-n^{-1} \bmod 2^{\bits/2}$ | odd $n < 2^{\bits/2 - 1}$ | single-word products only |
| [Barrett](algorithms/arithmetic/barrett#barrett) | `n_mod_barrett`, `NMOD_RED` | $\lfloor B/n \rfloor$, or derived from `nmod_t` ([](algorithms/arithmetic/barrett#nmod)) | reduces one word | 1 high + 1 low product per reduction |
| [Lemire–Kaser–Kurz](algorithms/arithmetic/barrett#lemire) | `n_mod_lemire` | $\lceil B/n \rceil$ | $n, x < 2^{\bits/2}$ | 1 low + 1 high product per reduction |
| Floating point | `n_mulmod_precomp`, all `_precomp` functions | `double` $1/n$ | most functions need $n < 2^{53}$ | double-precision operations and conversions |

The `nmod_t` type, initialised by `nmod_init`, bundles $n$ with its
Möller–Granlund inverse and normalisation shift. It is the representation
used by `nmod`, `nmod_vec`, `nmod_poly`, `nmod_mat` and the default `gr`
context for $\Z/n\Z$. The reduction macros `NMOD_RED`, `NMOD_RED2` and
`NMOD_RED3` reduce one, two or three words with it, and `nmod_addmul`
saves a reduction by forming $s + ab$ as a two-word integer before
reducing.

Montgomery arithmetic replaces $x$ by $xB \bmod n$, which makes a
reduction cheaper but requires conversions on entry and exit. It pays off
for long chains of multiplications on few values, such as exponentiation
and primality tests. Its *fast* (lazy) variants keep residues in $[0, 2n)$
and save the final comparison, for $n < 2^{\bits - 2}$; the *half*
variants work with $R = 2^{\bits/2}$ and need only single-word products.

## Many values with the same modulus {#vectors}

- **Fixed multiplier.** `_nmod_vec_scalar_mul_nmod` precomputes data for
  the scalar and uses [Shoup](algorithms/arithmetic/shoup) multiplication
  when $n < 2^{\bits - 1}$, Montgomery multiplication for odd full-word
  $n$, and the generic method otherwise, depending on a tuned length
  threshold (see [](nmod_vec/scalar)).
- **Sums of products.** A dot product of length $\ell$ is less than
  $\ell (n-1)^2$, which fits in at most three words. `_nmod_vec_dot`
  accumulates the exact sum and reduces once at the end, after choosing
  the number of words and the kernel with `_nmod_vec_dot_params`. See
  [](algorithms/linear-algebra/dot-products). The `_nmod_vec_nored`
  functions apply the same idea to row updates in Gaussian elimination.
- **Many inverses.** `_nmod_vec_invert` replaces $k$ modular inversions by
  one inversion and $3(k-1)$ multiplications
  ([](algorithms/arithmetic/batch-inversion)).
- **Exact arithmetic on SIMD lanes.** The `vec4d_mulmod` family of
  `machine_vectors` computes products modulo $n < 2^{50}$ on 4 or 8
  integers stored in doubles at once, using a fused multiply-add to
  recover the low part of each product
  ([](algorithms/simd/fp-modular-arithmetic)). This is the arithmetic of
  the small-prime FFT ([](algorithms/fft/small-prime-fft)), and through it of large polynomial and integer
  multiplication. `vec4n_addmod` adds integer lanes modulo $n$.
- **Matrices.** `nmod_mat_mul_blas` converts entries to doubles, multiplies
  with `flint_dgemm` while the exact result stays below $2^{53}$, and
  reduces only at the end, using several moduli and Chinese remaindering
  for large $n$ ([](algorithms/linear-algebra/gemm)). The other
  multiplication methods and how `nmod_mat_mul` chooses between them are
  compared in [](algorithms/linear-algebra/matrix-multiplication).

## Which one to use {#advice}

- **Scalar code:** use an `nmod_t` and the `nmod` functions. They accept
  any modulus, and `nmod_mul`, `nmod_addmul` and the `NMOD_RED` macros
  cover most needs. In code that only has $n$ and its inverse, the
  `_preinv` functions of `ulong_extras` do the same.
- **One operand fixed:** precompute `n_mulmod_precomp_shoup(a, n)` and use
  `n_mulmod_shoup`, if `NMOD_CAN_USE_SHOUP(mod)` holds. For whole vectors,
  call `_nmod_vec_scalar_mul_nmod` or `_nmod_vec_scalar_addmul_nmod`, which
  make this choice.
- **Long chains of multiplications on few values** (exponentiation,
  primality tests, Pollard rho): convert to Montgomery form with an
  `nmod_redc_ctx_t` for odd $n$, and use the `fast` functions when
  `nmod_redc_can_use_fast` allows. `nmod_pow_ui` does this itself for odd
  $n$ and exponents of at least $2^{11}$.
- **Sums of many products:** use `_nmod_vec_dot` or `NMOD_VEC_DOT` with
  parameters from `_nmod_vec_dot_params`, rather than reducing every term.
- **Bulk arithmetic on many independent values with $n < 2^{50}$:** the
  floating-point vector arithmetic of `machine_vectors` is the fastest on
  processors with wide FMA units; it is low level and requires care with
  input ranges.
- **Reducing single words that are not products** (hashing, sieving):
  Barrett reduction is cheaper than a two-word reduction, and Lemire's
  method is cheaper still for half-word $n$ and $x$.
- The `_precomp` (floating-point, scalar) functions of `ulong_extras`
  predate the others; the `_preinv` functions are generally faster.
