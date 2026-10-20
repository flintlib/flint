---
title: machine_vectors.h
header: machine_vectors.h
description: Short SIMD vectors of doubles, floats and words, exact modular arithmetic in floating point, and FLINT's own sgemm and dgemm.
preprocess:
  - -mavx2 -mfma
  - -march=skylake-avx512
  - -DFLINT_MACHINE_VECTORS_FORCE_GENERIC
  - -DFLINT_MACHINE_VECTORS_FORCE_GENERIC -DFLINT_MACHINE_VECTORS_STRICT_C
  - -D__ARM_NEON -D__aarch64__ -I{stubs}
category: Foundations
order: 5
requires: 64
---

This header defines short vectors of `double`, `float` and `ulong` values
with inline operations that compile to the SIMD instructions of the
target: AVX2 and AVX-512 on x86-64, NEON on 64-bit ARM, and portable C
elsewhere. On top of these it provides exact modular arithmetic on
integers stored in doubles, which is the arithmetic of the small-prime FFT
(`fft_small`), and the matrix multiplication kernels `flint_sgemm` and
`flint_dgemm`, which some of the matrix multiplication methods of
`nmod_mat` and `fmpz_mat` use
([](algorithms/linear-algebra/gemm)).

The vector operations are internal building blocks. Each backend provides
the operations its users in FLINT need, and a missing operation is not
emulated, so the available set differs between backends. Code that must
build everywhere should use only operations marked as available on all
backends; the guide [](machine_vectors/portable) lists them, with the
few cases where their results differ. [](machine_vectors/gemm) describes
the matrix multiplication functions in detail.

::: summary algorithms/simd/fp-modular-arithmetic :::

::: summary algorithms/linear-algebra/gemm :::

## Backends {#backends}

The header chooses one backend per translation unit from the target
macros of the compiler:

| Backend | Selected when | Representation |
|---|---|---|
| AVX2 | `__AVX2__` is defined and `ulong` has 64 bits | `__m128d`, `__m256d`, `__m256i`, ...; 8-lane types are pairs |
| AVX-512 | in addition `__AVX512F__` | the AVX2 types, plus `vec8dz`, `vec16fz` etc. in single zmm registers |
| NEON | `__ARM_NEON` or `_M_ARM64`, and 64-bit `ulong` | `float64x2_t`, `uint64x2_t`, `float32x4_t`; wider types are pairs |
| GNU vectors | otherwise, with GCC or clang | `__attribute__((vector_size))` types |
| Strict ISO C | otherwise, with other compilers | structs of arrays |

The last two are the *generic* backends.

- **Inline only.** Every vector function is a forced-inline static
  function, so the backend is fixed by the compiler flags of the
  translation unit that includes the header, which need not be those
  FLINT was built with. `flint_sgemm` and `flint_dgemm` are compiled into
  the library and use the backend of FLINT's own build.
- **Word size.** The AVX2 and NEON backends need 64-bit `ulong`, because
  their integer lanes are `ulong` and some intrinsics exist only in 64-bit
  mode. A 32-bit build uses a generic backend, in which `vec4n` has
  32-bit lanes.
- **Strict ISO C.** In this tier, whether an operation becomes a SIMD
  instruction is up to the compiler's SLP vectoriser, so its performance
  varies between compilers. It exists so that the header works on
  compilers without GNU vector extensions; GCC and clang use the GNU
  vector tier.

::: macro
FLINT_MACHINE_VECTORS_FORCE_GENERIC
FLINT_MACHINE_VECTORS_STRICT_C

If defined before the header is included, select a generic backend even
on a target with AVX2 or NEON. `FLINT_MACHINE_VECTORS_STRICT_C` selects
the strict ISO C tier, also with GCC or clang. Both are intended for
testing and profiling the portable code paths.
:::

::: macro
FLINT_MACHINE_VECTORS_GENERIC
FLINT_MACHINE_VECTORS_GNU_VECTOR_EXTENSIONS

Defined by the header. `FLINT_MACHINE_VECTORS_GENERIC` is defined when a
generic backend is forced or `FLINT_BITS` is not 64. A 64-bit target
without AVX2 or NEON uses a generic backend without defining it, so a
generic backend is in use if this macro is defined or if neither
`__AVX2__` nor NEON is. `FLINT_MACHINE_VECTORS_GNU_VECTOR_EXTENSIONS` is
defined when the GNU vector tier is in use.
:::

## Availability {#availability}

Each function below states on which backends it is available. In these
statements "AVX2" includes builds with AVX-512, and a statement such as
"AVX2 and NEON" means that the generic backends lack the operation. The
following table summarises them.

| Operations | AVX2 | NEON | Generic |
|---|---|---|---|
| Loads, stores, `set`, `zero`, `add`, `sub`, `mul`, `fmadd`, `fnmadd` on `vec1d`, `vec4d`, `vec8d`, `vec1f`, `vec4f`, `vec8f`, and `floor` on these types except `vec8d` | yes | yes | yes |
| The same on `vec2d` | `floor` only | yes | yes |
| The same on `vec16f` | no | no | yes |
| The same on `vec8dz` and `vec16fz` | AVX-512 | no | no |
| `vec4d`: `one`, `set_d4`, `get_index`, `neg`, `abs`, `div`, `min`, `max`, `sqrt`, `round`, `fmsub`, `_fused`, comparisons, `blendv`, masks and bits, `unpacklo`, `unpackhi`, permutations, gathers | yes | yes | yes |
| `vec4n`: loads, stores, `set_n`, `zero`, `add`, `sub`, bitwise and/or, shifts, `shift_lanes_up` | yes | yes | yes |
| Modular arithmetic, `addmod`, conversions, `vec4d_fnmsub`, the fused unpack-and-permute operations, `vec4n` unpacks and `permute2`, the transpose macros, `vec4n_store_unaligned_halves`, `vec8n`, the remaining `vec1d` and `vec8d` operations | yes | yes | no |
| `round_mul`, `reduce_2n_to_n`, `vec4d_reduce_0n_to_pmhn`, `vec4d_addsub`, `vec4d_half`, `vec4d_same`, `vec4d_same_mod`, `vec8d_same`, `vec4n_mul`, `vec4n_horizontal_sum`, `vec4n_permute_3_2_1_0`, `vec8n_bit_shift_right`, printing | yes | no | no |
| `vec4n_set_n4` | yes | no | yes |
| Other `vec2d` operations, `vec2n`, `vec1n`, `vec8d_one`, `vec8n_add`, `vec8n_sub` | no | yes | no |
| `vec8d_floor` | no | no | yes |

The generic backends implement the operations that `flint_sgemm` and
`flint_dgemm` use, the `vec4d` operations that `dfloat` uses and the
`vec4n` operations that `mp_real` uses. Everything that `fft_small` needs
requires AVX2 or NEON.

## Compiler flags and warnings {#compiling}

- **Compiler flags.** `configure` chooses `-march` from the detected
  CPU. `--enable-avx2` adds `-mavx2 -mfma`. `--enable-avx512` implies
  `--enable-avx2` and adds `-mavx512f`, and `-mavx512dq` if the compiler
  accepts it.
- **FMA.** The AVX2 versions of `fmadd`, `fmsub`, `fnmadd` and `fnmsub`
  use FMA instructions. Without `__FMA__` (for instance with `-mavx2`
  alone) the header still compiles, but code calling these functions
  does not. The `_fused` variants fall back to the C
  library `fma` in that case. With MSVC, `/arch:AVX2` is assumed to
  include FMA.
- **Contraction.** The generic backends write `fmadd`, `fnmadd` and
  `vec4d_fmsub` as `a * b + c`, which becomes a single FMA instruction
  only if the compiler contracts floating-point expressions. GCC does not
  in a strict ISO mode such as `-std=c11`, which is how FLINT is built,
  and clang does by default. Performance-critical code should request
  contraction, as `machine_vectors/gemm.c` does with
  `#pragma GCC optimize("fp-contract=fast")`. The AVX2, AVX-512 and NEON
  backends use fused instructions. The `_fused` variants round once on
  every backend.
- **ABI warnings.** On x86-64 without AVX, GCC and clang warn
  (`-Wpsabi`) when a function returns or takes a 32-byte generic vector
  such as `vec4d`, because the calling convention for such types changed.
  The header silences this for its own definitions, but the warning can
  still appear in the calling code, where it is harmless as long as the
  functions are inlined. Compile with `-Wno-psabi` or use a pragma, as in
  the examples on [](machine_vectors/portable).

## Conventions {#conventions}

- Lanes are numbered from 0. A load puts the element at the lowest
  address in lane 0, and `vec4d_set_d4(a0, a1, a2, a3)` puts `a0` in
  lane 0. In the descriptions, $a_i$ denotes lane $i$ of `a`.
- Vectors are passed and returned by value. Since some platforms pass
  vector types by value incorrectly, helper functions that take vectors
  should be declared `FLINT_FORCE_INLINE`, like those of the header.
- In the names of modular functions, `0n` stands for $[0, n)$, `2n` for
  $[0, 2n)$, `pm1n` for $[-n, n]$, `pm1no` for the open interval
  $(-n, n)$ and `pmhn` for $[-n/2, n/2]$. For example
  `vec4d_reduce_pm1no_to_0n` maps $(-n, n)$ to $[0, n)$.
- The modular functions and `round` assume that the floating-point
  rounding mode is round to nearest, the default.

## Types {#types}

::: type
vec1d
vec2d
vec4d
vec8d

Vectors of 1, 2, 4 and 8 doubles. `vec1d` is `double`. On AVX2, `vec4d`
is `__m256d` and `vec8d` is a struct of two `vec4d` named `e1` and `e2`;
`vec2d` is `__m128d` but only has `vec2d_floor`. On NEON, `vec2d` is
`float64x2_t`, `vec4d` is a struct of two `vec2d` and `vec8d` a struct of
two `vec4d`. On the generic backends `vec2d` and `vec4d` are GNU vector
types or structs with an array member `v`, and `vec8d` is a struct of two
`vec4d`. Code that accesses lanes directly is not portable; use
`vec4d_get_index` or a store.
:::

::: type
vec1f
vec4f
vec8f
vec16f

Vectors of 1, 4, 8 and 16 floats. `vec1f` is `float`. On AVX2, `vec4f`
and `vec8f` are `__m128` and `__m256`. On NEON `vec4f` is `float32x4_t`
and `vec8f` a pair. `vec16f` is a pair of `vec8f` everywhere, and has
operations only on the generic backends.
:::

::: type
vec1n
vec2n
vec4n
vec8n

Vectors of 1, 2, 4 and 8 `ulong` values. `vec1n` is `ulong`. On AVX2
`vec2n` and `vec4n` are `__m128i` and `__m256i`, and `vec8n` is a pair.
On NEON `vec2n` is `uint64x2_t` and the wider types are pairs. The
generic backends only have `vec4n`.
:::

::: type
vec8dz
vec16dz
vec16fz
vec32fz
vec8nz

AVX-512F only: 8 doubles (`__m512d`), 16 doubles (two `vec8dz`), 16
floats (`__m512`), 32 floats (two `vec16fz`) and 8 `ulong` values
(`__m512i`). The `z` suffix distinguishes them from `vec8d` and friends,
which remain pairs of AVX2 registers. Only `vec8dz` and `vec16fz` have
operations.
:::

::: type
vec4di

Four table indices, as produced by `vec4d_to_index` and consumed by
`vec4d_gather`. Available on all backends.
:::

## Loads, stores and gathers {#memory}

Plain and `_aligned` loads and stores are the same operation. On AVX2
they require the address to be aligned to the size of the register: 16
bytes for `vec4f`, 32 bytes for `vec4d`, `vec8f`, `vec8d` and `vec4n`,
and 64 bytes for `vec8dz` and `vec16fz`. The `_unaligned` forms accept
any address. NEON and the generic backends never require alignment.
Arrays for the aligned forms can be declared with `ALIGN_STRUCT` or
`_Alignas`.

::: macro
ALIGN_STRUCT(x)

Expands to the attribute that aligns a variable or type to `x` bytes:
`__declspec(align(x))` with MSVC and `__attribute__((aligned(x)))`
otherwise.
:::

::: function
vec1d vec1d_load(const double * a)
vec2d vec2d_load(const double * a)
vec4d vec4d_load(const double * a)
vec8d vec8d_load(const double * a)
vec8dz vec8dz_load(const double * a)
vec1d vec1d_load_aligned(const double * a)
vec2d vec2d_load_aligned(const double * a)
vec4d vec4d_load_aligned(const double * a)
vec8d vec8d_load_aligned(const double * a)
vec8dz vec8dz_load_aligned(const double * a)
vec1d vec1d_load_unaligned(const double * a)
vec2d vec2d_load_unaligned(const double * a)
vec4d vec4d_load_unaligned(const double * a)
vec8d vec8d_load_unaligned(const double * a)
vec8dz vec8dz_load_unaligned(const double * a)

Loads 1, 2, 4 or 8 consecutive doubles. The `vec2d` versions are
available on NEON and the generic backends, the `vec8dz` versions with
AVX-512 only, and the others on all backends.
:::

::: function
void vec1d_store(double * z, vec1d a)
void vec2d_store(double * z, vec2d a)
void vec4d_store(double * z, vec4d a)
void vec8d_store(double * z, vec8d a)
void vec8dz_store(double * z, vec8dz a)
void vec1d_store_aligned(double * z, vec1d a)
void vec2d_store_aligned(double * z, vec2d a)
void vec4d_store_aligned(double * z, vec4d a)
void vec8d_store_aligned(double * z, vec8d a)
void vec8dz_store_aligned(double * z, vec8dz a)
void vec1d_store_unaligned(double * z, vec1d a)
void vec2d_store_unaligned(double * z, vec2d a)
void vec4d_store_unaligned(double * z, vec4d a)
void vec8d_store_unaligned(double * z, vec8d a)
void vec8dz_store_unaligned(double * z, vec8dz a)

Stores the lanes of `a` at consecutive addresses starting at `z`.
Availability as for the loads.
:::

::: function
vec1f vec1f_load(const float * a)
vec4f vec4f_load(const float * a)
vec8f vec8f_load(const float * a)
vec16f vec16f_load(const float * a)
vec16fz vec16fz_load(const float * a)
vec1f vec1f_load_aligned(const float * a)
vec4f vec4f_load_aligned(const float * a)
vec8f vec8f_load_aligned(const float * a)
vec16f vec16f_load_aligned(const float * a)
vec16fz vec16fz_load_aligned(const float * a)
vec1f vec1f_load_unaligned(const float * a)
vec4f vec4f_load_unaligned(const float * a)
vec8f vec8f_load_unaligned(const float * a)
vec16f vec16f_load_unaligned(const float * a)
vec16fz vec16fz_load_unaligned(const float * a)
void vec1f_store(float * z, vec1f a)
void vec4f_store(float * z, vec4f a)
void vec8f_store(float * z, vec8f a)
void vec16f_store(float * z, vec16f a)
void vec16fz_store(float * z, vec16fz a)
void vec1f_store_aligned(float * z, vec1f a)
void vec4f_store_aligned(float * z, vec4f a)
void vec8f_store_aligned(float * z, vec8f a)
void vec16f_store_aligned(float * z, vec16f a)
void vec16fz_store_aligned(float * z, vec16fz a)
void vec1f_store_unaligned(float * z, vec1f a)
void vec4f_store_unaligned(float * z, vec4f a)
void vec8f_store_unaligned(float * z, vec8f a)
void vec16f_store_unaligned(float * z, vec16f a)
void vec16fz_store_unaligned(float * z, vec16fz a)

Loads and stores of floats. The `vec16f` versions are available on the
generic backends only, the `vec16fz` versions with AVX-512 only, and the
others on all backends.
:::

::: function
vec4n vec4n_load_unaligned(const ulong * a)
vec4n vec4n_load_aligned(const ulong * a)
void vec4n_store_unaligned(ulong * z, vec4n a)
void vec4n_store_aligned(ulong * z, vec4n a)

Loads and stores of four words. All backends.
:::

::: function
void vec4n_store_unaligned_halves(ulong * z, vec4n a)

Same as `vec4n_store_unaligned`, but done as two 16-byte stores on AVX2.
AVX2 and NEON.
:::

::: function
vec8n vec8n_load_unaligned(const ulong * a)

Loads eight words. AVX2 and NEON.
:::

::: function
vec4di vec4d_to_index(vec4d a)
vec4d vec4d_gather(const double * t, vec4di i)

`vec4d_to_index` converts lanes holding integers in $[0, 2^{31})$ to
table indices, and `vec4d_gather` returns
$\{t[i_0], t[i_1], t[i_2], t[i_3]\}$. All backends.
:::

## Constructing vectors and reading lanes {#lanes}

::: function
vec1d vec1d_set_d(double a)
vec2d vec2d_set_d(double a)
vec4d vec4d_set_d(double a)
vec8d vec8d_set_d(double a)
vec8dz vec8dz_set_d(double a)
vec1f vec1f_set_f(float a)
vec4f vec4f_set_f(float a)
vec8f vec8f_set_f(float a)
vec16f vec16f_set_f(float a)
vec16fz vec16fz_set_f(float a)

Returns a vector with every lane equal to `a`. Availability by type as
for the loads.
:::

::: function
vec1d vec1d_zero(void)
vec2d vec2d_zero(void)
vec4d vec4d_zero(void)
vec8d vec8d_zero(void)
vec8dz vec8dz_zero(void)
vec1f vec1f_zero(void)
vec4f vec4f_zero(void)
vec8f vec8f_zero(void)
vec16f vec16f_zero(void)
vec16fz vec16fz_zero(void)

Returns a vector of zeros. Availability by type as for the loads.
:::

::: function
vec1d vec1d_one(void)
vec4d vec4d_one(void)

Returns a vector of ones. `vec4d_one` is available on all backends,
`vec1d_one` on AVX2 and NEON.
:::

::: function
vec4d vec4d_set_d4(double a0, double a1, double a2, double a3)
vec8d vec8d_set_d8(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7)

Returns the vector with lanes $a_0, a_1, \ldots$ in this order.
`vec4d_set_d4` is available on all backends, `vec8d_set_d8` on AVX2 and
NEON.
:::

::: function
double vec4d_get_index(vec4d a, const int i)
double vec8d_get_index(vec8d a, int i)

Returns lane $i$ of `a`, for $0 \le i < 4$ and $0 \le i < 8$
respectively. `vec4d_get_index` is available on all backends,
`vec8d_get_index` on AVX2 and NEON.
:::

::: function
void vec4d_print(vec4d a)
void vec4n_print(vec4n a)

Prints `a` to standard output. `vec4d_print` prints the lanes in order
as `{a0, a1, a2, a3}` with `%f`. `vec4n_print` prints them in hexadecimal
from lane 3 down to lane 0, as `[hi a3_a2_a1_a0 lo]`. AVX2 only.
:::

## Floating-point arithmetic {#arith}

Unless stated otherwise the operations act lane by lane and round as the
corresponding IEEE 754 operation.

::: function
vec1d vec1d_add(vec1d a, vec1d b)
vec2d vec2d_add(vec2d a, vec2d b)
vec4d vec4d_add(vec4d a, vec4d b)
vec8d vec8d_add(vec8d a, vec8d b)
vec8dz vec8dz_add(vec8dz a, vec8dz b)
vec1d vec1d_sub(vec1d a, vec1d b)
vec2d vec2d_sub(vec2d a, vec2d b)
vec4d vec4d_sub(vec4d a, vec4d b)
vec8d vec8d_sub(vec8d a, vec8d b)
vec8dz vec8dz_sub(vec8dz a, vec8dz b)
vec1d vec1d_mul(vec1d a, vec1d b)
vec2d vec2d_mul(vec2d a, vec2d b)
vec4d vec4d_mul(vec4d a, vec4d b)
vec8d vec8d_mul(vec8d a, vec8d b)
vec8dz vec8dz_mul(vec8dz a, vec8dz b)

Returns $a + b$, $a - b$ and $ab$. Availability by type as for the
loads.
:::

::: function
vec1f vec1f_add(vec1f a, vec1f b)
vec4f vec4f_add(vec4f a, vec4f b)
vec8f vec8f_add(vec8f a, vec8f b)
vec16f vec16f_add(vec16f a, vec16f b)
vec16fz vec16fz_add(vec16fz a, vec16fz b)
vec1f vec1f_sub(vec1f a, vec1f b)
vec4f vec4f_sub(vec4f a, vec4f b)
vec8f vec8f_sub(vec8f a, vec8f b)
vec16f vec16f_sub(vec16f a, vec16f b)
vec16fz vec16fz_sub(vec16fz a, vec16fz b)
vec1f vec1f_mul(vec1f a, vec1f b)
vec4f vec4f_mul(vec4f a, vec4f b)
vec8f vec8f_mul(vec8f a, vec8f b)
vec16f vec16f_mul(vec16f a, vec16f b)
vec16fz vec16fz_mul(vec16fz a, vec16fz b)

Returns $a + b$, $a - b$ and $ab$ in single precision. Availability by
type as for the loads.
:::

::: function
vec1d vec1d_div(vec1d a, vec1d b)
vec4d vec4d_div(vec4d a, vec4d b)
vec8d vec8d_div(vec8d a, vec8d b)

Returns $a / b$. `vec4d_div` is available on all backends, the others on
AVX2 and NEON.
:::

::: function
vec1d vec1d_neg(vec1d a)
vec4d vec4d_neg(vec4d a)
vec8d vec8d_neg(vec8d a)
vec1d vec1d_abs(vec1d a)
vec4d vec4d_abs(vec4d a)

Returns $-a$ and $|a|$, which flip and clear the sign bit (also of zeros
and NaNs). The `vec4d` versions are available on all backends, the
others on AVX2 and NEON.
:::

::: function
vec1d vec1d_min(vec1d a, vec1d b)
vec4d vec4d_min(vec4d a, vec4d b)
vec8d vec8d_min(vec8d a, vec8d b)
vec1d vec1d_max(vec1d a, vec1d b)
vec4d vec4d_max(vec4d a, vec4d b)
vec8d vec8d_max(vec8d a, vec8d b)

Returns the smaller and the larger of $a$ and $b$. When a lane is
unordered (a NaN is involved), AVX2 and the generic backends return the
lane of `b`, NEON returns a NaN, and `vec1d_min` and `vec1d_max`, which
are `fmin` and `fmax`, return the operand that is not a NaN. The `vec4d`
versions are available on all backends, the others on AVX2 and NEON.
:::

::: function
vec1d vec1d_half(vec1d a)
vec4d vec4d_half(vec4d a)

Returns $a/2$, computed as a multiplication by $0.5$. `vec1d_half` is
available on AVX2 and NEON, `vec4d_half` on AVX2 only.
:::

::: function
vec1d vec1d_addsub(vec1d a, vec1d b)
vec4d vec4d_addsub(vec4d a, vec4d b)

Subtracts in even lanes and adds in odd lanes:
`vec4d_addsub` returns $\{a_0 - b_0, a_1 + b_1, a_2 - b_2, a_3 + b_3\}$
and `vec1d_addsub` returns $a - b$. `vec1d_addsub` is available on AVX2
and NEON, `vec4d_addsub` on AVX2 only.
:::

::: function
vec4d vec4d_sqrt(vec4d a)

Returns the correctly rounded square root. All backends.
:::

::: function
vec1d vec1d_round(vec1d a)
vec4d vec4d_round(vec4d a)
vec8d vec8d_round(vec8d a)

Rounds to an integer in the current rounding mode: to the nearest
integer, ties to even, by default. On NEON, `vec4d_round` and
`vec8d_round` always round to nearest with ties to even. `vec4d_round`
is available on all backends, the others on AVX2 and NEON.
:::

::: function
vec1d vec1d_floor(vec1d a)
vec2d vec2d_floor(vec2d a)
vec4d vec4d_floor(vec4d a)
vec8d vec8d_floor(vec8d a)
vec8dz vec8dz_floor(vec8dz a)
vec1f vec1f_floor(vec1f a)
vec4f vec4f_floor(vec4f a)
vec8f vec8f_floor(vec8f a)
vec16f vec16f_floor(vec16f a)
vec16fz vec16fz_floor(vec16fz a)

Returns $\lfloor a \rfloor$. `vec8d_floor` and `vec16f_floor` are
available on the generic backends only, the `z` versions with AVX-512
only, and the others on all backends (`vec2d_floor` also on AVX2).
:::

## Fused multiply-add {#fma}

::: function
vec1d vec1d_fmadd(vec1d a, vec1d b, vec1d c)
vec2d vec2d_fmadd(vec2d a, vec2d b, vec2d c)
vec4d vec4d_fmadd(vec4d a, vec4d b, vec4d c)
vec8d vec8d_fmadd(vec8d a, vec8d b, vec8d c)
vec8dz vec8dz_fmadd(vec8dz a, vec8dz b, vec8dz c)
vec1d vec1d_fnmadd(vec1d a, vec1d b, vec1d c)
vec2d vec2d_fnmadd(vec2d a, vec2d b, vec2d c)
vec4d vec4d_fnmadd(vec4d a, vec4d b, vec4d c)
vec8d vec8d_fnmadd(vec8d a, vec8d b, vec8d c)
vec8dz vec8dz_fnmadd(vec8dz a, vec8dz b, vec8dz c)

Returns $ab + c$ and $-ab + c$. On AVX2, AVX-512 and NEON, and for
`vec1d` and `vec1f` on every backend, the result is rounded once. On the generic
backends the wider types compute `a * b + c` and `c - a * b`, which round
twice unless the compiler contracts them. Availability by type as for
the loads.
:::

::: function
vec1f vec1f_fmadd(vec1f a, vec1f b, vec1f c)
vec4f vec4f_fmadd(vec4f a, vec4f b, vec4f c)
vec8f vec8f_fmadd(vec8f a, vec8f b, vec8f c)
vec16f vec16f_fmadd(vec16f a, vec16f b, vec16f c)
vec16fz vec16fz_fmadd(vec16fz a, vec16fz b, vec16fz c)
vec1f vec1f_fnmadd(vec1f a, vec1f b, vec1f c)
vec4f vec4f_fnmadd(vec4f a, vec4f b, vec4f c)
vec8f vec8f_fnmadd(vec8f a, vec8f b, vec8f c)
vec16f vec16f_fnmadd(vec16f a, vec16f b, vec16f c)
vec16fz vec16fz_fnmadd(vec16fz a, vec16fz b, vec16fz c)

Single-precision versions, with the same rounding behaviour.
Availability by type as for the loads.
:::

::: function
vec1d vec1d_fmsub(vec1d a, vec1d b, vec1d c)
vec4d vec4d_fmsub(vec4d a, vec4d b, vec4d c)
vec8d vec8d_fmsub(vec8d a, vec8d b, vec8d c)
vec1d vec1d_fnmsub(vec1d a, vec1d b, vec1d c)
vec4d vec4d_fnmsub(vec4d a, vec4d b, vec4d c)
vec8d vec8d_fnmsub(vec8d a, vec8d b, vec8d c)

Returns $ab - c$ and $-ab - c$, rounded once except for `vec4d_fmsub`
on the generic backends. NEON computes $ab - c$ as $-(c - ab)$, which
differs from `vec4d_fmsub_fused` in the sign of an exact zero.
`vec4d_fmsub` is available on all backends, the others on AVX2 and NEON.
:::

::: function
vec4d vec4d_fmadd_fused(vec4d a, vec4d b, vec4d c)
vec4d vec4d_fmsub_fused(vec4d a, vec4d b, vec4d c)
vec4d vec4d_fnmadd_fused(vec4d a, vec4d b, vec4d c)

Returns $ab + c$, $ab - c$ and $-ab + c$ rounded once on every backend,
like `fma(a, b, c)`, `fma(a, b, -c)` and `fma(-a, b, c)`, including the
sign of an exact zero. On AVX2 without FMA instructions they call `fma`
lane by lane. All backends.
:::

## Comparisons, masks and bits {#masks}

A mask has all bits set in the lanes where a condition holds and all
bits clear elsewhere. The `vec4d` operations in this section behave
identically on all backends.

::: function
vec4d vec4d_cmp_eq(vec4d a, vec4d b)
vec4d vec4d_cmp_ne(vec4d a, vec4d b)
vec4d vec4d_cmp_lt(vec4d a, vec4d b)
vec4d vec4d_cmp_le(vec4d a, vec4d b)
vec4d vec4d_cmp_gt(vec4d a, vec4d b)
vec4d vec4d_cmp_ge(vec4d a, vec4d b)
vec4d vec4d_cmp_nlt(vec4d a, vec4d b)
vec4d vec4d_cmp_nle(vec4d a, vec4d b)
vec4d vec4d_cmp_ngt(vec4d a, vec4d b)
vec4d vec4d_cmp_nge(vec4d a, vec4d b)

Returns the mask of the lanes where $a_i = b_i$, $a_i \ne b_i$,
$a_i < b_i$, and so on. The predicates `eq`, `lt`, `le`, `gt` and `ge`
are false when a lane of `a` or `b` is a NaN, and their negations `ne`,
`nlt`, `nle`, `ngt` and `nge` are true. All backends.

::: details
A lane of the mask has all 64 bits set where the predicate holds and all
bits clear elsewhere. As a double such a lane is a NaN or $+0$, so the
mask is meant to be passed to `vec4d_blendv`, `vec4d_movemask` or
bitwise operations such as `vec4d_bit_and`, not used in arithmetic.
:::
:::

::: function
vec1d vec1d_blendv(vec1d a, vec1d b, vec1d c)
vec4d vec4d_blendv(vec4d a, vec4d b, vec4d c)
vec8d vec8d_blendv(vec8d a, vec8d b, vec8d c)

Returns, in each lane, $b_i$ if the sign bit of $c_i$ is set and $a_i$
otherwise, so a mask from a comparison selects `b` where the comparison
holds. `vec1d_blendv` instead tests `c >= 0`, which differs for
$c = -0$ (selects `a`) and for NaNs (selects `b`). `vec4d_blendv` is
available on all backends, the others on AVX2 and NEON.
:::

::: function
int vec4d_movemask(vec4d a)

Returns the sign bits of the lanes as an integer, lane $i$ giving bit
$i$. All backends.
:::

::: function
int vec1d_same(double a, double b)
int vec4d_same(vec4d a, vec4d b)
int vec8d_same(vec8d a, vec8d b)

Returns whether every lane of `a` compares equal to the same lane of `b`
with `==` (so $+0 = -0$, and a NaN is never equal). `vec1d_same` is
available on AVX2 and NEON, the others on AVX2 only.
:::

::: function
vec4d vec4d_bit_and(vec4d a, vec4d b)
vec4d vec4d_bit_or(vec4d a, vec4d b)
vec4d vec4d_bit_xor(vec4d a, vec4d b)
vec4d vec4d_bit_andnot(vec4d a, vec4d b)

Bitwise operations on the 64-bit patterns of the lanes;
`vec4d_bit_andnot` computes `(~a) & b`. All backends.
:::

::: function
vec4d vec4d_bits_set_u64(uint64_t a)
vec4d vec4d_bits_add(vec4d a, vec4d b)
vec4d vec4d_bits_sub(vec4d a, vec4d b)
vec4d vec4d_bits_shl(vec4d a, int n)
vec4d vec4d_bits_shr(vec4d a, int n)

Operations on the lanes viewed as `uint64_t` bit patterns: a broadcast,
addition and subtraction modulo $2^{64}$, and logical shifts by
$0 \le n < 64$. All backends.
:::

## Permutations {#permutations}

::: function
vec4d vec4d_unpacklo(vec4d a, vec4d b)
vec4d vec4d_unpackhi(vec4d a, vec4d b)

Returns $\{a_0, b_0, a_2, b_2\}$ and $\{a_1, b_1, a_3, b_3\}$. All
backends.
:::

::: function
vec4d vec4d_permute_0_1_1_0(vec4d a)
vec4d vec4d_permute_0_2_0_2(vec4d a)
vec4d vec4d_permute_0_2_1_3(vec4d a)
vec4d vec4d_permute_1_3_1_3(vec4d a)
vec4d vec4d_permute_3_1_2_0(vec4d a)
vec4d vec4d_permute_3_2_1_0(vec4d a)

`vec4d_permute_i_j_k_l` returns $\{a_i, a_j, a_k, a_l\}$. Only these six
permutations exist. All backends.
:::

::: function
vec4d vec4d_permute2_0_2(vec4d a, vec4d b)
vec4d vec4d_permute2_1_3(vec4d a, vec4d b)

Viewing `a` and `b` as four 128-bit halves $a_0a_1$, $a_2a_3$, $b_0b_1$,
$b_2b_3$ numbered 0 to 3, `vec4d_permute2_i_j` returns halves $i$ and
$j$: $\{a_0, a_1, b_0, b_1\}$ and $\{a_2, a_3, b_2, b_3\}$. All
backends.
:::

::: function
vec4d vec4d_unpack_lo_permute_0_2_1_3(vec4d u, vec4d v)
vec4d vec4d_unpack_hi_permute_0_2_1_3(vec4d u, vec4d v)
vec4d vec4d_unpacklo_permute_3_1_2_0(vec4d u, vec4d v)
vec4d vec4d_unpackhi_permute_3_1_2_0(vec4d u, vec4d v)

An unpack followed by a permutation, fused where the backend allows:
$\{u_0, u_2, v_0, v_2\}$, $\{u_1, u_3, v_1, v_3\}$,
$\{v_2, v_0, u_2, u_0\}$ and $\{v_3, v_1, u_3, u_1\}$. AVX2 and NEON.
:::

::: macro
VEC4D_TRANSPOSE(z0, z1, z2, z3, a0, a1, a2, a3)

Sets the `vec4d` variables `z0`, ..., `z3` to the rows of the transpose
of the $4 \times 4$ matrix with rows `a0`, ..., `a3`. The outputs may be
the inputs. AVX2 and NEON.
:::

::: function
vec4n vec4n_unpacklo(vec4n a, vec4n b)
vec4n vec4n_unpackhi(vec4n a, vec4n b)
vec4n vec4n_permute2_0_2(vec4n a, vec4n b)
vec4n vec4n_permute2_1_3(vec4n a, vec4n b)

The `vec4n` versions of the `vec4d` operations of the same names. AVX2
and NEON.
:::

::: function
vec4n vec4n_permute_3_2_1_0(vec4n a)

Returns $\{a_3, a_2, a_1, a_0\}$. AVX2 only.
:::

::: macro
VEC4N_TRANSPOSE(z0, z1, z2, z3, a0, a1, a2, a3)

The `vec4n` version of `VEC4D_TRANSPOSE`. AVX2 and NEON.
:::

## Integer vectors {#integer}

Arithmetic on `vec4n` lanes is modulo $2^{\bits}$.

::: function
vec4n vec4n_set_n(ulong a)
vec8n vec8n_set_n(ulong a)
vec4n vec4n_zero(void)

Returns a vector with all lanes equal to `a`, or to zero.
`vec8n_set_n` is available on AVX2 and NEON, the others on all
backends.
:::

::: function
vec4n vec4n_set_n4(ulong a0, ulong a1, ulong a2, ulong a3)

Returns $\{a_0, a_1, a_2, a_3\}$. AVX2 and the generic backends; not
NEON.
:::

::: function
vec4n vec4n_add(vec4n a, vec4n b)
vec4n vec4n_sub(vec4n a, vec4n b)

Returns $a + b$ and $a - b$ lane by lane, modulo $2^{\bits}$. All
backends.
:::

::: function
vec4n vec4n_bit_and(vec4n a, vec4n b)
vec4n vec4n_bit_or(vec4n a, vec4n b)
vec8n vec8n_bit_and(vec8n a, vec8n b)

Bitwise and and or. `vec8n_bit_and` is available on AVX2 and NEON, the
others on all backends.
:::

::: function
vec4n vec4n_bit_shift_left(vec4n a, ulong b)
vec4n vec4n_bit_shift_right(vec4n a, ulong b)
vec8n vec8n_bit_shift_right(vec8n a, ulong b)

Logical shifts of every lane by the same count $0 \le b < \bits$.
`vec8n_bit_shift_right` is available on AVX2 only, the others on all
backends.
:::

::: function
vec4n vec4n_bit_shift_right_32(vec4n a)
vec8n vec8n_bit_shift_right_32(vec8n a)

Logical shift right by 32 bits. These are macros on AVX2 and functions
on NEON.
:::

::: function
vec4n vec4n_shift_lanes_up(vec4n a, vec4n b)

Returns $\{b_3, a_0, a_1, a_2\}$: the lanes of `a` moved up by one, with
the top lane of `b` entering at the bottom, as when passing a carry from
one vector to the next. All backends.
:::

::: function
vec4n vec4n_mul(vec4n u, vec4n v)

Returns the 64-bit products of the low 32 bits of each lane,
$(u_i \bmod 2^{32})(v_i \bmod 2^{32})$. AVX2 only.
:::

::: function
ulong vec4n_horizontal_sum(vec4n a)

Returns $a_0 + a_1 + a_2 + a_3$ modulo $2^{64}$. AVX2 only.
:::

::: function
vec4n vec4n_addmod(vec4n a, vec4n b, vec4n n)
vec8n vec8n_addmod(vec8n a, vec8n b, vec8n n)

Returns $(a + b) \bmod n$ in $[0, n)$ for $a, b < n$, with no condition
on $n$. AVX2 and NEON.
:::

::: function
vec4n vec4n_addmod_limited(vec4n a, vec4n b, vec4n n)
vec8n vec8n_addmod_limited(vec8n a, vec8n b, vec8n n)

The same for $n < 2^{63}$, which allows a cheaper test. AVX2 and NEON.
:::

## Conversions {#conversions}

::: function
vec4n vec4d_convert_limited_vec4n(vec4d a)
vec4d vec4n_convert_limited_vec4d(vec4n a)
vec8d vec8n_convert_limited_vec8d(vec8n a)

Converts between `double` lanes and `ulong` lanes holding integers in
$[0, 2^{52})$, which are represented exactly in both. Outside this range
the results are unspecified: most of these conversions work on the bit
pattern of $2^{52}$. AVX2 and NEON.
:::

## Modular arithmetic in floating point {#modular}

These functions compute modulo a positive integer $n$ on integers stored
exactly in doubles, with `ninv` the double-precision value of `1.0/n`.
They form the arithmetic of `fft_small` and are available on AVX2 and
NEON only. [](algorithms/simd/fp-modular-arithmetic) states and proves
the ranges below. None of the conditions are checked, and none of the
functions require $n$ to be odd or prime.

::: function
vec1d vec1d_mulmod(vec1d a, vec1d b, vec1d n, vec1d ninv)
vec4d vec4d_mulmod(vec4d a, vec4d b, vec4d n, vec4d ninv)
vec8d vec8d_mulmod(vec8d a, vec8d b, vec8d n, vec8d ninv)
vec1d vec1d_nmulmod(vec1d a, vec1d b, vec1d n, vec1d ninv)
vec4d vec4d_nmulmod(vec4d a, vec4d b, vec4d n, vec4d ninv)
vec8d vec8d_nmulmod(vec8d a, vec8d b, vec8d n, vec8d ninv)

Returns an integer $r \equiv ab$ (`mulmod`) or $r \equiv -ab$
(`nmulmod`) modulo $n$, for integers $a, b$ with $|ab| < 2^{103}$ and
$1 \le n < 2^{50}$, given `ninv = 1.0 / n`. The
result satisfies $|r| \le n/2 + \eta|ab|$ with
$\eta = (1 + 2^{-53})^3 - 1 < 3.0001 \cdot 2^{-53}$. In particular, if
$n < 2^{49}$ or `fft_small_mulmod_satisfies_bounds(n)` holds, then
$|ab| < 2n^2$ gives $|r| < n$ and $|ab| < 4n^2$ gives $|r| < 3n/2$.
Larger products still give a correct representative, only a larger one.
AVX2 and NEON.
:::

::: function
vec1d vec1d_mulmod_fast(vec1d a, vec1d b, vec1d n, vec1d ninv)
vec4d vec4d_mulmod_fast(vec4d a, vec4d b, vec4d n, vec4d ninv)
vec8d vec8d_mulmod_fast(vec8d a, vec8d b, vec8d n, vec8d ninv)
vec1d vec1d_nmulmod_fast(vec1d a, vec1d b, vec1d n, vec1d ninv)
vec4d vec4d_nmulmod_fast(vec4d a, vec4d b, vec4d n, vec4d ninv)
vec8d vec8d_nmulmod_fast(vec8d a, vec8d b, vec8d n, vec8d ninv)

Like `mulmod` and `nmulmod`, with the quotient rounded by `round_mul`.
Requires $|ab| < 2n^2$ and $n < 2^{50}$. The result then satisfies
$|r| \le n/2 + \eta|ab|$ as for `mulmod`, and also
$|r| \le n/2 + 2.0001 \cdot 2^{-53}|ab|$ on the magic-constant path, so
$|r| < n$ if $n < 2^{49}$ or `fft_small_mulmod_satisfies_bounds(n)`
holds. With the magic constant, a product beyond $2^{51} n$ in absolute
value gives a wrong result, possibly not an integer. With native
rounding, and always on NEON, these functions are `mulmod` and `nmulmod`.
AVX2 and NEON.
:::

::: function
vec1d vec1d_round_mul(vec1d a, vec1d b)
vec4d vec4d_round_mul(vec4d a, vec4d b)

Returns $ab$ rounded to an integer. With the magic constant (see
`FLINT_FFT_SMALL_ROUND_USE_NATIVE`) this is the nearest integer to the
exact product, with ties to even, and requires $|ab| \le 2^{51}$; a
larger product gives a wrong value or a non-integer. With native rounding
it is `round(mul(a, b))`, which rounds twice and has no size condition.
AVX2 only.

::: details
The magic-constant path computes `fmadd(a, b, R) - R` with
$R = 1.5 \cdot 2^{52}$. When $|ab| \le 2^{51}$, the exact sum $ab + R$
lies in $[2^{52}, 2^{53}]$, where consecutive doubles are 1 apart, so the
single rounding of the FMA rounds $ab$ to the nearest integer. This
avoids `vroundpd`, which is slow on Intel cores.
:::
:::

::: macro
FLINT_FFT_SMALL_ROUND_USE_NATIVE
FLINT_ROUND_MUL_USE_ROUND
FLINT_ROUND_MUL_USE_MAGIC

`FLINT_FFT_SMALL_ROUND_USE_NATIVE` is set in the `flint-mparam.h` of
each target: 0 (magic constant) for most Intel targets, 1 (native
rounding) for AMD, ARM and generic targets. The header defaults it to 1.
Defining `FLINT_ROUND_MUL_USE_ROUND` or `FLINT_ROUND_MUL_USE_MAGIC`
forces either choice, for testing.
:::

::: function
vec1d vec1d_reduce_to_pm1n(vec1d a, vec1d n, vec1d ninv)
vec4d vec4d_reduce_to_pm1n(vec4d a, vec4d n, vec4d ninv)
vec8d vec8d_reduce_to_pm1n(vec8d a, vec8d n, vec8d ninv)
vec1d vec1d_reduce_to_pm1no(vec1d a, vec1d n, vec1d ninv)
vec4d vec4d_reduce_to_pm1no(vec4d a, vec4d n, vec4d ninv)
vec8d vec8d_reduce_to_pm1no(vec8d a, vec8d n, vec8d ninv)

Returns an integer $r = a - qn \equiv a \pmod n$, with $q$ an
approximation of $a/n$ rounded to an integer, for an integer
$|a| < 2^{53}$ and an integer $n \ge 4$. The two functions are
identical. The result satisfies $|r| < n/2 + 2$, hence $|r| < n$. On
AVX2 the quotient is computed with `round_mul`, which with the magic
constant needs $|a| / n \le 2^{51}$; this is why $n \ge 4$ is required.
AVX2 and NEON.
:::

::: function
vec1d vec1d_reduce_to_0n(vec1d a, vec1d n, vec1d ninv)
vec4d vec4d_reduce_to_0n(vec4d a, vec4d n, vec4d ninv)
vec8d vec8d_reduce_to_0n(vec8d a, vec8d n, vec8d ninv)

Returns $a \bmod n$ in $[0, n)$ for an integer $|a| < 2^{53}$ and an
integer $n \ge 4$, by `reduce_to_pm1no` followed by
`reduce_pm1no_to_0n`. On AVX2 with the magic constant, the `vec4d` and
`vec8d` versions return $n$ for the input $-0.0$. AVX2 and NEON.
:::

::: function
int vec1d_same_mod(vec1d a, vec1d b, vec1d n, vec1d ninv)
int vec4d_same_mod(vec4d a, vec4d b, vec4d n, vec4d ninv)

Returns whether $a_i \equiv b_i \pmod n$ in every lane, by comparing the
results of `reduce_to_0n`, with the same conditions. `vec1d_same_mod` is
available on AVX2 and NEON, `vec4d_same_mod` on AVX2 only.
:::

::: function
vec1d vec1d_reduce_pm1no_to_0n(vec1d a, vec1d n)
vec4d vec4d_reduce_pm1no_to_0n(vec4d a, vec4d n)
vec8d vec8d_reduce_pm1no_to_0n(vec8d a, vec8d n)

Returns $a \bmod n$ in $[0, n)$ for $a \in (-n, n)$, by adding $n$ to
negative lanes. The `vec4d` and `vec8d` versions on AVX2 test the sign
bit, so they map $-0.0$ to $n$; the others return $-0.0$. AVX2 and
NEON.
:::

::: function
vec1d vec1d_reduce_2n_to_n(vec1d a, vec1d n)
vec4d vec4d_reduce_2n_to_n(vec4d a, vec4d n)
vec8d vec8d_reduce_2n_to_n(vec8d a, vec8d n)

Returns $a \bmod n$ in $[0, n)$ for $a \in [0, 2n)$. AVX2 only.
:::

::: function
vec1d vec1d_reduce_0n_to_pmhn(vec1d a, vec1d n)
vec4d vec4d_reduce_0n_to_pmhn(vec4d a, vec4d n)

Returns $a$ or $a - n$, whichever lies in $(-n/2, n/2]$, for
$a \in [0, n]$. `vec1d_reduce_0n_to_pmhn` is available on AVX2 and NEON,
`vec4d_reduce_0n_to_pmhn` on AVX2 only.
:::

::: function
vec1d vec1d_reduce_pm1n_to_pmhn(vec1d a, vec1d n)
vec4d vec4d_reduce_pm1n_to_pmhn(vec4d a, vec4d n)
vec8d vec8d_reduce_pm1n_to_pmhn(vec8d a, vec8d n)

Returns $a$, $a - n$ or $a + n$, whichever lies in $[-n/2, n/2]$, for
$|a| \le 3n/2$ (in particular for $a \in [-n, n]$). AVX2 and NEON.
:::

## Matrix multiplication {#gemm}

These functions compute a matrix product in single or double precision.
They are always available: when FLINT is built with BLAS they call it by
default, and otherwise FLINT's own kernels are used. They are building
blocks for exact linear algebra over $\Z$ and $\Z/n\Z$, for example in
`nmod_mat_mul_blas` and `fmpz_mat_mul_blas`.

All of them set $C = AB$ for row-major matrices, where $C$ is
$m \times n$, $A$ is $m \times k$ and $B$ is $k \times n$, and
`lda`, `ldb` and `ldc` are the distances between the starts of
consecutive rows. $C$ is overwritten and must not overlap $A$ or $B$.
[](machine_vectors/gemm) compares this interface with
`cblas_dgemm` and discusses threading and exactness.

::: function
void flint_sgemm(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
void flint_dgemm(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

Sets $C = AB$, with the BLAS or the FLINT implementation according to
`flint_gemm_use_blas`.
:::

::: function
void flint_sgemm_fallback(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
void flint_dgemm_fallback(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

Sets $C = AB$ with FLINT's own kernels, which are always available,
handle any dimensions and leading dimensions without padding or
alignment, and use up to `flint_get_num_threads()` threads of FLINT's
thread pool.
:::

::: function
void flint_sgemm_blas(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
void flint_dgemm_blas(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

Sets $C = AB$ with `cblas_sgemm` or `cblas_dgemm`. Throws an exception
if FLINT was built without BLAS.
:::

::: variable
int flint_gemm_use_blas

Selects the implementation used by `flint_sgemm` and `flint_dgemm`:
initialised to 1 if FLINT was built with BLAS and 0 otherwise. It may be
changed at run time, for example to compare the two implementations;
setting it to 1 without BLAS makes every call throw.
:::

## NEON-only operations {#neon}

The NEON backend builds `vec4d` and `vec8d` from `vec2d` and exposes the
`vec2d`, `vec2n` and `vec1n` operations it uses. They have the semantics
of the `vec4d`, `vec4n` and `vec1d` operations of the same names.

::: function
vec2d vec2d_set_d2(double a0, double a1)
vec2d vec2d_one(void)
vec8d vec8d_one(void)
vec4d vec4d_set_vec2d2(vec2d a, vec2d b)
vec8d vec8d_set_vec4d2(vec4d a, vec4d b)
double vec2d_get_index(vec2d a, int i)
int vec2d_same(vec2d a, vec2d b)

Constructors and lane access, NEON only. `vec4d_set_vec2d2` and `vec8d_set_vec4d2`
concatenate two halves, `a` giving the low lanes.
:::

::: function
vec2d vec2d_neg(vec2d a)
vec2d vec2d_abs(vec2d a)
vec2d vec2d_min(vec2d a, vec2d b)
vec2d vec2d_max(vec2d a, vec2d b)
vec2d vec2d_half(vec2d a)
vec2d vec2d_div(vec2d a, vec2d b)
vec2d vec2d_sqrt(vec2d a)
vec2d vec2d_round(vec2d a)
vec2d vec2d_fmsub(vec2d a, vec2d b, vec2d c)
vec2d vec2d_fnmsub(vec2d a, vec2d b, vec2d c)
vec2d vec2d_blendv(vec2d a, vec2d b, vec2d c)

Arithmetic on two doubles, NEON only. `vec2d_fmsub` and `vec2d_fnmsub` are computed
as negated fused operations, which the header notes may be slower.
:::

::: function
vec2d vec2d_cmp_eq(vec2d a, vec2d b)
vec2d vec2d_cmp_ne(vec2d a, vec2d b)
vec2d vec2d_cmp_lt(vec2d a, vec2d b)
vec2d vec2d_cmp_le(vec2d a, vec2d b)
vec2d vec2d_cmp_gt(vec2d a, vec2d b)
vec2d vec2d_cmp_ge(vec2d a, vec2d b)
vec2d vec2d_cmp_nlt(vec2d a, vec2d b)
vec2d vec2d_cmp_nle(vec2d a, vec2d b)
vec2d vec2d_cmp_ngt(vec2d a, vec2d b)
vec2d vec2d_cmp_nge(vec2d a, vec2d b)
vec2d vec2d_bit_and(vec2d a, vec2d b)
vec2d vec2d_bit_or(vec2d a, vec2d b)
vec2d vec2d_bit_xor(vec2d a, vec2d b)
vec2d vec2d_bit_andnot(vec2d a, vec2d b)
int vec2d_movemask(vec2d a)
vec2d vec2d_bits_set_u64(uint64_t a)
vec2d vec2d_bits_add(vec2d a, vec2d b)
vec2d vec2d_bits_sub(vec2d a, vec2d b)
vec2d vec2d_bits_shl(vec2d a, int n)
vec2d vec2d_bits_shr(vec2d a, int n)

Comparisons, masks and bit operations on two doubles, with the
semantics of the `vec4d` versions. NEON only.
:::

::: function
vec2d vec2d_unpacklo(vec2d a, vec2d b)
vec2d vec2d_unpackhi(vec2d a, vec2d b)
vec8d vec8d_unpacklo(vec8d a, vec8d b)
vec8d vec8d_unpackhi(vec8d a, vec8d b)

`vec2d_unpacklo` returns $\{a_0, b_0\}$ and `vec2d_unpackhi` returns
$\{a_1, b_1\}$; the `vec8d` versions act on each `vec4d` half. NEON
only.
:::

::: function
vec2d vec2d_mulmod(vec2d a, vec2d b, vec2d n, vec2d ninv)
vec2d vec2d_nmulmod(vec2d a, vec2d b, vec2d n, vec2d ninv)
vec2d vec2d_mulmod_fast(vec2d a, vec2d b, vec2d n, vec2d ninv)
vec2d vec2d_nmulmod_fast(vec2d a, vec2d b, vec2d n, vec2d ninv)
vec2d vec2d_reduce_to_pm1n(vec2d a, vec2d n, vec2d ninv)
vec2d vec2d_reduce_to_pm1no(vec2d a, vec2d n, vec2d ninv)
vec2d vec2d_reduce_to_0n(vec2d a, vec2d n, vec2d ninv)
vec2d vec2d_reduce_pm1no_to_0n(vec2d a, vec2d n)
vec2d vec2d_reduce_0n_to_pmhn(vec2d a, vec2d n)
vec2d vec2d_reduce_pm1n_to_pmhn(vec2d a, vec2d n)
int vec2d_same_mod(vec2d a, vec2d b, vec2d n, vec2d ninv)

Modular arithmetic on two doubles, with native rounding and the
semantics of the `vec4d` versions. NEON only.
:::

::: function
vec2n vec2n_load_unaligned(const ulong * a)
void vec2n_store_unaligned(ulong * z, vec2n a)
void vec1n_store_unaligned(ulong * z, vec1n a)
vec2n vec2n_set_n(ulong a)
vec2n vec2n_add(vec2n a, vec2n b)
vec2n vec2n_sub(vec2n a, vec2n b)
vec8n vec8n_add(vec8n a, vec8n b)
vec8n vec8n_sub(vec8n a, vec8n b)
vec2n vec2n_bit_and(vec2n a, vec2n b)
vec2n vec2n_bit_or(vec2n a, vec2n b)
vec2n vec2n_bit_shift_left(vec2n a, ulong b)
vec2n vec2n_bit_shift_right(vec2n a, ulong b)
vec2n vec2n_bit_shift_right_32(vec2n a)
vec2n vec2n_unpacklo(vec2n a, vec2n b)
vec2n vec2n_unpackhi(vec2n a, vec2n b)
vec2n vec2n_addmod(vec2n a, vec2n b, vec2n n)
vec2n vec2n_addmod_limited(vec2n a, vec2n b, vec2n n)
vec1n vec1n_addmod(vec1n a, vec1n b, vec1n n)

Operations on one, two and eight words, with the semantics of the
`vec4n` operations of the same names. NEON only.
:::

::: function
vec1n vec1d_convert_limited_vec1n(vec1d a)
vec2n vec2d_convert_limited_vec2n(vec2d a)
vec2d vec2n_convert_limited_vec2d(vec2n a)

Conversions for integers in $[0, 2^{52})$. Conversions to words round
non-integers to nearest for `vec2d` and truncate for `vec1d`. NEON
only.
:::
