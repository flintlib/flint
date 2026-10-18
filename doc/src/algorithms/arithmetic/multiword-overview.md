---
title: Arithmetic modulo multiword integers
kind: topic
summary: FLINT represents $\Z/n\Z$ for $n \ge 2^{64}$ in two ways, with `fmpz` residues in `fmpz_mod` and with fixed-length limb arrays in the generic ring `mpn_mod`. This page compares the representations, the reductions they use, and how they relate to `nmod` and to the generic rings.
order: 50
---

For a modulus $n$ that fits in a word, FLINT stores residues as `ulong`
values and reduces with the precomputed one-word inverses of `nmod`
([](algorithms/arithmetic/choosing)). For larger moduli there are two
representations:

| | `nmod` | `mpn_mod` | `fmpz_mod` |
|---|---|---|---|
| Moduli | $1 \le n < 2^{64}$ | $2^{64} \le n < 2^{1024}$ (2 to 16 limbs) | any $n \ge 1$ |
| Residue | `ulong` | $\ell$ limbs, fixed length | `fmpz` (inline below $2^{62}$, else `mpz`) |
| Memory per residue | 1 word | $\ell$ words, no allocation | 1 word plus a heap `mpz` |
| Modulus data | `nmod_t` by value | `gr_ctx_t` | `fmpz_mod_ctx_t` |
| Vectors, polynomials, matrices | `nmod_vec`, `nmod_poly`, `nmod_mat` | `gr_poly_t`, `gr_mat_t` over the ring | `fmpz_mod_vec`, `fmpz_mod_poly`, `fmpz_mod_mat` |
| Generic ring | `gr_ctx_init_nmod` | `gr_ctx_init_mpn_mod` | `gr_ctx_init_fmpz_mod` |

(For 32-bit machines, read $2^{32}$ and $2^{512}$.)

## fmpz_mod {#fmpz-mod}

`fmpz_mod` stores residues as ordinary `fmpz` integers in $[0, n)$ and
works for every modulus. The context `fmpz_mod_ctx_t` holds $n$ and
function pointers for addition, subtraction and multiplication that
`fmpz_mod_ctx_init` chooses by the size of $n$:

| Modulus | Multiplication | Reduction |
|---|---|---|
| $n < 2^{64}$ | `_fmpz_mod_mul1` | `nmod_mul` ([Möller–Granlund](algorithms/arithmetic/moller-granlund)) |
| $n = 2^{64}$ | `_fmpz_mod_mul2s` | wrapping word product |
| $2^{64} < n < 2^{128}$ | `_fmpz_mod_mul2` | [two-limb Barrett](algorithms/arithmetic/multiword-mulmod#two-limb), at most two corrections |
| up to 20 limbs | `_fmpz_mod_mulN` | `fmpz_mul` and `fmpz_mod` (FLINT's general division, `flint_mpz_mod`) |
| more than 20 limbs | `_fmpz_mod_mulN` | `fmpz_mul` and `fmpz_fdiv_r_preinvn` (Barrett with a precomputed inverse) |

Additions and subtractions use one conditional correction in the same
size classes. Powers use `fmpz_powm` (GMP's `mpz_powm`, which works in
Montgomery representation); when `fft_small` is available, moduli of
70000 bits and more use sliding window exponentiation with
`fmpz_mod_mul`, whose products then use the FFT. Because a
residue above $2^{62}$ is a pointer to an `mpz`, every operation goes
through `fmpz` arithmetic with its size checks and possible allocation,
which dominates the cost for moduli of a few limbs.

The vector and matrix modules reduce sums once: `_fmpz_mod_vec_dot`
computes the integer dot product with `_fmpz_vec_dot`, and
`fmpz_mod_mat_mul` multiplies with `fmpz_mat_mul` (which uses
multimodular and BLAS-based methods) and reduces the entries. Most
`fmpz_mod_mat` and `fmpz_mod_poly` algorithms beyond that are wrappers
around the generic `gr_mat` and `gr_poly` algorithms over the
`fmpz_mod` ring (`_gr_ctx_init_fmpz_mod_from_ref` builds a generic
context that refers to an existing `fmpz_mod_ctx_t`).

## mpn_mod {#mpn-mod}

`mpn_mod` stores a residue as exactly $\ell$ limbs, where $\ell$ is the
number of limbs of $n$. Residues have no header, own no memory and have
a size known from the context, so vectors, polynomials and matrices are
contiguous arrays and temporaries live on the stack. Arithmetic is done
directly on the limbs with `mpn` functions and FLINT's fixed-size
multiplication code:

- additions and subtractions use one conditional correction, with
  inline two-limb code;
- multiplication reduces the $2\ell$-limb product with a precomputed
  $\ell$-limb inverse of the normalised modulus
  ([](algorithms/arithmetic/multiword-mulmod#preinvn)); multiplication by
  a fixed scalar uses Shoup's method or the matrix method;
- sums of products (dot products, polynomial products and divisions,
  Gaussian elimination, matrix products) are accumulated in $2\ell - 1$
  to $2\ell + 1$ limbs and reduced once
  ([](algorithms/arithmetic/multiword-delayed-reduction)).

The limit of 16 limbs keeps temporaries small enough for the stack and
covers the sizes for which the module is tuned. `mpn_mod` exists only as
a generic ring (`gr_ctx_init_mpn_mod`): there are no `mpn_mod_poly` or
`mpn_mod_mat` types, and polynomials and matrices are `gr_poly_t` and
`gr_mat_t`. The ring overrides the vector, polynomial and matrix methods
that the generic algorithms call, so `gr_poly_gcd`, `gr_mat_lu`,
`gr_poly_factor` and the like run with the specialised kernels.

## Comparison {#comparison}

Both representations implement the same generic ring interface, so
generic code can be run over either, and the difference is speed. For a
single multiplication, `fmpz_mod_mul` with a two-limb modulus is a
competitive inline Barrett reduction, but beyond two limbs it calls a
general division (`fmpz_mod`, or the Barrett division of
`fmpz_fdiv_r_preinvn` above 20 limbs), and every operation pays for the
`fmpz` layer. The
fixed-length representation wins in the basecase algorithms (Euclid's
algorithm and Gaussian elimination with delayed reduction), allocates no
temporaries in divide-and-conquer algorithms such as half-gcd and block
LU decomposition, and makes additions cheap enough that Karatsuba
multiplication and Waksman's matrix multiplication pay off. The
reported speedups of `mpn_mod` over the `fmpz_mod` ring are between 1x
and 3x for `gr_poly_gcd` and `gr_mat_nonsingular_solve` with moduli of
128 to 1024 bits [@Joh2025gr].

For long polynomials, both rely on fast multiplication: `fmpz_mod_poly`
through `fmpz_poly` (Kronecker substitution and the small-prime FFT),
`mpn_mod` through its own Kronecker substitution and small-prime FFT
code. The FFT path of `mpn_mod` and its transformed representation of
polynomials are available for moduli up to about 190 bits.

## Which to use {#advice}

- For one-word moduli, use `nmod` (or the `nmod` generic ring).
- For moduli of 2 to 16 limbs and computations with polynomials or
  matrices, use generic `gr_poly` and `gr_mat` code over
  `gr_ctx_init_mpn_mod`.
- For larger moduli, for scalar code that must work for every modulus
  with one type, or when the `fmpz_mod_poly` and `fmpz_mod_mat`
  interfaces (factorisation, Howell forms, ...) are needed, use
  `fmpz_mod`.

Code that handles all sizes, such as FLINT's own number-theoretic
algorithms, chooses the ring at run time: `nmod` for one limb, `mpn_mod`
when `gr_ctx_init_mpn_mod` succeeds, and the `fmpz_mod` ring otherwise.

## Use in FLINT {#use}

- The APRCL primality test (`aprcl_is_prime`) computes in
  $\Z[\zeta_{p^k}]/(n)$ with `mpn_mod` coefficients when $n$ has 2 to 16
  limbs (file `aprcl/unity_zp_mpn.c`); this is the bottleneck of the
  test.
- The ECPP primality prover does its elliptic curve point arithmetic over
  `mpn_mod` for moduli of 2 to 16 limbs and over `fmpz_mod` beyond (file
  `ecpp/point_gr.c`).
- `fmpz_lucas_chain` evaluates Lucas sequences over `nmod`, `mpn_mod` or
  `fmpz_mod`.
- `fmpz_poly_hensel_lift_tree_recursive` lifts factorisations over
  `nmod`, `mpn_mod` (up to 12 limbs, where it beats `fmpz_mod` at the
  lengths that occur) or `fmpz_mod`.
- `fmpz_mod_poly`, `fmpz_mod_mat` and `fq` with a multiword prime use
  the `fmpz_mod` representation.
