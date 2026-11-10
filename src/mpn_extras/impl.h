/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_EXTRAS_IMPL_H
#define MPN_EXTRAS_IMPL_H

#include "longlong.h"
#include "nmod.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "crt_helpers.h"
#include "mpn_extras.h"

/* maximum number of primes batched into a single-limb leaf */
#define FLINT_MPN_CRT_MAX_LEAF_PRIMES 8

/* slack limbs on the slots of the ping-pong buffers: the CRT chunk values
   have up to 2 extra limbs (unreduced sum of up to 2^64 * b terms) and
   the tree adds at most one more bit per level */
#define FLINT_MPN_CRT_SLOT_EXTRA 3

/* maximum limb size for the fixed-length (single chunk) paths */
#define FLINT_MPN_CRT_FIXED_MAX 8

/* product size from which the small-value shortcut is set up, and the
   size of the product Q of leading primes it reconstructs modulo */
#define FLINT_MPN_CRT_SMALL_CHECK_LIMBS 16
#define FLINT_MPN_CRT_SMALL_CHECK_LIMBS_GENERIC 64
#define FLINT_MPN_CRT_SMALL_CHECK_Q_LIMBS 4

#define MAC(h, l, a, b)                 \
do {                                    \
    ulong p1, p0;                       \
    umul_ppmm(p1, p0, a, b);            \
    add_ssaaaa(h, l, h, l, p1, p0);     \
} while (0)

/*
    Tail reductions for the dot products. For nonfullword moduli
    (n < 2^63), Barrett reduction with the two-limb precomputed quotient
    (qhi, qlo) = floor(2^128 / n), which takes any two limbs; for fullword
    moduli, Granlund-Möller. Three limbs are first folded into two using
    r64 = 2^64 mod n and r128 = 2^128 mod n, before the test on n (the fold
    of n_lll_mod_barrett and NMOD_RED3_FOLD_FULLWORD; done first, it is a
    little faster here).
*/

/* arbitrary (hi, lo) mod n */
FLINT_FORCE_INLINE ulong
_nmod_red2_any(ulong hi, ulong lo, const flint_mpn_crt_leaf_struct * L)
{
    ulong u;
    nmod_t mod = L->mod;

    if (mod.norm != 0)
        return n_ll_mod_barrett(hi, lo, mod.n, L->qhi, L->qlo);

    /* hi < 2^64 < 2 n */
    if (hi >= mod.n)
        hi -= mod.n;
    NMOD_RED2_FULLWORD(u, hi, lo, mod);
    return u;
}

/* arbitrary (y2, y1, y0) mod n */
FLINT_FORCE_INLINE ulong
_nmod_red3_any(ulong y2, ulong y1, ulong y0, const flint_mpn_crt_leaf_struct * L)
{
    ulong h, l, t1, t0, r;
    nmod_t mod = L->mod;

    umul_ppmm(h, l, y2, L->r128);
    umul_ppmm(t1, t0, y1, L->r64);
    add_ssaaaa(h, l, h, l, t1, t0);
    add_ssaaaa(h, l, h, l, UWORD(0), y0);

    if (mod.norm != 0)
        return n_ll_mod_barrett(h, l, mod.n, L->qhi, L->qlo);

    /* h < 2 n for fullword n (see NMOD_RED3_FOLD_FULLWORD) */
    if (h >= mod.n)
        h -= mod.n;
    NMOD_RED2_FULLWORD(r, h, l, mod);
    return r;
}

/* 1-limb reduction with Barrett precomputation (n = 1 supported with npre = UWORD_MAX) */
FLINT_FORCE_INLINE ulong
_n_mod_barrett_any(ulong x, ulong n, ulong npre)
{
    return n_mod_barrett(x, n, npre);
}

/* whether to use threads for n nodes of the given limb size */
#define FLINT_MPN_CRT_PARALLEL_LIMBS 2000

FLINT_FORCE_INLINE int
_crt_use_threads(slong n, slong limbs)
{
    return n >= 2 && limbs >= FLINT_MPN_CRT_PARALLEL_LIMBS && flint_get_num_threads() > 1;
}

/* layout of the workspace */
typedef struct
{
    nn_ptr bufA, bufB;      /* level value buffers (ping-pong) */
    nn_ptr * ptrA, * ptrB;  /* per node: pointer to value */
    slong * lenA, * lenB;   /* per node: length of value */
    nn_ptr scratch;         /* division / multiplication scratch */
    slong scratch_limbs;
}
_crt_work_struct;

FLINT_FORCE_INLINE void
_crt_work_init(_crt_work_struct * W, const flint_mpn_crt_t C, nn_ptr tmp)
{
    W->bufA = tmp;
    W->bufB = W->bufA + C->work_level_limbs;
    W->ptrA = (nn_ptr *) (W->bufB + C->work_level_limbs);
    W->ptrB = W->ptrA + C->work_max_count;
    W->lenA = (slong *) (W->ptrB + C->work_max_count);
    W->lenB = W->lenA + C->work_max_count;
    W->scratch = (nn_ptr) (W->lenB + C->work_max_count);
    W->scratch_limbs = C->tmp_limbs - (W->scratch - tmp);
}

/* gcd and half-gcd (hgcd.c, gcd.c, gcdext.c, hgcd_fft.c) ********************/

/* Cutoffs in limbs, with the meaning of GMP's HGCD_THRESHOLD,
   HGCD_APPR_THRESHOLD, HGCD_REDUCE_THRESHOLD, GCD_DC_THRESHOLD,
   GCDEXT_DC_THRESHOLD and MATRIX22_STRASSEN_THRESHOLD; flint-mparam.h may
   override the defaults. */
#ifndef FLINT_MPN_HGCD_CUTOFF
# define FLINT_MPN_HGCD_CUTOFF 160
#endif
#ifndef FLINT_MPN_HGCD_APPR_CUTOFF
# define FLINT_MPN_HGCD_APPR_CUTOFF 60
#endif
#ifndef FLINT_MPN_HGCD_REDUCE_CUTOFF
# define FLINT_MPN_HGCD_REDUCE_CUTOFF 2681
#endif
/* The first step of the half-gcd based gcd computes the half-gcd of the
   top n - p = n/3 limbs, and that of the extended gcd of n/2 and then
   2n/3 limbs, so that these algorithms pay off from about where those
   half-gcds are past the half-gcd cutoff (as measured on Intel and AMD
   machines with tune-gcd, within the flat optimum). */
#ifndef FLINT_MPN_GCD_DC_CUTOFF
# define FLINT_MPN_GCD_DC_CUTOFF (3 * FLINT_MPN_HGCD_CUTOFF)
#endif
#ifndef FLINT_MPN_GCDEXT_DC_CUTOFF
# define FLINT_MPN_GCDEXT_DC_CUTOFF (3 * FLINT_MPN_HGCD_CUTOFF / 2)
#endif
#ifndef FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF
# define FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF 16
#endif

/* The cofactor t = (g - s a) / b of b, given g = gcd(a, b) and a cofactor s
   of a with g = s a + t b, where a, b > 0 and sn is the signed size of s.
   Writes t to tp (room for an limbs) and its signed size to *tn. */
void _flint_mpn_gcdext_cofactor(mp_ptr tp, mp_size_t * tn,
    mp_srcptr gp, mp_size_t gn, mp_srcptr sp, mp_size_t sn,
    mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn);

/* Transform-space (fft_small) versions of the half-gcd matrix operations.
   Each returns 0 (doing nothing) when the operands are too small for
   fft_small to pay off, in which case the caller falls back to ordinary
   products. */
#if FLINT_HAVE_FFT_SMALL
int _flint_mpn_hgcd_matrix_mul_fft(flint_mpn_hgcd_matrix_struct * M,
    const flint_mpn_hgcd_matrix_struct * M1);
mp_size_t _flint_mpn_hgcd_matrix_adjust_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p);
int _flint_mpn_hgcd_matrix_apply_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_size_t mn[2][2], mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t nn);
mp_size_t _flint_mpn_hgcd_matrix_adjust_mul_fft(flint_mpn_hgcd_matrix_struct * M,
    const flint_mpn_hgcd_matrix_struct * M1,
    mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p);
mp_size_t _flint_mpn_hgcd_mul_matrix_vector_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n);
#endif

#endif
