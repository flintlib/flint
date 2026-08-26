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
    (n < 2^63) we fold the top limbs using r64 = 2^64 mod n and
    r128 = 2^128 mod n and then apply Barrett reduction with the two-limb
    precomputed quotient (qhi, qlo) = floor(2^128 / n), as in fft_small.
    For fullword moduli we use Granlund-Möller after the folding.
*/

/* 2 -> 1 limb mod, n < 2^(FLINT_BITS-1), x < ~2^(FLINT_BITS) n */
FLINT_FORCE_INLINE ulong
_n_ll_rem_l_nonfullword(ulong xhi, ulong xlo, ulong n, ulong qhi, ulong qlo)
{
    ulong c2, c1, c0;
    FLINT_MPN_MUL_3P2X2(c2, c1, c0, qhi, qlo, xhi, xlo);
    (void) c1;
    (void) c0;
    xlo -= c2 * n;
    if (xlo >= n)
        xlo -= n;
    return xlo;
}

/* arbitrary (hi, lo) mod n */
FLINT_FORCE_INLINE ulong
_nmod_red2_any(ulong hi, ulong lo, const flint_mpn_crt_leaf_struct * L)
{
    ulong h, l, u;
    nmod_t mod = L->mod;
    umul_ppmm(h, l, hi, L->r64);
    add_ssaaaa(h, l, h, l, 0, lo);

    if (mod.norm != 0)
    {
        return _n_ll_rem_l_nonfullword(h, l, mod.n, L->qhi, L->qlo);
    }
    else
    {
        if (h >= mod.n)
            h -= mod.n;
        NMOD_RED2_FULLWORD(u, h, l, mod);
        return u;
    }
}

/* (y2, y1, y0) mod n with y2 < n */
FLINT_FORCE_INLINE ulong
_nmod_red3_any(ulong y2, ulong y1, ulong y0, const flint_mpn_crt_leaf_struct * L)
{
    ulong c1, c0, t1, t0, xhi, xlo, cy;
    nmod_t mod = L->mod;

    umul_ppmm(t1, t0, y2, L->r128);
    umul_ppmm(c1, c0, y1, L->r64);

    if (mod.norm != 0)
    {
        /* no carry: sum < 2^64 (n + n) + 2^64 <= 2^128 */
        add_ssaaaa(xhi, xlo, t1, t0, c1, c0);
        add_ssaaaa(xhi, xlo, xhi, xlo, 0, y0);
        return _n_ll_rem_l_nonfullword(xhi, xlo, mod.n, L->qhi, L->qlo);
    }
    else
    {
        add_sssaaaaaa(cy, xhi, xlo, 0, t1, t0, 0, c1, c0);
        add_sssaaaaaa(cy, xhi, xlo, cy, xhi, xlo, 0, 0, y0);
        while (cy)
            add_sssaaaaaa(cy, xhi, xlo, 0, xhi, xlo, 0, 0, L->r128);
        if (xhi >= mod.n)
            xhi -= mod.n;
        NMOD_RED2_FULLWORD(xlo, xhi, xlo, mod);
        return xlo;
    }
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

#endif
