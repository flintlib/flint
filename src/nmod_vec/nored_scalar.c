/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "nmod_vec.h"


#if FLINT_BITS == 64 && defined(__AVX2__)

#include "machine_vectors.h"

void
_nmod_vec_nored_scalar_addmul_halflimb_avx2(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    vec4n vc = _mm256_set1_epi32((unsigned int) c);
    slong i = 0;

    /* vec4n_mul does a 32x32 -> 64 bit mul */

    for ( ; i + 3 < len; i += 4)
    {
        vec4n va = vec4n_load_unaligned(res + i);
        vec4n vb = vec4n_load_unaligned(vec + i);
        vec4n_store_unaligned(res + i, vec4n_add(va, vec4n_mul(vb, vc)));
    }

    for ( ; i < len; i++)
        res[i] += vec[i] * c;
}

#endif



#define ADD2(r, hi, lo) add_ssaaaa(*((r)+1), *((r)), *((r)+1), *(r), hi, lo);
#define ADD3(r, hi, lo) add_sssaaaaaa(*((r)+2), *((r)+1), *(r), *((r)+2), *((r)+1), *(r), 0, hi, lo);


static void
_nmod_vec_nored_ll_scalar_addmul_halflimb_unroll4(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    slong i;

    for (i = 0; i + 3 < len; i += 4)
    {
        ADD2(res + 2 * (i + 0), 0, vec[i + 0] * c);
        ADD2(res + 2 * (i + 1), 0, vec[i + 1] * c);
        ADD2(res + 2 * (i + 2), 0, vec[i + 2] * c);
        ADD2(res + 2 * (i + 3), 0, vec[i + 3] * c);
    }

    for ( ; i < len; i++)
        ADD2(res + 2 * i, 0, vec[i] * c);
}

void
_nmod_vec_nored_ll_scalar_addmul_halflimb(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    slong i;

    if (len >= 12)
    {
        _nmod_vec_nored_ll_scalar_addmul_halflimb_unroll4(res, vec, len, c);
        return;
    }

    for (i = 0; i < len; i++)
        ADD2(res + 2 * i, 0, vec[i] * c);
}

static void
_nmod_vec_nored_ll_scalar_addmul_1(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        ulong hi, lo;
        umul_ppmm(hi, lo, vec[i], c);
        ADD2(res + 2 * i, hi, lo);
    }
}


static void
_nmod_vec_nored_ll_scalar_addmul_unroll4(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    slong i;

    for (i = 0; i + 3 < len; i += 4)
    {
        ulong l0, l1, l2, l3, h0, h1, h2, h3;

        umul_ppmm(h0, l0, vec[i + 0], c);
        ADD2(res + 2 * (i + 0), h0, l0);
        umul_ppmm(h1, l1, vec[i + 1], c);
        ADD2(res + 2 * (i + 1), h1, l1);
        umul_ppmm(h2, l2, vec[i + 2], c);
        ADD2(res + 2 * (i + 2), h2, l2);
        umul_ppmm(h3, l3, vec[i + 3], c);
        ADD2(res + 2 * (i + 3), h3, l3);
    }

    for ( ; i < len; i++)
    {
        ulong hi, lo;
        umul_ppmm(hi, lo, vec[i], c);
        ADD2(res + 2 * i, hi, lo);
    }
}

void
_nmod_vec_nored_ll_scalar_addmul(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    if (len < 24)
        _nmod_vec_nored_ll_scalar_addmul_1(res, vec, len, c);
    else
        _nmod_vec_nored_ll_scalar_addmul_unroll4(res, vec, len, c);
}


static void
_nmod_vec_nored_lll_scalar_addmul_unroll4(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    slong i;

    for (i = 0; i + 3 < len; i += 4)
    {
        ulong l0, l1, l2, l3, h0, h1, h2, h3;

        umul_ppmm(h0, l0, vec[i + 0], c);
        umul_ppmm(h1, l1, vec[i + 1], c);
        umul_ppmm(h2, l2, vec[i + 2], c);
        umul_ppmm(h3, l3, vec[i + 3], c);

        ADD3(res + 3 * (i + 0), h0, l0);
        ADD3(res + 3 * (i + 1), h1, l1);
        ADD3(res + 3 * (i + 2), h2, l2);
        ADD3(res + 3 * (i + 3), h3, l3);
    }

    for ( ; i < len; i++)
    {
        ulong hi, lo;
        umul_ppmm(hi, lo, vec[i], c);
        ADD3(res + 3 * i, hi, lo);
    }
}

static void
_nmod_vec_nored_lll_scalar_addmul_1(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        ulong hi, lo;
        umul_ppmm(hi, lo, vec[i], c);
        ADD3(res + 3 * i, hi, lo);
    }
}

void
_nmod_vec_nored_lll_scalar_addmul(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
    if (len < 24)
        _nmod_vec_nored_lll_scalar_addmul_1(res, vec, len, c);
    else
        _nmod_vec_nored_lll_scalar_addmul_unroll4(res, vec, len, c);
}

