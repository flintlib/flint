/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"
#include "crt_helpers.h"

/* transpose a block */
void _convert_block(
    ulong* Xs,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    ulong np,
    ulong I)
{
    for (ulong l = 0; l < np; l++)
    {
        vec4d p = vec4d_set_d(Rffts[l].p);
        vec4d pinv = vec4d_set_d(Rffts[l].pinv);
        double* x = sd_fft_ctx_blk_index(d + l*dstride, I);
        ulong j = 0; do {
            vec4d x0, x1, x2, x3;
            vec4n y0, y1, y2, y3;
            x0 = vec4d_load(x + j + 0*VEC_SZ);
            x1 = vec4d_load(x + j + 1*VEC_SZ);
            x2 = vec4d_load(x + j + 2*VEC_SZ);
            x3 = vec4d_load(x + j + 3*VEC_SZ);
            x0 = vec4d_reduce_to_0n(x0, p, pinv);
            x1 = vec4d_reduce_to_0n(x1, p, pinv);
            x2 = vec4d_reduce_to_0n(x2, p, pinv);
            x3 = vec4d_reduce_to_0n(x3, p, pinv);
            y0 = vec4d_convert_limited_vec4n(x0);
            y1 = vec4d_convert_limited_vec4n(x1);
            y2 = vec4d_convert_limited_vec4n(x2);
            y3 = vec4d_convert_limited_vec4n(x3);
            vec4n_store_unaligned(Xs + l*BLK_SZ + j + 0*VEC_SZ, y0);
            vec4n_store_unaligned(Xs + l*BLK_SZ + j + 1*VEC_SZ, y1);
            vec4n_store_unaligned(Xs + l*BLK_SZ + j + 2*VEC_SZ, y2);
            vec4n_store_unaligned(Xs + l*BLK_SZ + j + 3*VEC_SZ, y3);
        } while (j += 4*VEC_SZ, j < BLK_SZ);
        FLINT_ASSERT(j == BLK_SZ);
    }
}

ulong flint_mpn_nbits(const ulong* a, ulong an)
{
    while (an > 0 && a[an-1] == 0)
        an--;

    if (an == 0)
        return 0;

    return FLINT_BITS*(an - 1) + n_nbits(a[an-1]);
}

/* cmp(a, b*2^e), a does not have to be normalized */
int flint_mpn_cmp_ui_2exp(const ulong* a, ulong an, ulong b, ulong e)
{
    ulong q = e/FLINT_BITS;
    ulong r = e%FLINT_BITS;
    ulong x, b0, b1;

    while (an > 0 && a[an-1] == 0)
        an--;

    if (an == 0)
        return b != 0;

    // b*2^e = (b*2^r       )*2^(64*q)
    //       = (b0 + b1*2^64)*2^(64*q)
    if (r == 0)
    {
        b0 = b;
        b1 = 0;
    }
    else
    {
        b0 = b << r;
        b1 = b >> (FLINT_BITS - r);
    }

    //      check words [q+2,infty)
    // then check words [q+1, 64*q+128) against b1
    // then check words [q, q+1) against b0
    // then check words [0, q)

    if (an > q + 2)
        return 1;

    x = (q+1 < an) ? a[q+1] : 0;
    if (x != b1)
        return x > b1 ? 1 : -1;

    x = (q < an) ? a[q] : 0;
    if (x != b0)
        return x > b0 ? 1 : -1;

    q = n_min(q, an);
    while (q > 0)
    {
        q--;
        if (a[q] != 0)
            return 1;
    }

    return 0;
}


unsigned char flint_mpn_add_inplace_c(ulong* z, ulong zn, ulong* a, ulong an, unsigned char cf)
{
    FLINT_ASSERT(zn >= an);
    FLINT_ASSERT(an > 0);

    ulong i = 0;
    do {
        cf = _addcarry_ulong(cf, z[i], a[i], &z[i]);
    } while (i++, i < an);

    while (i < zn && cf != 0)
    {
        cf = _addcarry_ulong(cf, z[i], 0, &z[i]);
        i++;
    }

    return cf;
}
