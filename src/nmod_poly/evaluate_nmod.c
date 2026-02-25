/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod.h"

/*--------------------------------------*/
/* special value: 1 (one) and -1 (mone) */
/*--------------------------------------*/

/* general */
static ulong _nmod_poly_evaluate_one(nn_srcptr poly, slong len, ulong modn)
{
    if (len == 0)
        return 0;

    if (len <= 3)
    {
        ulong val = poly[0];
        for (slong m = 1; m < len; m++)
            val = n_addmod(val, poly[m], modn);
        return val;
    }

    ulong val0, val1, val2, val3;

    val0 = poly[0];
    val1 = poly[1];
    val2 = poly[2];
    val3 = poly[3];

    slong m = 4;

    for ( ; m+3 < len; m += 4)
    {
        val0 = n_addmod(val0, poly[m+0], modn);
        val1 = n_addmod(val1, poly[m+1], modn);
        val2 = n_addmod(val2, poly[m+2], modn);
        val3 = n_addmod(val3, poly[m+3], modn);
    }

    /* gather results */
    val2 = n_addmod(val2, val3, modn);
    val0 = n_addmod(val0, val1, modn);
    val0 = n_addmod(val0, val2, modn);

    /* last few terms */
    for ( ; m < len; m++)
        val0 = n_addmod(val0, poly[m], modn);

    return val0;
}

/* requires nbits(modn) <= FLINT_BITS - 1 */
static ulong _nmod_poly_evaluate_one1(nn_srcptr poly, slong len, ulong modn)
{
    if (len == 0)
        return 0;

    if (len <= 3)
    {
        ulong val = poly[0];
        for (slong m = 1; m < len; m++)
            val = n_addmod(val, poly[m], modn);
        return val;
    }

    ulong val0, val1, val2, val3;

    val0 = poly[0];
    val1 = poly[1];
    val2 = poly[2];
    val3 = poly[3];

    slong m = 4;

    for ( ; m+3 < len; m += 4)
    {
        val0 = val0 + poly[m+0];
        val1 = val1 + poly[m+1];
        val2 = val2 + poly[m+2];
        val3 = val3 + poly[m+3];
        if (val0 >= modn)
            val0 -= modn;
        if (val1 >= modn)
            val1 -= modn;
        if (val2 >= modn)
            val2 -= modn;
        if (val3 >= modn)
            val3 -= modn;
    }

    /* handle 2 more terms if they exist */
    if (m+1 < len)
    {
        val0 = val0 + poly[m+0];
        val1 = val1 + poly[m+1];
        if (val0 >= modn)
            val0 -= modn;
        if (val1 >= modn)
            val1 -= modn;
        m += 2;
    }

    /* remains 0 or 1 term */
    if (m < len)
    {
        val0 = val0 + poly[m+0];
        if (val0 >= modn)
            val0 -= modn;
    }

    /* gather results */
    val2 = n_addmod(val2, val3, modn);
    val0 = n_addmod(val0, val1, modn);
    val0 = n_addmod(val0, val2, modn);

    return val0;
}

/* requires nbits(modn) <= FLINT_BITS - 2 */
static ulong _nmod_poly_evaluate_one2(nn_srcptr poly, slong len, ulong modn)
{
    if (len == 0)
        return 0;

    if (len <= 7)
    {
        ulong val = poly[0];
        for (slong m = 1; m < len; m++)
            val = n_addmod(val, poly[m], modn);
        return val;
    }

    ulong val0, val1, val2, val3;

    /* initial values, in [0, 2*modn) */
    val0 = poly[0] + poly[4];
    val1 = poly[1] + poly[5];
    val2 = poly[2] + poly[6];
    val3 = poly[3] + poly[7];

    slong m = 8;
    for ( ; m+7 < len; m += 8)
    {
        /* add new values */
        val0 += poly[m+0] + poly[m+4];
        val1 += poly[m+1] + poly[m+5];
        val2 += poly[m+2] + poly[m+6];
        val3 += poly[m+3] + poly[m+7];

        /* reduce to [0, 2*modn) */
        if (val0 >= 2*modn)
            val0 -= 2*modn;
        if (val1 >= 2*modn)
            val1 -= 2*modn;
        if (val2 >= 2*modn)
            val2 -= 2*modn;
        if (val3 >= 2*modn)
            val3 -= 2*modn;
    }

    /* handle 4 more terms if they exist */
    if (m+3 < len)
    {
        val0 += poly[m+0];
        val1 += poly[m+1];
        val2 += poly[m+2];
        val3 += poly[m+3];
        m += 4;
    }

    /* handle 2 more terms if they exist */
    if (m+1 < len)
    {
        val0 += poly[m+0];
        val1 += poly[m+1];
        m += 2;
    }

    /* remains 0 or 1 term */
    if (m < len)
        val2 += poly[m+0];

    /* reduce to [0, 2*modn) */
    if (val0 >= 2*modn)
        val0 -= 2*modn;
    if (val1 >= 2*modn)
        val1 -= 2*modn;
    if (val2 >= 2*modn)
        val2 -= 2*modn;
    if (val3 >= 2*modn)
        val3 -= 2*modn;

    /* gather results */
    val0 = val0 + val1;
    if (val0 >= 2*modn)
        val0 -= 2*modn;
    val2 = val2 + val3;
    if (val2 >= 2*modn)
        val2 -= 2*modn;
    val0 = val0 + val2;
    if (val0 >= 2*modn)
        val0 -= 2*modn;
    if (val0 >= modn)
        val0 -= modn;

    return val0;
}

/* general */
static ulong _nmod_poly_evaluate_mone(nn_srcptr poly, slong len, ulong modn)
{
    if (len == 0)
        return 0;

    if (len == 1)
        return poly[0];

    if (len == 2)
        return n_submod(poly[0], poly[1], modn);

    if (len == 3)
        return n_addmod(n_submod(poly[0], poly[1], modn),
                        poly[2], modn);

    ulong val0, val1, val2, val3;

    val0 = poly[0];
    val1 = poly[1];
    val2 = poly[2];
    val3 = poly[3];

    slong m = 4;

    for ( ; m+3 < len; m += 4)
    {
        val0 = n_addmod(val0, poly[m+0], modn);
        val1 = n_addmod(val1, poly[m+1], modn);
        val2 = n_addmod(val2, poly[m+2], modn);
        val3 = n_addmod(val3, poly[m+3], modn);
    }

    /* handle 2 more terms if they exist */
    if (m+1 < len)
    {
        val0 = n_addmod(val0, poly[m+0], modn);
        val1 = n_addmod(val1, poly[m+1], modn);
        m += 2;
    }

    /* remains 0 or 1 term */
    if (m < len)
        val0 = n_addmod(val0, poly[m+0], modn);

    /* gather results */
    val2 = n_submod(val2, val3, modn);
    val0 = n_submod(val0, val1, modn);
    val0 = n_addmod(val0, val2, modn);

    return val0;
}

/* requires nbits(modn) <= FLINT_BITS - 1 */
static ulong _nmod_poly_evaluate_mone1(nn_srcptr poly, slong len, ulong modn)
{
    if (len == 0)
        return 0;

    if (len == 1)
        return poly[0];

    if (len == 2)
        return n_submod(poly[0], poly[1], modn);

    if (len == 3)
        return n_addmod(n_submod(poly[0], poly[1], modn),
                        poly[2], modn);

    ulong val0, val1, val2, val3;

    val0 = poly[0];
    val1 = poly[1];
    val2 = poly[2];
    val3 = poly[3];

    slong m = 4;

    for ( ; m+3 < len; m += 4)
    {
        val0 = val0 + poly[m+0];
        val1 = val1 + poly[m+1];
        val2 = val2 + poly[m+2];
        val3 = val3 + poly[m+3];
        if (val0 >= modn)
            val0 -= modn;
        if (val1 >= modn)
            val1 -= modn;
        if (val2 >= modn)
            val2 -= modn;
        if (val3 >= modn)
            val3 -= modn;
    }

    /* handle 2 more terms if they exist */
    if (m+1 < len)
    {
        val0 = val0 + poly[m+0];
        val1 = val1 + poly[m+1];
        if (val0 >= modn)
            val0 -= modn;
        if (val1 >= modn)
            val1 -= modn;
        m += 2;
    }

    /* remains 0 or 1 term */
    if (m < len)
    {
        val0 = val0 + poly[m+0];
        if (val0 >= modn)
            val0 -= modn;
    }

    /* gather results */
    val2 = n_submod(val2, val3, modn);
    val0 = n_submod(val0, val1, modn);
    val0 = n_addmod(val0, val2, modn);

    return val0;
}

/* requires nbits(modn) <= FLINT_BITS - 2 */
static ulong _nmod_poly_evaluate_mone2(nn_srcptr poly, slong len, ulong modn)
{
    if (len == 0)
        return 0;

    if (len == 1)
        return poly[0];

    if (len == 2)
        return n_submod(poly[0], poly[1], modn);

    if (len == 3)
        return n_addmod(n_submod(poly[0], poly[1], modn),
                        poly[2], modn);

    if (len <= 7)
    {
        ulong val0, val1;

        val0 = poly[0] + poly[2];
        val1 = poly[1] + poly[3];

        /* handle 2 more terms if they exist */
        slong m = 4;
        if (m+1 < len)
        {
            val0 += poly[m+0];
            val1 += poly[m+1];
            m += 2;
        }

        if (m < len)
            val0 += poly[m];

        /* remains 0 or 1 term */

        if (val0 >= 2*modn)
            val0 -= 2*modn;
        if (val1 >= 2*modn)
            val1 -= 2*modn;
        val0 = val0 + 2*modn - val1;
        if (val0 >= 2*modn)
            val0 -= 2*modn;
        if (val0 >= modn)
            val0 -= modn;

        return val0;
    }

    ulong val0, val1, val2, val3;

    /* initial values, in [0, 2*modn) */
    val0 = poly[0] + poly[4];
    val1 = poly[1] + poly[5];
    val2 = poly[2] + poly[6];
    val3 = poly[3] + poly[7];

    slong m = 8;
    for ( ; m+7 < len; m += 8)
    {
        /* add new values */
        val0 += poly[m+0] + poly[m+4];
        val1 += poly[m+1] + poly[m+5];
        val2 += poly[m+2] + poly[m+6];
        val3 += poly[m+3] + poly[m+7];

        /* reduce to [0, 2*modn) */
        if (val0 >= 2*modn)
            val0 -= 2*modn;
        if (val1 >= 2*modn)
            val1 -= 2*modn;
        if (val2 >= 2*modn)
            val2 -= 2*modn;
        if (val3 >= 2*modn)
            val3 -= 2*modn;
    }

    /* handle 4 more terms if they exist */
    if (m+3 < len)
    {
        val0 += poly[m+0];
        val1 += poly[m+1];
        val2 += poly[m+2];
        val3 += poly[m+3];
        m += 4;
    }

    /* handle 2 more terms if they exist */
    if (m+1 < len)
    {
        val0 += poly[m+0];
        val1 += poly[m+1];
        m += 2;
    }

    /* remains 0 or 1 term */
    if (m < len)
        val2 += poly[m+0];

    /* reduce to [0, 2*modn) */
    if (val0 >= 2*modn)
        val0 -= 2*modn;
    if (val1 >= 2*modn)
        val1 -= 2*modn;
    if (val2 >= 2*modn)
        val2 -= 2*modn;
    if (val3 >= 2*modn)
        val3 -= 2*modn;

    /* gather results */
    val0 = val0 + val2;
    if (val0 >= 2*modn)
        val0 -= 2*modn;
    val1 = val1 + val3;
    if (val1 >= 2*modn)
        val1 -= 2*modn;
    val0 = val0 + 2*modn - val1;
    if (val0 >= 2*modn)
        val0 -= 2*modn;
    if (val0 >= modn)
        val0 -= modn;

    return val0;
}

/*-----------*/
/* general c */
/*-----------*/

ulong _nmod_poly_evaluate_nmod(nn_srcptr poly, slong len, ulong c, nmod_t mod)
{
    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    if (len <= 11)
    {
        slong m = len - 1;
        ulong val = poly[m];
        m -= 1;
        for ( ; m >= 0; m--)
        {
            val = nmod_mul(val, c, mod);
            val = n_addmod(val, poly[m], mod.n);
        }
        return val;
    }

    ulong c4, val0, val1, val2, val3;
    slong m;

    c4 = nmod_mul(c, c, mod);
    c4 = nmod_mul(c4, c4, mod);

    m = len - 1;
    val0 = poly[m-3];
    val1 = poly[m-2];
    val2 = poly[m-1];
    val3 = poly[m-0];
    m -= 4;

    for ( ; m-3 >= 0; m -= 4)
    {
        val0 = nmod_mul(val0, c4, mod);
        val0 = n_addmod(val0, poly[m-3], mod.n);
        val1 = nmod_mul(val1, c4, mod);
        val1 = n_addmod(val1, poly[m-2], mod.n);
        val2 = nmod_mul(val2, c4, mod);
        val2 = n_addmod(val2, poly[m-1], mod.n);
        val3 = nmod_mul(val3, c4, mod);
        val3 = n_addmod(val3, poly[m-0], mod.n);
    }

    /* gather results with Horner */
    val2 = n_addmod(val2, nmod_mul(val3, c, mod), mod.n);
    val1 = n_addmod(val1, nmod_mul(val2, c, mod), mod.n);
    val0 = n_addmod(val0, nmod_mul(val1, c, mod), mod.n);

    /* last few terms */
    for ( ; m >= 0; m--)
    {
        val0 = nmod_mul(val0, c, mod);
        val0 = n_addmod(val0, poly[m], mod.n);
    }

    return val0;
}

ulong _nmod_poly_evaluate_nmod_precomp(nn_srcptr poly, slong len, ulong c, ulong c_precomp, ulong modn)
{
    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    if (len <= 11)
    {
        slong m = len - 1;
        ulong val = poly[m];
        m -= 1;
        for ( ; m >= 0; m--)
        {
            val = n_mulmod_shoup(c, val, c_precomp, modn);
            val = n_addmod(val, poly[m], modn);
        }
        return val;
    }

    ulong val0, val1, val2, val3;
    slong m;

    /* precomputations for c**4 */
    ulong c4, c4_precomp;
    val0 = n_mulmod_precomp_shoup_rem_from_quo(c_precomp, modn);
    n_mulmod_and_precomp_shoup(&c4, &c4_precomp, c, c, c_precomp, val0, c_precomp, modn);
    val0 = n_mulmod_precomp_shoup_rem_from_quo(c4_precomp, modn);
    n_mulmod_and_precomp_shoup(&c4, &c4_precomp, c4, c4, c4_precomp, val0, c4_precomp, modn);

    m = len - 1;
    val0 = poly[m-3];
    val1 = poly[m-2];
    val2 = poly[m-1];
    val3 = poly[m-0];
    m -= 4;

    for ( ; m-3 >= 0; m -= 4)
    {
        val0 = n_mulmod_shoup(c4, val0, c4_precomp, modn);
        val0 = n_addmod(val0, poly[m-3], modn);
        val1 = n_mulmod_shoup(c4, val1, c4_precomp, modn);
        val1 = n_addmod(val1, poly[m-2], modn);
        val2 = n_mulmod_shoup(c4, val2, c4_precomp, modn);
        val2 = n_addmod(val2, poly[m-1], modn);
        val3 = n_mulmod_shoup(c4, val3, c4_precomp, modn);
        val3 = n_addmod(val3, poly[m-0], modn);
    }

    /* gather results with Horner */
    val3 = n_mulmod_shoup(c, val3, c_precomp, modn);
    val2 = n_addmod(val2, val3, modn);
    val2 = n_mulmod_shoup(c, val2, c_precomp, modn);
    val1 = n_addmod(val1, val2, modn);
    val1 = n_mulmod_shoup(c, val1, c_precomp, modn);
    val0 = n_addmod(val0, val1, modn);

    /* last few terms */
    for ( ; m >= 0; m--)
    {
        val0 = n_mulmod_shoup(c, val0, c_precomp, modn);
        val0 = n_addmod(val0, poly[m], modn);
    }

    return val0;
}

ulong _nmod_poly_evaluate_nmod_precomp_lazy(nn_srcptr poly, slong len, ulong c, ulong c_precomp, ulong modn)
{
    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    if (len <= 15)
    {
        ulong p_hi, p_lo;
        slong m = len - 1;
        ulong val = poly[m];
        m -= 1;
        for ( ; m >= 0; m--)
        {
            umul_ppmm(p_hi, p_lo, c_precomp, val);
            val = poly[m] + c * val - p_hi * modn;
        }
        return val;
    }

    ulong p_hi, p_lo;
    ulong val0, val1, val2, val3;
    slong m;

    /* precomputations for c**4 */
    ulong c4, c4_precomp;
    val0 = n_mulmod_precomp_shoup_rem_from_quo(c_precomp, modn);
    n_mulmod_and_precomp_shoup(&c4, &c4_precomp, c, c, c_precomp, val0, c_precomp, modn);
    val0 = n_mulmod_precomp_shoup_rem_from_quo(c4_precomp, modn);
    n_mulmod_and_precomp_shoup(&c4, &c4_precomp, c4, c4, c4_precomp, val0, c4_precomp, modn);

    m = len - 1;
    val0 = poly[m-3];
    val1 = poly[m-2];
    val2 = poly[m-1];
    val3 = poly[m-0];
    m -= 4;

    for ( ; m-3 >= 0; m -= 4)
    {
        umul_ppmm(p_hi, p_lo, c4_precomp, val0);
        val0 = poly[m-3] + c4 * val0 - p_hi * modn;
        umul_ppmm(p_hi, p_lo, c4_precomp, val1);
        val1 = poly[m-2] + c4 * val1 - p_hi * modn;
        umul_ppmm(p_hi, p_lo, c4_precomp, val2);
        val2 = poly[m-1] + c4 * val2 - p_hi * modn;
        umul_ppmm(p_hi, p_lo, c4_precomp, val3);
        val3 = poly[m-0] + c4 * val3 - p_hi * modn;
    }

    /* each val is in [0, max(poly) + 2*modn) */
    /*  --> reduce to [0, max(poly)] */
    if (val2 >= 2*modn)
        val2 -= 2*modn;
    else if (val2 >= modn)
        val2 -= modn;
    if (val1 >= 2*modn)
        val1 -= 2*modn;
    else if (val1 >= modn)
        val1 -= modn;
    if (val0 >= 2*modn)
        val0 -= 2*modn;
    else if (val0 >= modn)
        val0 -= modn;
    /* gather results with Horner */
    umul_ppmm(p_hi, p_lo, c_precomp, val3);
    val2 = val2 + c * val3 - p_hi * modn;
    umul_ppmm(p_hi, p_lo, c_precomp, val2);
    val1 = val1 + c * val2 - p_hi * modn;
    umul_ppmm(p_hi, p_lo, c_precomp, val1);
    val0 = val0 + c * val1 - p_hi * modn;

    /* last few terms */
    for ( ; m >= 0; m--)
    {
        umul_ppmm(p_hi, p_lo, c_precomp, val0);
        val0 = poly[m] + c * val0 - p_hi * modn;
    }

    return val0;
}

/*----------------*/
/* main interface */
/*----------------*/

ulong nmod_poly_evaluate_nmod(const nmod_poly_t poly, ulong c)
{
    if (poly->length == 0)
        return 0;

    if (poly->length == 1 || c == 0)
        return poly->coeffs[0];

    if (c == 1)
    {
        if (poly->mod.norm >= 2)
            return _nmod_poly_evaluate_one2(poly->coeffs, poly->length, poly->mod.n);
        else if (poly->mod.norm == 1)
            return _nmod_poly_evaluate_one1(poly->coeffs, poly->length, poly->mod.n);
        else
            return _nmod_poly_evaluate_one(poly->coeffs, poly->length, poly->mod.n);
    }

    if (c == poly->mod.n - 1)
    {
        if (poly->mod.norm >= 2)
            return _nmod_poly_evaluate_mone2(poly->coeffs, poly->length, poly->mod.n);
        else if (poly->mod.norm == 1)
            return _nmod_poly_evaluate_mone1(poly->coeffs, poly->length, poly->mod.n);
        else
            return _nmod_poly_evaluate_mone(poly->coeffs, poly->length, poly->mod.n);
    }

    // if degree below the n_mulmod_shoup threshold
    // or modulus forbids n_mulmod_shoup usage, use nmod_mul
#if FLINT_MULMOD_SHOUP_THRESHOLD <= 2
    if (poly->mod.norm == 0)  // here poly->length >= threshold
#else
    if ((poly->length < FLINT_MULMOD_SHOUP_THRESHOLD)
           || (poly->mod.norm == 0))
#endif
    {
        return _nmod_poly_evaluate_nmod(poly->coeffs, poly->length, c, poly->mod);
    }

    // if 3*mod.n - 1 <= 2**FLINT_BITS, use n_mulmod_shoup, lazy variant
    else
    {
        const ulong modn = poly->mod.n;

#if FLINT_BITS == 64
        if (modn <= UWORD(6148914691236517205))
#else // FLINT_BITS == 32
        if (modn <= UWORD(1431655765))
#endif
        {
            const ulong c_precomp = n_mulmod_precomp_shoup(c, modn);
            ulong val = _nmod_poly_evaluate_nmod_precomp_lazy(poly->coeffs, poly->length, c, c_precomp, modn);
            // correct excess
            if (val >= 2*modn)
                val -= 2*modn;
            else if (val >= modn)
                val -= modn;
            return val;
        }

        // use n_mulmod_shoup, non-lazy variant
        else
        {
            const ulong c_precomp = n_mulmod_precomp_shoup(c, modn);
            return _nmod_poly_evaluate_nmod_precomp(poly->coeffs, poly->length, c, c_precomp, modn);
        }
    }
}
