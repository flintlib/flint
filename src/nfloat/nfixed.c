/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "mpn_extras.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_generic.h"
#include "acf.h"
#include "acb.h"
#include "nfloat.h"
#include "gr_special.h"
#include "fmpz_mat.h"

/* For printing */
#include "arf.h"

/* based on output of p-nfixed_mat_mul.c */
static int nfixed_mat_mul_use_waksman(slong n, slong nlimbs)
{
    if (nlimbs <= 8) return 0;
    if (nlimbs <= 11) return (n >= 6);
    if (nlimbs == 12) return (n >= 5);
    if (nlimbs <= 21) return (n >= 4);
    if (nlimbs <= 46) return (n >= 3);
    return (n >= 2);
}

/* based on output of p-nfixed_mat_mul.c */
static slong nfixed_mat_mul_strassen_cutoff(slong n, int parity, slong nlimbs)
{
    if (nlimbs <= 3)
        return parity ? 57 : 50;
    else
        return parity ? 37 : 26;
}

/* Arithmetic on fixed-point numbers in (-1,1) */
/* x[0] stores the sign bit, x[1], ..., x[n] store the absolute value */

void
_nfixed_print(nn_srcptr x, slong nlimbs, slong exp)
{
    arf_t t;
    arf_init(t);
    _arf_set_mpn_fixed(t, x + 1, nlimbs, nlimbs, x[0], nlimbs * FLINT_BITS, ARF_RND_DOWN);
    arf_mul_2exp_si(t, t, exp);
    arf_printd(t, nlimbs * FLINT_BITS / 3.321928 + 1);
    arf_clear(t);
}

#define DEF_NFIXED_ADD(n) \
FLINT_FORCE_INLINE \
void _nfixed_add_ ## n(nn_ptr res, nn_srcptr a, nn_srcptr b) \
{ \
    int asgn, bsgn; \
    asgn = a[0]; \
    bsgn = b[0]; \
 \
    if (asgn == bsgn) \
    { \
        res[0] = asgn; \
        NN_ADD_ ## n(res + 1, a + 1, b + 1); \
    } \
    else \
    { \
        res[0] = asgn ^ flint_mpn_signed_sub_ ## n(res + 1, a + 1, b + 1); \
    } \
}

#define DEF_NFIXED_SUB(n) \
FLINT_FORCE_INLINE \
void _nfixed_sub_ ## n(nn_ptr res, nn_srcptr a, nn_srcptr b) \
{ \
    int asgn, bsgn; \
    asgn = a[0]; \
    bsgn = b[0]; \
 \
    if (asgn != bsgn) \
    { \
        res[0] = asgn; \
        NN_ADD_ ## n(res + 1, a + 1, b + 1); \
    } \
    else \
    { \
        res[0] = asgn ^ flint_mpn_signed_sub_ ## n(res + 1, a + 1, b + 1); \
    } \
}

DEF_NFIXED_ADD(2)
DEF_NFIXED_ADD(3)
DEF_NFIXED_ADD(4)
DEF_NFIXED_ADD(5)
DEF_NFIXED_ADD(6)
DEF_NFIXED_ADD(7)
DEF_NFIXED_ADD(8)

DEF_NFIXED_SUB(2)
DEF_NFIXED_SUB(3)
DEF_NFIXED_SUB(4)
DEF_NFIXED_SUB(5)
DEF_NFIXED_SUB(6)
DEF_NFIXED_SUB(7)
DEF_NFIXED_SUB(8)

FLINT_FORCE_INLINE
void _nfixed_add(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];

    if (asgn == bsgn)
    {
        res[0] = asgn;
        mpn_add_n(res + 1, a + 1, b + 1, nlimbs);
    }
    else
    {
        res[0] = asgn ^ flint_mpn_signed_sub_n(res + 1, a + 1, b + 1, nlimbs);
    }
}

FLINT_FORCE_INLINE
void _nfixed_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];

    if (asgn != bsgn)
    {
        res[0] = asgn;
        mpn_add_n(res + 1, a + 1, b + 1, nlimbs);
    }
    else
    {
        res[0] = asgn ^ flint_mpn_signed_sub_n(res + 1, a + 1, b + 1, nlimbs);
    }
}

void _nfixed_vec_add(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs)
{
    slong i;

    if (nlimbs == 2)
    {
        for (i = 0; i < len; i++)
            _nfixed_add_2(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 3)
    {
        for (i = 0; i < len; i++)
            _nfixed_add_3(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 4)
    {
        for (i = 0; i < len; i++)
            _nfixed_add_4(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 5)
    {
        for (i = 0; i < len; i++)
            _nfixed_add_5(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 6)
    {
        for (i = 0; i < len; i++)
            _nfixed_add_6(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 7)
    {
        for (i = 0; i < len; i++)
            _nfixed_add_7(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 8)
    {
        for (i = 0; i < len; i++)
            _nfixed_add_8(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else
    {
        for (i = 0; i < len; i++)
            _nfixed_add(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1), nlimbs);
    }
}

void _nfixed_vec_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs)
{
    slong i;

    if (nlimbs == 2)
    {
        for (i = 0; i < len; i++)
            _nfixed_sub_2(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 3)
    {
        for (i = 0; i < len; i++)
            _nfixed_sub_3(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 4)
    {
        for (i = 0; i < len; i++)
            _nfixed_sub_4(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 5)
    {
        for (i = 0; i < len; i++)
            _nfixed_sub_5(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 6)
    {
        for (i = 0; i < len; i++)
            _nfixed_sub_6(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 7)
    {
        for (i = 0; i < len; i++)
            _nfixed_sub_7(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else if (nlimbs == 8)
    {
        for (i = 0; i < len; i++)
            _nfixed_sub_8(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1));
    }
    else
    {
        for (i = 0; i < len; i++)
            _nfixed_sub(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1), nlimbs);
    }
}

FLINT_FORCE_INLINE
void nfixed_mul(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];
    res[0] = asgn ^ bsgn;
    flint_mpn_mulhigh_n(res + 1, a + 1, b + 1, nlimbs);
}

#if 0
FLINT_FORCE_INLINE
void nfixed_sqr(nn_ptr res, nn_srcptr a, slong nlimbs)
{
    res[0] = 0;
    flint_mpn_sqrhigh(res + 1, a + 1, nlimbs);
}
#endif

FLINT_FORCE_INLINE
void nfixed_div2(nn_ptr res, nn_srcptr a, slong nlimbs)
{
    res[0] = a[0];
    mpn_rshift(res + 1, a + 1, nlimbs, 1);
}

void _nfixed_dot_2(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
{
    slong j;

    ulong as, bs, a0, a1, b0, b1, s0, s1, t0, t1, u0, u1, hi, lo;

    s0 = s1 = t0 = t1 = u0 = u1 = 0;

        /*
                    s0    s1
                  |a1 b1----|
                    u0    u1
            |a1 b0----|
                    t0    t1
            |a0 b1----|
        */

    for (j = 0; j < len; j++)
    {
        as = x[j * xstride];
        a0 = x[j * xstride + 1];
        a1 = x[j * xstride + 2];
        bs = y[j * ystride];
        b0 = y[j * ystride + 1];
        b1 = y[j * ystride + 2];

        if (as == bs)
        {
            umul_ppmm(hi, lo, a1, b1);
            add_ssaaaa(s1, s0, s1, s0, hi, lo);
            umul_ppmm(hi, lo, a0, b1);
            add_ssaaaa(t1, t0, t1, t0, 0, hi);
            umul_ppmm(hi, lo, a1, b0);
            add_ssaaaa(u1, u0, u1, u0, 0, hi);
        }
        else
        {
            umul_ppmm(hi, lo, a1, b1);
            sub_ddmmss(s1, s0, s1, s0, hi, lo);
            umul_ppmm(hi, lo, a0, b1);
            sub_ddmmss(t1, t0, t1, t0, 0, hi);
            umul_ppmm(hi, lo, a1, b0);
            sub_ddmmss(u1, u0, u1, u0, 0, hi);
        }
    }

    add_ssaaaa(s1, s0, s1, s0, t1, t0);
    add_ssaaaa(s1, s0, s1, s0, u1, u0);

    if ((slong) s1 < WORD(0))
    {
        sub_ddmmss(s1, s0, 0, 0, s1, s0);
        res[0] = 1;
    }
    else
    {
        res[0] = 0;
    }

    res[1] = s0;
    res[2] = s1;
}

void _nfixed_dot_3(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
{
    slong j;

    ulong as, bs, a0, a1, a2, b0, b1, b2, s0, s1, s2, hi, lo;
    ulong u0, u1, v0, v1;

    s0 = s1 = s2 = 0;
    u0 = u1 = v0 = v1 = 0;

    for (j = 0; j < len; j++)
    {
        as = x[j * xstride];
        a0 = x[j * xstride + 1];
        a1 = x[j * xstride + 2];
        a2 = x[j * xstride + 3];
        bs = y[j * ystride];
        b0 = y[j * ystride + 1];
        b1 = y[j * ystride + 2];
        b2 = y[j * ystride + 3];

        /*
                      |a2 b2----|

                |a1 b2----|
                |a2 b1----|

           |a0 b2----|
           |a1 b1----|
           |a2 b0----|
        */

        if (as == bs)
        {
            umul_ppmm(hi, lo, a0, b2);
            add_ssaaaa(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a1, b1);
            add_ssaaaa(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a2, b0);
            add_ssaaaa(u1, u0, u1, u0, 0, hi);

            umul_ppmm(hi, lo, a2, b2);
            add_ssaaaa(s2, s1, s2, s1, hi, lo);

            umul_ppmm(hi, lo, a1, b2);
            add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
            umul_ppmm(hi, lo, a2, b1);
            add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, hi, lo);
        }
        else
        {
            umul_ppmm(hi, lo, a0, b2);
            sub_ddmmss(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a1, b1);
            sub_ddmmss(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a2, b0);
            sub_ddmmss(u1, u0, u1, u0, 0, hi);

            umul_ppmm(hi, lo, a2, b2);
            sub_ddmmss(s2, s1, s2, s1, hi, lo);

            umul_ppmm(hi, lo, a1, b2);
            sub_dddmmmsss(s2, s1, s0, s2, s1, s0, 0, hi, lo);
            umul_ppmm(hi, lo, a2, b1);
            sub_dddmmmsss(s2, s1, s0, s2, s1, s0, 0, hi, lo);
        }
    }

    if ((slong) u1 < WORD(0))
    {
        sub_ddmmss(u1, u0, 0, 0, u1, u0);
        sub_dddmmmsss(s2, s1, s0, s2, s1, s0, 0, u1, u0);
    }
    else
    {
        add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, u1, u0);
    }

    if ((slong) s2 < WORD(0))
    {
        sub_dddmmmsss(s2, s1, s0, 0, 0, 0, s2, s1, s0);
        res[0] = 1;
    }
    else
    {
        res[0] = 0;
    }


    res[1] = s0;
    res[2] = s1;
    res[3] = s2;
}

void _nfixed_dot_4(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
{
    slong j;

    /*
                 s1     s2    s3
                      |a3 b3----|

                |a2 b3----|
                |a3 b2----|

             t0   t1   t2

           |a1 b3----|
           |a2 b2----|
           |a3 b1----|

             u0   u1

      |a0 b3----|
      |a1 b2----|
      |a2 b1----|
      |a3 b0----|

    */

    ulong as, a0, a1, a2, a3;
    ulong bs, b0, b1, b2, b3;
    ulong s0, s1, s2, s3, t0, t1, t2, u0, u1;
    ulong hi, lo;

    s0 = s1 = s2 = s3 = t0 = t1 = t2 = u0 = u1 = 0;

    for (j = 0; j < len; j++)
    {
        as = x[j * xstride];
        a0 = x[j * xstride + 1];
        a1 = x[j * xstride + 2];
        a2 = x[j * xstride + 3];
        a3 = x[j * xstride + 4];

        bs = y[j * ystride];
        b0 = y[j * ystride + 1];
        b1 = y[j * ystride + 2];
        b2 = y[j * ystride + 3];
        b3 = y[j * ystride + 4];

        if (as == bs)
        {
            umul_ppmm(hi, lo, a3, b3);
            add_ssaaaa(s3, s2, s3, s2, hi, lo);

            umul_ppmm(hi, lo, a2, b3);
            add_sssaaaaaa(s3, s2, s1, s3, s2, s1, 0, hi, lo);
            umul_ppmm(hi, lo, a3, b2);
            add_sssaaaaaa(s3, s2, s1, s3, s2, s1, 0, hi, lo);

            umul_ppmm(hi, lo, a1, b3);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, hi, lo);
            umul_ppmm(hi, lo, a2, b2);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, hi, lo);
            umul_ppmm(hi, lo, a3, b1);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, hi, lo);

            umul_ppmm(hi, lo, a0, b3);
            add_ssaaaa(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a1, b2);
            add_ssaaaa(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a2, b1);
            add_ssaaaa(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a3, b0);
            add_ssaaaa(u1, u0, u1, u0, 0, hi);
        }
        else
        {
            umul_ppmm(hi, lo, a3, b3);
            sub_ddmmss(s3, s2, s3, s2, hi, lo);

            umul_ppmm(hi, lo, a2, b3);
            sub_dddmmmsss(s3, s2, s1, s3, s2, s1, 0, hi, lo);
            umul_ppmm(hi, lo, a3, b2);
            sub_dddmmmsss(s3, s2, s1, s3, s2, s1, 0, hi, lo);

            umul_ppmm(hi, lo, a1, b3);
            sub_dddmmmsss(t2, t1, t0, t2, t1, t0, 0, hi, lo);
            umul_ppmm(hi, lo, a2, b2);
            sub_dddmmmsss(t2, t1, t0, t2, t1, t0, 0, hi, lo);
            umul_ppmm(hi, lo, a3, b1);
            sub_dddmmmsss(t2, t1, t0, t2, t1, t0, 0, hi, lo);

            umul_ppmm(hi, lo, a0, b3);
            sub_ddmmss(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a1, b2);
            sub_ddmmss(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a2, b1);
            sub_ddmmss(u1, u0, u1, u0, 0, hi);
            umul_ppmm(hi, lo, a3, b0);
            sub_ddmmss(u1, u0, u1, u0, 0, hi);
        }
    }

    if ((slong) u1 < WORD(0))
    {
        sub_ddmmss(u1, u0, 0, 0, u1, u0);
        sub_dddmmmsss(t2, t1, s0, t2, t1, s0, 0, u1, u0);
    }
    else
    {
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, u1, u0);
    }

    if ((slong) t2 < WORD(0))
    {
        sub_dddmmmsss(t2, t1, t0, 0, 0, 0, t2, t1, t0);
        sub_ddddmmmmssss(s3, s2, s1, s0, s3, s2, s1, s0, 0, t2, t1, t0);
    }
    else
    {
        add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, 0, t2, t1, t0);
    }

    if ((slong) s3 < WORD(0))
    {
        sub_ddddmmmmssss(s3, s2, s1, s0, 0, 0, 0, 0, s3, s2, s1, s0);
        res[0] = 1;
    }
    else
    {
        res[0] = 0;
    }

    res[1] = s0;
    res[2] = s1;
    res[3] = s2;
    res[4] = s3;
}

void _nfixed_dot_5(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
{
    slong j;
    slong nlimbs = 5;

    ulong tmp[6];
    ulong spos[6] = { 0, 0, 0, 0, 0, 0 };
    ulong sneg[6] = { 0, 0, 0, 0, 0, 0 };

    if (x[0] == y[0])
        flint_mpn_mulhigh_n(spos + 1, x + 1, y + 1, nlimbs);
    else
        flint_mpn_mulhigh_n(sneg + 1, x + 1, y + 1, nlimbs);

    for (j = 1; j < len; j++)
    {
        nfixed_mul(tmp, x + j * xstride, y + j * ystride, nlimbs);

        if (tmp[0] == 0)
            NN_ADD_5(spos + 1, spos + 1, tmp + 1);
        else
            NN_ADD_5(sneg + 1, sneg + 1, tmp + 1);
    }

    _nfixed_sub_5(res, spos, sneg);
}

void _nfixed_dot_6(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
{
    slong j;
    slong nlimbs = 6;

    ulong tmp[7];
    ulong spos[7] = { 0, 0, 0, 0, 0, 0, 0 };
    ulong sneg[7] = { 0, 0, 0, 0, 0, 0, 0 };

    if (x[0] == y[0])
        flint_mpn_mulhigh_n(spos + 1, x + 1, y + 1, nlimbs);
    else
        flint_mpn_mulhigh_n(sneg + 1, x + 1, y + 1, nlimbs);

    for (j = 1; j < len; j++)
    {
        nfixed_mul(tmp, x + j * xstride, y + j * ystride, nlimbs);

        if (tmp[0] == 0)
            NN_ADD_6(spos + 1, spos + 1, tmp + 1);
        else
            NN_ADD_6(sneg + 1, sneg + 1, tmp + 1);
    }

    _nfixed_sub_6(res, spos, sneg);
}

void _nfixed_dot_7(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
{
    slong j;
    slong nlimbs = 7;

    ulong tmp[8];
    ulong spos[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    ulong sneg[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

    if (x[0] == y[0])
        flint_mpn_mulhigh_n(spos + 1, x + 1, y + 1, nlimbs);
    else
        flint_mpn_mulhigh_n(sneg + 1, x + 1, y + 1, nlimbs);

    for (j = 1; j < len; j++)
    {
        nfixed_mul(tmp, x + j * xstride, y + j * ystride, nlimbs);

        if (tmp[0] == 0)
            NN_ADD_7(spos + 1, spos + 1, tmp + 1);
        else
            NN_ADD_7(sneg + 1, sneg + 1, tmp + 1);
    }

    _nfixed_sub_7(res, spos, sneg);
}

void _nfixed_dot_8(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
{
    slong j;
    slong nlimbs = 8;

    ulong tmp[9];
    ulong spos[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    ulong sneg[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    if (x[0] == y[0])
        flint_mpn_mulhigh_n(spos + 1, x + 1, y + 1, nlimbs);
    else
        flint_mpn_mulhigh_n(sneg + 1, x + 1, y + 1, nlimbs);

    for (j = 1; j < len; j++)
    {
        nfixed_mul(tmp, x + j * xstride, y + j * ystride, nlimbs);

        if (tmp[0] == 0)
            NN_ADD_8(spos + 1, spos + 1, tmp + 1);
        else
            NN_ADD_8(sneg + 1, sneg + 1, tmp + 1);
    }

    _nfixed_sub_8(res, spos, sneg);
}

/* A is (m x n), B is (n x p), C is (m x p) */
void
_nfixed_mat_mul_classical(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
{
    slong i, j, k;
    nn_ptr t;

#define A_ENTRY(i, j) ((A) + ((i) * n + (j)) * (nlimbs + 1))
#define B_ENTRY(i, j) ((B) + ((i) * p + (j)) * (nlimbs + 1))
#define C_ENTRY(i, j) ((C) + ((i) * p + (j)) * (nlimbs + 1))

    if (n == 1)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                nfixed_mul(C_ENTRY(i, j), A_ENTRY(i, 0), B_ENTRY(0, j), nlimbs);
        return;
    }

    if (nlimbs == 2)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_2(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), (nlimbs + 1) * p, n);
    }
    else if (nlimbs == 3)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_3(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), (nlimbs + 1) * p, n);
    }
    else if (nlimbs == 4)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_4(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), (nlimbs + 1) * p, n);
    }
    else if (nlimbs == 5)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_5(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), (nlimbs + 1) * p, n);
    }
    else if (nlimbs == 6)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_6(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), (nlimbs + 1) * p, n);
    }
    else if (nlimbs == 7)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_7(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), (nlimbs + 1) * p, n);
    }
    else if (nlimbs == 8)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_8(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), (nlimbs + 1) * p, n);
    }
    else
    {
        TMP_INIT;
        TMP_START;

        t = TMP_ALLOC((nlimbs + 1) * sizeof(ulong));

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < p; j++)
            {
                nfixed_mul(C_ENTRY(i, j), A_ENTRY(i, 0), B_ENTRY(0, j), nlimbs);

                for (k = 1; k < n; k++)
                {
                    nfixed_mul(t, A_ENTRY(i, k), B_ENTRY(k, j), nlimbs);
                    _nfixed_add(C_ENTRY(i, j), C_ENTRY(i, j), t, nlimbs);
                }
            }
        }

        TMP_END;
    }

#undef A_ENTRY
#undef B_ENTRY
#undef C_ENTRY
}

/* todo: optimize */
/* A is (m x n), B is (n x p), C is (m x p) */
void
_nfixed_mat_mul_classical_precise(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
{
    slong i;
    nn_ptr t, tA, tB, tC, u;
    nn_srcptr s;

    t = flint_malloc(((m * n) + (n * p) + (m * p)) * (nlimbs + 2) * sizeof(ulong));
    tA = t;
    tB = tA + (m * n) * (nlimbs + 2);
    tC = tB + (n * p) * (nlimbs + 2);

    for (i = 0; i < m * n; i++)
    {
        s = A + i * (nlimbs + 1);
        u = tA + i* (nlimbs + 2);
        flint_mpn_copyi(u + 2, s + 1, nlimbs);
        u[0] = s[0];
        u[1] = 0;
    }

    for (i = 0; i < n * p; i++)
    {
        s = B + i * (nlimbs + 1);
        u = tB + i * (nlimbs + 2);
        flint_mpn_copyi(u + 2, s + 1, nlimbs);
        u[0] = s[0];
        u[1] = 0;
    }

    _nfixed_mat_mul_classical(tC, tA, tB, m, n, p, nlimbs + 1);

    for (i = 0; i < m * p; i++)
    {
        s = tC + i * (nlimbs + 2);
        u = C + i * (nlimbs + 1);
        flint_mpn_copyi(u + 1, s + 2, nlimbs);
        u[0] = s[0];
    }

    flint_free(t);
}

/* compute c += (a1 + b1) * (a2 + b2) */
/* val0, val1, val2 are scratch space */
FLINT_FORCE_INLINE void
addmul_addadd(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2, slong nlimbs)
{
    _nfixed_add(val1, a1, b1, nlimbs);
    _nfixed_add(val2, a2, b2, nlimbs);
    nfixed_mul(val0, val1, val2, nlimbs);
    _nfixed_add(c, c, val0, nlimbs);
}

/* compute c += (a1 - b1) * (a2 - b2) */
/* val0, val1, val2 are scratch space */
FLINT_FORCE_INLINE void
addmul_subsub(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2, slong nlimbs)
{
    _nfixed_sub(val1, a1, b1, nlimbs);
    _nfixed_sub(val2, a2, b2, nlimbs);
    nfixed_mul(val0, val1, val2, nlimbs);
    _nfixed_add(c, c, val0, nlimbs);
}

/*
    Inlining speeds up Waksman multiplication with small nlimbs.
    Further speedups are possible by reordering the loops so that
    the O(n^3) part is done using dot products.
    However, these tricks currently do not suffice to beat
    classical multiplications in the relevant ranges, so we
    do not bother here.
*/

#define WAKSMAN_WANT_INLINING 0

#if WAKSMAN_WANT_INLINING

FLINT_FORCE_INLINE void
addmul_addadd_4(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2)
{
    _nfixed_add_4(val1, a1, b1);
    _nfixed_add_4(val2, a2, b2);
    nfixed_mul(val0, val1, val2, 4);
    _nfixed_add_4(c, c, val0);
}

FLINT_FORCE_INLINE void
addmul_addadd_5(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2)
{
    _nfixed_add_5(val1, a1, b1);
    _nfixed_add_5(val2, a2, b2);
    nfixed_mul(val0, val1, val2, 5);
    _nfixed_add_5(c, c, val0);
}

FLINT_FORCE_INLINE void
addmul_addadd_6(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2)
{
    _nfixed_add_6(val1, a1, b1);
    _nfixed_add_6(val2, a2, b2);
    nfixed_mul(val0, val1, val2, 6);
    _nfixed_add_6(c, c, val0);
}

FLINT_FORCE_INLINE void
addmul_addadd_7(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2)
{
    _nfixed_add_7(val1, a1, b1);
    _nfixed_add_7(val2, a2, b2);
    nfixed_mul(val0, val1, val2, 7);
    _nfixed_add_7(c, c, val0);
}

FLINT_FORCE_INLINE void
addmul_addadd_8(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2)
{
    _nfixed_add_8(val1, a1, b1);
    _nfixed_add_8(val2, a2, b2);
    nfixed_mul(val0, val1, val2, 8);
    _nfixed_add_8(c, c, val0);
}

#endif

void
_nfixed_mat_mul_waksman2(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs, slong Cstride, slong Astride, slong Bstride)
{
    slong l, j, k;

    slong np = n >> 1;

    nn_ptr Ctmp = flint_calloc((nlimbs + 1) * ((p + m) + 5), sizeof(ulong));

    /* remaining temp space */
    nn_ptr Crow = Ctmp;                     /* Crow has p entries */
    nn_ptr Ccol = Crow + (nlimbs + 1) * p;  /* Ccol has m entries */
    nn_ptr val0 = Ccol + (nlimbs + 1) * m;  /* val0 has room for 2 sums */
    nn_ptr val1 = val0 + (nlimbs + 1) * 2;  /* val1 has room for 1 sum   */
    nn_ptr val2 = val1 + (nlimbs + 1);      /* val2 has room for 1 sum   */
    nn_ptr crow = val2 + (nlimbs + 1);      /* crow has room for 1 sum   */

#define A_ENTRY(i, j) ((A) + (i) * Astride + (j) * (nlimbs + 1))
#define B_ENTRY(i, j) ((B) + (i) * Bstride + (j) * (nlimbs + 1))
#define C_ENTRY(i, j) ((C) + (i) * Cstride + (j) * (nlimbs + 1))

#define Crow_ENTRY(ii) (Crow + (ii) * (nlimbs + 1))
#define Ccol_ENTRY(ii) (Ccol + (ii) * (nlimbs + 1))

    /* todo: zero only where needed */
    for (j = 0; j < m; j++)
        flint_mpn_zero(C_ENTRY(j, 0), p * (nlimbs + 1));

    for (j = 1; j <= np; j++)
    {
        slong j2 = (j << 1) - 1;

        for (k = 0; k < p; k++)
        {
            addmul_addadd(val0, val1, val2, C_ENTRY(0, k), A_ENTRY(0, j2-1), B_ENTRY(j2, k), A_ENTRY(0, j2), B_ENTRY(j2-1, k), nlimbs);
            addmul_subsub(val0, val1, val2, Crow_ENTRY(k), A_ENTRY(0, j2-1), B_ENTRY(j2, k), A_ENTRY(0, j2), B_ENTRY(j2-1, k), nlimbs);
        }

        for (l = 1; l < m; l++)
        {
            addmul_addadd(val0, val1, val2, C_ENTRY(l, 0), A_ENTRY(l, j2-1), B_ENTRY(j2, 0), A_ENTRY(l, j2), B_ENTRY(j2-1, 0), nlimbs);
            addmul_subsub(val0, val1, val2, Ccol_ENTRY(l), A_ENTRY(l, j2-1), B_ENTRY(j2, 0), A_ENTRY(l, j2), B_ENTRY(j2-1, 0), nlimbs);
        }

#if WAKSMAN_WANT_INLINING
        if (nlimbs == 5)
        {
            for (k = 1; k < p; k++)
                for (l = 1; l < m; l++)
                    addmul_addadd_5(val0, val1, val2, C_ENTRY(l, k), A_ENTRY(l, j2-1), B_ENTRY(j2, k), A_ENTRY(l, j2), B_ENTRY(j2-1, k));
        }
        else if (nlimbs == 6)
        {
            for (k = 1; k < p; k++)
                for (l = 1; l < m; l++)
                    addmul_addadd_6(val0, val1, val2, C_ENTRY(l, k), A_ENTRY(l, j2-1), B_ENTRY(j2, k), A_ENTRY(l, j2), B_ENTRY(j2-1, k));
        }
        else if (nlimbs == 7)
        {
            for (k = 1; k < p; k++)
                for (l = 1; l < m; l++)
                    addmul_addadd_7(val0, val1, val2, C_ENTRY(l, k), A_ENTRY(l, j2-1), B_ENTRY(j2, k), A_ENTRY(l, j2), B_ENTRY(j2-1, k));
        }
        else if (nlimbs == 8)
        {
            for (k = 1; k < p; k++)
                for (l = 1; l < m; l++)
                    addmul_addadd_8(val0, val1, val2, C_ENTRY(l, k), A_ENTRY(l, j2-1), B_ENTRY(j2, k), A_ENTRY(l, j2), B_ENTRY(j2-1, k));
        }
        else
#endif
        {
            for (k = 1; k < p; k++)
                for (l = 1; l < m; l++)
                    addmul_addadd(val0, val1, val2, C_ENTRY(l, k), A_ENTRY(l, j2-1), B_ENTRY(j2, k), A_ENTRY(l, j2), B_ENTRY(j2-1, k), nlimbs);
        }
    }

    for (l = 1; l < m; l++)
    {
        _nfixed_add(val1, Ccol_ENTRY(l), C_ENTRY(l, 0), nlimbs);
        nfixed_div2(Ccol_ENTRY(l), val1, nlimbs);
        _nfixed_sub(C_ENTRY(l, 0), C_ENTRY(l, 0), Ccol_ENTRY(l), nlimbs);
    }

    _nfixed_add(val1, Crow, C_ENTRY(0, 0), nlimbs);
    nfixed_div2(val0, val1, nlimbs);
    _nfixed_sub(C_ENTRY(0, 0), C_ENTRY(0, 0), val0, nlimbs);

    for (k = 1; k < p; k++)
    {
        _nfixed_add(crow, Crow_ENTRY(k), C_ENTRY(0, k), nlimbs);
        nfixed_div2(val1, crow, nlimbs);
        _nfixed_sub(C_ENTRY(0, k), C_ENTRY(0, k), val1, nlimbs);
        _nfixed_sub(crow, val1, val0, nlimbs);

        for (l = 1; l < m; l++)
        {
            _nfixed_sub(val2, C_ENTRY(l, k), crow, nlimbs);
            _nfixed_sub(C_ENTRY(l, k), val2, Ccol_ENTRY(l), nlimbs);
        }
    }

    if ((n & 1) == 1)
    {
        for (l = 0; l < m; l++)
        {
            for (k = 0; k < p; k++)
            {
                nfixed_mul(val0, A_ENTRY(l, n-1), B_ENTRY(n-1, k), nlimbs);
                _nfixed_add(C_ENTRY(l, k), C_ENTRY(l, k), val0, nlimbs);
            }
        }
    }

    flint_free(Ctmp);

#undef A_ENTRY
#undef B_ENTRY
#undef C_ENTRY
}

void
_nfixed_mat_mul_waksman(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
{
    _nfixed_mat_mul_waksman2(C, A, B, m, n, p, nlimbs, p * (nlimbs + 1), n * (nlimbs + 1), p * (nlimbs + 1));
}


typedef struct
{
    nn_ptr start;
    slong r;
    slong c;
    slong row_stride;
}
_nfixed_mat_struct;

typedef _nfixed_mat_struct _nfixed_mat_t[1];

static void
_nfixed_mat_init(_nfixed_mat_t A, slong r, slong c, slong nlimbs)
{
    A->start = flint_malloc((nlimbs + 1) * (r * c) * sizeof(ulong));
    A->r = r;
    A->c = c;
    A->row_stride = c * (nlimbs + 1);
}

static void
_nfixed_mat_clear(_nfixed_mat_t A, slong nlimbs)
{
    flint_free(A->start);
}

static void
_nfixed_mat_window_init(_nfixed_mat_t A, const _nfixed_mat_t mat, slong r1, slong c1, slong r2, slong c2, slong nlimbs)
{
    A->start = mat->start + (r1 * mat->row_stride) + c1 * (nlimbs + 1);
    A->r = r2 - r1;
    A->c = c2 - c1;
    A->row_stride = mat->row_stride;
}

static void
_nfixed_mat_window_clear(_nfixed_mat_t A, slong nlimbs)
{
}

/*
static void
nfixed_mat_print(nn_ptr A, slong ar, slong ac, slong nlimbs)
{
    slong i, j;
    flint_printf("{%wd y %wd : [", ar, ac);
    for (i = 0; i < ar; i++)
        for (j = 0; j < ac; j++)
        {
            _nfixed_print(A + i * ac * (nlimbs + 1) + j * (nlimbs + 1), nlimbs, 0);
            flint_printf(", ");
        }

    flint_printf("]}\n");
}

static void
_nfixed_mat_print(_nfixed_mat_t A, slong nlimbs)
{
    slong i, j;
    flint_printf("{%wd x %wd : [", A->r, A->c);
    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
        {
            _nfixed_print(A->start + i * A->row_stride + j * (nlimbs + 1), nlimbs, 0);
            flint_printf(", ");
        }

    flint_printf("]}\n");
}
*/

static void
_nfixed_mat_add(_nfixed_mat_t C, const _nfixed_mat_t A, const _nfixed_mat_t B, slong nlimbs)
{
    nn_srcptr Aptr, Bptr;
    nn_ptr Cptr;

    Aptr = A->start;
    Bptr = B->start;
    Cptr = C->start;

    slong Astride = A->row_stride;
    slong Bstride = B->row_stride;
    slong Cstride = C->row_stride;

    slong i, r = A->r, c = A->c;

    for (i = 0; i < r; i++)
        _nfixed_vec_add(Cptr + i * Cstride, Aptr + i * Astride, Bptr + i * Bstride, c, nlimbs);
}

static void
_nfixed_mat_sub(_nfixed_mat_t C, const _nfixed_mat_t A, const _nfixed_mat_t B, slong nlimbs)
{
    nn_srcptr Aptr, Bptr;
    nn_ptr Cptr;

    Aptr = A->start;
    Bptr = B->start;
    Cptr = C->start;

    slong Astride = A->row_stride;
    slong Bstride = B->row_stride;
    slong Cstride = C->row_stride;

    slong i, r = A->r, c = A->c;

    for (i = 0; i < r; i++)
        _nfixed_vec_sub(Cptr + i * Cstride, Aptr + i * Astride, Bptr + i * Bstride, c, nlimbs);
}

void
_nfixed_mat_mul_waksman3(_nfixed_mat_t C, const _nfixed_mat_t A, const _nfixed_mat_t B, slong nlimbs)
{
    nn_srcptr Aptr, Bptr;
    nn_ptr Cptr;

    Aptr = A->start;
    Bptr = B->start;
    Cptr = C->start;

    slong Astride = A->row_stride;
    slong Bstride = B->row_stride;
    slong Cstride = C->row_stride;

    slong m = A->r;
    slong n = A->c;
    slong p = B->c;

    _nfixed_mat_mul_waksman2(Cptr, Aptr, Bptr, m, n, p, nlimbs, Cstride, Astride, Bstride);
}

static void
_nfixed_mat_mul_classical2(_nfixed_mat_t C, const _nfixed_mat_t A, const _nfixed_mat_t B, slong nlimbs)
{
    nn_srcptr Aptr, Bptr;
    nn_ptr Cptr;

    Aptr = A->start;
    Bptr = B->start;
    Cptr = C->start;

    slong Astride = A->row_stride;
    slong Bstride = B->row_stride;
    slong Cstride = C->row_stride;

    slong m = A->r;
    slong n = A->c;
    slong p = B->c;

    slong i, j, k;
    nn_ptr t;

#define A_ENTRY(i, j) ((Aptr) + (i) * Astride + (j) * (nlimbs + 1))
#define B_ENTRY(i, j) ((Bptr) + (i) * Bstride + (j) * (nlimbs + 1))
#define C_ENTRY(i, j) ((Cptr) + (i) * Cstride + (j) * (nlimbs + 1))

    if (n == 1)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                nfixed_mul(C_ENTRY(i, j), A_ENTRY(i, 0), B_ENTRY(0, j), nlimbs);
        return;
    }

    if (nlimbs == 2)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_2(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), Bstride, n);
    }
    else if (nlimbs == 3)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_3(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), Bstride, n);
    }
    else if (nlimbs == 4)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_4(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), Bstride, n);
    }
    else if (nlimbs == 5)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_5(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), Bstride, n);
    }
    else if (nlimbs == 6)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_6(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), Bstride, n);
    }
    else if (nlimbs == 7)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_7(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), Bstride, n);
    }
    else if (nlimbs == 8)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfixed_dot_8(C_ENTRY(i, j), A_ENTRY(i, 0), nlimbs + 1, B_ENTRY(0, j), Bstride, n);
    }
    else
    {
        TMP_INIT;
        TMP_START;

        t = TMP_ALLOC((nlimbs + 1) * sizeof(ulong));

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < p; j++)
            {
                nfixed_mul(C_ENTRY(i, j), A_ENTRY(i, 0), B_ENTRY(0, j), nlimbs);

                for (k = 1; k < n; k++)
                {
                    nfixed_mul(t, A_ENTRY(i, k), B_ENTRY(k, j), nlimbs);
                    _nfixed_add(C_ENTRY(i, j), C_ENTRY(i, j), t, nlimbs);
                }
            }
        }

        TMP_END;
    }

#undef A_ENTRY
#undef B_ENTRY
#undef C_ENTRY
}


static void
_nfixed_mat_mul_strassen2(_nfixed_mat_t C, const _nfixed_mat_t A, const _nfixed_mat_t B, slong cutoff, slong nlimbs)
{
    slong ar, ac, bc, nn;
    slong anr, anc, bnr, bnc;

    _nfixed_mat_t A11, A12, A21, A22;
    _nfixed_mat_t B11, B12, B21, B22;
    _nfixed_mat_t C11, C12, C21, C22;
    _nfixed_mat_t X1, X2;

    ar = A->r;
    ac = A->c;
    bc = B->c;

    nn = FLINT_MIN(ar, ac);
    nn = FLINT_MIN(nn, bc);

    /* Important: if the cutoff handling changes, _nfixed_mat_mul_bound_strassen must change too. */
    if (cutoff < 0)
        cutoff = nfixed_mat_mul_strassen_cutoff(nn, ac & 1, nlimbs);
    else
        cutoff = FLINT_MAX(cutoff, 2);

    if (nn < cutoff)
    {
        if (nfixed_mat_mul_use_waksman(nn, nlimbs))
            _nfixed_mat_mul_waksman3(C, A, B, nlimbs);
        else
            _nfixed_mat_mul_classical2(C, A, B, nlimbs);
        return;
    }

    anr = ar / 2;
    anc = ac / 2;
    bnr = anc;
    bnc = bc / 2;

    _nfixed_mat_window_init(A11, A, 0, 0, anr, anc, nlimbs);
    _nfixed_mat_window_init(A12, A, 0, anc, anr, 2 * anc, nlimbs);
    _nfixed_mat_window_init(A21, A, anr, 0, 2 * anr, anc, nlimbs);
    _nfixed_mat_window_init(A22, A, anr, anc, 2 * anr, 2 * anc, nlimbs);

    _nfixed_mat_window_init(B11, B, 0, 0, bnr, bnc, nlimbs);
    _nfixed_mat_window_init(B12, B, 0, bnc, bnr, 2 * bnc, nlimbs);
    _nfixed_mat_window_init(B21, B, bnr, 0, 2 * bnr, bnc, nlimbs);
    _nfixed_mat_window_init(B22, B, bnr, bnc, 2 * bnr, 2 * bnc, nlimbs);

    _nfixed_mat_window_init(C11, C, 0, 0, anr, bnc, nlimbs);
    _nfixed_mat_window_init(C12, C, 0, bnc, anr, 2 * bnc, nlimbs);
    _nfixed_mat_window_init(C21, C, anr, 0, 2 * anr, bnc, nlimbs);
    _nfixed_mat_window_init(C22, C, anr, bnc, 2 * anr, 2 * bnc, nlimbs);

    _nfixed_mat_init(X1, anr, FLINT_MAX(bnc, anc), nlimbs);
    _nfixed_mat_init(X2, anc, bnc, nlimbs);

    X1->c = anc;

    _nfixed_mat_add(X1, A22, A12, nlimbs);
    _nfixed_mat_add(X2, B22, B12, nlimbs);
    _nfixed_mat_mul_strassen2(C21, X1, X2, cutoff, nlimbs);

    _nfixed_mat_sub(X1, A22, A21, nlimbs);
    _nfixed_mat_sub(X2, B22, B21, nlimbs);
    _nfixed_mat_mul_strassen2(C22, X1, X2, cutoff, nlimbs);

    _nfixed_mat_add(X1, X1, A12, nlimbs);
    _nfixed_mat_add(X2, X2, B12, nlimbs);
    _nfixed_mat_mul_strassen2(C11, X1, X2, cutoff, nlimbs);

    _nfixed_mat_sub(X1, X1, A11, nlimbs);
    _nfixed_mat_mul_strassen2(C12, X1, B12, cutoff, nlimbs);

    X1->c = bnc;
    _nfixed_mat_mul_strassen2(X1, A12, B21, cutoff, nlimbs);
    _nfixed_mat_add(C11, C11, X1, nlimbs);
    _nfixed_mat_add(C12, C12, C22, nlimbs);
    _nfixed_mat_sub(C12, C11, C12, nlimbs);
    _nfixed_mat_sub(C11, C21, C11, nlimbs);
    _nfixed_mat_sub(X2, X2, B11, nlimbs);
    _nfixed_mat_mul_strassen2(C21, A21, X2, cutoff, nlimbs);

    _nfixed_mat_clear(X2, nlimbs);

    _nfixed_mat_sub(C21, C11, C21, nlimbs);
    _nfixed_mat_add(C22, C22, C11, nlimbs);
    _nfixed_mat_mul_strassen2(C11, A11, B11, cutoff, nlimbs);

    _nfixed_mat_add(C11, X1, C11, nlimbs);

    X1->c = FLINT_MAX(bnc, anc);
    _nfixed_mat_clear(X1, nlimbs);

    _nfixed_mat_window_clear(A11, nlimbs);
    _nfixed_mat_window_clear(A12, nlimbs);
    _nfixed_mat_window_clear(A21, nlimbs);
    _nfixed_mat_window_clear(A22, nlimbs);

    _nfixed_mat_window_clear(B11, nlimbs);
    _nfixed_mat_window_clear(B12, nlimbs);
    _nfixed_mat_window_clear(B21, nlimbs);
    _nfixed_mat_window_clear(B22, nlimbs);

    _nfixed_mat_window_clear(C11, nlimbs);
    _nfixed_mat_window_clear(C12, nlimbs);
    _nfixed_mat_window_clear(C21, nlimbs);
    _nfixed_mat_window_clear(C22, nlimbs);

    if (bc > 2 * bnc)
    {
        _nfixed_mat_t Bc, Cc;
        _nfixed_mat_window_init(Bc, B, 0, 2 * bnc, ac, bc, nlimbs);
        _nfixed_mat_window_init(Cc, C, 0, 2 * bnc, ar, bc, nlimbs);

        _nfixed_mat_mul_classical2(Cc, A, Bc, nlimbs);
        _nfixed_mat_window_clear(Bc, nlimbs);
        _nfixed_mat_window_clear(Cc, nlimbs);
    }

    if (ar > 2 * anr)
    {
        _nfixed_mat_t Ar, Bc, Cr;
        _nfixed_mat_window_init(Ar, A, 2 * anr, 0, ar, ac, nlimbs);
        _nfixed_mat_window_init(Bc, B, 0, 0, ac, 2 * bnc, nlimbs);
        _nfixed_mat_window_init(Cr, C, 2 * anr, 0, ar, 2 * bnc, nlimbs);

        _nfixed_mat_mul_classical2(Cr, Ar, Bc, nlimbs);
        _nfixed_mat_window_clear(Ar, nlimbs);
        _nfixed_mat_window_clear(Bc, nlimbs);
        _nfixed_mat_window_clear(Cr, nlimbs);
    }

    if (ac > 2 * anc)
    {
        _nfixed_mat_t Ac, Br, Cb, tmp;
        slong mt, nt;

        _nfixed_mat_window_init(Ac, A, 0, 2 * anc, 2 * anr, ac, nlimbs);
        _nfixed_mat_window_init(Br, B, 2 * bnr, 0, ac, 2 * bnc, nlimbs);
        _nfixed_mat_window_init(Cb, C, 0, 0, 2 * anr, 2 * bnc, nlimbs);

        mt = Ac->r;
        nt = Br->c;

        /* todo: faster */
        _nfixed_mat_init(tmp, mt, nt, nlimbs);
        _nfixed_mat_mul_classical2(tmp, Ac, Br, nlimbs);
        _nfixed_mat_add(Cb, Cb, tmp, nlimbs);
        _nfixed_mat_clear(tmp, nlimbs);
        _nfixed_mat_window_clear(Ac, nlimbs);
        _nfixed_mat_window_clear(Br, nlimbs);
        _nfixed_mat_window_clear(Cb, nlimbs);
    }
}

void
_nfixed_mat_mul_strassen(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong cutoff, slong nlimbs)
{
    _nfixed_mat_t CC, AA, BB;

    AA->start = (nn_ptr) A;
    AA->r = m;
    AA->c = n;
    AA->row_stride = n * (nlimbs + 1);

    BB->start = (nn_ptr) B;
    BB->r = n;
    BB->c = p;
    BB->row_stride = p * (nlimbs + 1);

    CC->start = C;
    CC->r = m;
    CC->c = p;
    CC->row_stride = p * (nlimbs + 1);

    _nfixed_mat_mul_strassen2(CC, AA, BB, cutoff, nlimbs);
}

void
_nfixed_mat_mul(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
{
    slong d, cutoff;

    d = FLINT_MIN(m, n);
    d = FLINT_MIN(d, p);

    /* Important: if the cutoff handling changes, _nfixed_mat_mul_bound must change too. */
    if (d > 10)
    {
        cutoff = nfixed_mat_mul_strassen_cutoff(d, n & 1, nlimbs);

        if (n >= cutoff)
        {
            _nfixed_mat_mul_strassen(C, A, B, m, n, p, -1, nlimbs);
            return;
        }
    }

    if (nfixed_mat_mul_use_waksman(d, nlimbs))
        _nfixed_mat_mul_waksman(C, A, B, m, n, p, nlimbs);
    else
        _nfixed_mat_mul_classical(C, A, B, m, n, p, nlimbs);
}

/*
  Given an m x n x p matrix multiplication with inputs bounded
  entrywise by A, B and nlimbs precision:

  - Set *bound* to a bound for the entries in all intermediate
    variables of the computation. This should be < 1 to
    ensure correctness.
  - Set *error* to a bound for the output error, measured in ulp.

   The caller can assume that the bound is a nondecreasing function
   of A and B.

   IMPORTANT: when changing the algorithm in _nfixed_mat_mul, this
   must be changed to correspond.
*/

void
_nfixed_mat_mul_bound_classical(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
{
    double fixup = 1.0 + 1e-6;

    /* Error bound (in ulp) for naive scalar multiplication, and for dot product */
    double error_mul = (2 * nlimbs - 1);
    double error_dot = n * error_mul;

    *bound = (n * A * B) * fixup;
    *error = error_dot * fixup;
}

void
_nfixed_mat_mul_bound_waksman(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
{
    double fixup = 1.0 + 1e-6;

    /* Error bound (in ulp) for naive scalar multiplication */
    double error_mul = (2 * nlimbs - 1);

    *bound = FLINT_MAX(A + B, 6 * (n / 2) * (A + B) * (A + B) + A * B) * fixup;
    *error = ((6 * (n / 2) + 1) * error_mul + 5) * fixup;
}

void
_nfixed_mat_mul_bound_strassen(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong cutoff, slong nlimbs)
{
    double fixup = 1.0 + 1e-6;
    slong d;

    /* Error bound (in ulp) for naive scalar multiplication, and for dot product */
    double error_mul = (2 * nlimbs - 1);
    double error_dot = n * error_mul;

    d = FLINT_MIN(m, n);
    d = FLINT_MIN(d, p);

    if (cutoff < 0)
        cutoff = nfixed_mat_mul_strassen_cutoff(d, n, nlimbs);
    else
        cutoff = FLINT_MAX(cutoff, 2);

    if (d < cutoff)
    {
        if (nfixed_mat_mul_use_waksman(d, nlimbs))
            _nfixed_mat_mul_bound_waksman(bound, error, m, n, p, A, B, nlimbs);
        else
            _nfixed_mat_mul_bound_classical(bound, error, m, n, p, A, B, nlimbs);
        return;
    }

    slong mm, nn, pp;
    double bound_transformed_A, bound_transformed_B;
    double bound_everything, bound_recursive, error_recursive;
    double bound_recombination, error_recombination;
    double ulp;

    /* Bound for entries of transformed block matrices */
    /* S1 = A22 + A12   <= 2A
       S2 = A22 - A21   <= 2A
       S3 = S2 + A12    <= 3A
       S4 = S3 - A11    <= 3A, and similarly Ti for B */
    bound_transformed_A = 3.0 * A;
    bound_transformed_B = 3.0 * B;
    bound_everything = FLINT_MAX(bound_transformed_A, bound_transformed_B);

    /* Bound intermediate entries and errors for recursive multiplications. */
    mm = m / 2;
    nn = n / 2;
    pp = p / 2;
    _nfixed_mat_mul_bound_strassen(&bound_recursive, &error_recursive, mm, nn, pp, bound_transformed_A, bound_transformed_B, cutoff, nlimbs);
    bound_everything = FLINT_MAX(bound_everything, bound_recursive);

    /* Bound for recombinations. We don't use bound_recursive here,
       because this can be a huge overestimate if the basecase
       is nonclassical multiplication. Instead, we use the
       theoretical bounds for the subproducts and add the
       bound u = error_recursive (in the event the recursive
       multiplications were not rounded down). */
    /* P1 = S1 T1      <= 4 nn A B + u
       P2 = S2 T2      <= 4 nn A B + u
       P3 = S3 T3      <= 9 nn A B + u
       P4 = A11 B11    <= nn A B + u
       P5 = A12 B21    <= nn A B + u
       P6 = S4 B12     <= 3 nn A B + u
       P7 = A21 T4     <= 3 nn A B + u
       U1 = P3 + P5    <= 10 nn A B + 2u
       U2 = P1 - U1    <= 14 nn A B + 3u
       U3 = U1 - P2    <= 14 nn A B + 3u
       C11 = P4 + P5   <= 2 nn A B + 2u
       C12 = U3 - P6   <= 17 nn A B + 4u
       C21 = U2 - P7   <= 17 nn A B + 4u
       C22 = P2 + U2   <= 18 nn A B + 4u
    */

    ulp = ldexp(1.0, FLINT_MAX(-128, -nlimbs * FLINT_BITS));

    error_recombination = 4 * error_recursive;
    bound_recombination = 18 * nn * A * B + error_recombination * ulp;

    /* Bound for border corrections when m, n, and/or p is odd.
       Todo: could be added conditionally. Assumes border
       corrections use classical multiplication. */
    bound_recombination += A * B;
    bound_recombination = FLINT_MAX(bound_recombination, n * A * B);
    error_recombination += error_mul;
    error_recombination = FLINT_MAX(error_recombination, error_dot);

    bound_everything = FLINT_MAX(bound_everything, bound_recombination);

    *bound = bound_everything * fixup;
    *error = error_recombination * fixup;
}

void
_nfixed_mat_mul_bound(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
{
    slong d, cutoff;

    d = FLINT_MIN(m, n);
    d = FLINT_MIN(d, p);

    if (d > 10)
    {
        cutoff = nfixed_mat_mul_strassen_cutoff(d, n & 1, nlimbs);

        if (n >= cutoff)
        {
            _nfixed_mat_mul_bound_strassen(bound, error, m, n, p, A, B, -1, nlimbs);
            return;
        }
    }

    if (nfixed_mat_mul_use_waksman(d, nlimbs))
        _nfixed_mat_mul_bound_waksman(bound, error, m, n, p, A, B, nlimbs);
    else
        _nfixed_mat_mul_bound_classical(bound, error, m, n, p, A, B, nlimbs);
}

/* Karatsuba formula */
/* (A C - B D) + ((A + B)(C + D) - A C - B D) i */
void
_nfixed_complex_mat_mul_bound(double * bound, double * error, slong m, slong n, slong p, double A, double B, double C, double D, slong nlimbs)
{
    double rbound, rerror, fixup = 1.0 + 1e-6;

    _nfixed_mat_mul_bound(&rbound, &rerror, m, n, p, A + B, C + D, nlimbs);

    (*bound) = FLINT_MAX(3.0 * rbound, FLINT_MAX(A + B, C + D)) * fixup;
    (*error) = rerror * (3.0 * fixup);
}
