/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

/* currently defined in poly_mullow_karatsuba.c */
void _mpn_dot_rev_2x2_3(mp_ptr s, mp_srcptr a, mp_srcptr b, slong len);
void _mpn_dot_rev_2x2_4(mp_ptr s, mp_srcptr a, mp_srcptr b, slong len);
void _mpn_dot_rev_2x2_5(mp_ptr s, mp_srcptr a, mp_srcptr b, slong len);
void _mpn_dot_rev_3x3_5(mp_ptr s, mp_srcptr a, mp_srcptr b, slong len);
void _mpn_dot_rev_nxn_2n(mp_ptr res, mp_srcptr a, mp_srcptr b, slong len, mp_size_t nlimbs);
void _mpn_dot_rev_nxn_2nm1(mp_ptr res, mp_srcptr a, mp_srcptr b, slong len, mp_size_t nlimbs);
void _mpn_dot_rev_nxn_2np1(mp_ptr res, mp_srcptr a, mp_srcptr b, slong len, mp_size_t nlimbs);

static void
mpn_mod_set_mpn2(mp_ptr res, mp_srcptr s, mp_size_t l, gr_ctx_t ctx)
{
    MPN_NORM(s, l);
    mpn_mod_set_mpn(res, s, l, ctx);
}

int
_mpn_mod_poly_mullow_classical(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
{
    slong i, top1, top2;
    mp_size_t nlimbs, slimbs;
    mp_bitcnt_t sbits;
    mp_limb_t s[2 * MPN_MOD_MAX_LIMBS + 1];
    int squaring;

    if (len == 1)
        return mpn_mod_mul(res, poly1, poly2, ctx);

    len1 = FLINT_MIN(len1, len);
    len2 = FLINT_MIN(len2, len);

    if (len1 == 1)
        return _mpn_mod_vec_mul_scalar(res, poly2, len, poly1, ctx);

    if (len2 == 1)
        return _mpn_mod_vec_mul_scalar(res, poly1, len, poly2, ctx);

    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    mpn_mod_mul(res, poly1, poly2, ctx);

    if (len == len1 + len2 - 1)
    {
        mpn_mod_mul(res + (len1 + len2 - 2) * nlimbs, poly1 + (len1 - 1) * nlimbs, poly2 + (len2 - 1) * nlimbs, ctx);
        len--;
    }

    if (len <= 1)
        return GR_SUCCESS;

    squaring = (poly1 == poly2) && (len1 == len2);
    sbits = MPN_MOD_CTX_MODULUS_BITS(ctx);
    sbits = 2 * sbits + FLINT_BIT_COUNT(FLINT_MIN(len1, len2));
    slimbs = (sbits + FLINT_BITS - 1) / FLINT_BITS;

    FLINT_ASSERT((slimbs == 2 * nlimbs) || (slimbs == 2 * nlimbs + 1) || (slimbs == 2 * nlimbs - 1));

    if (squaring)
    {
        mp_limb_t t[2 * MPN_MOD_MAX_LIMBS];
        slong start, stop;

        if (slimbs == 2 * nlimbs - 1)
        {
            for (i = 1; i < FLINT_MIN(len, 2 * len1 - 2); i++)
            {
                start = FLINT_MAX(0, i - len1 + 1);
                stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                if (nlimbs == 2)
                    _mpn_dot_rev_2x2_3(s, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1);
                else if (nlimbs == 3)
                    _mpn_dot_rev_3x3_5(s, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1);
                else
                    _mpn_dot_rev_nxn_2nm1(s, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1, nlimbs);

                mpn_lshift(s, s, slimbs, 1);
                if (i % 2 == 0 && i / 2 < len1)
                {
                    flint_mpn_sqr(t, poly1 + (i / 2) * nlimbs, nlimbs);
                    mpn_add_n(s, s, t, slimbs);
                }

                mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
            }
        }
        else if (slimbs == 2 * nlimbs)
        {
            for (i = 1; i < FLINT_MIN(len, 2 * len1 - 2); i++)
            {
                start = FLINT_MAX(0, i - len1 + 1);
                stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                if (nlimbs == 2)
                    _mpn_dot_rev_2x2_4(s, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1);
                else
                    _mpn_dot_rev_nxn_2n(s, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1, nlimbs);

                mpn_lshift(s, s, slimbs, 1);
                if (i % 2 == 0 && i / 2 < len1)
                {
                    flint_mpn_sqr(t, poly1 + (i / 2) * nlimbs, nlimbs);
                    mpn_add_n(s, s, t, slimbs);
                }

                mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
            }
        }
        else
        {
            for (i = 1; i < FLINT_MIN(len, 2 * len1 - 2); i++)
            {
                start = FLINT_MAX(0, i - len1 + 1);
                stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

                if (nlimbs == 2)
                    _mpn_dot_rev_2x2_5(s, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1);
                else
                    _mpn_dot_rev_nxn_2np1(s, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1, nlimbs);
                mpn_lshift(s, s, slimbs, 1);
                if (i % 2 == 0 && i / 2 < len1)
                {
                    flint_mpn_sqr(t, poly1 + (i / 2) * nlimbs, nlimbs);
                    s[slimbs - 1] += mpn_add_n(s, s, t, slimbs - 1);
                }

                mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
            }
        }

        return GR_SUCCESS;
    }

    if (nlimbs == 2)
    {
        if (slimbs == 3)
        {
            for (i = 1; i < len; i++)
            {
                top1 = FLINT_MIN(len1 - 1, i);
                top2 = FLINT_MIN(len2 - 1, i);
                _mpn_dot_rev_2x2_3(s, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1);
                mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
            }
        }
        else if (slimbs == 4)
        {
            for (i = 1; i < len; i++)
            {
                top1 = FLINT_MIN(len1 - 1, i);
                top2 = FLINT_MIN(len2 - 1, i);
                _mpn_dot_rev_2x2_4(s, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1);
                mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
            }
        }
        else
        {
            for (i = 1; i < len; i++)
            {
                top1 = FLINT_MIN(len1 - 1, i);
                top2 = FLINT_MIN(len2 - 1, i);
                _mpn_dot_rev_2x2_5(s, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1);
                mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
            }
        }
    }
    else if (nlimbs == 3 && slimbs == 5)
    {
        for (i = 1; i < len; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            _mpn_dot_rev_3x3_5(s, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1);
            mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
        }
    }
    else if (slimbs == 2 * nlimbs)
    {
        for (i = 1; i < len; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            _mpn_dot_rev_nxn_2n(s, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1, nlimbs);
            mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
        }
    }
    else if (slimbs == 2 * nlimbs + 1)
    {
        for (i = 1; i < len; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            _mpn_dot_rev_nxn_2np1(s, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1, nlimbs);
            mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
        }
    }
    else if (slimbs == 2 * nlimbs - 1)
    {
        for (i = 1; i < len; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            _mpn_dot_rev_nxn_2np1(s, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1, nlimbs);
            mpn_mod_set_mpn2(res + i * nlimbs, s, slimbs, ctx);
        }
    }

    return GR_SUCCESS;
}
