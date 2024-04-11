/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

void _mpn_dot_rev_2x2_3(mp_ptr s, mp_srcptr a, mp_srcptr b, slong len)
{
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t p2, p1, p0;
    mp_limb_t s2, s1, s0;
    mp_limb_t u2, u1;
    mp_limb_t v2;
    slong k;

    s2 = s1 = s0 = 0;
    u2 = u1 = 0;
    v2 = 0;

    for (k = 0; k < len; k++)
    {
        A0 = a[2 * k + 0];
        A1 = a[2 * k + 1];
        B0 = b[2 * (len - 1 - k) + 0];
        B1 = b[2 * (len - 1 - k) + 1];

        umul_ppmm(p2, p1, A1, B0);
        add_ssaaaa(u2, u1, u2, u1, p2, p1);

        p2 = A1 * B1;
        umul_ppmm(p1, p0, A0, B0);
        add_sssaaaaaa(s2, s1, s0, s2, s1, s0, p2, p1, p0);

        umul_ppmm(p2, p1, A0, B1);
        add_ssaaaa(v2, u1, v2, u1, p2, p1);
    }

    u2 = u2 + v2;
    add_ssaaaa(s2, s1, s2, s1, u2, u1);

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
}

void _mpn_dot_rev_2x2_4(mp_ptr s, mp_srcptr a, mp_srcptr b, slong len)
{
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t p3, p2, p1, p0;
    mp_limb_t s3, s2, s1, s0;
    mp_limb_t u3, u2, u1;
    mp_limb_t v3, v2;
    slong k;

    s3 = s2 = s1 = s0 = 0;
    u3 = u2 = u1 = 0;
    v3 = v2 = 0;

    for (k = 0; k < len; k++)
    {
        A0 = a[2 * k + 0];
        A1 = a[2 * k + 1];
        B0 = b[2 * (len - 1 - k) + 0];
        B1 = b[2 * (len - 1 - k) + 1];

        umul_ppmm(p2, p1, A1, B0);
        add_sssaaaaaa(u3, u2, u1, u3, u2, u1, UWORD(0), p2, p1);

        umul_ppmm(p3, p2, A1, B1);
        umul_ppmm(p1, p0, A0, B0);
        add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, p3, p2, p1, p0);

        umul_ppmm(p2, p1, A0, B1);
        add_sssaaaaaa(v3, v2, u1, v3, v2, u1, UWORD(0), p2, p1);
    }

    add_ssaaaa(u3, u2, u3, u2, v3, v2);
    add_sssaaaaaa(s3, s2, s1, s3, s2, s1, u3, u2, u1);

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
}

void _mpn_dot_rev_2x2_5(mp_ptr s, mp_srcptr a, mp_srcptr b, slong len)
{
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t p3, p2, p1, p0;
    mp_limb_t s4, s3, s2, s1, s0;
    mp_limb_t u2, u1;
    mp_limb_t v3, v2;
    slong k;

    s4 = s3 = s2 = s1 = s0 = 0;
    u2 = u1 = 0;
    v3 = v2 = 0;

    for (k = 0; k < len; k++)
    {
        A0 = a[2 * k + 0];
        A1 = a[2 * k + 1];
        B0 = b[2 * (len - 1 - k) + 0];
        B1 = b[2 * (len - 1 - k) + 1];

        umul_ppmm(p2, p1, A1, B0);
        add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);
        umul_ppmm(p1, p0, B0, A0);
        add_sssaaaaaa(u2, u1, s0, u2, u1, s0, UWORD(0), p1, p0);
        umul_ppmm(p2, p1, B1, A0);
        add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);
        umul_ppmm(p3, p2, B1, A1);
        add_sssaaaaaa(s4, v3, v2, s4, v3, v2, UWORD(0), p3, p2);
    }

    /* s3 is small, so this doesn't overflow */
    add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), u2, u1);
    add_sssaaaaaa(s4, s3, s2, s4, s3, s2, UWORD(0), v3, v2);

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
    s[4] = s4;
}

void _mpn_dot_rev_3x3_5(mp_ptr s, mp_srcptr a, mp_srcptr b, slong len)
{
    mp_limb_t A0, A1, A2, B0, B1, B2;
    mp_limb_t p4, p3, p2, p1, p0;
    mp_limb_t s4, s3, s2, s1, s0;
    mp_limb_t u2, u1;
    mp_limb_t v3, v2;
    slong k;

    s4 = s3 = s2 = s1 = s0 = 0;
    u2 = u1 = 0;
    v3 = v2 = 0;

    for (k = 0; k < len; k++)
    {
        A0 = a[3 * k + 0];
        A1 = a[3 * k + 1];
        A2 = a[3 * k + 2];
        B0 = b[3 * (len - 1 - k) + 0];
        B1 = b[3 * (len - 1 - k) + 1];
        B2 = b[3 * (len - 1 - k) + 2];

        umul_ppmm(p2, p1, A1, B0);
        add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);
        umul_ppmm(p1, p0, B0, A0);
        add_sssaaaaaa(u2, u1, s0, u2, u1, s0, UWORD(0), p1, p0);
        umul_ppmm(p2, p1, B1, A0);
        add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);
        umul_ppmm(p3, p2, B1, A1);
        add_sssaaaaaa(s4, v3, v2, s4, v3, v2, UWORD(0), p3, p2);
        umul_ppmm(p3, p2, A2, B0);
        add_sssaaaaaa(s4, v3, v2, s4, v3, v2, UWORD(0), p3, p2);
        umul_ppmm(p3, p2, B2, A0);
        add_sssaaaaaa(s4, v3, v2, s4, v3, v2, UWORD(0), p3, p2);
        umul_ppmm(p4, p3, A2, B1);
        add_ssaaaa(s4, v3, s4, v3, p4, p3);
        umul_ppmm(p4, p3, B2, A1);
        add_ssaaaa(s4, v3, s4, v3, p4, p3);
        s4 += B2 * A2;
    }

    /* s3 is small, so this doesn't overflow */
    add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), u2, u1);
    add_sssaaaaaa(s4, s3, s2, s4, s3, s2, UWORD(0), v3, v2);

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
    s[4] = s4;
}

void
_mpn_dot_rev_nxn_2n(mp_ptr res, mp_srcptr a, mp_srcptr b, slong len, mp_size_t nlimbs)
{
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS + 3];
    mp_size_t slimbs = 2 * nlimbs;
    slong j;

    flint_mpn_mul_n(res, a, b + (len - 1) * nlimbs, nlimbs);

    for (j = 1; j < len; j++)
    {
        flint_mpn_mul_n(t, a + j * nlimbs, b + (len - 1 - j) * nlimbs, nlimbs);
        mpn_add_n(res, res, t, slimbs);
    }
}

void
_mpn_dot_rev_nxn_2nm1(mp_ptr res, mp_srcptr a, mp_srcptr b, slong len, mp_size_t nlimbs)
{
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS + 3];
    mp_size_t slimbs = 2 * nlimbs - 1;
    slong j;

    flint_mpn_mul_n(t, a, b + (len - 1) * nlimbs, nlimbs);
    flint_mpn_copyi(res, t, slimbs);

    for (j = 1; j < len; j++)
    {
        flint_mpn_mul_n(t, a + j * nlimbs, b + (len - 1 - j) * nlimbs, nlimbs);
        mpn_add_n(res, res, t, slimbs);
    }
}

void
_mpn_dot_rev_nxn_2np1(mp_ptr res, mp_srcptr a, mp_srcptr b, slong len, mp_size_t nlimbs)
{
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS + 3];
    mp_size_t slimbs = 2 * nlimbs + 1;
    slong j;

    flint_mpn_mul_n(res, a, b + (len - 1) * nlimbs, nlimbs);
    res[slimbs - 1] = 0;

    for (j = 1; j < len; j++)
    {
        flint_mpn_mul_n(t, a + j * nlimbs, b + (len - 1 - j) * nlimbs, nlimbs);
        res[slimbs - 1] += mpn_add_n(res, res, t, 2 * nlimbs);
    }
}


/* Sets {res,len*slimbs} viewed as a polynomial with length-slimbs coefficients
   to the product of {poly1,len1*nlimbs} and {poly2,len2*nlimbs} viewed as
   polynomials with length-nlimbs coefficients. */
static void
_mpn_poly_mullow_classical(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, mp_size_t nlimbs, mp_size_t slimbs)
{
    slong i, top1, top2;

    FLINT_ASSERT((slimbs == 2 * nlimbs) || (slimbs == 2 * nlimbs + 1) || (slimbs == 2 * nlimbs - 1));

    if (nlimbs == 2)
    {
        if (slimbs == 3)
        {
            for (i = 0; i < len; i++)
            {
                top1 = FLINT_MIN(len1 - 1, i);
                top2 = FLINT_MIN(len2 - 1, i);
                _mpn_dot_rev_2x2_3(res + i * slimbs, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1);
            }
        }
        else if (slimbs == 4)
        {
            for (i = 0; i < len; i++)
            {
                top1 = FLINT_MIN(len1 - 1, i);
                top2 = FLINT_MIN(len2 - 1, i);
                _mpn_dot_rev_2x2_4(res + i * slimbs, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1);
            }
        }
        else
        {
            for (i = 0; i < len; i++)
            {
                top1 = FLINT_MIN(len1 - 1, i);
                top2 = FLINT_MIN(len2 - 1, i);
                _mpn_dot_rev_2x2_5(res + i * slimbs, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1);
            }
        }
    }
    else if (nlimbs == 3 && slimbs == 5)
    {
        for (i = 0; i < len; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            _mpn_dot_rev_3x3_5(res + i * slimbs, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1);
        }
    }
    else if (slimbs == 2 * nlimbs - 1)
    {
        for (i = 0; i < len; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            _mpn_dot_rev_nxn_2nm1(res + i * slimbs, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1, nlimbs);
        }
    }
    else if (slimbs == 2 * nlimbs)
    {
        for (i = 0; i < len; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            _mpn_dot_rev_nxn_2n(res + i * slimbs, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1, nlimbs);
        }
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            _mpn_dot_rev_nxn_2np1(res + i * slimbs, poly1 + (i - top2) * nlimbs, poly2 + (i - top1) * nlimbs, top1 + top2 - i + 1, nlimbs);
        }
    }
}

static void
_mpn_poly_sqrlow_classical(mp_ptr res, mp_srcptr poly1, slong len1, slong len, mp_size_t nlimbs, mp_size_t slimbs)
{
    slong i, start, stop;
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS + 3];
    mp_ptr rp;

    FLINT_ASSERT((slimbs == 2 * nlimbs) || (slimbs == 2 * nlimbs + 1) || (slimbs == 2 * nlimbs - 1));

    if (slimbs == 2 * nlimbs - 1)
    {
        flint_mpn_sqr(t, poly1, nlimbs);
        flint_mpn_copyi(res, t, slimbs);

        for (i = 1; i < FLINT_MIN(len, 2 * len1 - 2); i++)
        {
            rp = res + i * slimbs;
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            if (nlimbs == 2)
                _mpn_dot_rev_2x2_3(rp, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1);
            else if (nlimbs == 3)
                _mpn_dot_rev_3x3_5(rp, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1);
            else
                _mpn_dot_rev_nxn_2nm1(rp, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1, nlimbs);

            mpn_lshift(rp, rp, slimbs, 1);
            if (i % 2 == 0 && i / 2 < len1)
            {
                flint_mpn_sqr(t, poly1 + (i / 2) * nlimbs, nlimbs);
                mpn_add_n(rp, rp, t, slimbs);
            }
        }

        if (len1 >= 2 && len >= 2 * len1 - 1)
        {
            flint_mpn_sqr(t, poly1 + (len1 - 1) * nlimbs, nlimbs);
            flint_mpn_copyi(res + (2 * len1 - 2) * slimbs, t, slimbs);
        }
    }
    else if (slimbs == 2 * nlimbs)
    {
        flint_mpn_sqr(res, poly1, nlimbs);

        for (i = 1; i < FLINT_MIN(len, 2 * len1 - 2); i++)
        {
            rp = res + i * slimbs;
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            if (nlimbs == 2)
                _mpn_dot_rev_2x2_4(rp, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1);
            else
                _mpn_dot_rev_nxn_2n(rp, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1, nlimbs);

            mpn_lshift(rp, rp, slimbs, 1);
            if (i % 2 == 0 && i / 2 < len1)
            {
                flint_mpn_sqr(t, poly1 + (i / 2) * nlimbs, nlimbs);
                mpn_add_n(rp, rp, t, slimbs);
            }
        }

        if (len1 >= 2 && len >= 2 * len1 - 1)
            flint_mpn_sqr(res + (2 * len1 - 2) * slimbs, poly1 + (len1 - 1) * nlimbs, nlimbs);
    }
    else
    {
        flint_mpn_sqr(res, poly1, nlimbs);
        res[slimbs - 1] = 0;

        for (i = 1; i < FLINT_MIN(len, 2 * len1 - 2); i++)
        {
            rp = res + i * slimbs;
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);
            if (nlimbs == 2)
                _mpn_dot_rev_2x2_5(rp, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1);
            else
                _mpn_dot_rev_nxn_2np1(rp, poly1 + start * nlimbs, poly1 + (i - stop) * nlimbs, stop - start + 1, nlimbs);
            mpn_lshift(rp, rp, slimbs, 1);
            if (i % 2 == 0 && i / 2 < len1)
            {
                flint_mpn_sqr(t, poly1 + (i / 2) * nlimbs, nlimbs);
                rp[slimbs - 1] += mpn_add_n(rp, rp, t, slimbs - 1);
            }
        }

        if (len1 >= 2 && len >= 2 * len1 - 1)
        {
            flint_mpn_sqr(res + (2 * len1 - 2) * slimbs, poly1 + (len1 - 1) * nlimbs, nlimbs);
            res[(2 * len1 - 1) * slimbs - 1] = 0;
        }
    }
}

FLINT_FORCE_INLINE void
_mpn_poly_add_n(mp_ptr res, mp_srcptr f, slong flen, mp_srcptr g, slong glen, mp_size_t nlimbs)
{
    slong m = FLINT_MIN(flen, glen);

    mpn_add_n(res, f, g, nlimbs * m);

    if (flen > glen)
        flint_mpn_copyi(res + nlimbs * m, f + nlimbs * m, nlimbs * (flen - m));
    if (glen > flen)
        flint_mpn_copyi(res + nlimbs * m, g + nlimbs * m, nlimbs * (glen - m));
}

/* inputs have nlimbs, output has nlimbs+1 */
FLINT_FORCE_INLINE void
_mpn_poly_add_n_carry(mp_ptr res, mp_srcptr f, slong flen, mp_srcptr g, slong glen, mp_size_t nlimbs)
{
    slong m = FLINT_MIN(flen, glen);
    slong i;
    mp_size_t nlimbs2 = nlimbs + 1;

    if (nlimbs == 2)
    {
        for (i = 0; i < m; i++)
            add_sssaaaaaa(res[i * nlimbs2 + 2], res[i * nlimbs2 + 1], res[i * nlimbs2],
                          0, f[i * nlimbs + 1], f[i * nlimbs],
                          0, g[i * nlimbs + 1], g[i * nlimbs]);

        for (i = m; i < flen; i++)
        {
            flint_mpn_copyi(res + i * nlimbs2, f + i * nlimbs, nlimbs);
            res[(i + 1) * nlimbs2 - 1] = 0;
        }

        for (i = m; i < glen; i++)
        {
            flint_mpn_copyi(res + i * nlimbs2, g + i * nlimbs, nlimbs);
            res[(i + 1) * nlimbs2 - 1] = 0;
        }
    }
    else
    {
        for (i = 0; i < m; i++)
            res[(i + 1) * nlimbs2 - 1] = mpn_add_n(res + i * nlimbs2, f + i * nlimbs, g + i * nlimbs, nlimbs);

        for (i = m; i < flen; i++)
        {
            flint_mpn_copyi(res + i * nlimbs2, f + i * nlimbs, nlimbs);
            res[(i + 1) * nlimbs2 - 1] = 0;
        }

        for (i = m; i < glen; i++)
        {
            flint_mpn_copyi(res + i * nlimbs2, g + i * nlimbs, nlimbs);
            res[(i + 1) * nlimbs2 - 1] = 0;
        }
    }
}

/* Sets {res,len*slimbs} viewed as a polynomial with length-slimbs coefficients
   to the product of {f,flen*nlimbs} and {g,glen*nlimbs} viewed as
   polynomials with length-nlimbs coefficients.

   Switches to classical multiplication below *cutoff*.

   The parameter *norm* bounds the number of leading zero bits in the
   top limb of the inputs. When this reaches zero, we need to increase
   *nlimbs* internally (this can occur on at most one level of recursion).
   It is assumed that *slimbs* is large enough to hold not only the
   final results but also the intermediate coefficients in the
   recursive additions and multiplications.

   TODO: consider doing modular additions and subtractions instead
   of incrementing nlimbs.
*/
static void
_mpn_poly_mul_karatsuba(mp_ptr res, mp_srcptr f, slong flen, mp_srcptr g, slong glen, mp_size_t nlimbs, mp_size_t slimbs, slong cutoff, int norm)
{
    slong m, f1len, g1len, tlen, ulen, vlen, alloc;
    mp_ptr t, u, v;
    mp_srcptr f0, f1, g0, g1;
    int squaring = (f == g) && (flen == glen);
    TMP_INIT;

    if (FLINT_MIN(flen, glen) < FLINT_MAX(cutoff, 2))
    {
        if (squaring)
            _mpn_poly_sqrlow_classical(res, f, flen, flen + glen - 1, nlimbs, slimbs);
        else
            _mpn_poly_mullow_classical(res, f, flen, g, glen, flen + glen - 1, nlimbs, slimbs);
        return;
    }

    /* split at X = x^m */
    /* res = f0 g0 + (f0 g1 + f1 g0) X + f1 g1 X^2
           = f0 g0 + ((f0 + f1) (g0 + g1) - f0 g0 - f1 g1) X + f1 g1 X^2 */
    m = (FLINT_MIN(flen, glen) + 1) / 2;

    f0 = f;
    g0 = g;
    f1 = f + m * nlimbs;
    g1 = g + m * nlimbs;
    f1len = flen - m;
    g1len = glen - m;

    /* low part: res[0, ..., 2m-2] = f0 g0 */
    _mpn_poly_mul_karatsuba(res, f, m, g, m, nlimbs, slimbs, cutoff, norm);

    /* res[2m-1] = 0 */
    flint_mpn_zero(res + (2 * m - 1) * slimbs, slimbs);

    /* high part: res[2m, ..., flen+glen-2] = f1 g1 */
    _mpn_poly_mul_karatsuba(res + (2 * m) * slimbs, f1, f1len, g1, g1len, nlimbs, slimbs, cutoff, norm);

    /* temporary space for the middle part */
    tlen = FLINT_MAX(m, f1len);
    ulen = FLINT_MAX(m, g1len);
    vlen = tlen + ulen - 1;
    alloc = tlen * (nlimbs + 1) + ulen * (nlimbs + 1) + vlen * slimbs;

    TMP_START;
    t = TMP_ALLOC(alloc * sizeof(mp_limb_t));
    u = t + tlen * (nlimbs + 1);
    v = u + ulen * (nlimbs + 1);

    /* f0 + f1 or g0 + g1 can overflow the limb boundary, so we repack them
       with one extra limb. */
    if (norm == 0)
    {
        /* t = f0 + f1 */
        _mpn_poly_add_n_carry(t, f0, m, f1, f1len, nlimbs);

        if (!squaring)
        {
            /* u = g0 + g1 */
            _mpn_poly_add_n_carry(u, g0, m, g1, g1len, nlimbs);
            /* v = (f0 + f1) (g0 + g1) */
            _mpn_poly_mul_karatsuba(v, t, tlen, u, ulen, nlimbs + 1, slimbs, cutoff, -1);
        }
        else 
        {
            _mpn_poly_mul_karatsuba(v, t, tlen, t, tlen, nlimbs + 1, slimbs, cutoff, -1);
        }
    }
    else
    {
        /* t = f0 + f1 */
        _mpn_poly_add_n(t, f0, m, f1, f1len, nlimbs);

        if (!squaring)
        {
            /* u = g0 + g1 */
            _mpn_poly_add_n(u, g0, m, g1, g1len, nlimbs);
            /* v = (f0 + f1) (g0 + g1) */
            _mpn_poly_mul_karatsuba(v, t, tlen, u, ulen, nlimbs, slimbs, cutoff, norm - 1);
        }
        else
        {
            _mpn_poly_mul_karatsuba(v, t, tlen, t, tlen, nlimbs, slimbs, cutoff, norm - 1);
        }
    }

    /* v -= f0 g0 */
    mpn_sub_n(v, v, res, (2 * m - 1) * slimbs);
    /* v -= f1 g1 */
    mpn_sub_n(v, v, res + 2 * m * slimbs, (f1len + g1len - 1) * slimbs);

    mpn_add_n(res + m * slimbs, res + m * slimbs, v, vlen * slimbs);

    TMP_END;
}

int
_mpn_mod_poly_mullow_karatsuba(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, slong cutoff, gr_ctx_t ctx)
{
    mp_ptr t;
    slong i, l;
    mp_size_t nlimbs, slimbs;
    mp_bitcnt_t sbits;
    int norm;
    TMP_INIT;
    TMP_START;

    norm = MPN_MOD_CTX_NORM(ctx);

    len1 = FLINT_MIN(len1, len);
    len2 = FLINT_MIN(len2, len);

    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    sbits = MPN_MOD_CTX_MODULUS_BITS(ctx);

    if (cutoff == -1)
    {
        if (poly1 == poly2 && len1 == len2)
        {
            if (sbits <= 128)
                cutoff = 32;
            else if (sbits <= 154)
                cutoff = 24;
            else if (sbits <= 448)
                cutoff = 12;
            else if (sbits <= 600)
                cutoff = 6;
            else
                cutoff = 4;
        }
        else
        {
            if (sbits <= 155)
                cutoff = 12;
            else if (sbits <= 320)
                cutoff = 8;
            else if (sbits <= 560)
                cutoff = 6;
            else
                cutoff = 4;
        }

        cutoff = FLINT_MIN(cutoff, len1);
        cutoff = FLINT_MIN(cutoff, len2);
    }

    /* b = FLINT_BIT_COUNT(FLINT_MIN(len1, len2)) larger than the bound for
       classical multiplication, since we need padding for up to b
       recursions */
    sbits = 2 * sbits + 2 * FLINT_BIT_COUNT(FLINT_MIN(len1, len2));
    slimbs = (sbits + FLINT_BITS - 1) / FLINT_BITS;

    t = TMP_ALLOC(sizeof(mp_limb_t) * slimbs * (len1 + len2 - 1));

    _mpn_poly_mul_karatsuba(t, poly1, len1, poly2, len2, nlimbs, slimbs, cutoff, norm);

    /* We support len < len1 + len2 - 1 but currently only save time
       on reductions. TODO: fix this. */
    for (i = 0; i < len; i++)
    {
        l = slimbs;
        MPN_NORM(t + i * slimbs, l);
        mpn_mod_set_mpn(res + i * nlimbs, t + i * slimbs, l, ctx);
    }

    TMP_END;
    return GR_SUCCESS;
}
