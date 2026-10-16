/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"

/* {r, s} -= ({w, wn} << shift) mod 2^(FLINT_BITS s). Needs s limbs of scratch space t. */
static void
_mpn_sub_shifted(nn_ptr r, slong s, nn_srcptr w, slong wn, ulong shift, nn_ptr t)
{
    slong q = shift / FLINT_BITS;
    slong rb = shift % FLINT_BITS;
    slong m, m2, j;

    if (q >= s)
        return;

    m = s - q;
    m2 = FLINT_MIN(m, wn);

    if (rb == 0)
    {
        for (j = 0; j < m2; j++)
            t[j] = w[j];
        for (j = m2; j < m; j++)
            t[j] = 0;
    }
    else
    {
        ulong cy = mpn_lshift(t, w, m2, rb);
        if (m2 < m)
        {
            t[m2] = cy;
            for (j = m2 + 1; j < m; j++)
                t[j] = 0;
        }
    }

    mpn_sub_n(r + q, r + q, t, m);
}

/* {r, s} += (c << shift) mod 2^(FLINT_BITS s) */
static void
_mpn_add_1_shifted(nn_ptr r, slong s, ulong c, ulong shift)
{
    slong q = shift / FLINT_BITS;
    slong rb = shift % FLINT_BITS;

    if (q >= s)
        return;

    if (rb == 0)
    {
        mpn_add_1(r + q, r + q, s - q, c);
    }
    else
    {
        ulong lo = c << rb;
        ulong hi = c >> (FLINT_BITS - rb);
        ulong cy;
        add_ssaaaa(cy, r[q], 0, r[q], 0, lo);
        if (q + 1 < s)
            mpn_add_1(r + q + 1, r + q + 1, s - q - 1, hi + cy);
    }
}

/* Sliding window sum with n + 1 limbs. */
FLINT_FORCE_INLINE void
_window_add(nn_ptr w, nn_srcptr x, slong n)
{
    w[n] += mpn_add_n(w, w, x, n);
}

FLINT_FORCE_INLINE void
_window_sub(nn_ptr w, nn_srcptr x, slong n)
{
    w[n] -= mpn_sub_n(w, w, x, n);
}

FLINT_FORCE_INLINE void
_set_output(fmpz * res, nn_srcptr s, slong slimbs, int method)
{
    if (method == FLINT_MPN_POLY_MUL_UNSIGNED)
    {
        if (slimbs == 1)
            fmpz_set_ui(res, s[0]);
        else if (slimbs == 2)
            fmpz_set_uiui(res, s[1], s[0]);
        else
            fmpz_set_ui_array(res, s, slimbs);
    }
    else
    {
        if (slimbs == 1)
            fmpz_set_si(res, s[0]);
        else if (slimbs == 2)
            fmpz_set_signed_uiui(res, s[1], s[0]);
        else
            fmpz_set_signed_ui_array(res, s, slimbs);
    }
}

void
_fmpz_poly_mpn_set_outputs(fmpz * res, nn_ptr U, slong s, slong nlo, slong nhi, int method,
    nn_srcptr a, slong len1, slong n1, nn_srcptr b, slong len2, slong n2,
    slong ub1, slong ub2, int sgn1, int sgn2, int squaring)
{
    slong i, cnt, top1, top2, jlo, jhi;
    nn_ptr wa, wb, t, sbuf;
    TMP_INIT;

    if (method != FLINT_MPN_POLY_MUL_BIAS)
    {
        for (i = nlo; i < nhi; i++)
            _set_output(res + i - nlo, U + (i - nlo) * s, s, method);
        return;
    }

    TMP_START;
    wa = TMP_ALLOC(sizeof(ulong) * (n1 + 1 + n2 + 1 + s));
    wb = wa + n1 + 1;
    t = wb + n2 + 1;

    /* initial window sums for i = nlo */
    flint_mpn_zero(wa, n1 + 1);
    flint_mpn_zero(wb, n2 + 1);

    if (sgn2 || squaring)
    {
        jlo = FLINT_MAX(0, nlo - len2 + 1);
        jhi = FLINT_MIN(len1 - 1, nlo);
        for (i = jlo; i <= jhi; i++)
            _window_add(wa, a + i * n1, n1);
    }

    if (sgn1 && !squaring)
    {
        jlo = FLINT_MAX(0, nlo - len1 + 1);
        jhi = FLINT_MIN(len2 - 1, nlo);
        for (i = jlo; i <= jhi; i++)
            _window_add(wb, b + i * n2, n2);
    }

    for (i = nlo; i < nhi; i++)
    {
        top1 = FLINT_MIN(len1 - 1, i);
        top2 = FLINT_MIN(len2 - 1, i);
        cnt = top1 + top2 - i + 1;
        sbuf = U + (i - nlo) * s;

        /* r = U - 2^ub2 Wa - 2^ub1 Wb + 2^(ub1+ub2) cnt */
        if (squaring)
        {
            _mpn_sub_shifted(sbuf, s, wa, n1 + 1, ub1 + 1, t);
            _mpn_add_1_shifted(sbuf, s, cnt, 2 * ub1);
        }
        else
        {
            if (sgn2)
                _mpn_sub_shifted(sbuf, s, wa, n1 + 1, ub2, t);
            if (sgn1)
                _mpn_sub_shifted(sbuf, s, wb, n2 + 1, ub1, t);
            if (sgn1 && sgn2)
                _mpn_add_1_shifted(sbuf, s, cnt, ub1 + ub2);
        }

        /* slide the windows to i + 1 */
        if (sgn2 || squaring)
        {
            if (i + 1 <= len1 - 1)
                _window_add(wa, a + (i + 1) * n1, n1);
            if (i + 1 - len2 >= 0)
                _window_sub(wa, a + (i + 1 - len2) * n1, n1);
        }

        if (sgn1 && !squaring)
        {
            if (i + 1 <= len2 - 1)
                _window_add(wb, b + (i + 1) * n2, n2);
            if (i + 1 - len1 >= 0)
                _window_sub(wb, b + (i + 1 - len1) * n2, n2);
        }

        _set_output(res + i - nlo, sbuf, s, method);
    }

    TMP_END;
}
