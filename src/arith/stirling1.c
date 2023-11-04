/*
    Copyright (C) 2010, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "arith.h"

/* compute single coefficient in polynomial product */
static void
_fmpz_poly_mulmid_single(fmpz_t res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong i)
{
    slong j, top1, top2;

    top1 = FLINT_MIN(len1 - 1, i);
    top2 = FLINT_MIN(len2 - 1, i);

    fmpz_mul(res, poly1 + i - top2, poly2 + top2);

    for (j = 1; j < top1 + top2 - i + 1; j++)
        fmpz_addmul(res, poly1 + i - top2 + j, poly2 + top2 - j);
}

#define MAX_BASECASE 16

/*
which == 1: expand (a+x)...((b-1)+x) truncated to length len
which == 2: expand (1+ax)...(1+(b-1)x) truncated to length len

final: only extract final coefficient
*/
static void
stirling_1u_ogf_bsplit(fmpz * res, ulong a, ulong b, slong len, int which, int final)
{
    ulong c, n, cbc;

    len = FLINT_MIN(len, b - a + 1);

    /* (c+x)^n has coefficients bounded by max(c,n)^n */
    n = b - a;
    c = FLINT_MAX(b - 1, n);
    cbc = FLINT_BIT_COUNT(c);

    if (n == 1 || (len <= MAX_BASECASE && n * cbc <= FLINT_BITS))
    {
        mp_limb_t v[MAX_BASECASE];
        slong i, j;

        if (which == 1)
        {
            v[0] = a;
            v[1] = 1;

            /* multiply by ((a+i) + x) */
            for (i = 1; i < n; i++)
            {
                if (i + 1 < len)
                    v[i + 1] = 1;
                for (j = FLINT_MIN(i, len - 1); j >= 1; j--)
                    v[j] = (a + i) * v[j] + v[j - 1];
                v[0] *= (a + i);
            }
        }
        else
        {
            v[0] = 1;
            v[1] = a;

            /* multiply by (1 + (a+i) x) */
            for (i = 1; i < n; i++)
            {
                if (i + 1 < len)
                    v[i + 1] = v[i] * (a + i);

                for (j = FLINT_MIN(i, len - 1); j >= 1; j--)
                    v[j] = v[j] + (a + i) * v[j - 1];
            }
        }

        if (final)
            fmpz_set_ui(res, v[len - 1]);
        else
            for (i = 0; i < len; i++)
                fmpz_set_ui(res + i, v[i]);
    }
    else
    {
        fmpz *L, *R;
        slong len1, len2;
        slong m = a + (b - a) / 2;

        len1 = FLINT_MIN(m - a + 1, len);
        len2 = FLINT_MIN(b - m + 1, len);

        L = _fmpz_vec_init(len1 + len2);
        R = L + len1;

        stirling_1u_ogf_bsplit(L, a, m, len, which, 0);
        /* Remark: R can be computed from L using a Taylor shift. In theory,
           this improves the complexity from O(M(n) log n) to O(M(n)).
           In practice, the Taylor shift is rarely much faster than computing
           R from scratch, and sometimes slower, so we do not bother. */
        stirling_1u_ogf_bsplit(R, m, b, len, which, 0);

        if (final)
            _fmpz_poly_mulmid_single(res, L, len1, R, len2, len - 1);
        else
            _fmpz_poly_mullow(res, R, len2, L, len1, FLINT_MIN(len, len1 + len2 - 1));

        _fmpz_vec_clear(L, len1 + len2);
    }
}

static void
stirling_1u_egf(fmpz_t res, ulong n, ulong k)
{
    fmpz * num;
    fmpz * rnum;
    fmpz_t den;
    fmpz_t rden;
    slong i, len;

    if (k >= n || k == 0)
    {
        fmpz_set_ui(res, n == k);
        return;
    }

    len = n - k + 1;

    num = _fmpz_vec_init(len + 1);
    rnum = _fmpz_vec_init(len);
    fmpz_init(den);
    fmpz_init(rden);

    fmpz_one(den);
    for (i = 0; i < len; i++)
        fmpz_one(num + i);

    _fmpq_poly_integral(num, den, num, den, len + 1);
    for (i = 0; i < len; i++)
        fmpz_swap(num + i, num + i + 1);

    _fmpq_poly_pow_trunc(rnum, rden, num, den, len, k, len);

    fmpz_set_ui(num, k);
    fmpz_add_ui(num, num, 1);
    fmpz_rfac_ui(num, num, n - k);

    fmpz_mul(num, num, rnum + n - k);
    fmpz_divexact(res, num, rden);

    _fmpz_vec_clear(num, len + 1);
    _fmpz_vec_clear(rnum, len);
    fmpz_clear(den);
    fmpz_clear(rden);
}


void
arith_stirling_number_1u(fmpz_t res, ulong n, ulong k)
{
    if (k >= n || k == 0)
    {
        fmpz_set_ui(res, n == k);
    }
    else if (k == 1)
    {
        fmpz_fac_ui(res, n - 1);
    }
    else if (n > 140 && k > 0.87 * n)
    {
        stirling_1u_egf(res, n, k);
    }
    else
    {
        if (k < n / 2)
            stirling_1u_ogf_bsplit(res, 1, n, k, 1, 1);
        else
            stirling_1u_ogf_bsplit(res, 1, n, n - k + 1, 2, 1);
    }
}

void
arith_stirling_number_1u_vec(fmpz * res, ulong n, slong klen)
{
    slong k, len;

    if (klen <= 0)
        return;

    len = FLINT_MIN(klen - 1, n - 1);

    if (n >= 1 && len >= 1)
        stirling_1u_ogf_bsplit(res + 1, 1, n, len, 1, 0);

    fmpz_set_ui(res + 0, n == 0);
    for (k = n; k < klen; k++)
        fmpz_set_ui(res + k, n == k);
}

void
arith_stirling_number_1(fmpz_t s, ulong n, ulong k)
{
    arith_stirling_number_1u(s, n, k);
    if ((n + k) % 2)
        fmpz_neg(s, s);
}

void
arith_stirling_number_1_vec(fmpz * row, ulong n, slong klen)
{
    slong k;

    arith_stirling_number_1u_vec(row, n, klen);

    for (k = (n + 1) % 2; k < klen; k += 2)
        fmpz_neg(row + k, row + k);
}

