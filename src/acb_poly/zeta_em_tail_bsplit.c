/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "bernoulli.h"

static void
bsplit(acb_ptr P, acb_ptr T, const acb_t s, const acb_t Na,
    slong a, slong b, int cont, slong len, slong prec)
{
    slong plen = FLINT_MIN(2 * (b - a) + 1, len);

    if (b - a == 1)
    {
        slong j = a;

        fmpz_t t;
        fmpz_init(t);

        if (j == 0)
        {
            acb_set(P, s);
            if (len > 1) acb_one(P + 1);
            if (len > 2) acb_zero(P + 2);
        }
        else
        {
            /* P = ((s+x)+2j-1)(s+x)+2j) */

            /* constant term = s^2 + (4j-1)s + (4j^2-2j)*/
            acb_mul(P, s, s, prec);
            acb_addmul_ui(P, s, 4*j-1, prec);
            fmpz_set_ui(t, 2*j);
            fmpz_mul_ui(t, t, 2*j-1);
            acb_add_fmpz(P, P, t, prec);

            /* x term = (2s+4j-1) */
            if (len > 1)
            {
                acb_mul_ui(P + 1, s, 2, prec);
                acb_add_ui(P + 1, P + 1, 4*j-1, prec);
            }

            /* x^2 term = 1 */
            if (len > 2)
                acb_one(P + 2);

        }

        /* P = P / ((2j+1)*(2j+2)) */
        fmpz_set_ui(t, 2*j+1);
        fmpz_mul_ui(t, t, 2*j+2);
        _acb_vec_scalar_div_fmpz(P, P, plen, t, prec);

        /* P = P / Na  or  P / Na^2 */
        if (j == 0)
            acb_set(T, Na);
        else
            acb_mul(T, Na, Na, prec);
        _acb_vec_scalar_div(P, P, plen, T, prec);

        /* T = P * B_{2j+2} */
        _acb_vec_scalar_mul_fmpz(T, P, plen, fmpq_numref(bernoulli_cache + 2 * j + 2), prec);
        _acb_vec_scalar_div_fmpz(T, T, plen, fmpq_denref(bernoulli_cache + 2 * j + 2), prec);

        fmpz_clear(t);
    }
    else
    {
        slong m, len1, len2, alloc;
        acb_ptr P1, T1, P2, T2;

        m = a + (b - a) / 2;

        len1 = FLINT_MIN(2 * (m - a) + 1, len);
        len2 = FLINT_MIN(2 * (b - m) + 1, len);

        alloc = 2 * len1 + 2 * len2;

        P1 = _acb_vec_init(alloc);
        T1 = P1 + len1;
        P2 = T1 + len1;
        T2 = P2 + len2;

        bsplit(P1, T1, s, Na, a, m, 1, len, prec);
        bsplit(P2, T2, s, Na, m, b, 1, len, prec);

        /* P = P1 * P2 */
        if (cont)
            _acb_poly_mullow(P, P2, len2, P1, len1, plen, prec);

        /* T = T1 + P1 * T2 */
        _acb_poly_mullow(T, T2, len2, P1, len1, plen, prec);
        _acb_vec_add(T, T, T1, len1, prec);

        _acb_vec_clear(P1, alloc);
    }
}

void
_acb_poly_zeta_em_tail_bsplit(acb_ptr z, const acb_t s, const acb_t Na, acb_srcptr Nasx, slong M, slong len, slong prec)
{
    acb_ptr P, T;

    if (M < 1)
    {
        _acb_vec_zero(z, len);
        return;
    }

    BERNOULLI_ENSURE_CACHED(2 * M);

    P = _acb_vec_init(len);
    T = _acb_vec_init(len);

    bsplit(P, T, s, Na, 0, M, 0, len, prec);

    _acb_poly_mullow(z, T, len, Nasx, len, len, prec);

    _acb_vec_clear(P, len);
    _acb_vec_clear(T, len);
}

