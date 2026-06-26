/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2014, 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "padic_radix.h"
#include "gr.h"

slong
_padic_radix_exp_bound(slong v, slong N, ulong p)
{
    ulong p1 = p - 1;
    ulong Nhi, Nlo, vhi, vlo, D, q = 0, r = 0;
    int use_fast = 0;

    umul_ppmm(vhi, vlo, (ulong) v, p1);          /* v*(p-1) */

    if (vhi == 0 && vlo > 1)                      /* D = v*(p-1) - 1 >= 1 */
    {
        D = vlo - 1;
        umul_ppmm(Nhi, Nlo, (ulong) N, p1);       /* N*(p-1) */

        if (Nhi == 0)
        {
            ulong A = Nlo - 1;                    /* N>=1, p>=2 => Nlo>=1 */
            q = A / D;
            r = A - q * D;
            use_fast = 1;
        }
        else
        {
            ulong Ahi = Nhi - (Nlo == 0);         /* (Nhi:Nlo) - 1 */
            ulong Alo = Nlo - 1;
            if (Ahi < D)                          /* ensures the quotient is a limb */
            {
                udiv_qrnnd(q, r, Ahi, Alo, D);
                use_fast = 1;
            }
        }
    }

    if (use_fast)
        return FLINT_MAX((slong) (q + (r != 0)), 1);

    {
        fmpz_t num, den, qq;
        slong n;

        fmpz_init(num);
        fmpz_init(den);
        fmpz_init(qq);

        fmpz_set_ui(qq, p - 1);
        fmpz_mul_si(num, qq, N);
        fmpz_sub_ui(num, num, 1);
        fmpz_mul_si(den, qq, v);
        fmpz_sub_ui(den, den, 1);
        fmpz_cdiv_q(qq, num, den);
        n = fmpz_get_si(qq);

        fmpz_clear(num);
        fmpz_clear(den);
        fmpz_clear(qq);

        return FLINT_MAX(n, 1);
    }
}

void
_padic_radix_exp(radix_integer_t rop, const radix_integer_t u,
    slong v, slong N, const radix_t radix)
{
    slong cutoff;
    ulong pbits;

    pbits = NMOD_BITS(radix->b);

    if (pbits <= 21)
        cutoff = 200;
    else if (pbits <= 32)
        cutoff = 150;
    else
        cutoff = 80;

    if (N < cutoff)
        _padic_radix_exp_rectangular(rop, u, v, N, radix);
    else
        _padic_radix_exp_balanced(rop, u, v, N, radix);
}

static int
_padic_radix_exp_wrapper(padic_radix_t res, const padic_radix_t x, int algorithm, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    ulong p = GR_PADIC_RADIX_CTX(ctx)->p;
    slong vx = x->v, Nx = x->N;

    slong prec_abs = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong prec_rel = PADIC_RADIX_CTX_PREC_REL(ctx);
    slong prec = FLINT_MIN(prec_rel, prec_abs);

    /* 0            -> 1 */
    /* 0 + O(x)     -> unable (p >= 3) */
    /* 0 + O(x^2)   -> unable (p == 2) */
    /* 0 + O(x^3)   -> 1 + O(x^5) */

    if (radix_integer_is_zero(&x->u, radix))
    {
        if (Nx == PADIC_RADIX_EXACT)
            return padic_radix_one(res, ctx);

        if (Nx < ((p == 2) ? 2 : 1) || prec == PADIC_RADIX_PREC_INF)
            return GR_UNABLE;

        radix_integer_one(&res->u, radix);
        res->v = 0;
        res->N = Nx;
        return _padic_radix_finalize(res, ctx);
    }

    if (vx < ((p == 2) ? 2 : 1))
        return GR_DOMAIN;

    if (prec == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;

    if (Nx != PADIC_RADIX_EXACT)
        prec = FLINT_MIN(prec, Nx);

    if (algorithm == 0)
        _padic_radix_exp(&res->u, &x->u, vx, prec, radix);
    else if (algorithm == 1)
        _padic_radix_exp_rectangular(&res->u, &x->u, vx, prec, radix);
    else
        _padic_radix_exp_balanced(&res->u, &x->u, vx, prec, radix);

    res->v = 0;
    res->N = prec;
    return _padic_radix_finalize(res, ctx);
}

int
padic_radix_exp(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    return _padic_radix_exp_wrapper(res, x, 0, ctx);
}

int
padic_radix_exp_rectangular(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    return _padic_radix_exp_wrapper(res, x, 1, ctx);
}

int
padic_radix_exp_balanced(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    return _padic_radix_exp_wrapper(res, x, 2, ctx);
}

