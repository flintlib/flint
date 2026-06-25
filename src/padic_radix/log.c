/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "padic_radix.h"
#include "gr.h"

/*
    Smallest n such that for all i >= n, i v - ord_p(i) >= N, so that the tail
    of  sum_{i>=1} x^i / i  with ord_p(x) >= v is below p^N.  Word-sized p (a
    radix digit always fits a word); port of the corresponding branch of
    _padic_log_bound.
*/
slong
_padic_radix_log_bound(slong v, slong N, ulong p)
{
    slong b, c;

    c = N - (slong) n_flog((ulong) v, p);
    c = FLINT_MAX(c, 1);
    b = ((c + (slong) n_clog((ulong) c, p) + 1) + (v - 1)) / v;

    while (--b >= 2)
    {
        slong t = b * v - (slong) n_clog((ulong) b, p);

        if (t < N)
            return b + 1;
    }

    return 2;
}

void
_padic_radix_log(radix_integer_t rop, const radix_integer_t y, slong N,
    const radix_t radix)
{
    slong cutoff;
    ulong pbits = NMOD_BITS(radix->b);

    if (pbits <= 5)
        cutoff = 350;
    else if (pbits <= 9)
        cutoff = 280;
    else if (pbits <= 21)
        cutoff = 240;
    else if (pbits <= 32)
        cutoff = 160;
    else if (pbits <= 46)
        cutoff = 120;
    else
        cutoff = 80;

    if (N < cutoff)
        _padic_radix_log_rectangular(rop, y, N, radix);
    else
        _padic_radix_log_balanced(rop, y, N, radix);
}

static int
_padic_radix_log_wrapper(padic_radix_t res, const padic_radix_t x, int algorithm, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    ulong p = GR_PADIC_RADIX_CTX(ctx)->p;
    slong thr = (p == 2) ? 2 : 1;
    slong Nx = x->N;
    slong prec_abs = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong prec_rel = PADIC_RADIX_CTX_PREC_REL(ctx);
    slong prec = FLINT_MIN(prec_rel, prec_abs);
    radix_integer_t y, one;
    slong vy;

    /*
        op = u p^v with u a unit.  log converges only for a 1-unit, i.e.
        op == 1 (mod p), which forces v == 0.  Any v != 0 (op divisible by p,
        or non-integral) is outside the domain; handling it here also avoids
        materialising u p^v for a large valuation.
    */
    if (radix_integer_is_zero(&x->u, radix))
    {
        /*
            op == 0 (mod p^Nx): the valuation is only known to be >= Nx.  If
            Nx >= 1 the valuation is positive, so op is not a unit and log is
            undefined (this includes exact zero, Nx == EXACT).  If Nx <= 0
            (e.g. O(p^0) or O(p^{-1})) op may still be a 1-unit, so we cannot
            decide whether the series converges.
        */
        if (Nx <= 0)
            return GR_UNABLE;
        return GR_DOMAIN;
    }

    if (x->v != 0)
        return GR_DOMAIN;

    /* log(1) = 0 exactly (op == 1 to its stored precision). */
    if (radix_integer_is_one(&x->u, radix))
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = Nx;
        return _padic_radix_finalize(res, ctx);
    }

    /* Resolve the working absolute precision. */
    if (prec == PADIC_RADIX_PREC_INF)
    {
        if (Nx == PADIC_RADIX_EXACT)
            return GR_UNABLE;
        prec = Nx;
    }
    else if (Nx != PADIC_RADIX_EXACT)
    {
        prec = FLINT_MIN(prec, Nx);
    }

    if (prec > PADIC_RADIX_ERR_MAX)
        prec = PADIC_RADIX_ERR_MAX;

    if (prec <= 0)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = prec;
        return _padic_radix_finalize(res, ctx);
    }

    /* y = 1 - op = 1 - u  (mod p^prec); op = u since v == 0 */
    radix_integer_init(y, radix);
    radix_integer_init(one, radix);
    radix_integer_one(one, radix);

    radix_integer_mod_digits(y, &x->u, prec, radix);      /* u mod p^prec */
    radix_integer_sub(y, one, y, radix);
    radix_integer_mod_digits(y, y, prec, radix);          /* nonnegative residue */

    if (radix_integer_is_zero(y, radix))                  /* op == 1 (mod p^prec) */
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = prec;
        radix_integer_clear(y, radix);
        radix_integer_clear(one, radix);
        return _padic_radix_finalize(res, ctx);
    }

    vy = radix_integer_valuation_digits(y, radix);

    if (vy < thr)                                         /* does not converge */
    {
        radix_integer_clear(y, radix);
        radix_integer_clear(one, radix);
        return GR_DOMAIN;
    }

    if (vy >= prec)                                       /* log == 0 (mod p^prec) */
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = prec;
        radix_integer_clear(y, radix);
        radix_integer_clear(one, radix);
        return _padic_radix_finalize(res, ctx);
    }

    if (algorithm == 0)
        _padic_radix_log(&res->u, y, prec, radix);
    else if (algorithm == 1)
        _padic_radix_log_rectangular(&res->u, y, prec, radix);
    else
        _padic_radix_log_balanced(&res->u, y, prec, radix);

    radix_integer_clear(y, radix);
    radix_integer_clear(one, radix);

    res->v = 0;
    res->N = prec;
    return _padic_radix_finalize(res, ctx);
}

int
padic_radix_log(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    return _padic_radix_log_wrapper(res, x, 0, ctx);
}

int
padic_radix_log_balanced(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    return _padic_radix_log_wrapper(res, x, 2, ctx);
}

int
padic_radix_log_rectangular(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    return _padic_radix_log_wrapper(res, x, 1, ctx);
}
