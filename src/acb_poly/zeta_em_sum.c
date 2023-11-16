/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_dirichlet.h"

/* res = src * (c + x) */
void _acb_poly_mullow_cpx(acb_ptr res, acb_srcptr src, slong len, const acb_t c, slong trunc, slong prec)
{
    slong i;

    if (len < trunc)
        acb_set(res + len, src + len - 1);

    for (i = len - 1; i > 0; i--)
    {
        acb_mul(res + i, src + i, c, prec);
        acb_add(res + i, res + i, src + i - 1, prec);
    }

    acb_mul(res, src, c, prec);
}

/* todo: don't hardcode this */
#define SIEVE_ALLOC_LIMIT 4e9   /* 4 GB */

void
_acb_poly_zeta_em_sum(acb_ptr z, const acb_t s, const acb_t a, int deflate, ulong N, ulong M, slong d, slong prec)
{
    acb_ptr t, u, v, sum;
    acb_t Na, one;
    slong i;

    t = _acb_vec_init(d + 1);
    u = _acb_vec_init(d);
    v = _acb_vec_init(d);
    sum = _acb_vec_init(d);
    acb_init(Na);
    acb_init(one);

    prec += 2 * (FLINT_BIT_COUNT(N) + FLINT_BIT_COUNT(d));
    acb_one(one);

    /* sum 1/(k+a)^(s+x) */
    if (acb_is_one(a) && d <= 2 && _acb_vec_estimate_allocated_bytes(d * N / 6, prec) < SIEVE_ALLOC_LIMIT)
        acb_dirichlet_powsum_sieved(sum, s, N, d, prec);
    else if (acb_is_one(a) && d <= 4) /* todo: also better for slightly larger d, if N and prec large enough */
        acb_dirichlet_powsum_smooth(sum, s, N, d, prec);
    else if (N > 50 && flint_get_num_threads() > 1)
        _acb_poly_powsum_series_naive_threaded(sum, s, a, one, N, d, prec);
    else
        _acb_poly_powsum_series_naive(sum, s, a, one, N, d, prec);

    /* t = 1/(N+a)^(s+x); we might need one extra term for deflation */
    acb_add_ui(Na, a, N, prec);
    _acb_poly_acb_invpow_cpx(t, Na, s, d + 1, prec);

    /* sum += (N+a) * 1/((s+x)-1) * t */
    if (!deflate)
    {
        /* u = (N+a)^(1-(s+x)) */
        acb_sub_ui(v, s, 1, prec);
        _acb_poly_acb_invpow_cpx(u, Na, v, d, prec);

        /* divide by 1/((s-1) + x) */
        acb_sub_ui(v, s, 1, prec);
        acb_div(u, u, v, prec);

        for (i = 1; i < d; i++)
        {
            acb_sub(u + i, u + i, u + i - 1, prec);
            acb_div(u + i, u + i, v, prec);
        }

        _acb_vec_add(sum, sum, u, d, prec);
    }
    /* sum += ((N+a)^(1-(s+x)) - 1) / ((s+x) - 1) */
    else
    {
        /* at s = 1, this becomes (N*t - 1)/x, i.e. just remove one coeff  */
        if (acb_is_one(s))
        {
            for (i = 0; i < d; i++)
                acb_mul(u + i, t + i + 1, Na, prec);
            _acb_vec_add(sum, sum, u, d, prec);
        }
        else
        {
            /* TODO: this is numerically unstable for large derivatives,
                and divides by zero if s contains 1. We want a good
                way to evaluate the power series ((N+a)^y - 1) / y where y has
                nonzero constant term, without doing a division.
                How is this best done? */

            _acb_vec_scalar_mul(t, t, d, Na, prec);
            acb_sub_ui(t + 0, t + 0, 1, prec);
            acb_sub_ui(u + 0, s, 1, prec);
            acb_inv(u + 0, u + 0, prec);
            for (i = 1; i < d; i++)
                acb_mul(u + i, u + i - 1, u + 0, prec);
            for (i = 1; i < d; i += 2)
                acb_neg(u + i, u + i);
            _acb_poly_mullow(v, u, d, t, d, d, prec);
            _acb_vec_add(sum, sum, v, d, prec);
            _acb_poly_acb_invpow_cpx(t, Na, s, d, prec);
        }
    }

    /* sum += u = 1/2 * t */
    _acb_vec_scalar_mul_2exp_si(u, t, d, -WORD(1));
    _acb_vec_add(sum, sum, u, d, prec);

    /* Euler-Maclaurin formula tail */
    if (d < 5 || d < M / 10)
        _acb_poly_zeta_em_tail_naive(u, s, Na, t, M, d, prec);
    else
        _acb_poly_zeta_em_tail_bsplit(u, s, Na, t, M, d, prec);

    _acb_vec_add(z, sum, u, d, prec);

    _acb_vec_clear(t, d + 1);
    _acb_vec_clear(u, d);
    _acb_vec_clear(v, d);
    _acb_vec_clear(sum, d);
    acb_clear(Na);
    acb_clear(one);
}
