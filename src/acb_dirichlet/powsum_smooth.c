/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "acb_poly.h"

/* bound number of 5-smooth numbers up to n; see http://oeis.org/A188425 */
/* this does not need to be tight */
static slong
smooth_bound(ulong n)
{
    if (n <= 256) return 52;
    if (n <= 65536) return 284;
    if (n <= 16777216) return 836;
    return 13283;  /* ok up to 2^64 */
}

static ulong smul(ulong x, ulong y)
{
    ulong hi, lo;
    umul_ppmm(hi, lo, x, y);
    if (hi)
        return UWORD_MAX;
    else
        return lo;
}

static slong _index(const ulong * v, ulong c)
{
    slong i;
    for (i = 0; ; i++)
        if (v[i] == c)
            return i;
}

void
acb_dirichlet_powsum_smooth(acb_ptr res, const acb_t s, ulong N, slong d, slong prec)
{
    ulong * smooth;     /* numbers 2^i 3^j 5^k <= N */
    slong num_smooth;   /* the number of such numbers */
    acb_ptr sums;       /* partial sum for each smooth prefactor */
    acb_ptr powers;     /* (3^j 5^k)^(-s) */
    acb_ptr t;          /* temporary */
    arb_t log_n;
    ulong x2, x3, x5, m, n, nprev;
    slong i2, i3, i5, i, j, iu;
    int critical_line, integer;

    if (N <= 1)
    {
        acb_set_ui(res, N);
        _acb_vec_zero(res + 1, d - 1);
        return;
    }

    if (N >= UWORD_MAX - 2)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    critical_line = arb_is_exact(acb_realref(s)) &&
        (arf_cmp_2exp_si(arb_midref(acb_realref(s)), -1) == 0);

    integer = arb_is_zero(acb_imagref(s)) && arb_is_int(acb_realref(s));

    /* generate the smooth numbers */
    smooth = flint_malloc(smooth_bound(N) * sizeof(ulong));
    smooth[0] = 1;
    num_smooth = 1;
    x2 = 2; x3 = 3; x5 = 5;
    i2 = i3 = i5 = 0;

    while ((m = FLINT_MIN(FLINT_MIN(x2, x3), x5)) <= N)
    {
        smooth[num_smooth++] = m;
        if (m == x2) x2 = smul(2, smooth[++i2]);
        if (m == x3) x3 = smul(3, smooth[++i3]);
        if (m == x5) x5 = smul(5, smooth[++i5]);
    }

    sums = _acb_vec_init(num_smooth * d);
    powers = _acb_vec_init(num_smooth * d);
    t = _acb_vec_init(d);
    arb_init(log_n);

    arb_zero(log_n);
    nprev = 1;

    /* add 1^-s */
    for (i = 0; i < num_smooth; i++)
        acb_one(sums + i * d);

    /* compute all the non-smooth index terms (bulk of the work) */
    for (n = 7; n <= N; n += 2)
    {
        if ((n % 3 != 0) && (n % 5 != 0))
        {
            acb_dirichlet_powsum_term(t, log_n, &nprev, s, n, integer, critical_line, d, prec);
            _acb_vec_add(sums, sums, t, d, prec);

            for (i = 1; i < num_smooth && (smooth[i] <= (N / n)); i++)
                _acb_vec_add(sums + i * d, sums + i * d, t, d, prec);
        }
    }

    /* compute 2^(-s) and powers (3^j 3^k)^(-s) */
    arb_zero(log_n);
    nprev = 1;

    for (i = 1; i < num_smooth; i++)
    {
        n = smooth[i];

        if (n == 2)
        {
            acb_dirichlet_powsum_term(powers + i * d, log_n, &nprev, s,
                n, integer, critical_line, d, prec);
        }
        else if (n % 2 != 0)
        {
            if (n <= 5)
            {
                acb_dirichlet_powsum_term(powers + i * d, log_n, &nprev, s,
                    n, integer, critical_line, d, prec);
            }
            else if (n % 3 == 0)
            {
                i3 = _index(smooth, 3);
                iu = _index(smooth, n / 3);
                _acb_poly_mullow(powers + i * d,
                    powers + i3 * d, d, powers + iu * d, d, d, prec);
            }
            else
            {
                i5 = _index(smooth, 5);
                iu = _index(smooth, n / 5);
                _acb_poly_mullow(powers + i * d,
                    powers + i5 * d, d, powers + iu * d, d, d, prec);
            }
        }
    }

    /* merge the sums into the power-of-two sums */
    for (i = 0; i < num_smooth; i++)
    {
        ulong u, v;

        m = smooth[i];
        u = m;
        v = 0;

        while ((u & 1) == 0)
        {
            u >>= 1;
            v++;
        }

        if ((UWORD(1) << v) != m)
        {
            j = _index(smooth, UWORD(1) << v);
            iu = _index(smooth, u);

            if (u == 1)
            {
                _acb_vec_add(sums + j * d, sums + j * d, sums + i * d, d, prec);
            }
            else
            {
                _acb_poly_mullow(t, sums + i * d, d, powers + iu * d, d, d, prec);
                _acb_vec_add(sums + j * d, sums + j * d, t, d, prec);
            }
        }
    }

    /* finally evaluate with respect to powers of 2 using horner */
    _acb_vec_zero(res, d);
    i2 = _index(smooth, 2);

    for (i = num_smooth - 1; i >= 0; i--)
    {
        n = smooth[i];

        if ((n & (n - 1)) == 0)
        {
            _acb_poly_mullow(t, powers + i2 * d, d, res, d, d, prec);
            _acb_vec_add(res, sums + i * d, t, d, prec);
        }
    }

    _acb_vec_clear(sums, num_smooth * d);
    _acb_vec_clear(powers, num_smooth * d);
    _acb_vec_clear(t, d);
    arb_clear(log_n);
    flint_free(smooth);
}

