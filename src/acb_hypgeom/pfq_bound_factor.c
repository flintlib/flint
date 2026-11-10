/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "acb.h"
#include "acb_hypgeom.h"

/*
    Let T(k) denote the terms of the series and r(k) = |T(k+1) / T(k)|
    = |z| prod |a_i + k| / prod |b_i + k|, where
    |a_i + k| / |b_i + k| <= 1 + |a_i - b_i| / |b_i + k|. We compute C
    with sum_{k >= n} |T(k)| <= C |T(n)|.

    If Re(b_i + n) > 0 for all i, the lower bounds for |b_i + k| are
    increasing in k, so r(k) <= D = r(n) for k >= n and C = 1 / (1 - D).

    Otherwise, some b_i is close to (or past) a negative integer -m_i with
    m_i >= n, the terms can grow again near k = m_i, and there is no
    uniform D < 1. We split [n, inf) into ranges on each of which every
    |b_i + k| is monotone (or which consist of the single index m_i),
    bound r(k) by D_j on the j-th range [k_j, k_{j+1}), and sum

        sum_{k >= n} |T(k)| / |T(n)| <= sum_j P_j S_j,

    where P_j = prod_{l < j} D_l^(k_{l+1} - k_l) bounds |T(k_j)| / |T(n)|
    and S_j bounds sum_{t < k_{j+1} - k_j} D_j^t. The ranges are refined
    geometrically on both sides of each m_i, so that only O(log m_i) of
    them are needed.
*/

/* lower bound for the distance from b to the integers (at least |Im(b)|) */
static void
_lower_bound_dist_int(mag_t u, const acb_t b)
{
    const arf_struct * x = arb_midref(acb_realref(b));

    arb_get_mag_lower(u, acb_imagref(b));

    if (arf_is_finite(x) && arf_cmpabs_2exp_si(x, FLINT_BITS - 2) < 0)
    {
        mag_t v;
        arb_t t;
        fmpz_t m;

        mag_init(v);
        arb_init(t);
        fmpz_init(m);

        /* Re(b) - round(Re(b)), exactly */
        arf_get_fmpz(m, x, ARF_RND_NEAR);
        arf_sub_fmpz(arb_midref(t), x, m, ARF_PREC_EXACT, ARF_RND_DOWN);
        mag_set(arb_radref(t), arb_radref(acb_realref(b)));
        arb_get_mag_lower(v, t);
        mag_max(u, u, v);

        mag_clear(v);
        arb_clear(t);
        fmpz_clear(m);
    }
}

/* lower bound for |b + k|, uniformly for integers k in [k0, k1), where
   k1 = 0 means infinity */
static void
_lower_bound_range(mag_t u, const acb_t b, ulong k0, ulong k1)
{
    acb_t w;

    acb_init(w);

    acb_add_ui(w, b, k0, MAG_BITS);

    if (arb_is_positive(acb_realref(w)))
    {
        /* increasing on the range */
        acb_get_mag_lower(u, w);
    }
    else
    {
        if (k1 != 0)
            acb_add_ui(w, b, k1 - 1, MAG_BITS);

        if (k1 != 0 && arb_is_negative(acb_realref(w)))
            acb_get_mag_lower(u, w);   /* decreasing on the range */
        else
            _lower_bound_dist_int(u, b);
    }

    acb_clear(w);
}

/* uniform bound: for k >= n, |b + k| >= |b + n| if Re(b + n) > 0; otherwise
   |b + k| >= u = max(1/2, |Im(b)|) for all k except at most one index,
   where |b + k| >= d = the distance from b to the integers; the ratio at
   that index is larger by at most a factor u / d */
static void
_bound_factor_uniform(mag_t C, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const mag_t zmag, ulong n)
{
    mag_t t, u, d, extra;
    acb_t w;
    slong i;

    mag_init(t);
    mag_init(u);
    mag_init(d);
    mag_init(extra);
    acb_init(w);

    mag_set(C, zmag);
    mag_one(extra);

    for (i = 0; i < q; i++)
    {
        acb_add_ui(w, b + i, n, MAG_BITS);

        if (arb_is_positive(acb_realref(w)))
        {
            acb_get_mag_lower(u, w);
        }
        else
        {
            arb_get_mag_lower(u, acb_imagref(b + i));
            mag_one(t);
            mag_mul_2exp_si(t, t, -1);
            mag_max(u, u, t);

            _lower_bound_dist_int(d, b + i);

            if (mag_cmp(d, u) < 0)
            {
                mag_div(t, u, d);
                mag_mul(extra, extra, t);
            }
        }

        if (i < p)
        {
            acb_sub(w, a + i, b + i, MAG_BITS);
            acb_get_mag(t, w);
            mag_div(t, t, u);
            mag_one(u);
            mag_add(t, t, u);
            mag_mul(C, C, t);
        }
        else
        {
            mag_div(C, C, u);
        }
    }

    mag_one(t);
    mag_sub_lower(u, t, C);

    if (mag_is_zero(u))
        mag_inf(C);
    else
    {
        mag_div(C, t, u);
        mag_mul(C, C, extra);
    }

    mag_clear(t);
    mag_clear(u);
    mag_clear(d);
    mag_clear(extra);
    acb_clear(w);
}

/* bound for the term ratio on [k0, k1), where k1 = 0 means infinity */
static void
_bound_factor_ranges_D(mag_t D, acb_srcptr a, slong p, acb_srcptr b,
    slong q, const mag_t zmag, ulong k0, ulong k1)
{
    mag_t t, u;
    acb_t w;
    slong i;

    mag_init(t);
    mag_init(u);
    acb_init(w);

    mag_set(D, zmag);

    for (i = 0; i < q; i++)
    {
        _lower_bound_range(u, b + i, k0, k1);

        if (i < p)
        {
            acb_sub(w, a + i, b + i, MAG_BITS);
            acb_get_mag(t, w);
            mag_div(t, t, u);
            mag_one(u);
            mag_add(t, t, u);
            mag_mul(D, D, t);
        }
        else
        {
            mag_div(D, D, u);
        }
    }

    mag_clear(t);
    mag_clear(u);
    acb_clear(w);
}

static int
_ulong_cmp(const void * a, const void * b)
{
    ulong x = *((const ulong *) a);
    ulong y = *((const ulong *) b);
    return (x > y) - (x < y);
}

#define MAX_POINT (UWORD(1) << (FLINT_BITS - 1))

static void
_bound_factor_ranges(mag_t C, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, ulong n)
{
    ulong * pts;
    slong i, j, npts, alloc;
    mag_t D, P, S, t, u, zmag;
    acb_t w;
    fmpz_t m;

    alloc = 1 + q * 2 * (FLINT_BITS + 2);
    pts = flint_malloc(alloc * sizeof(ulong));
    npts = 0;
    pts[npts++] = n;

    fmpz_init(m);

    /* breakpoints around each m_i = round(-Re(b_i)) >= n */
    for (i = 0; i < q; i++)
    {
        const arf_struct * x = arb_midref(acb_realref(b + i));
        ulong mi, s;

        if (!arf_is_finite(x) || arf_sgn(x) >= 0
            || arf_cmpabs_2exp_si(x, FLINT_BITS - 2) >= 0)
            continue;

        arf_get_fmpz(m, x, ARF_RND_NEAR);
        fmpz_neg(m, m);
        mi = fmpz_get_ui(m);

        if (mi < n)
            continue;

        for (s = 1; s <= mi - n; s *= 2)
            if (mi - s > n)
                pts[npts++] = mi - s;

        pts[npts++] = mi;

        /* (coarser far out; these are rarely needed) */
        for (s = 1; s < MAX_POINT - mi; )
        {
            pts[npts++] = mi + s;

            if (s <= mi)
                s *= 2;
            else if (s < (MAX_POINT - mi) / 16)
                s *= 16;
            else
                break;
        }
    }

    fmpz_clear(m);

    qsort(pts, npts, sizeof(ulong), _ulong_cmp);

    for (i = j = 0; i < npts; i++)
        if (j == 0 || pts[i] != pts[j - 1])
            pts[j++] = pts[i];
    npts = j;

    mag_init(D);
    mag_init(P);
    mag_init(S);
    mag_init(t);
    mag_init(u);
    mag_init(zmag);
    acb_init(w);

    acb_get_mag(zmag, z);
    mag_one(P);
    mag_zero(C);

    for (j = 0; j < npts && mag_is_finite(C); j++)
    {
        ulong k0 = pts[j];
        ulong k1 = (j + 1 < npts) ? pts[j + 1] : 0;

        /* past all the poles, the ratio bound on [k0, inf) is decreasing:
           finish here if it is < 1 */
        if (k1 != 0)
        {
            for (i = 0; i < q; i++)
            {
                acb_add_ui(w, b + i, k0, MAG_BITS);
                if (!arb_is_positive(acb_realref(w)))
                    break;
            }

            if (i == q)
            {
                _bound_factor_ranges_D(D, a, p, b, q, zmag, k0, 0);
                mag_one(t);
                mag_sub_lower(u, t, D);

                if (!mag_is_zero(u))
                {
                    mag_div(S, t, u);
                    mag_mul(S, S, P);
                    mag_add(C, C, S);
                    break;
                }
            }
        }

        /* ratio bound on [k0, k1) */
        _bound_factor_ranges_D(D, a, p, b, q, zmag, k0, k1);

        /* bound for the sum of D^t over the range */
        mag_one(t);
        mag_sub_lower(u, t, D);

        if (!mag_is_zero(u))
        {
            mag_div(S, t, u);
        }
        else if (k1 != 0 && mag_is_finite(D))
        {
            mag_pow_ui(S, D, k1 - k0 - 1);
            mag_mul_ui(S, S, k1 - k0);
        }
        else
        {
            mag_inf(C);
            break;
        }

        mag_mul(S, S, P);
        mag_add(C, C, S);

        if (k1 != 0)
        {
            mag_pow_ui(t, D, k1 - k0);
            mag_mul(P, P, t);
        }
    }

    flint_free(pts);
    mag_clear(D);
    mag_clear(P);
    mag_clear(S);
    mag_clear(t);
    mag_clear(u);
    mag_clear(zmag);
    acb_clear(w);
}

void
acb_hypgeom_pfq_bound_factor(mag_t C,
    acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, ulong n)
{
    mag_t t, u;
    acb_t w;
    slong i;
    int increasing = 1;

    if (p > q)
    {
        mag_inf(C);
        return;
    }

    mag_init(t);
    mag_init(u);
    acb_init(w);

    acb_get_mag(C, z);

    for (i = 0; i < q && increasing; i++)
    {
        acb_add_ui(w, b + i, n, MAG_BITS);

        if (arb_is_positive(acb_realref(w)))
        {
            acb_get_mag_lower(u, w);

            if (i < p)
            {
                acb_sub(w, a + i, b + i, MAG_BITS);
                acb_get_mag(t, w);
                mag_div(t, t, u);
                mag_one(u);
                mag_add(t, t, u);
                mag_mul(C, C, t);
            }
            else
            {
                mag_div(C, C, u);
            }
        }
        else
        {
            increasing = 0;
        }
    }

    if (increasing)
    {
        mag_one(t);
        mag_sub_lower(u, t, C);

        if (mag_is_zero(u))
            mag_inf(C);
        else
            mag_div(C, t, u);
    }
    else if (n < MAX_POINT && mag_is_finite(C))
    {
        /* the uniform bound is cheaper; use it if it is good */
        acb_get_mag(t, z);
        _bound_factor_uniform(C, a, p, b, q, t, n);

        if (!mag_is_finite(C) || mag_cmp_2exp_si(C, 4) > 0)
        {
            _bound_factor_ranges(u, a, p, b, q, z, n);
            mag_min(C, C, u);
        }
    }
    else
    {
        mag_inf(C);
    }

    mag_clear(t);
    mag_clear(u);
    acb_clear(w);
}
