/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "hypgeom.h"

slong
hypgeom_root_bound(const mag_t z, int r)
{
    if (r == 0)
    {
        return 0;
    }
    else
    {
        arf_t t;
        slong v;
        arf_init(t);
        arf_set_mag(t, z);
        arf_root(t, t, r, MAG_BITS, ARF_RND_UP);
        arf_add_ui(t, t, 1, MAG_BITS, ARF_RND_UP);
        v = arf_get_si(t, ARF_RND_UP);
        arf_clear(t);
        return v;
    }
}

/*
Given T(K), compute bound for T(n) z^n.

We need to multiply by

z^n * 1/rf(K+1,m)^r * (rf(K+1,m)/rf(K+1-A,m)) * (rf(K+1-B,m)/rf(K+1-2B,m))

where m = n - K. This is equal to

z^n *

(K+A)! (K-2B)! (K-B+m)!
-----------------------    * ((K+m)! / K!)^(1-r)
(K-B)! (K-A+m)! (K-2B+m)!
*/
void
hypgeom_term_bound(mag_t Tn, const mag_t TK, slong K, slong A, slong B, int r, const mag_t z, slong n)
{
    mag_t t, u, num;
    slong m;

    mag_init(t);
    mag_init(u);
    mag_init(num);

    m = n - K;

    if (m < 0)
    {
        flint_throw(FLINT_ERROR, "hypgeom term bound\n");
    }

    /* TK * z^n */
    mag_pow_ui(t, z, n);
    mag_mul(num, TK, t);

    /* numerator: (K+A)! (K-2B)! (K-B+m)! */
    mag_fac_ui(t, K+A);
    mag_mul(num, num, t);

    mag_fac_ui(t, K-2*B);
    mag_mul(num, num, t);

    mag_fac_ui(t, K-B+m);
    mag_mul(num, num, t);

    /* denominator: (K-B)! (K-A+m)! (K-2B+m)! */
    mag_rfac_ui(t, K-B);
    mag_mul(num, num, t);

    mag_rfac_ui(t, K-A+m);
    mag_mul(num, num, t);

    mag_rfac_ui(t, K-2*B+m);
    mag_mul(num, num, t);

    /* ((K+m)! / K!)^(1-r) */
    if (r == 0)
    {
        mag_fac_ui(t, K+m);
        mag_mul(num, num, t);

        mag_rfac_ui(t, K);
        mag_mul(num, num, t);
    }
    else if (r != 1)
    {
        mag_fac_ui(t, K);
        mag_rfac_ui(u, K+m);
        mag_mul(t, t, u);

        mag_pow_ui(t, t, r-1);
        mag_mul(num, num, t);
    }

    mag_set(Tn, num);

    mag_clear(t);
    mag_clear(u);
    mag_clear(num);
}

slong
hypgeom_bound(mag_t error, int r,
    slong A, slong B, slong K, const mag_t TK, const mag_t z, slong tol_2exp)
{
    mag_t Tn, t, u, one, tol, num, den;
    slong n, m;

    mag_init(Tn);
    mag_init(t);
    mag_init(u);
    mag_init(one);
    mag_init(tol);
    mag_init(num);
    mag_init(den);

    mag_one(one);
    mag_set_ui_2exp_si(tol, UWORD(1), -tol_2exp);

    /* approximate number of needed terms */
    n = hypgeom_estimate_terms(z, r, tol_2exp);

    /* required for 1 + O(1/k) part to be decreasing */
    n = FLINT_MAX(n, K + 1);

    /* required for z^k / (k!)^r to be decreasing */
    m = hypgeom_root_bound(z, r);
    n = FLINT_MAX(n, m);

    /*  We now have |R(k)| <= G(k) where G(k) is monotonically decreasing,
        and can bound the tail using a geometric series as soon
        as soon as G(k) < 1. */

    /* bound T(n-1) */
    hypgeom_term_bound(Tn, TK, K, A, B, r, z, n-1);

    while (1)
    {
        /* bound R(n) */
        mag_mul_ui(num, z, n);
        mag_mul_ui(num, num, n - B);

        mag_set_ui_lower(den, n - A);
        mag_mul_ui_lower(den, den, n - 2*B);

        if (r != 0)
        {
            mag_set_ui_lower(u, n);
            mag_pow_ui_lower(u, u, r);
            mag_mul_lower(den, den, u);
        }

        mag_div(t, num, den);

        /* multiply bound for T(n-1) by bound for R(n) to bound T(n) */
        mag_mul(Tn, Tn, t);

        /* geometric series termination check */
        /* u = max(1-t, 0), rounding down [lower bound] */
        mag_sub_lower(u, one, t);

        if (!mag_is_zero(u))
        {
            mag_div(u, Tn, u);

            if (mag_cmp(u, tol) < 0)
            {
                mag_set(error, u);
                break;
            }
        }

        /* move on to next term */
        n++;
    }

    mag_clear(Tn);
    mag_clear(t);
    mag_clear(u);
    mag_clear(one);
    mag_clear(tol);
    mag_clear(num);
    mag_clear(den);

    return n;
}

