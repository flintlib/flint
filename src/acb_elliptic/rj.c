/*
    Copyright (C) 2017, 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"
#include "acb_calc.h"

#ifdef __GNUC__
# define pow __builtin_pow
#else
# include <math.h>
#endif

static const unsigned short den_ratio_tab[512] = {
    1,1,14,3,44,13,10,17,152,1,46,1,12,29,62,1,
    16,37,2,41,172,1,94,7,8,53,2,1,236,61,2,1,
    2144,1,142,73,20,1,158,3,664,1,2,89,4,1,2,97,
    16,101,206,1,428,109,2,113,8,1,2,11,4,1,254,1,
    8384,1,2,137,556,1,2,1,8,149,302,1,4,157,2,1,
    2608,1,334,13,4,173,2,1,1432,181,2,1,4,1,382,193,
    32,197,398,1,4,1,2,1,1688,1,2,1,4,1,446,1,
    3632,229,2,233,4,1,478,241,24,1,2,1,1004,1,2,257,
    128,1,526,1,4,269,542,1,8,277,2,281,1132,1,2,17,
    16,293,2,1,4,1,2,1,2456,1,622,313,4,317,2,1,
    32,1,2,1,1324,1,2,337,8,1,14,1,1388,349,2,353,
    16,1,718,19,4,1,734,1,8,373,10,1,1516,1,766,1,
    64,389,2,1,4,397,2,401,8,1,2,409,4,1,2,1,
    6704,421,2,1,4,1,862,433,8,1,878,1,1772,1,2,449,
    32,1,2,457,4,461,926,1,3736,1,2,1,4,1,958,1,
    16,1,974,1,1964,1,2,1,3992,1,1006,1,4,509,2,1,
    256,1,2,521,2092,1,2,23,8,1,2,1,4,541,2,1,
    8752,1,2,1,4,557,2,1,4504,1,2,569,2284,1,2,577,
    32,1,2,1,2348,1,2,593,8,1,1198,601,4,1,1214,1,
    16,613,2,617,2476,1,2,1,8,1,1262,1,4,1,2,641,
    41152,1,1294,1,4,653,2,1,5272,661,2,1,4,1,2,673,
    16,677,2,1,2732,1,2,1,5528,1,2,1,4,701,2,1,
    32,709,2,1,4,1,1438,1,8,1,1454,3,4,733,2,1,
    11824,1,1486,1,4,1,1502,1,8,757,2,761,4,1,2,769,
    128,773,2,1,4,1,2,1,6296,1,2,1,4,797,2,1,
    16,1,2,809,3244,1,2,1,8,821,1646,1,3308,829,2,1,
    32,1,1678,29,4,1,2,1,8,853,2,857,3436,1,1726,1,
    16,1,2,1,4,877,2,881,7064,1,1774,1,4,1,2,1,
    64,1,2,1,3628,1,1822,1,8,1,1838,1,4,1,2,929,
    16,1,2,937,4,941,2,1,7576,1,2,953,4,1,2,31,
    32,1,1934,1,3884,1,2,977,8,1,1966,1,4,1,1982,1,
    16,997,2,1,4,1,2,1009,8,1013,2,1,4076,1021,2,1
};

static inline slong rj_fdiv(slong x, slong y)
{
    if (x < 0)
        return -1;
    else
        return x / y;
}

void
acb_elliptic_rj_taylor_sum(acb_t res, const acb_t E2, const acb_t E3,
                const acb_t E4, const acb_t E5, slong nterms, slong prec)
{
    slong m2, m3, m4, m5, m2start, m3start, m4start, m5start, NMAX, N, M;
    slong m2dim, m3dim;
    acb_t s2, s3, s4, s5;
    acb_ptr powtab;
    fmpz_t c2, c3, c4, c5, den, t;

    acb_init(s2);
    acb_init(s3);
    acb_init(s4);
    acb_init(s5);
    fmpz_init(c2);
    fmpz_init(c3);
    fmpz_init(c4);
    fmpz_init(c5);
    fmpz_init(den);
    fmpz_init(t);

    NMAX = nterms - 1;

    m2dim = NMAX / 2 + 1;
    m3dim = NMAX / 3 + 1;

    powtab = _acb_vec_init(m2dim * m3dim);

    /* Compute universal denominator */
    fmpz_one(den);
    for (N = 1; N <= NMAX; N++)
        fmpz_mul_ui(den, den, den_ratio_tab[N]);

    /* Precompute powers of E2 and E3 */
    for (m2 = 0; m2 <= NMAX / 2; m2++)
    {
        for (m3 = 0; m3 <= rj_fdiv(NMAX - 2 * m2, 3); m3++)
        {
            slong i, j, k;

            i = m3 * m2dim + m2;

            if (m2 <= 1 && m3 <= 1)
            {
                if (m2 == 0 && m3 == 0)
                    acb_one(powtab + i);
                else if (m2 == 0 && m3 == 1)
                    acb_set(powtab + i, E3);
                else if (m2 == 1 && m3 == 0)
                    acb_set(powtab + i, E2);
                else
                    acb_mul(powtab + i, E2, E3, prec);
            }
            else
            {
                j = (m3 / 2) * m2dim + (m2 / 2);
                k = (m3 - (m3 / 2)) * m2dim + (m2 - (m2 / 2));
                acb_mul(powtab + i, powtab + j, powtab + k, prec);
            }
        }
    }

    acb_zero(s5);
    m5start = NMAX / 5;
    fmpz_mul_ui(c5, den, 3);
    for (m5 = 0; m5 < m5start; m5++)
    {
        fmpz_mul_ui(c5, c5, 2 * m5 + 1);
        fmpz_divexact_ui(c5, c5, 2 * m5 + 2);
    }

    for (m5 = m5start; m5 >= 0; m5--)
    {
        acb_zero(s4);
        m4start = rj_fdiv(NMAX - 5 * m5, 4);
        if (m5 != m5start)
        {
            fmpz_mul_ui(c5, c5, 2 * m5 + 2);
            fmpz_divexact_ui(c5, c5, 2 * m5 + 1);
        }
        fmpz_set(c4, c5);
        for (m4 = 0; m4 < m4start; m4++)
        {
            fmpz_mul_ui(c4, c4, 2 * m5 + 2 * m4 + 1);
            fmpz_divexact_ui(c4, c4, 2 * m4 + 2);
        }

        for (m4 = m4start; m4 >= 0; m4--)
        {
            acb_zero(s3);
            m3start = rj_fdiv(NMAX - 5 * m5 - 4 * m4, 3);
            if (m4 != m4start)
            {
                fmpz_mul_ui(c4, c4, 2 * m4 + 2);
                fmpz_divexact_ui(c4, c4, 2 * m5 + 2 * m4 + 1);
            }
            fmpz_set(c3, c4);

            for (m3 = 0; m3 <= m3start; m3++)
            {
                m2start = rj_fdiv(NMAX - 5 * m5 - 4 * m4 - 3 * m3, 2);
                fmpz_set(c2, c3);
                for (m2 = 0; m2 <= m2start; m2++)
                {
                    M = m5 + m4 + m3 + m2;
                    N = 5 * m5 + 4 * m4 + 3 * m3 + 2 * m2;

                    if (N > NMAX)
                        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

                    fmpz_divexact_ui(t, c2, 2 * N + 3);
                    if ((M + N) % 2)
                        fmpz_neg(t, t);

                    acb_addmul_fmpz(s3, powtab + m3 * m2dim + m2, t, prec);

                    if (m2 < m2start)
                    {
                        fmpz_mul_ui(c2, c2, 2 * m5 + 2 * m4 + 2 * m3 + 2 * m2 + 1);
                        fmpz_divexact_ui(c2, c2, 2 * m2 + 2);
                    }
                }

                if (m3 < m3start)
                {
                    fmpz_mul_ui(c3, c3, 2 * m5 + 2 * m4 + 2 * m3 + 1);
                    fmpz_divexact_ui(c3, c3, 2 * m3 + 2);
                }
            }

            /* Horner with respect to E4. */
            acb_mul(s4, s4, E4, prec);
            acb_add(s4, s4, s3, prec);
        }

        /* Horner with respect to E5. */
        acb_mul(s5, s5, E5, prec);
        acb_add(s5, s5, s4, prec);
    }

    acb_div_fmpz(res, s5, den, prec);

    _acb_vec_clear(powtab, m2dim * m3dim);
    acb_clear(s2);
    acb_clear(s3);
    acb_clear(s4);
    acb_clear(s5);
    fmpz_clear(c2);
    fmpz_clear(c3);
    fmpz_clear(c4);
    fmpz_clear(c5);
    fmpz_clear(den);
    fmpz_clear(t);
}

void
acb_elliptic_rj_carlson(acb_t res, const acb_t x, const acb_t y,
            const acb_t z, const acb_t p, int flags, slong prec)
{
    acb_t xx, yy, zz, pp, sx, sy, sz, sp, t, d, delta, S;
    acb_t A, AA, X, Y, Z, P, E2, E3, E4, E5;
    mag_t err, err2, prev_err;
    slong k, wp, accx, accy, accz, accp, order;
    int rd, real;

/*
    printf("RJ carlson   ");
    acb_printn(x, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(y, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(z, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(p, 6, ARB_STR_NO_RADIUS); printf("  ");
    printf("\n");
*/

    if (!acb_is_finite(x) || !acb_is_finite(y) || !acb_is_finite(z) ||
        !acb_is_finite(p))
    {
        acb_indeterminate(res);
        return;
    }

    if ((acb_contains_zero(x) + acb_contains_zero(y) + acb_contains_zero(z) > 1)
        || acb_contains_zero(p))
    {
        acb_indeterminate(res);
        return;
    }

    /* Special case computing R_D(x,y,z) */
    rd = (z == p) || acb_eq(z, p);

    acb_init(xx); acb_init(yy); acb_init(zz); acb_init(pp);
    acb_init(sx); acb_init(sy); acb_init(sz); acb_init(sp);
    acb_init(S); acb_init(A); acb_init(AA);
    acb_init(X); acb_init(Y); acb_init(Z); acb_init(P);
    acb_init(E2); acb_init(E3); acb_init(E4); acb_init(E5);
    acb_init(t); acb_init(d); acb_init(delta);
    mag_init(err);
    mag_init(err2);
    mag_init(prev_err);

    acb_set(xx, x);
    acb_set(yy, y);
    acb_set(zz, z);
    acb_set(pp, p);
    acb_zero(S);

    real = acb_is_real(x) && acb_is_real(y) && acb_is_real(z) && acb_is_real(p) &&
           arb_is_nonnegative(acb_realref(x)) && arb_is_nonnegative(acb_realref(y)) &&
           arb_is_nonnegative(acb_realref(z)) && arb_is_nonnegative(acb_realref(p));

    order = 5; /* will be set later */

    /* First guess precision based on the inputs. */
    /* This does not account for mixing. */
    accx = acb_rel_accuracy_bits(xx);
    accy = acb_rel_accuracy_bits(yy);
    accz = acb_rel_accuracy_bits(zz);
    accp = acb_rel_accuracy_bits(pp);
    accx = FLINT_MAX(accx, accy);
    accx = FLINT_MAX(accx, accz);
    accx = FLINT_MAX(accx, accp);
    if (accx < prec - 20)
        prec = FLINT_MAX(2, accx + 20);
    wp = prec + 4 + FLINT_BIT_COUNT(prec);

    if (!rd)
    {
        acb_mul_2exp_si(A, p, 1);
        acb_add(A, A, z, wp);
    }
    else
    {
        acb_mul_ui(A, z, 3, wp);
    }
    acb_add(A, A, x, wp);
    acb_add(A, A, y, wp);
    acb_div_ui(A, A, 5, wp);
    acb_set(AA, A);

    if (!rd)
    {
        acb_sub(delta, p, x, wp);
        acb_sub(t, p, y, wp);
        acb_mul(delta, delta, t, wp);
        acb_sub(t, p, z, wp);
        acb_mul(delta, delta, t, wp);
    }

    /* must do at least one iteration */
    for (k = 0; k < prec; k++)
    {
        acb_sqrt(sx, xx, wp);
        acb_sqrt(sy, yy, wp);
        acb_sqrt(sz, zz, wp);
        if (!rd) acb_sqrt(sp, pp, wp);

        acb_add(t, sy, sz, wp);
        acb_mul(t, t, sx, wp);
        acb_addmul(t, sy, sz, wp);

        acb_add(xx, xx, t, wp);
        acb_add(yy, yy, t, wp);
        acb_add(zz, zz, t, wp);
        if (!rd) acb_add(pp, pp, t, wp);
        acb_add(AA, AA, t, wp);

        acb_mul_2exp_si(xx, xx, -2);
        acb_mul_2exp_si(yy, yy, -2);
        acb_mul_2exp_si(zz, zz, -2);
        if (!rd) acb_mul_2exp_si(pp, pp, -2);
        acb_mul_2exp_si(AA, AA, -2);

        if (!rd)
        {
            /* d = (sp+sx)(sp+sy)(sp+sz) */
            /* e = 4^(-3k) delta / d^2 */
            /* S += 4^(-k) RC(1, 1+e) / d */
            acb_add(d, sp, sx, wp);
            acb_add(t, sp, sy, wp);
            acb_mul(d, d, t, wp);
            acb_add(t, sp, sz, wp);
            acb_mul(d, d, t, wp);

            /* E2 = e */
            acb_mul(E2, d, d, wp);
            acb_div(E2, delta, E2, wp);
            acb_mul_2exp_si(E2, E2, -6 * k);

            acb_elliptic_rc1(E4, E2, wp);
            acb_div(E4, E4, d, wp);
            acb_mul_2exp_si(E4, E4, -2 * k);

            acb_add(S, S, E4, wp);
        }
        else
        {
            acb_mul(t, sz, zz, wp);
            acb_mul_2exp_si(t, t, 2);
            acb_inv(t, t, wp);
            acb_mul_2exp_si(t, t, -2 * k);
            acb_mul_2exp_si(t, t, -1);

            acb_add(S, S, t, wp);
        }

        /* Improve precision estimate and set expansion order. */
        /* Should this done for other k also? */
        if (k == 0)
        {
            accx = acb_rel_accuracy_bits(xx);
            accy = acb_rel_accuracy_bits(yy);
            accz = acb_rel_accuracy_bits(zz);
            accp = acb_rel_accuracy_bits(pp);
            accx = FLINT_MAX(accx, accy);
            accx = FLINT_MAX(accx, accz);
            accx = FLINT_MAX(accx, accp);
            if (accx < prec - 20)
                prec = FLINT_MAX(2, accx + 20);
            wp = prec + 4 + FLINT_BIT_COUNT(prec);

            if (!rd)
                if (real)
                    order = 2.3 * pow(prec, 0.34);
                else
                    order = 2.5 * pow(prec, 0.35);
            else
                if (real)
                    order = 2.0 * pow(prec, 0.33);
                else
                    order = 2.2 * pow(prec, 0.33);
            order = FLINT_MIN(order, 500);
            order = FLINT_MAX(order, 5);
        }

        /* Close enough? */
        acb_sub(t, xx, yy, wp);
        acb_get_mag(err, t);
        acb_sub(t, xx, zz, wp);
        acb_get_mag(err2, t);
        mag_max(err, err, err2);
        if (!rd)
        {
            acb_sub(t, xx, pp, wp);
            acb_get_mag(err2, t);
            mag_max(err, err, err2);
        }
        acb_get_mag_lower(err2, xx);
        mag_div(err, err, err2);

        mag_pow_ui(err, err, order);

        if (mag_cmp_2exp_si(err, -prec) < 0 ||
                (k > 2 && mag_cmp(err, prev_err) > 0))
        {
            k++;
            break;
        }

        mag_set(prev_err, err);
    }

    /* X = (A-x)/(4^k AA) */
    /* Y = (A-y)/(4^k AA) */
    /* Z = (A-z)/(4^k AA) */
    /* P = (-X-Y-Z)/2 */
    acb_mul_2exp_si(t, AA, 2 * k);
    acb_inv(t, t, prec);
    acb_sub(X, A, x, prec);
    acb_mul(X, X, t, prec);
    acb_sub(Y, A, y, prec);
    acb_mul(Y, Y, t, prec);
    acb_sub(Z, A, z, prec);
    acb_mul(Z, Z, t, prec);
    acb_add(P, X, Y, prec);
    acb_add(P, P, Z, prec);
    acb_neg(P, P);
    acb_mul_2exp_si(P, P, -1);

    /* todo: improve for R_D */
    /* E2 = XY + XZ + YZ - 3 P^2 */
    /* E3 = XYZ + 2 E2 P + 4 P^3 */
    /* E4 = (2 XYZ + E2 P + 3 P^3) P */
    /* E5 = XYZP^2 */

    acb_mul(t, P, P, prec); /* t = P^2 */

    acb_mul(E2, X, Y, prec);
    acb_mul(E3, E2, Z, prec);
    acb_mul_2exp_si(E4, E3, 1);
    acb_mul(E5, E3, t, prec);

    acb_add(sx, X, Y, prec);
    acb_addmul(E2, sx, Z, prec);
    acb_submul_ui(E2, t, 3, prec);

    acb_mul(sx, E2, P, prec);
    acb_add(E4, E4, sx, prec);
    acb_mul_2exp_si(sx, sx, 1);
    acb_add(E3, E3, sx, prec);

    acb_mul(t, t, P, prec); /* t = P^3 */
    acb_addmul_ui(E3, t, 4, prec);
    acb_addmul_ui(E4, t, 3, prec);
    acb_mul(E4, E4, P, prec);

    /* Error bound. */
    acb_get_mag(err, X);
    acb_get_mag(err2, Y);
    mag_max(err, err, err2);
    acb_get_mag(err2, Z);
    mag_max(err, err, err2);
    acb_get_mag(err2, P);
    mag_max(err, err, err2);
    mag_mul_ui(err, err, 9);
    mag_mul_2exp_si(err, err, -3);
    mag_geom_series(err, err, order);
    mag_mul_2exp_si(err, err, 1);

    acb_elliptic_rj_taylor_sum(sx, E2, E3, E4, E5, order, wp);

    if (acb_is_real(X) && acb_is_real(Y) && acb_is_real(Z))
        arb_add_error_mag(acb_realref(sx), err);
    else
        acb_add_error_mag(sx, err);

    acb_rsqrt(t, AA, wp);
    acb_div(t, t, AA, wp);
    acb_mul_2exp_si(t, t, -2 * k);
    acb_mul(t, t, sx, wp);

    acb_addmul_ui(t, S, 6, prec);

    acb_set(res, t);

    acb_clear(xx); acb_clear(yy); acb_clear(zz); acb_clear(pp);
    acb_clear(sx); acb_clear(sy); acb_clear(sz); acb_clear(sp);
    acb_clear(S); acb_clear(A); acb_clear(AA);
    acb_clear(X); acb_clear(Y); acb_clear(Z); acb_clear(P);
    acb_clear(E2); acb_clear(E3); acb_clear(E4); acb_clear(E5);
    acb_clear(t); acb_clear(d); acb_clear(delta);
    mag_clear(err);
    mag_clear(err2);
    mag_clear(prev_err);
}

static int
acb_eq_conj(const acb_t x, const acb_t y)
{
    int res;
    acb_t t;
    acb_init(t);
    acb_conj(t, y);
    res = acb_eq(x, t);
    acb_clear(t);
    return res;
}

/* todo: speed up this function
   -- early abort
   -- only compute single rsqrt or sqrt when evaluating precisely
      (but need to make sure branch is correct!)
*/
static int
RJ_integrand(acb_ptr res, const acb_t t, void * param, slong order, slong prec)
{
    acb_ptr x, y, z, p;
    acb_t xt, yt, zt, pt;
    int analytic, deflated;

    if (order > 1)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);  /* Would be needed for Taylor method. */

    x = ((acb_ptr) param);
    y = ((acb_ptr) param) + 1;
    z = ((acb_ptr) param) + 2;
    p = ((acb_ptr) param) + 3;
    deflated = acb_is_zero(x);

    analytic = (order != 0);

    acb_init(xt);
    acb_init(yt);
    acb_init(zt);
    acb_init(pt);

    /* if x = 0, change of variables t -> t^2 to remove singularity at 0 */
    if (deflated)
    {
        acb_sqr(xt, t, prec);

        acb_add(yt, y, xt, prec);
        acb_add(zt, z, xt, prec);
        acb_add(pt, p, xt, prec);

        if (acb_contains_zero(yt) || acb_contains_zero(zt) || acb_contains_zero(pt))
        {
            acb_indeterminate(res);
        }
        else
        {
            acb_rsqrt_analytic(yt, yt, analytic, prec);
            acb_rsqrt_analytic(zt, zt, analytic, prec);

            acb_mul(xt, yt, zt, prec);
            acb_div(xt, xt, pt, prec);
            acb_mul_2exp_si(xt, xt, 1);

            acb_set(res, xt);
        }
    }
    else
    {
        acb_add(xt, x, t, prec);
        acb_add(yt, y, t, prec);
        acb_add(zt, z, t, prec);
        acb_add(pt, p, t, prec);

        if (acb_contains_zero(xt) || acb_contains_zero(yt) || acb_contains_zero(zt) || acb_contains_zero(pt))
        {
            acb_indeterminate(res);
        }
        else
        {
            acb_rsqrt_analytic(xt, xt, analytic, prec);
            acb_rsqrt_analytic(yt, yt, analytic, prec);
            acb_rsqrt_analytic(zt, zt, analytic, prec);

            acb_mul(xt, xt, yt, prec);
            acb_mul(xt, xt, zt, prec);
            acb_div(xt, xt, pt, prec);

            acb_set(res, xt);
        }
    }

    acb_clear(xt);
    acb_clear(yt);
    acb_clear(zt);
    acb_clear(pt);

    return 0;
}

void
acb_elliptic_rj_integration(acb_t res, const acb_t x, const acb_t y,
            const acb_t z, const acb_t p, int flags, slong prec)
{
    acb_t a, b, N, I, J;
    arb_t A;
    acb_ptr xyzp;
    mag_t tol;
    int deflated;

/*
    printf("RJ integration:   ");
    acb_printn(x, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(y, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(z, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(p, 6, ARB_STR_NO_RADIUS); printf("  ");
    printf("\n");
*/

    acb_init(N);
    acb_init(a);
    acb_init(b);
    acb_init(I);
    acb_init(J);
    arb_init(A);
    xyzp = _acb_vec_init(4);
    mag_init(tol);

    /* compute shift that puts parameters in right half-plane */
    arb_min(A, acb_realref(x), acb_realref(y), prec);
    arb_min(A, A, acb_realref(z), prec);
    arb_min(A, A, acb_realref(p), prec);
    arb_neg(A, A);
    arb_one(acb_realref(a));
    arb_max(A, A, acb_realref(a), prec);
    arb_add_ui(A, A, 2, prec);

    arb_get_ubound_arf(arb_midref(A), A, prec);
    mag_zero(arb_radref(A));

    acb_set(xyzp, x);
    acb_set(xyzp + 1, y);
    acb_set(xyzp + 2, z);
    acb_set(xyzp + 3, p);

    /* If there is a zero among x, y, z, put it first. */
    if (acb_is_zero(y))
        acb_swap(xyzp, xyzp + 1);
    if (acb_is_zero(z))
        acb_swap(xyzp, xyzp + 2);

    deflated = acb_is_zero(xyzp);

    acb_set_arb(N, A);

    /* Path deformation to avoid 0 */
    if ((arb_is_nonnegative(acb_imagref(x)) || arb_is_positive(acb_realref(x))) &&
        (arb_is_nonnegative(acb_imagref(y)) || arb_is_positive(acb_realref(y))) &&
        (arb_is_nonnegative(acb_imagref(z)) || arb_is_positive(acb_realref(z))) &&
        (arb_is_nonnegative(acb_imagref(p)) || arb_is_positive(acb_realref(p))))
    {
        arb_set_si(acb_imagref(N), 1);
    }
    else if ((arb_is_negative(acb_imagref(x)) || arb_is_positive(acb_realref(x))) &&
            (arb_is_negative(acb_imagref(y)) || arb_is_positive(acb_realref(y))) &&
            (arb_is_negative(acb_imagref(z)) || arb_is_positive(acb_realref(z))) &&
            (arb_is_negative(acb_imagref(p)) || arb_is_positive(acb_realref(p))))
    {
        arb_set_si(acb_imagref(N), -1);
    }
    else
    {
        int i;

        arb_set_si(acb_imagref(N), 2);

        /* Go through the upper half-plane, but low enough that any
           parameter starting in the lower plane doesn't cross the
           branch cut */
        for (i = 0; i < 4; i++)
        {
            if (deflated && (i == 0))
                continue;

            if (arb_is_nonnegative(acb_imagref(xyzp + i)) ||
                arb_is_positive(acb_realref(xyzp + i)))
                continue;

            arb_zero(acb_realref(a)); /* use as tmp var */
            arb_get_abs_lbound_arf(arb_midref(acb_realref(a)), acb_imagref(xyzp + i), prec);
            arb_min(acb_imagref(N), acb_imagref(N), acb_realref(a), prec);
        }

        arb_mul_2exp_si(acb_imagref(N), acb_imagref(N), -1);
    }

    mag_one(tol);
    mag_mul_2exp_si(tol, tol, -prec);
    acb_zero(a);
    if (deflated)
        acb_sqrt(b, N, prec);
    else
        acb_set(b, N);

/*
    flint_printf("integrate ");
    flint_printf("%d\n", deflated);
    flint_printf("x = "); acb_printd(xyzp, 10); flint_printf("\n");
    flint_printf("y = "); acb_printd(xyzp + 1, 10); flint_printf("\n");
    flint_printf("z = "); acb_printd(xyzp + 2, 10); flint_printf("\n");
    flint_printf("p = "); acb_printd(xyzp + 3, 10); flint_printf("\n");
    flint_printf("a = "); acb_printd(a, 10); flint_printf("\n");
    flint_printf("b = "); acb_printd(b, 10); flint_printf("\n");
    flint_printf("N = "); acb_printd(N, 10); flint_printf("\n");
*/

    acb_calc_integrate(I, RJ_integrand, xyzp, a, b, prec, tol, NULL, prec);
    acb_mul_ui(I, I, 3, prec);
    acb_mul_2exp_si(I, I, -1);

/*
    flint_printf("I = "); acb_printd(I, 10); flint_printf("\n");
*/

    acb_add(xyzp, x, N, prec);
    acb_add(xyzp + 1, y, N, prec);
    acb_add(xyzp + 2, z, N, prec);
    acb_add(xyzp + 3, p, N, prec);

    acb_elliptic_rj_carlson(J, xyzp, xyzp + 1, xyzp + 2, xyzp + 3, 0, prec);

    acb_add(res, I, J, prec);

    acb_clear(N);
    acb_clear(a);
    acb_clear(b);
    acb_clear(I);
    acb_clear(J);
    arb_clear(A);
    _acb_vec_clear(xyzp, 4);
    mag_clear(tol);
}

void
acb_elliptic_rj(acb_t res, const acb_t x, const acb_t y,
            const acb_t z, const acb_t p, int flags, slong prec)
{
/*
    printf("RJ:   ");
    acb_printn(x, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(y, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(z, 6, ARB_STR_NO_RADIUS); printf("  ");
    acb_printn(p, 6, ARB_STR_NO_RADIUS); printf("  ");
    printf("\n");
*/

    if (!acb_is_finite(x) || !acb_is_finite(y) || !acb_is_finite(z) ||
        !acb_is_finite(p))
    {
        acb_indeterminate(res);
        return;
    }

    if ((acb_contains_zero(x) + acb_contains_zero(y) + acb_contains_zero(z) > 1)
        || acb_contains_zero(p))
    {
        acb_indeterminate(res);
        return;
    }

    /* Carlson's algorithm is correct in the degenerate case
       computing R_D */
    if (x == p || acb_eq(x, p))
    {
        acb_elliptic_rj_carlson(res, y, z, x, p, flags, prec);
        return;
    }

    if (y == p || acb_eq(y, p))
    {
        acb_elliptic_rj_carlson(res, x, z, y, p, flags, prec);
        return;
    }

    if (z == p || acb_eq(z, p))
    {
        acb_elliptic_rj_carlson(res, x, y, z, p, flags, prec);
        return;
    }

    /* Sufficient condition for correctness */
    if (arb_is_nonnegative(acb_realref(x)) &&
          arb_is_nonnegative(acb_realref(y)) &&
          arb_is_nonnegative(acb_realref(z)) &&
          arb_is_positive(acb_realref(p)))
    {
        acb_elliptic_rj_carlson(res, x, y, z, p, flags, prec);
        return;
    }

    /* Sufficient condition for correctness */
    if (acb_is_real(x) && acb_is_real(y) && acb_is_real(z) && acb_is_real(p))
    {
        acb_elliptic_rj_carlson(res, x, y, z, p, flags, prec);
        return;
    }

    /* Also a sufficient condition */
    if (arb_is_nonnegative(acb_realref(p)) || arb_is_nonzero(acb_imagref(p)))
    {
        if ((arb_is_zero(acb_imagref(x)) && arb_is_nonnegative(acb_realref(x)) && acb_eq_conj(y, z)) ||
            (arb_is_zero(acb_imagref(y)) && arb_is_nonnegative(acb_realref(y)) && acb_eq_conj(x, z)) ||
            (arb_is_zero(acb_imagref(z)) && arb_is_nonnegative(acb_realref(z)) && acb_eq_conj(x, y)))
        {
            acb_elliptic_rj_carlson(res, x, y, z, p, flags, prec);
            return;
        }
    }

    /* Fast abort for input straddling branch cuts */
    if ((arb_contains_zero(acb_imagref(x)) && !(arb_is_nonnegative(acb_imagref(x)) || arb_is_nonnegative(acb_realref(x)))) ||
        (arb_contains_zero(acb_imagref(y)) && !(arb_is_nonnegative(acb_imagref(y)) || arb_is_nonnegative(acb_realref(y)))) ||
        (arb_contains_zero(acb_imagref(z)) && !(arb_is_nonnegative(acb_imagref(z)) || arb_is_nonnegative(acb_realref(z)))) ||
        (arb_contains_zero(acb_imagref(p)) && !(arb_is_nonnegative(acb_imagref(p)) || arb_is_nonnegative(acb_realref(p)))))
    {
        acb_indeterminate(res);
        return;
    }

    /* Use integration as fallback */
    acb_elliptic_rj_integration(res, x, y, z, p, flags, prec);
}

