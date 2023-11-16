/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb_mat.h"
#include "arb_hypgeom.h"
#include "acb_dirichlet.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

#define VERBOSE 0

static double
log_gamma_upper_approx(double a, double z)
{
    if (a < z)
        return (a - 1) * log(z) - z;
    else
        return a * (log(a) - 1);
}

void
acb_dirichlet_root_number2(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    acb_dirichlet_root_number(res, G, chi, prec);

    if (dirichlet_char_is_real(G, chi))
        arb_zero(acb_imagref(res));
}

static void
_arf_trunc(arf_t x)
{
    if (arf_sgn(x) < 0)
        arf_ceil(x, x);
    else
        arf_floor(x, x);
}

static void
arb_extract_bits(arb_t t, const arb_t z, slong b)
{
    arb_mul_2exp_si(t, z, b);
    _arf_trunc(arb_midref(t));
    mag_zero(arb_radref(t));
    arb_mul_2exp_si(t, t, -b);
}

static void
acb_dirichlet_afe_tail_bound(mag_t res, const fmpq_t sd2, slong N, ulong q, int parity)
{
    mag_t pi_n2_q, t, u;
    fmpz_t sprime;

    mag_init(pi_n2_q);
    mag_init(t);
    mag_init(u);
    fmpz_init(sprime);

    /* pi_n2_q = pi * N^2 / q (lower bound) */
    mag_const_pi_lower(pi_n2_q);
    mag_mul_ui_lower(pi_n2_q, pi_n2_q, N);
    mag_mul_ui_lower(pi_n2_q, pi_n2_q, N);
    mag_set_ui(t, q);
    mag_div_lower(pi_n2_q, pi_n2_q, t);

    /* upper bound for sd2 */
    fmpz_cdiv_q(sprime, fmpq_numref(sd2), fmpq_denref(sd2));

    /* require pi_n2_q > s' */
    mag_set_fmpz(t, sprime);
    if (fmpz_sgn(sprime) > 0 && mag_cmp(pi_n2_q, t) <= 0)
    {
        mag_inf(res);
    }
    else
    {
        mag_expinv(res, pi_n2_q);

        mag_div_ui(res, res, N);
        if (!parity)
            mag_div_ui(res, res, N);

        /* (1 + q/pi) */
        mag_set_ui(t, q);
        mag_const_pi_lower(u);
        mag_div(t, t, u);
        mag_add_ui(t, t, 1);
        mag_mul(res, res, t);

        /* max(1, 2^s') */
        if (fmpz_sgn(sprime) > 0)
            mag_mul_2exp_fmpz(res, res, sprime);

        /* (pi/q)^(s'-1) */
        fmpz_sub_ui(sprime, sprime, 1);
        if (fmpz_sgn(sprime) >= 0)
        {
            mag_const_pi(t);
            mag_set_ui_lower(u, q);
            mag_div(t, t, u);
            mag_pow_fmpz(t, t, sprime);
        }
        else
        {
            mag_const_pi_lower(t);
            mag_set_ui(u, q);
            mag_div_lower(t, t, u);

            fmpz_neg(sprime, sprime);
            mag_pow_fmpz_lower(t, t, sprime);
            mag_inv(t, t);
        }

        mag_mul(res, res, t);
    }

    mag_clear(pi_n2_q);
    mag_clear(t);
    mag_clear(u);
    fmpz_clear(sprime);
}


void
acb_dirichlet_fmpq_sum_afe(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, const mag_t abs_tol, slong prec)
{
    slong NN, n, start_bits, bits, wp, wp2, gamma_cached_prec;
    mag_t AE, err, abs_tol_gamma;
    arb_t ns, t, u, v, z, z0, z1, x, x2, Ga, Gz1, Gz0, expmz0, z0_prevn, Gz0_prevn, expmz0_prevn;
    acb_t c;
    fmpq_t s2;
    int parity, gamma_singular;
    ulong q;

    double abs_tol_mag;
    double gammainc_mag, gamma_mag, ns_mag;
    double aa, zz;

#if VERBOSE
    double t1, t2, t3;
#endif

    mag_init(AE);
    mag_init(err);
    mag_init(abs_tol_gamma);
    arb_init(ns);
    arb_init(t);
    arb_init(u);
    arb_init(v);
    arb_init(z);
    arb_init(z0);
    arb_init(z1);
    arb_init(x);
    arb_init(x2);
    arb_init(Ga);
    arb_init(Gz0);
    arb_init(Gz1);
    arb_init(expmz0);
    arb_init(z0_prevn);
    arb_init(Gz0_prevn);
    arb_init(expmz0_prevn);
    acb_init(c);
    fmpq_init(s2);

    if (G == NULL)
    {
        parity = 0;
        q = 1;
    }
    else
    {
        parity = dirichlet_parity_char(G, chi);
        q = G->q;
    }

    acb_zero(res);

    /* Initial precision for gamma; may have to be increased later. */
    gamma_cached_prec = prec * 1.05 + 30;

    /* s2 = (s+parity)/2 */
    fmpq_add_ui(s2, s, parity);
    fmpq_div_2exp(s2, s2, 1);

    gamma_singular = (fmpz_is_one(fmpq_denref(s2)) && fmpz_sgn(fmpq_numref(s2)) <= 0);

    if (!gamma_singular)
        arb_gamma_fmpq(Ga, s2, gamma_cached_prec);

    for (n = 1; ; n += 1)
    {
#if VERBOSE
        printf("-----------------------------------------------------------\n");
        flint_printf("n = %wd   (s+parity)/2 = %f  z = %f   q = %wu\n", n, fmpq_get_d(s2), 3.1415926535897932 * n * n / q, q);
#endif

        acb_dirichlet_afe_tail_bound(err, s2, n, q, parity);

#if VERBOSE
        printf("  abs_tol = "); mag_printd(abs_tol, 5); printf("\n");
        printf("  truncation error = "); mag_printd(err, 5); printf("\n");
#endif

        if (mag_cmp(err, abs_tol) < 0)
        {
            if (G == NULL || dirichlet_char_is_real(G, chi))
                arb_add_error_mag(acb_realref(res), err);
            else
                acb_add_error_mag(res, err);

            break;
        }

        /* Compute local precision and tolerances. */

        abs_tol_mag = mag_get_d_log2_approx(abs_tol);
        aa = fmpq_get_d(s2);
        zz = 3.1415926535897932385 * n * n / q;

        /* Gamma((s+parity)/2, z)  (want lower bound, to estimate cancellation) */
        gammainc_mag = log_gamma_upper_approx(aa, zz) / log(2);
        /* n^-s */
        ns_mag = -fmpq_get_d(s) * log(n) / log(2);

        /* Want Gamma(a,z) n^-s with abs_tol --> want Gamma(a,z) with abs_tol * n^s */
        mag_set_ui_2exp_si(abs_tol_gamma, 1, abs_tol_mag - ns_mag);

        /* wp = Precision needed sans cancellation. */
        wp = gammainc_mag + ns_mag - abs_tol_mag + 5;
        wp = FLINT_MAX(wp, 30);

        /* wp2 = Precision needed with cancellation. */
        if (gamma_singular)
        {
            /* Max term is roughly n^(-s) * z^((s+parity)/2) * exp(z) */
            wp2 = -abs_tol_mag + ns_mag + aa * log(zz) + zz / log(2) + 5;
            wp2 = FLINT_MAX(wp2, 30);
        }
        else
        {
            /* Estimate of Gamma((s+parity)/2) */
            gamma_mag = ARF_EXP(arb_midref(Ga));
            wp2 = FLINT_MAX(gamma_mag, gammainc_mag) + ns_mag - abs_tol_mag + 5;
            wp2 = FLINT_MAX(wp2, 30);
        }

#if VERBOSE
        printf("  abs_tol_gamma = "); mag_printd(abs_tol_gamma, 5); printf("\n");
        printf("  gamma(a)      = "); arb_printd(Ga, 10); printf("\n");
        printf("  wp = %ld   wp2 = %ld\n", wp, wp2);
#endif

        if (G == NULL)
            acb_one(c);
        else
            acb_dirichlet_chi(c, G, chi, n, wp);

        if (acb_is_zero(c))
            continue;

        arb_const_pi(z, wp2);
        arb_mul_ui(z, z, n, wp2);
        arb_mul_ui(z, z, n, wp2);
        arb_div_ui(z, z, q, wp2);

        start_bits = 32;
        arb_extract_bits(z0, z, start_bits);

        /* Can we use the asymptotic series? */
        NN = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N(AE, s2, z0, abs_tol_gamma);

#if VERBOSE
        t1 = clock();
#endif

        /* For each new point, evaluate from scratch or use continuation? The
           former seems to be faster. */
        if (1)
        {
            if (NN != -1)
            {
                /* OK to use the asymptotic series. */
                _arb_hypgeom_gamma_upper_fmpq_inf_bsplit(Gz0, s2, z0, NN, wp);
                arb_add_error_mag(Gz0, AE);
#if VERBOSE
                flint_printf("  asymptotic series with N = %wd: ", NN); arb_printd(Gz0, 10); printf("\n");
#endif
            }
            else
            {
                /* Otherwise fallback to the series at 0. */
                if (gamma_singular)
                {
                    slong nn;

                    nn = *fmpq_numref(s2);

                    if (COEFF_IS_MPZ(nn))
                    {
                        arb_indeterminate(Gz0);
                    }
                    else
                    {
                        nn = -nn;

                        NN = _arb_hypgeom_gamma_upper_singular_si_choose_N(AE, nn, z0, abs_tol_gamma);
                        _arb_hypgeom_gamma_upper_singular_si_bsplit(Gz0, nn, z0, NN, wp2);
                        arb_add_error_mag(Gz0, AE);

#if VERBOSE
                        flint_printf("  singular series with N = %wd,  z0 = ", NN); arb_printd(z0, 10); printf("   ");
                        mag_printd(AE, 10); printf("\n");
#endif
                    }
                }
                else
                {
                    NN = _arb_hypgeom_gamma_lower_fmpq_0_choose_N(AE, s2, z0, abs_tol_gamma);
                    _arb_hypgeom_gamma_lower_fmpq_0_bsplit(Gz0, s2, z0, NN, wp2);
                    arb_add_error_mag(Gz0, AE);

#if VERBOSE
                    flint_printf("  lower series with N = %wd,  z0 = ", NN); arb_printd(z0, 10); printf("   ");
                    mag_printd(AE, 10); printf("\n");
#endif

                    while (mag_cmp(arb_radref(Ga), abs_tol_gamma) > 0)
                    {
                        gamma_cached_prec *= 2;
                        arb_gamma_fmpq(Ga, s2, gamma_cached_prec);
                    }
#if VERBOSE
                    flint_printf("  lower series with N = %wd: ", NN); arb_printd(Gz0, 10); printf("\n");
#endif
                    arb_sub(Gz0, Ga, Gz0, wp);
#if VERBOSE
                    flint_printf("  G(a) - lower series: "); arb_printd(Gz0, 10); printf("\n");
#endif
                }
            }
        }
        else
        {
            _arb_gamma_upper_fmpq_step_bsplit(Gz0, s2, z0_prevn, z0, Gz0_prevn, expmz0_prevn, abs_tol_gamma, wp);
        }

#if VERBOSE
        printf("  Gz0 = "); arb_printd(Gz0, 10); printf("\n");
#endif

        if (n == 1)
        {
            arb_neg(expmz0, z0);
            arb_exp(expmz0, expmz0, wp);
        }
        else
        {
            arb_sub(t, z0_prevn, z0, wp);
            arb_exp(t, t, wp);
            arb_mul(expmz0, expmz0_prevn, t, wp);
        }

        arb_set(z0_prevn, z0);
        arb_set(expmz0_prevn, expmz0);
        arb_set(Gz0_prevn, Gz0);

#if VERBOSE
        t2 = clock();
#endif

        /* Bit-burst steps */
        for (bits = start_bits * 2; bits < wp / 8; bits *= 2)
        {
            arb_extract_bits(z1, z, bits);
            _arb_gamma_upper_fmpq_step_bsplit(Gz1, s2, z0, z1, Gz0, expmz0, abs_tol_gamma, wp);
            arb_sub(t, z0, z1, wp);
            arb_exp(t, t, wp);
            arb_mul(expmz0, expmz0, t, wp);
            arb_set(Gz0, Gz1);
            arb_set(z0, z1);
        }

        /* Final step, including error bound */
        _arb_gamma_upper_fmpq_step_bsplit(Gz1, s2, z0, z, Gz0, expmz0, abs_tol_gamma, wp);
        arb_set(Gz0, Gz1);

#if VERBOSE
        printf("  Gz0 = "); arb_printd(Gz0, 10); printf("   tol  "); mag_printd(abs_tol_gamma, 5); printf("\n");
#endif

        /* Multiply by prefactor n^-s */
        arb_set_ui(ns, n);
        arb_pow_fmpq(ns, ns, s, wp);
        arb_div(Gz0, Gz0, ns, wp);

#if VERBOSE
        printf("  1/n^s = "); arb_printn(ns, 5, ARB_STR_NO_RADIUS); printf("\n");
        printf("  Gz0 * pre = "); arb_printd(Gz0, 10); printf("   tol  "); mag_printd(abs_tol, 5); printf("\n");
#endif
        acb_addmul_arb(res, c, Gz0, prec);

#if VERBOSE
        printf("  sum = "); acb_printd(res, 10); printf("\n");
        t3 = clock();
        printf("  time: %f, %f\n", (t2 - t1) / CLOCKS_PER_SEC, (t3 - t2) / CLOCKS_PER_SEC);
#endif
    }

    mag_clear(AE);
    mag_clear(err);
    mag_clear(abs_tol_gamma);
    arb_clear(ns);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    arb_clear(z);
    arb_clear(z0);
    arb_clear(z1);
    arb_clear(x);
    arb_clear(x2);
    arb_clear(Ga);
    arb_clear(Gz0);
    arb_clear(Gz1);
    arb_clear(expmz0);
    arb_clear(z0_prevn);
    arb_clear(Gz0_prevn);
    arb_clear(expmz0_prevn);
    acb_clear(c);
    fmpq_clear(s2);
}

#define PI 3.1415926535897932385
#define INV_LOG2 1.4426950408889634074;

/* max(pi/q,s/2)**(s/2-1) * exp(-max(pi/q,s/2))
static double
estimate_sum1_mag(double s, double q)
{
    return ((0.5 * s - 1) * log(FLINT_MAX(PI / q, 0.5 * s)) - FLINT_MAX(PI / q, 0.5 * s)) * INV_LOG2;
}
 */

void
acb_dirichlet_l_fmpq_afe(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    arb_t t;
    acb_t S1, S2, w;
    fmpq_t s2;
    mag_t tol1, tol2;
    double ds, m1, m2, m2pre;
    slong origprec, prec1, prec2;
    ulong q;
    int parity;

    /* Todo: implement decomposition for imprimitive characters. */
    if (G != NULL && !dirichlet_char_is_primitive(G, chi))
    {
        acb_indeterminate(res);
        return;
    }

    q = (G == NULL) ? 1 : G->q;
    parity = (G == NULL) ? 0 : dirichlet_parity_char(G, chi);

    /* Division by gamma((s+parity)/2) at a pole. */
    if (fmpz_is_one(fmpq_denref(s)))
    {
        const fmpz * n = fmpq_numref(s);

        if ((parity == 0 && fmpz_sgn(n) <= 0 && fmpz_is_even(n)) ||
            (parity == 1 && fmpz_sgn(n) < 0 && fmpz_is_odd(n)))
        {
            /* Special case for zeta */
            if (q == 1 && fmpz_is_zero(n))
                acb_set_d(res, -0.5);
            else
                acb_zero(res);
            return;
        }
    }

    origprec = prec;
    prec = prec * 1.001 + 2 * FLINT_BIT_COUNT(q);

    acb_init(S1);
    acb_init(S2);
    acb_init(w);
    arb_init(t);
    fmpq_init(s2);
    mag_init(tol1);
    mag_init(tol2);

    ds = fmpq_get_d(s);

    m1 = log_gamma_upper_approx(0.5 * (ds + parity), PI / q) * INV_LOG2;
    m2 = log_gamma_upper_approx(0.5 * (1.0 - ds + parity), PI / q) * INV_LOG2;

    m2pre = (ds - 0.5) * log(PI / q) * INV_LOG2;

    mag_one(tol1);
    mag_mul_2exp_si(tol1, tol1, FLINT_MAX(m1, m2 + m2pre) - prec);

    mag_mul_2exp_si(tol2, tol1, -m2pre);

    prec1 = prec - (FLINT_MAX(m1, m2 + m2pre) - m1);
    prec1 = FLINT_MAX(prec1, 32);

    prec2 = prec - (FLINT_MAX(m1, m2 + m2pre) - (m2 + m2pre));
    prec2 = FLINT_MAX(prec2, 32);

#if VERBOSE
    printf("mag1 = %ld   mag2 = %ld   mag2 + pre = %ld    prec, prec1, prec2 = %ld, %ld, %ld\n",
        (slong) m1, (slong) m2, (slong) (m2 + m2pre), prec, prec1, prec2);
    printf("tol1 = %ld   tol2 = %ld\n", MAG_EXP(tol1), MAG_EXP(tol2));

#endif

    acb_dirichlet_fmpq_sum_afe(S1, s, G, chi, tol1, prec1);

#if VERBOSE
    printf("=====================================================\n");
    printf("S1 = "); acb_printd(S1, 20); printf("  estimate = "); printf(" %g\n", ldexp(1.0, m1));
    printf("=====================================================\n");
#endif

    if (q == 1 && fmpz_is_one(fmpq_numref(s)) && fmpz_equal_ui(fmpq_denref(s), 2))
    {
        acb_mul_2exp_si(res, S1, 1);
    }
    else
    {
        /* rootnum (pi/q)^(s-1/2) sum(1-s) */
        if (fmpz_is_one(fmpq_numref(s)) && fmpz_equal_ui(fmpq_denref(s), 2))
        {
            acb_conj(S2, S1);
        }
        else
        {
            fmpq_sub_ui(s2, s, 1);
            fmpq_neg(s2, s2);
            acb_dirichlet_fmpq_sum_afe(S2, s2, G, chi, tol2, prec2);
            acb_conj(S2, S2);
        }

#if VERBOSE
        printf("=====================================================\n");
        printf("S1 = "); acb_printd(S1, 20); printf("  estimate = "); printf(" %g\n", ldexp(1.0, m1));
        printf("S2 = "); acb_printd(S2, 20); printf("  estimate = "); printf(" %g\n", ldexp(1.0, m2));
        printf("=====================================================\n");
#endif

        arb_const_pi(t, prec);
        arb_div_ui(t, t, q, prec);
        fmpq_set_si(s2, 1, 2);
        fmpq_sub(s2, s, s2);
        arb_pow_fmpq(t, t, s2, prec);
        acb_mul_arb(S2, S2, t, prec);

        if (q != 1)
        {
            acb_dirichlet_root_number2(w, G, chi, prec);
            acb_mul(S2, S2, w, prec);
        }

#if VERBOSE
        printf("S2 * prefactor = "); acb_printd(S2, 20); printf("  estimate = "); printf(" %g\n", ldexp(1.0, m2 + m2pre));
#endif

        acb_add(res, S1, S2, prec);
    }

    /* add pi^(s/2) / (s (s-1)) */
    if (q == 1)
    {
        arb_const_pi(t, prec);
        fmpq_div_2exp(s2, s, 1);
        arb_pow_fmpq(t, t, s2, prec);

        fmpq_sub_ui(s2, s, 1);
        fmpq_mul(s2, s2, s);
        arb_div_fmpz(t, t, fmpq_numref(s2), prec);
        arb_mul_fmpz(t, t, fmpq_denref(s2), prec);

        acb_add_arb(res, res, t, prec);
    }

    /* divide by gamma((s+parity)/2) */
    fmpq_add_ui(s2, s, parity);
    fmpq_div_2exp(s2, s2, 1);
    arb_gamma_fmpq(t, s2, prec);
    acb_div_arb(res, res, t, prec);

    acb_set_round(res, res, origprec);

    acb_clear(S1);
    acb_clear(S2);
    acb_clear(w);
    arb_clear(t);
    fmpq_clear(s2);
    mag_clear(tol1);
    mag_clear(tol2);
}
