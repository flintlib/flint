/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

    Inspired by code written for Sage by Jonathan Bober.

******************************************************************************/

#include <math.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "ulong_extras.h"
#include "arith.h"

#define DOUBLE_PREC 53
#define PI 3.141592653589793238462643
#define HRR_A (1.1143183348516376904 + 1e-12)  /* 44*pi^2/(225*sqrt(3)) */
#define HRR_B (0.0592384391754448833 + 1e-12)  /* pi*sqrt(2)/75 */
#define HRR_C (2.5650996603237281911 + 1e-12)  /* pi*sqrt(2/3) */


static double partitions_remainder_bound(double n, double terms)
{
    return HRR_A/sqrt(terms)
            + HRR_B*sqrt(terms/(n-1)) * sinh(HRR_C * sqrt(n)/terms);
}

/* Crude upper bound, sufficient to estimate the precision */
static double log_sinh(double x)
{
    if (x > 4)
        return x;
    else
        return log(x) + x*x*(1/6.);
}

static double partitions_remainder_bound_log2(double n, double N)
{
    double t1, t2;

    t1 = log(HRR_A) - 0.5*log(N);
    t2 = log(HRR_B) + 0.5*(log(N) - log(n-1)) + log_sinh(HRR_C * sqrt(n)/N);

    return (FLINT_MAX(t1, t2) + 1) * 1.4426950408889634074;
}

double cos_pi_pq(mp_limb_signed_t p, mp_limb_signed_t q)
{
    /* Force 0 <= p < q */
    p = FLINT_ABS(p);
    p %= (2 * q);
    if (p >= q)
        p = 2 * q - p;

    if (4 * p <= q)
        return cos(p * PI / q);
    else if (4 * p < 3 * q)
        return sin((q - 2*p) * PI / (2 * q));
    else
        return -cos((q - p) * PI / q);
}

void mpfr_cos_pi_pq(mpfr_t t, mp_limb_signed_t p, mp_limb_signed_t q)
{
    /* Force 0 <= p < q */
    p = FLINT_ABS(p);
    p %= (2 * q);
    if (p >= q)
        p = 2 * q - p;

    mpfr_const_pi(t, MPFR_RNDN);

    if (4 * p <= q)
    {
        mpfr_mul_si(t, t, p, MPFR_RNDN);
        mpfr_div_ui(t, t, q, MPFR_RNDN);
        mpfr_cos(t, t, MPFR_RNDN);
    }
    else if (4 * p < 3 * q)
    {
        mpfr_mul_si(t, t, q - 2*p, MPFR_RNDN);
        mpfr_div_ui(t, t, 2 * q, MPFR_RNDN);
        mpfr_sin(t, t, MPFR_RNDN);
    }
    else
    {
        mpfr_mul_si(t, t, q - p, MPFR_RNDN);
        mpfr_div_ui(t, t, q, MPFR_RNDN);
        mpfr_cos(t, t, MPFR_RNDN);
        mpfr_neg(t, t, MPFR_RNDN);
    }
}

void
eval_trig_prod(mpfr_t sum, trig_prod_t prod)
{
    int i;

    if (prod->prefactor == 0)
    {
        mpfr_set_ui(sum, 0UL, MPFR_RNDN);
        return;
    }

    if (mpfr_get_prec(sum) <= DOUBLE_PREC)
    {
        double s;
        s = prod->prefactor * sqrt((double)prod->sqrt_p/(double)prod->sqrt_q);
        for (i = 0; i < prod->n; i++)
            s *= cos_pi_pq(prod->cos_p[i], prod->cos_q[i]);
        mpfr_set_d(sum, s, MPFR_RNDN);
    }
    else
    {
        mp_limb_t v;
        mpfr_t t;

        mpfr_init2(t, mpfr_get_prec(sum));
        mpfr_set_si(sum, prod->prefactor, MPFR_RNDN);
        v = n_gcd(FLINT_MAX(prod->sqrt_p, prod->sqrt_q),
                  FLINT_MIN(prod->sqrt_p, prod->sqrt_q));
        prod->sqrt_p /= v;
        prod->sqrt_q /= v;

        if (prod->sqrt_p != 1)
        {
            mpfr_sqrt_ui(t, prod->sqrt_p, MPFR_RNDN);
            mpfr_mul(sum, sum, t, MPFR_RNDN);
        }

        if (prod->sqrt_q != 1)
        {
            mpfr_sqrt_ui(t, prod->sqrt_q, MPFR_RNDN);
            mpfr_div(sum, sum, t, MPFR_RNDN);
        }

        for (i = 0; i < prod->n; i++)
        {
            mpfr_cos_pi_pq(t, prod->cos_p[i], prod->cos_q[i]);
            mpfr_mul(sum, sum, t, MPFR_RNDN);
        }

        mpfr_clear(t);
    }
}


void
number_of_partitions_mpfr(mpfr_t x, ulong n)
{
    trig_prod_t prod;
    mpfr_t acc, C, t1, t2, t3, t4;
    double Cd;
    long k, N, prec, guard_bits, mag_bound;

    if (n <= 2)
    {
        mpfr_set_ui(x, FLINT_MAX(1, n), MPFR_RNDN);
        return;
    }

    /* Compute number of needed terms */
    for (N = 1; partitions_remainder_bound_log2(n, N) > 10; N += 10);
    for ( ; partitions_remainder_bound(n, N) > 0.25; N += 10);

    /* Compute initial precision */
    guard_bits = 2 * FLINT_BIT_COUNT(N) + 32;
    prec = partitions_remainder_bound_log2(n, 1) + guard_bits;
    prec = FLINT_MAX(prec, 53);

    mpfr_set_prec(x, prec);
    mpfr_init2(acc, prec);
    mpfr_init2(C, prec);
    mpfr_init2(t1, prec);
    mpfr_init2(t2, prec);
    mpfr_init2(t3, prec);
    mpfr_init2(t4, prec);

    mpfr_set_ui(x, 0, MPFR_RNDN);
    mpfr_set_ui(acc, 0, MPFR_RNDN);

    /* C = (pi/6)*sqrt(24*n-1) */
    mpfr_const_pi(t1, MPFR_RNDN);
    mpfr_sqrt_ui(t2, 24*n - 1, MPFR_RNDN);
    mpfr_mul(t1, t1, t2, MPFR_RNDN);
    mpfr_div_ui(C, t1, 6, MPFR_RNDN);
    Cd = mpfr_get_d(C, MPFR_RNDN);

    for (k = 1; k <= N; k++)
    {
        trig_prod_init(prod);
        dedekind_cosine_sum_factored(prod, k, n % k);

        if (prod->prefactor != 0)
        {
            mag_bound = partitions_remainder_bound_log2(n, k);
            guard_bits = (long) FLINT_BIT_COUNT(n) / 2 - 
                        ((long) FLINT_BIT_COUNT(k));
            guard_bits = FLINT_MAX(guard_bits, (long)(FLINT_BIT_COUNT(N)));
            guard_bits += 5;
            prec = mag_bound + guard_bits;
            prec = FLINT_MAX(prec, DOUBLE_PREC);

            mpfr_set_prec(t1, prec);
            mpfr_set_prec(t2, prec);
            mpfr_set_prec(t3, prec);
            mpfr_set_prec(t4, prec);

            /* Compute A_k(n) * sqrt(3/k) * 4 / (24*n-1) */
            prod->prefactor *= 4;
            prod->sqrt_p *= 3;
            prod->sqrt_q *= k;
            eval_trig_prod(t1, prod);
            mpfr_div_ui(t1, t1, 24*n - 1, MPFR_RNDN);

            /* Multiply by (cosh(z) - sinh(z)/z) where z = C / k*/
            if (prec <= DOUBLE_PREC)
            {
                double z = Cd / k;
                mpfr_mul_d(t1, t1, cosh(z) - sinh(z)/z, MPFR_RNDN);
            }
            else
            {
                mpfr_div_ui(t2, C, k, MPFR_RNDN);
                mpfr_sinh_cosh(t3, t4, t2, MPFR_RNDN);
                mpfr_div(t3, t3, t2, MPFR_RNDN);
                mpfr_sub(t2, t4, t3, MPFR_RNDN);
                mpfr_mul(t1, t1, t2, MPFR_RNDN);
            }

            /* Add to accumulator */
            mpfr_add(acc, acc, t1, MPFR_RNDN);
            if (mpfr_get_prec(acc) > 2 * prec + 32)
            {
                mpfr_add(x, x, acc, MPFR_RNDN);
                mpfr_set_prec(acc, prec + 32);
                mpfr_set_ui(acc, 0, MPFR_RNDN);
            }
        }
    }

    mpfr_add(x, x, acc, MPFR_RNDN);
    mpfr_rint(x, x, MPFR_RNDN);

    mpfr_clear(acc);
    mpfr_clear(C);
    mpfr_clear(t1);
    mpfr_clear(t2);
    mpfr_clear(t3);
    mpfr_clear(t4);
}