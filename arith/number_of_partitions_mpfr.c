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

#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "math.h"
#include "ulong_extras.h"
#include "arith.h"

/* Parameters for evaluation of the Hardy-Ramanujan-Rademacher formula.
   See the documentation section for additional information. */

#define HRR_GUARD_BITS 15
#define HRR_CUTOFF 0.4      /* Quit when error < cutoff. Should be < 0.5 */

#define ROUNDUP 1e-13
#define HRR_A (1.1143183348516376904 + ROUNDUP)  /* 44*pi^2/(225*sqrt(3)) */
#define HRR_B (0.0592384391754448833 + ROUNDUP)  /* pi*sqrt(2)/75 */
#define HRR_C (2.5650996603237281911 + ROUNDUP)  /* pi*sqrt(mpf(2)/3) */


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

void
number_of_partitions_mpfr(mpfr_t x, ulong n)
{
    mpfr_t t, u, v, w, C0, C1, C2, C3;
    double s, D0, D1, D2, D3;
    long k, prec;

    if (n <= 2)
    {
        mpfr_set_ui(x, n + (n == 0), MPFR_RNDN);
        return;
    }

    prec = partitions_remainder_bound_log2(n, 1) + HRR_GUARD_BITS;

    mpfr_set_prec(x, prec);
    mpfr_set_ui(x, 0UL, GMP_RNDN);

    mpfr_init2(t, prec);
    mpfr_init2(u, prec);
    mpfr_init2(v, prec);
    mpfr_init2(w, prec);
    mpfr_init2(C0, prec);
    mpfr_init2(C1, prec);
    mpfr_init2(C2, prec);
    mpfr_init2(C3, prec);

    /* C0 = sqrt(24*n-1) */
    mpfr_set_ui(C0, n, GMP_RNDN);
    mpfr_mul_ui(C0, C0, 24UL, GMP_RNDN);
    mpfr_sub_ui(C0, C0, 1UL, GMP_RNDN);

    /* C2 = 4*sqrt(3) / (24*n-1) */
    mpfr_sqrt_ui(C2, 48UL, GMP_RNDN);
    mpfr_div(C2, C2, C0, GMP_RNDN);

    mpfr_sqrt(C0, C0, GMP_RNDN);

    /* C1 = C0 * pi / 6 */
    mpfr_const_pi(C1, GMP_RNDN);
    mpfr_mul(C1, C1, C0, GMP_RNDN);
    mpfr_div_ui(C1, C1, 6UL, GMP_RNDN);

    /* C3 = C2 / C1 */
    mpfr_div(C3, C2, C1, GMP_RNDN);

    for (k = 1; ; k++)
    {
        prec = partitions_remainder_bound_log2(n, k) + HRR_GUARD_BITS;

        /* printf("%ld: %ld\n", k, prec); */

        mpfr_set_prec(t, prec);  /* raw! */
        mpfr_set_prec(u, prec);
        mpfr_set_prec(v, prec);
        mpfr_set_prec(w, prec);

        /* t = A(h,k) */
        dedekind_cosine_sum_mpfr_fast(t, k, n);

        if (!mpfr_zero_p(t))  /* XXX: is this what breaks */
        {
            /* sum += A(h,k) * sqrt(k)*(C2*cosh(C1/k)/k - C3*sinh(C1/k) */
            mpfr_div_ui(u, C1, k, GMP_RNDN);
            mpfr_sinh_cosh(w, v, u, GMP_RNDN);
            mpfr_mul(v, v, C2, GMP_RNDN);
            mpfr_div_ui(v, v, k, GMP_RNDN);
            mpfr_mul(w, w, C3, GMP_RNDN);
            mpfr_sub(v, v, w, GMP_RNDN);
            mpfr_sqrt_ui(u, k, GMP_RNDN);
            mpfr_mul(v, v, u, GMP_RNDN);
            mpfr_mul(t, t, v, GMP_RNDN);
            mpfr_add(x, x, t, GMP_RNDN);
        }

        if (prec < FLINT_D_BITS)
            break;
    }

    mpfr_set_prec(t, FLINT_D_BITS + HRR_GUARD_BITS);
    mpfr_set_ui(t, 0UL, GMP_RNDN);

    D0 = mpfr_get_d(C0, GMP_RNDN);
    D1 = mpfr_get_d(C1, GMP_RNDN);
    D2 = mpfr_get_d(C2, GMP_RNDN);
    D3 = mpfr_get_d(C3, GMP_RNDN);

    for (k++; ; k++)
    {
        s = dedekind_cosine_sum_d(k, n);
        s = s * sqrt(k) * (D2*cosh(D1/k)/k - D3*sinh(D1/k));
        mpfr_add_d(t, t, s, GMP_RNDN);

        if (partitions_remainder_bound(n, k) < HRR_CUTOFF)
            break;
    }

    mpfr_add(x, x, t, GMP_RNDN);
    mpfr_rint(x, x, GMP_RNDN);

    mpfr_clear(t);
    mpfr_clear(u);
    mpfr_clear(v);
    mpfr_clear(w);
    mpfr_clear(C0);
    mpfr_clear(C1);
    mpfr_clear(C2);
    mpfr_clear(C3);
}
