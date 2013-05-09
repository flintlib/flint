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
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "ulong_extras.h"
#include "arith.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "profiler.h"


#define DOUBLE_PREC 53
#define PI 3.141592653589793238462643
#define INV_LOG2 (1.44269504088896340735992468 + 1e-12)
#define HRR_A (1.1143183348516376904 + 1e-12)  /* 44*pi^2/(225*sqrt(3)) */
#define HRR_B (0.0592384391754448833 + 1e-12)  /* pi*sqrt(2)/75 */
#define HRR_C (2.5650996603237281911 + 1e-12)  /* pi*sqrt(2/3) */
#define HRR_D (1.2424533248940001551 + 1e-12)  /* log(2) + log(3)/2 */


#define PI_USE_CHUDNOVSKY 1
#define PI_CHUDNOVSKY_CUTOFF 1000000

#define VERBOSE 0


static double
partitions_remainder_bound(double n, double terms)
{
    return HRR_A/sqrt(terms)
            + HRR_B*sqrt(terms/(n-1)) * sinh(HRR_C * sqrt(n)/terms);
}

/* Crude upper bound, sufficient to estimate the precision */
static double
log_sinh(double x)
{
    if (x > 4)
        return x;
    else
        return log(x) + x*x*(1/6.);
}

static double
partitions_remainder_bound_log2(double n, double N)
{
    double t1, t2;

    t1 = log(HRR_A) - 0.5*log(N);
    t2 = log(HRR_B) + 0.5*(log(N) - log(n-1)) + log_sinh(HRR_C * sqrt(n)/N);

    return (FLINT_MAX(t1, t2) + 1) * INV_LOG2;
}

len_t
partitions_needed_terms(ulong n)
{
    len_t N;
    for (N = 1; partitions_remainder_bound_log2(n, N) > 10; N++);
    for ( ; partitions_remainder_bound(n, N) > (n > 1500 ? 0.25 : 1); N++);
    return N;
}

static double
partitions_term_bound(double n, double k)
{
    return ((PI*sqrt(24*n-1) / (6.0*k)) + HRR_D - log(24.0*n-1) + 0.5*log(k)) * INV_LOG2;
}

/* Bound number of prime factors in k */
static mp_limb_t primorial_tab[] = {
    1, 2, 6, 30, 210, 2310, 30030, 510510, 9699690, 223092870,
#if FLINT64
    6469693230UL, 200560490130UL, 7420738134810UL, 304250263527210UL,
    13082761331670030UL, 614889782588491410UL
#endif
};

static __inline__ int
bound_primes(ulong k)
{
    int i;

    for (i = 0; i < sizeof(primorial_tab) / sizeof(mp_limb_t); i++)
        if (k <= primorial_tab[i])
            return i;

    return i;
}


static __inline__ len_t
log2_ceil(double x)
{
    /* ceil(log2(n)) = bitcount(n-1);
       this is too large if x is a power of two */
    return FLINT_BIT_COUNT((len_t) x);
}

static len_t
partitions_prec_bound(ulong n, len_t k, len_t N)
{
    len_t prec;

    prec = partitions_term_bound(n, k);
    prec += log2_ceil(8 * N * (26 * (sqrt(n) / k) + 7 * bound_primes(k) + 22));

    return prec;
}

double
cos_pi_pq(mp_limb_signed_t p, mp_limb_signed_t q)
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

void
mpfr_sqrt_z(mpfr_t x, mpz_t z, mpfr_rnd_t rnd)
{
    if (mpz_fits_ulong_p(z))
        mpfr_sqrt_ui(x, mpz_get_ui(z), rnd);
    else
    {
        mpfr_set_z(x, z, rnd);
        mpfr_sqrt(x, x, rnd);
    }
}

void
mpfr_set_fmpz(mpfr_t c, const fmpz_t b)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_set_z(c, COEFF_TO_PTR(*b), MPFR_RNDN);
    else
        mpfr_set_si(c, *b, MPFR_RNDN);
}

void
mpfr_mul_fmpz(mpfr_t c, mpfr_srcptr a, const fmpz_t b)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_mul_z(c, a, COEFF_TO_PTR(*b), MPFR_RNDN);
    else
        mpfr_mul_si(c, a, *b, MPFR_RNDN);
}

void
mpfr_add_fmpz(mpfr_t c, mpfr_srcptr a, const fmpz_t b)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_add_z(c, a, COEFF_TO_PTR(*b), MPFR_RNDN);
    else
        mpfr_add_si(c, a, *b, MPFR_RNDN);
}


void
_fmpz_poly_evaluate_mpfr(mpfr_t res, const fmpz * f, len_t len,
                           const mpfr_t a)
{
    if (len == 0)
        mpfr_set_ui(res, 0, MPFR_RNDN);
    else if (len == 1)
        mpfr_set_fmpz(res, f);
    else
    {
        len_t i = len - 1;
        mpfr_t t;
        mpfr_init2(t, mpfr_get_prec(res));
        mpfr_set_fmpz(res, f + i);
        for (i = len - 2; i >= 0; i--)
        {
            mpfr_mul(t, res, a, MPFR_RNDN);
            mpfr_add_fmpz(res, t, f + i);
        }
        mpfr_clear(t);
    }
}

void
fmpz_poly_evaluate_mpfr(mpfr_t res, const fmpz_poly_t f, const mpfr_t a)
{
    if (res == a)
    {
        mpfr_t t;
        mpfr_init2(t, mpfr_get_prec(res));
        _fmpz_poly_evaluate_mpfr(t, f->coeffs, f->length, a);
        mpfr_swap(res, t);
        mpfr_clear(t);
    }
    else
    {
        _fmpz_poly_evaluate_mpfr(res, f->coeffs, f->length, a);
    }
}

void
findroot(mpfr_t x, fmpz_poly_t poly, double x0)
{
    len_t i, prec, initial_prec, target_prec, guard_bits;
    len_t precs[FLINT_BITS];
    fmpz_poly_t poly2;
    mpfr_t t, u, xn;

    initial_prec = 48;
    target_prec = mpfr_get_prec(x) + 32;

    mpfr_init2(t, 53);
    mpfr_init2(u, 53);
    mpfr_init2(xn, 53);
    mpfr_set_d(xn, x0, MPFR_RNDN);

    fmpz_poly_init(poly2);
    fmpz_poly_derivative(poly2, poly);
    guard_bits = fmpz_poly_max_bits(poly2);
    guard_bits = FLINT_ABS(guard_bits);

    for (i = 0, prec = target_prec; prec >= initial_prec; i++)
    {
        precs[i] = prec;
        prec = prec / 2 + 8;
    }

    for (i--; i >= 0; i--)
    {
        mpfr_set_prec(t, precs[i] + guard_bits);
        mpfr_set_prec(u, precs[i] + guard_bits);
        mpfr_prec_round(xn, precs[i], MPFR_RNDN);
        fmpz_poly_evaluate_mpfr(t, poly, xn);
        fmpz_poly_evaluate_mpfr(u, poly2, xn);
        mpfr_div(t, t, u, MPFR_RNDN);
        mpfr_sub(xn, xn, t, MPFR_RNDN);
    }

    mpfr_set(x, xn, MPFR_RNDN);

    fmpz_poly_clear(poly2);
    mpfr_clear(t);
    mpfr_clear(u);
    mpfr_clear(xn);
}

void cos_minpoly(fmpz_poly_t poly, len_t p, len_t q)
{
    if (p % 2 == 0)
        arith_cos_minpoly(poly, q);
    else
        arith_cos_minpoly(poly, 2 * q);
}

int use_newton(len_t prec, len_t q)
{
    if (q < 250 && prec > 400 + 4*q*q)
        return 1;
    return 0;
}

void mpfr_cos_pi_pq(mpfr_t t, mp_limb_signed_t p, mp_limb_signed_t q)
{
    /* Force 0 <= p < q */
    p = FLINT_ABS(p);
    p %= (2 * q);
    if (p >= q)
        p = 2 * q - p;

    if (use_newton(mpfr_get_prec(t), q))
    {
        fmpz_poly_t poly;
        len_t d;
        fmpz_poly_init(poly);
        d = n_gcd(q, p);
        q /= d;
        p /= d;
        cos_minpoly(poly, p, q);
        findroot(t, poly, cos(3.1415926535897932385 * p / q));
        fmpz_poly_clear(poly);
    }
    else
    {
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
        v = n_gcd_full(prod->sqrt_p, prod->sqrt_q);
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
sinh_cosh_divk_precomp(mpfr_t sh, mpfr_t ch, mpfr_t ex, len_t k)
{
    mpfr_t t;
    mpfr_root(ch, ex, k, MPFR_RNDN);
    /* The second term doesn't need full precision,
       but this doesn't affect performance that much... */
    mpfr_init2(t, mpfr_get_prec(ch));
    mpfr_ui_div(t, 1, ch, MPFR_RNDN);
    mpfr_sub(sh, ch, t, MPFR_RNDN);
    mpfr_add(ch, ch, t, MPFR_RNDN);
    mpfr_div_2exp(ch, ch, 1, MPFR_RNDN);
    mpfr_div_2exp(sh, sh, 1, MPFR_RNDN);
    mpfr_clear(t);
}


void
_arith_number_of_partitions_mpfr(mpfr_t x, ulong n, len_t N0, len_t N)
{
    trig_prod_t prod;
    mpfr_t acc, C, t1, t2, t3, t4, exp1;
    mpz_t n24;
    double Cd;
    len_t k, prec, guard_bits;
#if VERBOSE
    timeit_t t0;
#endif

    if (n <= 2)
    {
        mpfr_set_ui(x, FLINT_MAX(1, n), MPFR_RNDN);
        return;
    }

    /* Compute initial precision */
    guard_bits = 2 * FLINT_BIT_COUNT(N) + 32;
    prec = partitions_remainder_bound_log2(n, N0) + guard_bits;
    prec = FLINT_MAX(prec, DOUBLE_PREC);

    mpfr_set_prec(x, prec);
    mpfr_init2(acc, prec);
    mpfr_init2(C, prec);
    mpfr_init2(t1, prec);
    mpfr_init2(t2, prec);
    mpfr_init2(t3, prec);
    mpfr_init2(t4, prec);

    mpfr_set_ui(x, 0, MPFR_RNDN);
    mpfr_set_ui(acc, 0, MPFR_RNDN);

    mpz_init(n24);
    mpz_set_ui(n24, n);
    mpz_mul_ui(n24, n24, 24);
    mpz_sub_ui(n24, n24, 1);

#if VERBOSE
    timeit_start(t0);
#endif

    /* C = (pi/6)*sqrt(24*n-1) */

    if (PI_USE_CHUDNOVSKY && prec > PI_CHUDNOVSKY_CUTOFF)
        mpfr_pi_chudnovsky(t1, MPFR_RNDN);
    else
        mpfr_const_pi(t1, MPFR_RNDN);

    mpfr_sqrt_z(t2, n24, MPFR_RNDN);
    mpfr_mul(t1, t1, t2, MPFR_RNDN);
    mpfr_div_ui(C, t1, 6, MPFR_RNDN);
    Cd = mpfr_get_d(C, MPFR_RNDN);

    mpfr_init2(exp1, prec);
    mpfr_exp(exp1, C, prec);

#if VERBOSE
    timeit_stop(t0);
    printf("TERM 1: %ld ms\n", t0->cpu);
#endif

    for (k = N0; k <= N; k++)
    {
        trig_prod_init(prod);
        arith_hrr_expsum_factored(prod, k, n % k);

        if (prod->prefactor != 0)
        {
            if (prec > DOUBLE_PREC)
            {
                prec = partitions_prec_bound(n, k, N);

                mpfr_set_prec(t1, prec);
                mpfr_set_prec(t2, prec);
                mpfr_set_prec(t3, prec);
                mpfr_set_prec(t4, prec);
            }

            /* Compute A_k(n) * sqrt(3/k) * 4 / (24*n-1) */
            prod->prefactor *= 4;
            prod->sqrt_p *= 3;
            prod->sqrt_q *= k;
            eval_trig_prod(t1, prod);
            mpfr_div_z(t1, t1, n24, MPFR_RNDN);

            /* Multiply by (cosh(z) - sinh(z)/z) where z = C / k */
            if (prec <= DOUBLE_PREC)
            {
                double z = Cd / k;
                mpfr_mul_d(t1, t1, cosh(z) - sinh(z)/z, MPFR_RNDN);
            }
            else
            {
                mpfr_div_ui(t2, C, k, MPFR_RNDN);

                if (k < 35)
                    sinh_cosh_divk_precomp(t3, t4, exp1, k);
                else
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

    mpz_clear(n24);
    mpfr_clear(acc);
    mpfr_clear(exp1);
    mpfr_clear(C);
    mpfr_clear(t1);
    mpfr_clear(t2);
    mpfr_clear(t3);
    mpfr_clear(t4);
}

void
arith_number_of_partitions_mpfr(mpfr_t x, ulong n)
{
    _arith_number_of_partitions_mpfr(x, n, 1, partitions_needed_terms(n));
}
