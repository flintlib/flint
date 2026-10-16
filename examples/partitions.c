/*
    Copyright (C) 2011, 2013, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    FLINT program for demonstrating the integer partition function.

    usage: partitions n [-threads t] [-quiet] [-legacy] [-timing]

    n may be an expression such as 10^12 or 2^40+1.  By default p(n) is
    computed by mp_real_partitions_hrr; with -legacy, by the previous
    implementation of the Hardy-Ramanujan-Rademacher formula in arb ball
    arithmetic, which is kept below as an example of arb programming
    (the sum of the terms with tight precision estimates, the factored
    exponential sums, exp(C/k) as roots of exp(C) for the first terms,
    and a cascade of accumulators to avoid quadratic cost).
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/thread_support.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>
#include <flint/arb.h>
#include <flint/arith.h>
#include <flint/partitions.h>
#include <flint/gr.h>
#include <flint/profiler.h>

/* ---- the legacy arb implementation ------------------------------------ */


#define DOUBLE_CUTOFF 40
#define DOUBLE_ERR 1e-12

#define DOUBLE_PREC 53
#define MIN_PREC 20
#define PI 3.141592653589793238462643
#define INV_LOG2 (1.44269504088896340735992468 + 1e-12)
#define HRR_A (1.1143183348516376904 + 1e-12)  /* 44*pi^2/(225*sqrt(3)) */
#define HRR_B (0.0592384391754448833 + 1e-12)  /* pi*sqrt(2)/75 */
#define HRR_C (2.5650996603237281911 + 1e-12)  /* pi*sqrt(2/3) */
#define HRR_D (1.2424533248940001551 + 1e-12)  /* log(2) + log(3)/2 */

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

static slong
legacy_needed_terms(double n)
{
    slong N;
    for (N = 1; partitions_remainder_bound_log2(n, N) > 10; N++);
    for ( ; partitions_remainder_bound(n, N) > 0.4; N++);
    return N;
}

static double
partitions_term_bound(double n, double k)
{
    return ((PI*sqrt(24*n-1) / (6.0*k)) + HRR_D - log(24.0*n-1) + 0.5*log(k)) * INV_LOG2;
}

/* Bound number of prime factors in k */
static ulong primorial_tab[] = {
    1, 2, 6, 30, 210, 2310, 30030, 510510, 9699690, 223092870,
#if FLINT64
    UWORD(6469693230), UWORD(200560490130), UWORD(7420738134810),
    UWORD(304250263527210), UWORD(13082761331670030), UWORD(614889782588491410)
#endif
};

static inline int
bound_primes(ulong k)
{
    unsigned int i;

    for (i = 0; i < sizeof(primorial_tab) / sizeof(ulong); i++)
        if (k <= primorial_tab[i])
            return i;

    return i;
}


static inline slong
log2_ceil(double x)
{
    /* ceil(log2(n)) = bitcount(n-1);
       this is too large if x is a power of two */
    if (x < 0x1p62)
        return FLINT_BIT_COUNT((slong) x);
    /* the cast would overflow from 2^63 (the bound for small k
       reaches it near n = 2.4e17; on x86-64 the undefined conversion
       gave 64 bits instead of 72 at n = 1e20) */
    return (slong) ceil(log2(x) * (1.0 + 1e-15)) + 1;
}

static slong
partitions_prec_bound(double n, slong k, slong N)
{
    slong prec;

    prec = partitions_term_bound(n, k);
    prec += log2_ceil(8 * N * (26 * (sqrt(n) / k) + 7 * bound_primes(k) + 22));

    return prec;
}

static void
eval_trig_prod(arb_t sum, trig_prod_t prod, slong prec)
{
    int i;
    ulong v;
    arb_t t;

    if (prod->prefactor == 0)
    {
        arb_zero(sum);
        return;
    }

    arb_init(t);

    arb_set_si(sum, prod->prefactor);
    v = n_gcd(FLINT_MAX(prod->sqrt_p, prod->sqrt_q),
              FLINT_MIN(prod->sqrt_p, prod->sqrt_q));
    prod->sqrt_p /= v;
    prod->sqrt_q /= v;

    if (prod->sqrt_p != 1)
    {
        arb_sqrt_ui(t, prod->sqrt_p, prec);
        arb_mul(sum, sum, t, prec);
    }

    if (prod->sqrt_q != 1)
    {
        arb_rsqrt_ui(t, prod->sqrt_q, prec);
        arb_mul(sum, sum, t, prec);
    }

    for (i = 0; i < prod->n; i++)
    {
        fmpq_t pq;
        *fmpq_numref(pq) = prod->cos_p[i];
        *fmpq_denref(pq) = prod->cos_q[i];
        arb_cos_pi_fmpq(t, pq, prec);
        arb_mul(sum, sum, t, prec);
    }

    arb_clear(t);
}

static void
sinh_cosh_divk_precomp(arb_t sh, arb_t ch, const arb_t ex, slong k, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_round(t, ex, prec);
    arb_root_ui(ch, t, k, prec);
    /* The second term doesn't need full precision,
       but this doesn't affect performance that much... */
    arb_inv(t, ch, prec);
    arb_sub(sh, ch, t, prec);
    arb_add(ch, ch, t, prec);
    arb_mul_2exp_si(ch, ch, -1);
    arb_mul_2exp_si(sh, sh, -1);
    arb_clear(t);
}

static void
legacy_hrr_sum_arb_range(arb_t x, const fmpz_t n, const arb_t C, const arb_t exp1, const fmpz_t n24, slong start, slong N, slong step, slong prec, slong acc_prec, slong res_prec)
{
    arb_t acc, t1, t2, t3, t4;
    trig_prod_t prod;
    slong k;
    double nd;

    arb_init(acc);
    arb_init(t1);
    arb_init(t2);
    arb_init(t3);
    arb_init(t4);

    nd = fmpz_get_d(n);

    for (k = start; k <= N; k += step)
    {
        trig_prod_init(prod);
        arith_hrr_expsum_factored(prod, k, fmpz_fdiv_ui(n, k));

        if (prod->prefactor != 0)
        {
            if (prec > MIN_PREC)
                prec = partitions_prec_bound(nd, k, N);

            prod->prefactor *= 4;
            prod->sqrt_p *= 3;
            prod->sqrt_q *= k;

            /* Compute A_k(n) * sqrt(3/k) * 4 / (24*n-1) */
            eval_trig_prod(t1, prod, prec);
            arb_div_fmpz(t1, t1, n24, prec);

            /* Multiply by (cosh(z) - sinh(z)/z) where z = C / k */
            arb_set_round(t2, C, prec);
            arb_div_ui(t2, t2, k, prec);

            if (k < 35 && prec > 1000)
                sinh_cosh_divk_precomp(t3, t4, exp1, k, prec);
            else
                arb_sinh_cosh(t3, t4, t2, prec);

            arb_div(t3, t3, t2, prec);
            arb_sub(t2, t4, t3, prec);
            arb_mul(t1, t1, t2, prec);

            /* Add to accumulator */
            arb_add(acc, acc, t1, acc_prec);

            if (acc_prec > 2 * prec + 32)
            {
                arb_add(x, x, acc, res_prec);
                arb_zero(acc);
                acc_prec = prec + 32;
            }
        }
    }

    arb_add(x, x, acc, res_prec);

    arb_clear(acc);
    arb_clear(t1);
    arb_clear(t2);
    arb_clear(t3);
    arb_clear(t4);
}

typedef struct
{
    arb_ptr x;
    const fmpz * n;
    arb_srcptr C;
    arb_srcptr exp1;
    const fmpz * n24;
    slong N0;
    slong N;
    slong step;
    slong prec;
    slong acc_prec;
    slong res_prec;

} work_t;

static void
worker(slong i, work_t * work)
{
    legacy_hrr_sum_arb_range(work->x + i, work->n, work->C, work->exp1, work->n24, work->N0 + i, work->N, work->step, work->prec, work->acc_prec, work->res_prec);
}

static void
legacy_hrr_sum_arb(arb_t x, const fmpz_t n, slong N0, slong N)
{
    arb_t C, t, exp1;
    fmpz_t n24;
    slong prec, res_prec, acc_prec, guard_bits;
    slong num_threads;
    double nd;

    if (fmpz_cmp_ui(n, 2) <= 0)
    {
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }

    nd = fmpz_get_d(n);

    /* compute initial precision */
    guard_bits = 2 * FLINT_BIT_COUNT(N) + 32;
    prec = partitions_remainder_bound_log2(nd, N0) + guard_bits;
    prec = FLINT_MAX(prec, DOUBLE_PREC);
    res_prec = acc_prec = prec;


    arb_init(C);
    arb_init(exp1);
    fmpz_init(n24);

    arb_zero(x);

    /* n24 = 24n - 1 */
    fmpz_set(n24, n);
    fmpz_mul_ui(n24, n24, 24);
    fmpz_sub_ui(n24, n24, 1);

    /* C = (pi/6) sqrt(24n-1) */
    arb_const_pi(C, prec);

    arb_init(t);
    arb_sqrt_fmpz(t, n24, prec);
    arb_mul(C, C, t, prec);
    arb_div_ui(C, C, 6, prec);
    arb_clear(t);


    /* exp1 = exp(C) */
    arb_exp(exp1, C, prec);

    num_threads = flint_get_num_threads();

    if (num_threads == 1)
    {
        legacy_hrr_sum_arb_range(x, n, C, exp1, n24, N0, N, 1, prec, acc_prec, res_prec);
    }
    else
    {
        arb_ptr s;
        slong i;
        work_t work;

        num_threads = FLINT_MIN(num_threads, 8);
        num_threads = FLINT_MAX(num_threads, 1);

        s = _arb_vec_init(num_threads);

        work.x = s;
        work.n = n;
        work.C = C;
        work.exp1 = exp1;
        work.n24 = n24;
        work.N0 = N0;
        work.N = N;
        work.step = num_threads;
        work.prec = prec;
        work.acc_prec = acc_prec;
        work.res_prec = res_prec;

        flint_parallel_do((do_func_t) worker, &work, num_threads, -1, FLINT_PARALLEL_UNIFORM);

        for (i = 0; i < num_threads; i++)
            arb_add(x, x, s + i, prec);

        _arb_vec_clear(s, num_threads);
    }

    fmpz_clear(n24);
    arb_clear(exp1);
    arb_clear(C);
}

/* p(n) for n > 2 by the legacy implementation */
static void
legacy_partitions(fmpz_t p, const fmpz_t n)
{
    arb_t x;
    arf_t bound;
    slong N;

    arb_init(x);
    arf_init(bound);

    N = legacy_needed_terms(fmpz_get_d(n));
    legacy_hrr_sum_arb(x, n, 1, N);
    partitions_rademacher_bound(bound, n, N);
    arb_add_error_arf(x, bound);

    if (!arb_get_unique_fmpz(p, x))
        flint_throw(FLINT_ERROR, "not unique!\n%s\n", arb_get_str(x, 50, 0));

    arb_clear(x);
    arf_clear(bound);
}

/* ---- the program ---------------------------------------------------------- */

int
main(int argc, char * argv[])
{
    fmpz_t n, p;
    slong i;
    int quiet = 0, legacy = 0, have_n = 0;
    slong num_threads = 1;
    int timing = 0;

    fmpz_init(n);
    fmpz_init(p);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-quiet"))
            quiet = 1;
        else if (!strcmp(argv[i], "-legacy"))
            legacy = 1;
        else if (!strcmp(argv[i], "-timing"))
            timing = 1;
        else if (!strcmp(argv[i], "-threads") && i + 1 < argc)
            num_threads = atol(argv[++i]);
        else if (!strcmp(argv[i], "-legacy"))
            legacy = 1;
        else
        {
            /* allow expression input like "10^12" */
            gr_ctx_t ZZ;
            gr_ctx_init_fmpz(ZZ);
            if (gr_set_str(n, argv[i], ZZ) != GR_SUCCESS)
            {
                flint_printf("unable to parse integer: %s\n", argv[i]);
                return 1;
            }
            gr_ctx_clear(ZZ);
            have_n = 1;
        }
    }

    if (!have_n)
    {
        flint_printf("usage: partitions n [-threads t] [-quiet] [-legacy] [-timing]\n");
        return 1;
    }

    if (fmpz_sgn(n) < 0 || fmpz_bits(n) > 2 * FLINT_BITS)
    {
        flint_printf("n out of range\n");
        return 1;
    }

    flint_set_num_threads(num_threads);

    flint_printf("p(");
    fmpz_print(n);
    flint_printf(") = \n");

    if (!timing)
    {
        TIMEIT_ONCE_START;
        if (legacy && fmpz_cmp_ui(n, 3) >= 0)
            legacy_partitions(p, n);
        else
            partitions_fmpz_fmpz(p, n, 0);
        TIMEIT_ONCE_STOP;
    }
    else
    {
        TIMEIT_START;
        if (legacy && fmpz_cmp_ui(n, 3) >= 0)
            legacy_partitions(p, n);
        else
            partitions_fmpz_fmpz(p, n, 0);
        TIMEIT_STOP;
    }

    if (!quiet)
    {
        fmpz_print(p);
        flint_printf("\n");
    }
    else
    {
        flint_printf("(%wd digits)\n", (slong) fmpz_sizeinbase(p, 10));
    }

    print_memory_usage();

    fmpz_clear(n);
    fmpz_clear(p);
    flint_cleanup_master();
    return 0;
}
