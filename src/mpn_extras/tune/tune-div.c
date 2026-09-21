/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Tuning program for the division and square root cutoffs in flint-mparam.h:

        FLINT_MPN_DIV_DC_CUTOFF             schoolbook -> divide and conquer division
        FLINT_MPN_DIVAPPR_DC_CUTOFF         same, approximate quotient
        FLINT_MPN_TDIV_QR_NEWTON_CUTOFF     divide and conquer -> Newton division
        FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF   same, quotient >= 2 x divisor
        FLINT_MPN_TDIV_Q_DC_CUTOFF          schoolbook -> divide and conquer, quotient only
        FLINT_MPN_DIVAPPROX_SHORT_CUTOFF    -> short division (quotient length)
        FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF   flint_mpn_divapprox -> Newton
        FLINT_MPN_INV_NEWTON_CUTOFF         flint_mpn_inv -> Newton
        FLINT_MPN_INV_NEWTON_LONG_CUTOFF    same, quotient 4 x the length of x
        FLINT_MPN_DIVREM_1_HW_CUTOFF        hardware division chain -> mpn_divrem_1
        FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF   same, normalized divisor
        FLINT_MPN_DIV_2_HW_CUTOFF           two-limb divisors: hardware divisions
                                            -> 3/2 inverse
        FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF    same for 3- to 7-limb divisors, by
                                            quotient length
                                            (these four only with
                                            FLINT_PREINVERT_LIMB_USE_NATIVE)
        FLINT_MPN_DIV_2_GMP_CUTOFF          two-limb divisors: FLINT -> GMP's
                                            assembly mpn_divrem_2 (if available)
        FLINT_MPN_SQRTREM_NEWTON_CUTOFF     divide and conquer -> Newton square root
        FLINT_MPN_DC_BDIV_QR_CUTOFF         schoolbook -> divide and conquer Hensel
                                            division with remainder
        FLINT_MPN_DC_BDIV_Q_CUTOFF          same, quotient only
        FLINT_MPN_DIVEXACT_NEWTON_CUTOFF    divide and conquer -> Newton Hensel
                                            (exact) division
        FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF   same, quotient >= 4 x divisor

    The division and square root sources are compiled into this program with
    the cutoffs turned into variables, so that all crossovers are measured in
    one run without rebuilding FLINT. Build and run it from the FLINT build
    directory with

        make tune
        build/mpn_extras/tune/tune-div

    (or build just this program with
    make build/mpn_extras/tune/tune-div).

    Run on an idle machine; a full run takes about 30 seconds. The suggested
    values are printed at the end, and the tables above them show the
    time ratios (new method / old method, below 1 means the new method is
    faster) from which they are derived.
*/

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "flint.h"
#include "mpn_extras.h"

#if FLINT_BITS != 64
# error "the tuned code is only used on 64-bit machines"
#endif

/* cutoffs as variables */
static slong tune_div_dc = 16, tune_divappr_dc = 80, tune_tdiv_qr_newton = 1024,
    tune_tdiv_qr_newton_long = 640,
    tune_tdiv_q_dc = 100, tune_divapprox_newton = 500,
    tune_sqrtrem_newton = 4500, tune_dc_bdiv_qr = 24, tune_dc_bdiv_q = 50,
    tune_divexact_newton = 1024, tune_divexact_unbalanced = 64,
    tune_inv_newton = 256, tune_inv_newton_long = 256, tune_divapprox_short = 58, tune_divrem_1_hw = 24, tune_divrem_1_norm_hw = 14,
    tune_div_2_hw = 22, tune_div_small_hw_qn = 4, tune_div_2_gmp = 1000000;

#undef FLINT_MPN_DIV_DC_CUTOFF
#undef FLINT_MPN_TDIV_Q_DC_CUTOFF
#undef FLINT_MPN_DIVAPPR_DC_CUTOFF
#undef FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF
#undef FLINT_MPN_INV_NEWTON_CUTOFF
#undef FLINT_MPN_DIVAPPROX_SHORT_CUTOFF
#undef FLINT_MPN_INV_NEWTON_LONG_CUTOFF
#if !FLINT_PREINVERT_LIMB_USE_NATIVE
/* not tuned: keep the values of flint-mparam.h */
static const slong FLINT_MPN_DIVREM_1_HW_CUTOFF_DEFAULT = FLINT_MPN_DIVREM_1_HW_CUTOFF;
static const slong FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF_DEFAULT = FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF;
static const slong FLINT_MPN_DIV_2_HW_CUTOFF_DEFAULT = FLINT_MPN_DIV_2_HW_CUTOFF;
static const slong FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF_DEFAULT = FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF;
#endif
#undef FLINT_MPN_DIVREM_1_HW_CUTOFF
#undef FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF
#undef FLINT_MPN_DIV_2_HW_CUTOFF
#undef FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF
#undef FLINT_MPN_DIV_2_GMP_CUTOFF
#undef FLINT_MPN_SQRTREM_NEWTON_CUTOFF
#undef FLINT_MPN_DC_BDIV_QR_CUTOFF
#undef FLINT_MPN_DC_BDIV_Q_CUTOFF
#undef FLINT_MPN_DIVEXACT_NEWTON_CUTOFF
#undef FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF
#define FLINT_MPN_DIV_DC_CUTOFF tune_div_dc
#define FLINT_MPN_TDIV_Q_DC_CUTOFF tune_tdiv_q_dc
#define FLINT_MPN_DIVAPPR_DC_CUTOFF tune_divappr_dc
#define FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF tune_divapprox_newton
#define FLINT_MPN_INV_NEWTON_CUTOFF tune_inv_newton
#define FLINT_MPN_DIVAPPROX_SHORT_CUTOFF tune_divapprox_short
#define FLINT_MPN_INV_NEWTON_LONG_CUTOFF tune_inv_newton_long
#define FLINT_MPN_DIVREM_1_HW_CUTOFF tune_divrem_1_hw
#define FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF tune_divrem_1_norm_hw
#define FLINT_MPN_DIV_2_HW_CUTOFF tune_div_2_hw
#define FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF tune_div_small_hw_qn
#define FLINT_MPN_DIV_2_GMP_CUTOFF tune_div_2_gmp
#define FLINT_MPN_SQRTREM_NEWTON_CUTOFF tune_sqrtrem_newton
#define FLINT_MPN_DC_BDIV_QR_CUTOFF tune_dc_bdiv_qr
#define FLINT_MPN_DC_BDIV_Q_CUTOFF tune_dc_bdiv_q
#define FLINT_MPN_DIVEXACT_NEWTON_CUTOFF tune_divexact_newton
#define FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF tune_divexact_unbalanced

/* private copies of the exported functions */
#define _flint_mpn_tdiv_qr tune_mpn_tdiv_qr
#define _flint_mpn_tdiv_qr_small tune_mpn_tdiv_qr_small
#define _flint_mpn_divrem_basecase_preinv1 tune_mpn_divrem_basecase_preinv1
#define _flint_mpn_divapprox_basecase_preinv1 tune_mpn_divapprox_basecase_preinv1
#define _flint_mpn_div_basecase_preinv1 tune_mpn_div_basecase_preinv1
#define _flint_mpn_divrem_n_divconquer_preinv1 tune_mpn_divrem_n_divconquer_preinv1
#define _flint_mpn_divapprox_n_divconquer_preinv1 tune_mpn_divapprox_n_divconquer_preinv1
#define _flint_mpn_divrem_preinv1 tune_mpn_divrem_preinv1
#define _flint_mpn_divapprox_preinv1 tune_mpn_divapprox_preinv1
#define _flint_mpn_tdiv_qr_divconquer tune_mpn_tdiv_qr_divconquer
#define _flint_mpn_tdiv_q_divconquer tune_mpn_tdiv_q_divconquer
#define flint_mpn_divapprox tune_mpn_divapprox
#define flint_mpn_inv tune_mpn_inv
#define _flint_mpn_inv_basecase tune_mpn_inv_basecase
#define _flint_mpn_sqrtrem_newton tune_mpn_sqrtrem_newton
#define _flint_mpn_sqrtrem_gmp tune_mpn_sqrtrem_gmp
#define _flint_mpn_sqrtrem_divconquer tune_mpn_sqrtrem_divconquer
#define _flint_mpn_sqrtrem tune_mpn_sqrtrem
#define _flint_mpn_divexact_hensel tune_mpn_divexact_hensel
#define _flint_mpn_divexact tune_mpn_divexact
#define _flint_mpn_divisible_bdiv tune_mpn_divisible_bdiv

/* prototypes for the renamed functions */
void tune_mpn_tdiv_qr(mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
void tune_mpn_tdiv_qr_small(mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
mp_limb_t tune_mpn_divrem_basecase_preinv1(mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_limb_t);
mp_limb_t tune_mpn_divapprox_basecase_preinv1(mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_limb_t);
mp_limb_t tune_mpn_div_basecase_preinv1(mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_limb_t);
mp_limb_t tune_mpn_divrem_n_divconquer_preinv1(mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_limb_t, mp_ptr);
mp_limb_t tune_mpn_divapprox_n_divconquer_preinv1(mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_limb_t, mp_ptr);
mp_limb_t tune_mpn_divrem_preinv1(mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_limb_t, mp_ptr);
mp_limb_t tune_mpn_divapprox_preinv1(mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_limb_t, mp_ptr);
void tune_mpn_tdiv_qr_divconquer(mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
void tune_mpn_tdiv_q_divconquer(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
void tune_mpn_divapprox(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
void tune_mpn_inv(mp_ptr, mp_srcptr, mp_size_t, mp_size_t);
void tune_mpn_inv_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_size_t);
void tune_mpn_sqrtrem_newton(mp_ptr, mp_ptr, mp_srcptr, mp_size_t);
void tune_mpn_sqrtrem_gmp(mp_ptr, mp_ptr, mp_srcptr, mp_size_t);
mp_size_t tune_mpn_sqrtrem_divconquer(mp_ptr, mp_ptr, mp_srcptr, mp_size_t);
mp_size_t tune_mpn_sqrtrem(mp_ptr, mp_ptr, mp_srcptr, mp_size_t);
void tune_mpn_divexact_hensel(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
void tune_mpn_divexact(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
int tune_mpn_divisible_bdiv(mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, unsigned int);

#include "../div_qr.c"
#include "../sqrtrem.c"
#include "../divexact.c"
#include "../inv.c"

/* timing: minimum over 5 runs of the time per call */
static double
now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

/* time per call of expr1 and expr2, measured alternately (the minimum over
   7 rounds each), to reduce the effect of frequency changes and noise */
#define TIME2(expr1, res1, expr2, res2) \
    do { \
        double __b1 = 1e30, __b2 = 1e30, __t, __t0; \
        long __reps, __i; \
        int __r; \
        for (__reps = 1; ; __reps *= 2) \
        { \
            __t0 = now(); \
            for (__i = 0; __i < __reps; __i++) { expr1; } \
            __t = now() - __t0; \
            if (__t > 0.003) break; \
        } \
        for (__r = 0; __r < 7; __r++) \
        { \
            __t0 = now(); \
            for (__i = 0; __i < __reps; __i++) { expr1; } \
            __t = (now() - __t0) / __reps; \
            if (__t < __b1) __b1 = __t; \
            __t0 = now(); \
            for (__i = 0; __i < __reps; __i++) { expr2; } \
            __t = (now() - __t0) / __reps; \
            if (__t < __b2) __b2 = __t; \
        } \
        (res1) = __b1; \
        (res2) = __b2; \
    } while (0)

static flint_rand_t state;

/* random normalized divisor and dividend for a 2n / n division */
static void
balanced_operands(mp_ptr a, mp_ptr b, mp_size_t n)
{
    flint_mpn_urandomb(a, state, 2 * n * FLINT_BITS);
    flint_mpn_urandomb(b, state, n * FLINT_BITS);
    b[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    a[2 * n - 1] = b[n - 1] - 1;   /* quotient below B^n */
}

/* random odd divisor and dividend for a 2n / n Hensel division */
static void
odd_operands(mp_ptr a, mp_ptr b, mp_size_t n)
{
    flint_mpn_urandomb(a, state, 2 * n * FLINT_BITS);
    flint_mpn_urandomb(b, state, n * FLINT_BITS);
    b[0] |= 1;
    b[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    a[2 * n - 1] |= UWORD(1) << (FLINT_BITS - 1);
}

/*
    Crossover from a table of ratios (new / old): the first size from which
    the new method is faster and stays within 2% of the old one at all larger
    sizes, ignoring isolated outliers (a single size whose neighbours are
    both faster). Returns -1 if there is none.
*/
static slong
crossover(const slong * sizes, const double * ratio, int num)
{
    int i, j, ok;

    for (i = 0; i < num; i++)
    {
        if (ratio[i] >= 1.0)
            continue;
        ok = 1;
        for (j = i + 1; j < num; j++)
        {
            if (ratio[j] > 1.02
                && !(j + 1 < num && ratio[j - 1] < 1.0 && ratio[j + 1] < 1.0))
                ok = 0;
        }
        if (ok)
            return sizes[i];
    }
    return -1;
}

static void
print_row(slong size, double t_old, double t_new)
{
    flint_printf("  %6wd   %10.3e %10.3e   %6.3f\n", size, t_old, t_new, t_new / t_old);
}

#define MAXS 64

int main(int argc, char ** argv)
{
    slong sizes[MAXS], result;
    double ratio[MAXS], t1, t2, r1, r2;
    int num, k;
    mp_ptr a, b, q, r, w, t, s;
    mp_size_t n;

    (void) argc;
    (void) argv;

    flint_rand_init(state);

    a = flint_malloc(40000 * sizeof(mp_limb_t));
    b = flint_malloc(20000 * sizeof(mp_limb_t));
    q = flint_malloc(20000 * sizeof(mp_limb_t));
    r = flint_malloc(20000 * sizeof(mp_limb_t));
    w = flint_malloc(40000 * sizeof(mp_limb_t));
    t = flint_malloc(20000 * sizeof(mp_limb_t));
    s = flint_malloc(20000 * sizeof(mp_limb_t));

    /* 1. schoolbook vs one level of divide and conquer, exact division */
    flint_printf("FLINT_MPN_DIV_DC_CUTOFF: 2n / n exact division, schoolbook vs divide and conquer\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "school", "dc", "ratio");
    num = 0;
    for (n = 6; n <= 80; n += (n < 24) ? 2 : 4)
    {
        mp_limb_t dinv;
        balanced_operands(a, b, n);
        dinv = flint_mpn_preinv1(b[n - 1], b[n - 2]);
        TIME2((tune_div_dc = WORD_MAX, flint_mpn_copyi(w, a, 2 * n), tune_mpn_divrem_preinv1(q, w, 2 * n, b, n, dinv, t)), t1,
              (tune_div_dc = n, flint_mpn_copyi(w, a, 2 * n), tune_mpn_divrem_preinv1(q, w, 2 * n, b, n, dinv, t)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_div_dc = (result == -1) ? 80 : result;
    flint_printf("  -> %wd\n\n", tune_div_dc);

    /* 2. the same for the approximate quotient */
    flint_printf("FLINT_MPN_DIVAPPR_DC_CUTOFF: 2n / n approximate division, schoolbook vs divide and conquer\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "school", "dc", "ratio");
    num = 0;
    /* below FLINT_MPN_DIV_DC_CUTOFF both variants are the schoolbook
       division */
    for (n = FLINT_MAX(8, tune_div_dc); n <= 240; n += (n < 64) ? 4 : 16)
    {
        mp_limb_t dinv;
        balanced_operands(a, b, n);
        dinv = flint_mpn_preinv1(b[n - 1], b[n - 2]);
        TIME2((tune_divappr_dc = WORD_MAX, flint_mpn_copyi(w, a, 2 * n), tune_mpn_divapprox_preinv1(q, w, 2 * n, b, n, dinv, t)), t1,
              (tune_divappr_dc = n, flint_mpn_copyi(w, a, 2 * n), tune_mpn_divapprox_preinv1(q, w, 2 * n, b, n, dinv, t)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_divappr_dc = (result == -1) ? 240 : result;
    flint_printf("  -> %wd\n\n", tune_divappr_dc);

    /* 3. divide and conquer vs Newton, exact division, for qn = n
       (FLINT_MPN_TDIV_QR_NEWTON_CUTOFF) and qn = 2n
       (FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF, where Newton wins earlier) */
    {
        slong res[2];
        int shape;

        for (shape = 1; shape <= 2; shape++)
        {
            mp_size_t an;

            flint_printf("%s: %dn / n division, divide and conquer vs Newton\n",
                (shape == 1) ? "FLINT_MPN_TDIV_QR_NEWTON_CUTOFF" : "FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF", shape + 1);
            flint_printf("  %6s   %10s %10s   %6s\n", "n", "dc", "newton", "ratio");
            num = 0;
            for (n = 250; n <= 4000; n = n * 5 / 4)
            {
                balanced_operands(a, b, n);
                b[n - 1] >>= 3;   /* typical unnormalized divisor */
                an = 2 * n;
                if (shape == 2)
                {
                    an = 3 * n - 1;
                    flint_mpn_urandomb(a, state, an * FLINT_BITS);
                }
                TIME2(tune_mpn_tdiv_qr_divconquer(q, r, a, an, b, n), t1,
                      _flint_mpn_tdiv_qr_newton(q, r, a, an, b, n), t2);
                print_row(n, t1, t2);
                sizes[num] = n; ratio[num++] = t2 / t1;
            }
            result = crossover(sizes, ratio, num);
            res[shape - 1] = (result == -1) ? 4000 : result;
            flint_printf("  -> %wd\n\n", res[shape - 1]);
        }
        tune_tdiv_qr_newton = res[0];
        tune_tdiv_qr_newton_long = FLINT_MIN(res[0], res[1]);
    }

    /* 4. quotient only: schoolbook vs divide and conquer */
    flint_printf("FLINT_MPN_TDIV_Q_DC_CUTOFF: 2n / n quotient only, schoolbook vs divide and conquer\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "school", "dc", "ratio");
    num = 0;
    tune_divapprox_short = WORD_MAX;   /* tuned next */
    for (n = 16; n <= 600; n = n * 5 / 4)
    {
        balanced_operands(a, b, n);
        b[n - 1] >>= 3;
        TIME2((tune_tdiv_q_dc = WORD_MAX, tune_mpn_tdiv_q_divconquer(q, a, 2 * n, b, n)), t1,
              (tune_tdiv_q_dc = n, tune_mpn_tdiv_q_divconquer(q, a, 2 * n, b, n)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_tdiv_q_dc = (result == -1) ? 600 : result;
    flint_printf("  -> %wd\n\n", tune_tdiv_q_dc);

    /* 4b. short division vs the above, for the (2n + 1) / n quotient of a
       correctly rounded n-limb division (qn = n + 2), for normalized
       divisors (floating-point mantissas, where the previous algorithms
       need no shifts) and unnormalized ones; the larger crossover */
    {
        slong res[2];
        int norm;

        for (norm = 1; norm >= 0; norm--)
        {
            flint_printf("FLINT_MPN_DIVAPPROX_SHORT_CUTOFF: (2n + 1) / n quotient, previous algorithms vs short division (%s divisor)\n",
                norm ? "normalized" : "unnormalized");
            flint_printf("  %6s   %10s %10s   %6s\n", "qn", "previous", "short", "ratio");
            num = 0;
            for (n = 12; n <= 400; n = n * 9 / 8 + 1)
            {
                balanced_operands(a, b, n + 1);
                b[n - 1] = norm ? (b[n - 1] | (UWORD(1) << (FLINT_BITS - 1))) : (b[n - 1] >> 3);
                TIME2((tune_divapprox_short = WORD_MAX, tune_mpn_tdiv_q_divconquer(q, a, 2 * n + 1, b, n)), t1,
                      (tune_divapprox_short = n + 2, tune_mpn_tdiv_q_divconquer(q, a, 2 * n + 1, b, n)), t2);
                print_row(n + 2, t1, t2);
                sizes[num] = n + 2; ratio[num++] = t2 / t1;
            }
            result = crossover(sizes, ratio, num);
            res[norm] = (result == -1) ? 402 : result;
            flint_printf("  -> %wd\n\n", res[norm]);
        }
        tune_divapprox_short = FLINT_MAX(res[0], res[1]);
    }

    /* 5. flint_mpn_divapprox: divide and conquer vs Newton */
    flint_printf("FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF: 2n / n approximate quotient, divide and conquer vs Newton\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "dc", "newton", "ratio");
    num = 0;
    for (n = 100; n <= 3000; n = n * 5 / 4)
    {
        balanced_operands(a, b, n);
        TIME2((tune_divapprox_newton = WORD_MAX, tune_mpn_divapprox(q, a, 2 * n, b, n)), t1,
              (tune_divapprox_newton = 1, tune_mpn_divapprox(q, a, 2 * n, b, n)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_divapprox_newton = (result == -1) ? 3000 : result;
    flint_printf("  -> %wd\n\n", tune_divapprox_newton);

    /* 6b. flint_mpn_inv: division vs Newton, for floor(B^(2n) / x) with
       an n-limb x */
    flint_printf("FLINT_MPN_INV_NEWTON_CUTOFF: floor(B^2n / x), n-limb x, division vs Newton\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "div", "newton", "ratio");
    num = 0;
    for (n = 24; n <= 1000; n = n * 5 / 4)
    {
        flint_mpn_urandomb(b, state, n * FLINT_BITS);
        b[n - 1] |= 1;
        TIME2((tune_inv_newton = WORD_MAX, tune_mpn_inv(q, b, n, 2 * n)), t1,
              (tune_inv_newton = 1, tune_mpn_inv(q, b, n, 2 * n)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_inv_newton = (result == -1) ? 1000 : result;
    flint_printf("  -> %wd\n\n", tune_inv_newton);

    flint_printf("FLINT_MPN_INV_NEWTON_LONG_CUTOFF: floor(B^5n / x), n-limb x, division vs Newton\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "div", "newton", "ratio");
    num = 0;
    for (n = 24; n <= 600; n = n * 5 / 4)
    {
        flint_mpn_urandomb(b, state, n * FLINT_BITS);
        b[n - 1] |= 1;
        TIME2((tune_inv_newton_long = WORD_MAX, tune_mpn_inv(q, b, n, 5 * n)), t1,
              (tune_inv_newton_long = 1, tune_mpn_inv(q, b, n, 5 * n)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_inv_newton_long = (result == -1) ? 600 : result;
    flint_printf("  -> %wd\n\n", tune_inv_newton_long);

    /* 6c. one-limb divisors: chain of hardware divisions vs GMP's
       mpn_divrem_1, for unnormalized and normalized divisors */
#if FLINT_PREINVERT_LIMB_USE_NATIVE
    {
        slong res[2];
        int norm;

        for (norm = 0; norm <= 1; norm++)
        {
            flint_printf("%s: an-limb / 1-limb division, hardware chain vs mpn_divrem_1 (%s divisor)\n",
                norm ? "FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF" : "FLINT_MPN_DIVREM_1_HW_CUTOFF",
                norm ? "normalized" : "unnormalized");
            flint_printf("  %6s   %10s %10s   %6s\n", "an", "hw", "divrem_1", "ratio");
            num = 0;
            for (n = 2; n <= 64; n += (n < 32) ? 2 : 8)
            {
                mp_limb_t d[16];
                double u1, u2;
                int k;

                /* several divisors, as the hardware latency depends on
                   the operands */
                for (k = 0; k < 16; k++)
                {
                    d[k] = n_randlimb(state);
                    if (norm)
                        d[k] |= UWORD(1) << (FLINT_BITS - 1);
                    else
                        d[k] = (d[k] >> (1 + n_randint(state, FLINT_BITS - 2))) | 1;
                }
                flint_mpn_urandomb(a, state, n * FLINT_BITS);
                t1 = t2 = 0;
                for (k = 0; k < 16; k++)
                {
                    TIME2(_div_small_1_hw(q, a, n, d[k]), u1,
                          mpn_divrem_1(q, 0, a, n, d[k]), u2);
                    t1 += u1;
                    t2 += u2;
                }
                print_row(n, t1, t2);
                sizes[num] = n; ratio[num++] = t2 / t1;
            }
            result = crossover(sizes, ratio, num);
            res[norm] = (result == -1) ? 64 : result;
            flint_printf("  -> %wd\n\n", res[norm]);
        }
        tune_divrem_1_hw = res[0];
        tune_divrem_1_norm_hw = res[1];
    }

    /* 6d. two-limb divisors: hardware divisions vs the 3/2 inverse */
    flint_printf("FLINT_MPN_DIV_2_HW_CUTOFF: an-limb / 2-limb division, hardware divisions vs 3/2 inverse\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "an", "hw", "inverse", "ratio");
    num = 0;
    for (n = 3; n <= 64; n += (n < 24) ? 1 : 8)
    {
        mp_limb_t d[16][2];
        double u1, u2;
        int k;

        for (k = 0; k < 16; k++)
        {
            d[k][0] = n_randlimb(state);
            d[k][1] = n_randlimb(state) >> n_randint(state, FLINT_BITS);
            if (d[k][1] == 0)
                d[k][1] = 1;
        }
        flint_mpn_urandomb(a, state, n * FLINT_BITS);
        t1 = t2 = 0;
        for (k = 0; k < 16; k++)
        {
            TIME2((tune_div_2_hw = WORD_MAX, tdiv_qr_small(q, r, a, n, d[k], 2)), u1,
                  (tune_div_2_hw = 0, tdiv_qr_small(q, r, a, n, d[k], 2)), u2);
            t1 += u1;
            t2 += u2;
        }
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_div_2_hw = (result == -1) ? 64 : result;
    flint_printf("  -> %wd\n\n", tune_div_2_hw);

    /* 6e. 3- to 7-limb divisors: the same by quotient length (total time
       over the divisor lengths) */
    flint_printf("FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF: 3..7-limb divisors, qn-limb quotient, hardware divisions vs 3/2 inverse\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "qn", "hw", "inverse", "ratio");
    num = 0;
    for (n = 1; n <= 16; n++)
    {
        mp_limb_t d[16][7];
        double u1, u2;
        int k, bn;

        t1 = t2 = 0;
        for (bn = 3; bn <= 7; bn++)
        {
            for (k = 0; k < 16; k++)
            {
                flint_mpn_urandomb(d[k], state, bn * FLINT_BITS);
                d[k][bn - 1] >>= n_randint(state, FLINT_BITS);
                if (d[k][bn - 1] == 0)
                    d[k][bn - 1] = 1;
            }
            flint_mpn_urandomb(a, state, (n + bn - 1) * FLINT_BITS);
            for (k = 0; k < 16; k++)
            {
                TIME2((tune_div_small_hw_qn = WORD_MAX, tdiv_qr_small(q, r, a, n + bn - 1, d[k], bn)), u1,
                      (tune_div_small_hw_qn = 0, tdiv_qr_small(q, r, a, n + bn - 1, d[k], bn)), u2);
                t1 += u1;
                t2 += u2;
            }
        }
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_div_small_hw_qn = (result == -1) ? 17 : result;
    flint_printf("  -> %wd\n\n", tune_div_small_hw_qn);
#else
    flint_printf("FLINT_MPN_DIVREM_1_HW_CUTOFF, FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF, "
        "FLINT_MPN_DIV_2_HW_CUTOFF, FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF: unused (FLINT_PREINVERT_LIMB_USE_NATIVE is 0)\n\n");
    tune_divrem_1_hw = FLINT_MPN_DIVREM_1_HW_CUTOFF_DEFAULT;
    tune_divrem_1_norm_hw = FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF_DEFAULT;
    tune_div_2_hw = FLINT_MPN_DIV_2_HW_CUTOFF_DEFAULT;
    tune_div_small_hw_qn = FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF_DEFAULT;
#endif

    /* 6f. two-limb divisors: FLINT (the hardware chain below the cutoff
       just tuned, then the 3/2 inverse) vs GMP's mpn_divrem_2 */
#if FLINT_HAVE_NATIVE_mpn_divrem_2
    flint_printf("FLINT_MPN_DIV_2_GMP_CUTOFF: an-limb / 2-limb division, FLINT vs GMP's mpn_divrem_2\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "an", "flint", "gmp", "ratio");
    num = 0;
    {
        slong n0 = 3;
#if FLINT_PREINVERT_LIMB_USE_NATIVE
        n0 = FLINT_MAX(n0, tune_div_2_hw);
#endif
        for (n = n0; n <= 256; n += (n < n0 + 16) ? 2 : (n < 64 ? 8 : 32))
        {
            mp_limb_t d[16][2];
            double u1, u2;
            int k;

            for (k = 0; k < 16; k++)
            {
                d[k][0] = n_randlimb(state);
                d[k][1] = n_randlimb(state) >> n_randint(state, FLINT_BITS);
                if (d[k][1] == 0)
                    d[k][1] = 1;
            }
            flint_mpn_urandomb(a, state, n * FLINT_BITS);
            t1 = t2 = 0;
            for (k = 0; k < 16; k++)
            {
                TIME2((tune_div_2_gmp = WORD_MAX, tdiv_qr_small(q, r, a, n, d[k], 2)), u1,
                      (tune_div_2_gmp = 0, tdiv_qr_small(q, r, a, n, d[k], 2)), u2);
                t1 += u1;
                t2 += u2;
            }
            print_row(n, t1, t2);
            sizes[num] = n; ratio[num++] = t2 / t1;
        }
    }
    result = crossover(sizes, ratio, num);
    tune_div_2_gmp = (result == -1) ? 1000000 : result;
    flint_printf("  -> %wd\n\n", tune_div_2_gmp);
#else
    flint_printf("FLINT_MPN_DIV_2_GMP_CUTOFF: unused (no assembly mpn_divrem_2)\n\n");
#endif

    /* 7. square root: divide and conquer vs Newton (geometric mean of the
       ratios with and without remainder) */
    flint_printf("FLINT_MPN_SQRTREM_NEWTON_CUTOFF: an-limb square root, divide and conquer vs Newton\n");
    flint_printf("  %6s   %10s %10s   %6s   %10s %10s   %6s\n", "an", "dc_qr", "newton_qr", "ratio", "dc_q", "newton_q", "ratio");
    num = 0;
    for (n = 1000; n <= 12000; n = n * 6 / 5)
    {
        flint_mpn_urandomb(a, state, n * FLINT_BITS);
        a[n - 1] |= UWORD(1) << (FLINT_BITS - 4);
        TIME2(tune_mpn_sqrtrem_divconquer(s, r, a, n), t1,
              _flint_mpn_sqrtrem_newton(s, r, a, n), t2);
        r1 = t2 / t1;
        flint_printf("  %6wd   %10.3e %10.3e   %6.3f", n, t1, t2, r1);
        TIME2(tune_mpn_sqrtrem_divconquer(s, NULL, a, n), t1,
              _flint_mpn_sqrtrem_newton(s, NULL, a, n), t2);
        r2 = t2 / t1;
        flint_printf("   %10.3e %10.3e   %6.3f\n", t1, t2, r2);
        sizes[num] = n; ratio[num++] = sqrt(r1 * r2);
    }
    result = crossover(sizes, ratio, num);
    tune_sqrtrem_newton = (result == -1) ? 12000 : result;
    flint_printf("  -> %wd\n\n", tune_sqrtrem_newton);

    /* 8. Hensel division with remainder (as used by the divisibility
       test): schoolbook vs divide and conquer */
    flint_printf("FLINT_MPN_DC_BDIV_QR_CUTOFF: 2n / n Hensel remainder, schoolbook vs divide and conquer\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "school", "dc", "ratio");
    num = 0;
    for (n = 8; n <= 120; n += (n < 40) ? 4 : 8)
    {
        odd_operands(a, b, n);
        TIME2((tune_dc_bdiv_qr = WORD_MAX, tune_mpn_divisible_bdiv(a, 2 * n, b, n, 0)), t1,
              (tune_dc_bdiv_qr = n, tune_mpn_divisible_bdiv(a, 2 * n, b, n, 0)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_dc_bdiv_qr = (result == -1) ? 120 : result;
    flint_printf("  -> %wd\n\n", tune_dc_bdiv_qr);

    /* 9. Hensel division, quotient only: schoolbook vs divide and conquer */
    flint_printf("FLINT_MPN_DC_BDIV_Q_CUTOFF: 2n / n exact division, schoolbook vs divide and conquer Hensel\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "school", "dc", "ratio");
    num = 0;
    for (n = 12; n <= 200; n += (n < 60) ? 4 : 10)
    {
        odd_operands(a, b, n);
        TIME2((tune_dc_bdiv_q = WORD_MAX, bdiv_q_preinv1(q, a, n + 1, b, n, 0)), t1,
              (tune_dc_bdiv_q = n, bdiv_q_preinv1(q, a, n + 1, b, n, 0)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_dc_bdiv_q = (result == -1) ? 200 : result;
    flint_printf("  -> %wd\n\n", tune_dc_bdiv_q);

    /* 10. exact division: divide and conquer vs Newton Hensel division, for
       qn = n (FLINT_MPN_DIVEXACT_NEWTON_CUTOFF) and qn = 4n
       (FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF) */
    {
        slong res[2];
        int shape;

        for (shape = 1; shape <= 2; shape++)
        {
            mp_size_t qn;

            flint_printf("%s: (%dn + n) / n exact division, divide and conquer vs Newton Hensel\n",
                (shape == 1) ? "FLINT_MPN_DIVEXACT_NEWTON_CUTOFF" : "FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF",
                (shape == 1) ? 1 : 4);
            flint_printf("  %6s   %10s %10s   %6s\n", "n", "dc", "newton", "ratio");
            num = 0;
            for (n = (shape == 1) ? 200 : 16; n <= ((shape == 1) ? 3000 : 1000); n = n * 5 / 4)
            {
                qn = (shape == 1) ? n : 4 * n;
                flint_mpn_urandomb(a, state, (qn + n - 1) * FLINT_BITS);
                flint_mpn_urandomb(b, state, n * FLINT_BITS);
                b[0] |= 1;
                a[qn + n - 2] |= UWORD(1) << (FLINT_BITS - 1);
                TIME2(bdiv_q_preinv1(q, a, qn, b, n, 0), t1,
                      tune_mpn_divexact_hensel(q, a, qn + n - 1, b, n), t2);
                print_row(n, t1, t2);
                sizes[num] = n; ratio[num++] = t2 / t1;
            }
            result = crossover(sizes, ratio, num);
            res[shape - 1] = (result == -1) ? ((shape == 1) ? 3000 : 1000) : result;
            flint_printf("  -> %wd\n\n", res[shape - 1]);
        }
        tune_divexact_newton = res[0];
        tune_divexact_unbalanced = FLINT_MIN(res[0], res[1]);
    }

    /* sanity check of the tuned configuration against GMP */
    for (k = 0; k < 2000; k++)
    {
        mp_size_t bn = 2 + n_randint(state, 1500), an = bn + n_randint(state, 1500);
        flint_mpn_urandomb(a, state, an * FLINT_BITS);
        flint_mpn_urandomb(b, state, bn * FLINT_BITS);
        if (b[bn - 1] == 0)
            b[bn - 1] = 1;
        mpn_tdiv_qr(q, r, 0, a, an, b, bn);
        tune_mpn_tdiv_qr_divconquer(t, s, a, an, b, bn);
        if (mpn_cmp(q, t, an - bn + 1) != 0 || mpn_cmp(r, s, bn) != 0)
        {
            flint_printf("sanity check failed: an = %wd, bn = %wd\n", an, bn);
            return 1;
        }
        /* floor(a/b) <= t <= floor(a/b) + 1 */
        tune_mpn_divapprox(t, a, an, b, bn);
        flint_mpn_copyi(w, q, an - bn + 1);
        mpn_add_1(w, w, an - bn + 1, 1);
        if (mpn_cmp(t, q, an - bn + 1) < 0 || mpn_cmp(t, w, an - bn + 1) > 0)
        {
            flint_printf("sanity check failed (divapprox): an = %wd, bn = %wd\n", an, bn);
            return 1;
        }
        /* exact division and divisibility of a - r = q b */
        mpn_sub(w, a, an, r, bn);
        tune_mpn_divexact(t, w, an, b, bn);
        if (mpn_cmp(q, t, an - bn + 1) != 0 || !flint_mpn_divisible(w, an, b, bn))
        {
            flint_printf("sanity check failed (divexact): an = %wd, bn = %wd\n", an, bn);
            return 1;
        }
    }

    flint_printf("Suggested values for flint-mparam.h:\n\n");
    flint_printf("#define FLINT_MPN_DIV_DC_CUTOFF %wd\n", tune_div_dc);
    flint_printf("#define FLINT_MPN_DIVAPPR_DC_CUTOFF %wd\n", tune_divappr_dc);
    flint_printf("#define FLINT_MPN_TDIV_QR_NEWTON_CUTOFF %wd\n", tune_tdiv_qr_newton);
    flint_printf("#define FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF %wd\n", tune_tdiv_qr_newton_long);
    flint_printf("#define FLINT_MPN_TDIV_Q_DC_CUTOFF %wd\n", tune_tdiv_q_dc);
    flint_printf("#define FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF %wd\n", tune_divapprox_newton);
    flint_printf("#define FLINT_MPN_INV_NEWTON_CUTOFF %wd\n", tune_inv_newton);
    flint_printf("#define FLINT_MPN_INV_NEWTON_LONG_CUTOFF %wd\n", tune_inv_newton_long);
    flint_printf("#define FLINT_MPN_DIVAPPROX_SHORT_CUTOFF %wd\n", tune_divapprox_short);
    flint_printf("#define FLINT_MPN_DIVREM_1_HW_CUTOFF %wd\n", tune_divrem_1_hw);
    flint_printf("#define FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF %wd\n", tune_divrem_1_norm_hw);
    flint_printf("#define FLINT_MPN_DIV_2_HW_CUTOFF %wd\n", tune_div_2_hw);
    flint_printf("#define FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF %wd\n", tune_div_small_hw_qn);
    flint_printf("#define FLINT_MPN_DIV_2_GMP_CUTOFF %wd\n", tune_div_2_gmp);
    flint_printf("#define FLINT_MPN_SQRTREM_NEWTON_CUTOFF %wd\n", tune_sqrtrem_newton);
    flint_printf("#define FLINT_MPN_DC_BDIV_QR_CUTOFF %wd\n", tune_dc_bdiv_qr);
    flint_printf("#define FLINT_MPN_DC_BDIV_Q_CUTOFF %wd\n", tune_dc_bdiv_q);
    flint_printf("#define FLINT_MPN_DIVEXACT_NEWTON_CUTOFF %wd\n", tune_divexact_newton);
    flint_printf("#define FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF %wd\n", tune_divexact_unbalanced);

    flint_free(a); flint_free(b); flint_free(q); flint_free(r);
    flint_free(w); flint_free(t); flint_free(s);
    flint_rand_clear(state);
    return 0;
}
