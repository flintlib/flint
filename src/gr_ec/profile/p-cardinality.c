/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Counting the points of E(F_q), over several families of curves chosen
    to have different endomorphism rings, since that is what the algorithms
    with a special case might react to:

      generic   random a4, a6; ordinary, with a large CM discriminant
      CM -3     j = 0 with p = 1 mod 3; ordinary, CM by Z[zeta_3]
      CM -4     j = 1728 with p = 1 mod 4; ordinary, CM by Z[i]
      CM -7     j = -3375; ordinary, CM by the order of discriminant -7
      ss (-3)   j = 0 with p = 2 mod 3; supersingular, #E = p + 1
      ss (-4)   j = 1728 with p = 3 mod 4; supersingular, #E = p + 1

    Three algorithms: the O(q) walk over the field, baby-step giant-step at
    O(q^(1/4)) group operations, and Schoof, polynomial in log q.

    What each needs of the base ring differs, and that decides the columns:
    naive and BSGS need square roots, to count the roots in y and to find a
    random point at all, so they are only available where gr_sqrt is; on
    mpn_mod it currently returns GR_UNABLE for a general modulus, and only
    Schoof can run. Schoof needs no square root, but it does need the short
    model and residue characteristic above 3.

    Large moduli also need gr_ctx_set_is_field: point counting asks the
    base ring whether it is a field, and neither fmpz_mod nor mpn_mod
    proves primality on its own.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_ec.h"

#define TIME_US(dest, stmt) \
    do { \
        double _tc, _tw; \
        TIMEIT_START \
        stmt; \
        TIMEIT_STOP_VALUES(_tc, _tw); \
        (dest) = _tw * 1e6; \
        (void) _tc; \
    } while (0)

typedef enum { FAM_GENERIC, FAM_CM3, FAM_CM4, FAM_CM7, FAM_SS3, FAM_SS4, FAM_NUM } family_t;

static const char * family_name[FAM_NUM] =
    { "generic", "CM -3", "CM -4", "CM -7", "ss (-3)", "ss (-4)" };

/* a prime of the given size, congruent to r mod m when m > 1 */
static void
pick_prime(fmpz_t p, flint_rand_t state, slong bits, ulong m, ulong r)
{
    slong tries;

    for (tries = 0; tries < 10000; tries++)
    {
        fmpz_randprime(p, state, bits, 0);

        if (m <= 1 || fmpz_fdiv_ui(p, m) == r)
            return;
    }
}

/* what congruence the family needs of p */
static void
family_prime(fmpz_t p, flint_rand_t state, slong bits, family_t fam)
{
    switch (fam)
    {
        case FAM_CM3: pick_prime(p, state, bits, 3, 1); break;
        case FAM_SS3: pick_prime(p, state, bits, 3, 2); break;
        case FAM_CM4: pick_prime(p, state, bits, 4, 1); break;
        case FAM_SS4: pick_prime(p, state, bits, 4, 3); break;
        default:      pick_prime(p, state, bits, 1, 0); break;
    }
}

/* y^2 = x^3 + 3c x + 2c with c = j / (1728 - j) realises the j-invariant */
static int
curve_from_j(gr_ec_ctx_t E, gr_ctx_t R, slong j)
{
    gr_ptr c, a4, a6, t;
    int status = GR_SUCCESS;

    GR_TMP_INIT4(c, a4, a6, t, R);

    status |= gr_set_si(c, j, R);
    status |= gr_set_si(t, 1728 - j, R);
    status |= gr_div(c, c, t, R);
    status |= gr_mul_ui(a4, c, 3, R);
    status |= gr_mul_two(a6, c, R);

    if (status == GR_SUCCESS)
        status = gr_ec_ctx_init_short_weierstrass(E, R, a4, a6);

    GR_TMP_CLEAR4(c, a4, a6, t, R);

    return status;
}

static int
build_curve(gr_ec_ctx_t E, gr_ctx_t R, flint_rand_t state, family_t fam)
{
    gr_ptr a, z;
    int status = GR_UNABLE;
    slong tries;

    if (fam == FAM_CM7)
        return curve_from_j(E, R, -3375);

    GR_TMP_INIT2(a, z, R);

    for (tries = 0; tries < 40 && status != GR_SUCCESS; tries++)
    {
        if (gr_randtest_not_zero(a, state, R) != GR_SUCCESS
                || gr_zero(z, R) != GR_SUCCESS)
            break;

        if (fam == FAM_GENERIC)
        {
            gr_ptr b;
            GR_TMP_INIT(b, R);
            if (gr_randtest(b, state, R) == GR_SUCCESS)
                status = gr_ec_ctx_init_short_weierstrass(E, R, a, b);
            GR_TMP_CLEAR(b, R);
        }
        else if (fam == FAM_CM3 || fam == FAM_SS3)
            status = gr_ec_ctx_init_short_weierstrass(E, R, z, a);   /* j = 0 */
        else
            status = gr_ec_ctx_init_short_weierstrass(E, R, a, z);   /* j = 1728 */
    }

    GR_TMP_CLEAR2(a, z, R);

    return status;
}

/* median of a small sample, sorted in place */
static int
cmp_double(const void * a, const void * b)
{
    double x = *(const double *) a, y = *(const double *) b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static double
median(double * v, slong n)
{
    if (n == 0)
        return 0.0;
    qsort(v, n, sizeof(double), cmp_double);
    return v[n / 2];
}

#define MAX_CURVES 9

static void
run_size(flint_rand_t state, slong bits, int use_mpn, slong ncurves,
        int do_naive, int do_bsgs, int do_schoof)
{
    family_t fam;
    fmpz_t p;

    fmpz_init(p);

    flint_printf("--- %wd-bit p, %s ---\n", bits,
            use_mpn ? "mpn_mod" : (bits < 62 ? "nmod" : "fmpz_mod"));
    flint_printf("%-9s %7s %12s %12s %12s %12s   %s\n",
            "family", "curves", "naive", "bsgs", "schoof", "cm", "agree");

    for (fam = 0; fam < FAM_NUM; fam++)
    {
        double tn[MAX_CURVES], tb[MAX_CURVES], ts[MAX_CURVES], tc[MAX_CURVES];
        slong nn = 0, nb = 0, ns = 0, nc = 0, built = 0, i;
        int agree = 1;

        for (i = 0; i < ncurves && i < MAX_CURVES; i++)
        {
            gr_ctx_t R;
            gr_ec_ctx_t E;
            fmpz_t c_naive, c_bsgs, c_schoof, c_cm;
            int have_n = 0, have_b = 0, have_s = 0, have_c = 0;

            family_prime(p, state, bits, fam);

            if (use_mpn)
            {
                if (gr_ctx_init_mpn_mod(R, p) != GR_SUCCESS)
                    continue;
            }
            else if (bits < 62)
            {
                if (!fmpz_abs_fits_ui(p)
                        || gr_ctx_init_nmod(R, fmpz_get_ui(p)) != GR_SUCCESS)
                    continue;
            }
            else
                gr_ctx_init_fmpz_mod(R, p);

            /* neither fmpz_mod nor mpn_mod proves primality by itself */
            GR_IGNORE(gr_ctx_set_is_field(R, T_TRUE));

            if (build_curve(E, R, state, fam) != GR_SUCCESS)
            {
                gr_ctx_clear(R);
                continue;
            }

            built++;

            fmpz_init(c_naive); fmpz_init(c_bsgs); fmpz_init(c_schoof); fmpz_init(c_cm);

            if (do_naive)
            {
                if (gr_ec_ctx_cardinality_naive(c_naive, E) == GR_SUCCESS)
                {
                    have_n = 1;
                    TIME_US(tn[nn], { GR_IGNORE(gr_ec_ctx_cardinality_naive(c_naive, E)); });
                    nn++;
                }
            }

            if (do_bsgs)
            {
                if (gr_ec_ctx_cardinality_bsgs(c_bsgs, E) == GR_SUCCESS)
                {
                    have_b = 1;
                    TIME_US(tb[nb], { GR_IGNORE(gr_ec_ctx_cardinality_bsgs(c_bsgs, E)); });
                    nb++;
                }
            }

            if (do_schoof)
            {
                if (gr_ec_ctx_cardinality_schoof(c_schoof, E) == GR_SUCCESS)
                {
                    have_s = 1;
                    TIME_US(ts[ns], { GR_IGNORE(gr_ec_ctx_cardinality_schoof(c_schoof, E)); });
                    ns++;
                }
            }

            if (gr_ec_ctx_cardinality_cm(c_cm, E) == GR_SUCCESS)
            {
                have_c = 1;
                TIME_US(tc[nc], { GR_IGNORE(gr_ec_ctx_cardinality_cm(c_cm, E)); });
                nc++;
            }

            /* whatever ran must have produced the same number */
            if (have_n && have_b && !fmpz_equal(c_naive, c_bsgs)) agree = 0;
            if (have_n && have_s && !fmpz_equal(c_naive, c_schoof)) agree = 0;
            if (have_b && have_s && !fmpz_equal(c_bsgs, c_schoof)) agree = 0;
            if (have_c && have_n && !fmpz_equal(c_cm, c_naive)) agree = 0;
            if (have_c && have_b && !fmpz_equal(c_cm, c_bsgs)) agree = 0;
            if (have_c && have_s && !fmpz_equal(c_cm, c_schoof)) agree = 0;

            /* a supersingular family must come out at exactly p + 1 */
            if ((fam == FAM_SS3 || fam == FAM_SS4) && (have_b || have_s || have_n || have_c))
            {
                fmpz_t want;
                fmpz_init(want);
                fmpz_add_ui(want, p, 1);
                if (have_n && !fmpz_equal(c_naive, want)) agree = 0;
                if (have_b && !fmpz_equal(c_bsgs, want)) agree = 0;
                if (have_s && !fmpz_equal(c_schoof, want)) agree = 0;
                if (have_c && !fmpz_equal(c_cm, want)) agree = 0;
                fmpz_clear(want);
            }

            fmpz_clear(c_naive); fmpz_clear(c_bsgs); fmpz_clear(c_schoof); fmpz_clear(c_cm);
            gr_ec_ctx_clear(E);
            gr_ctx_clear(R);
        }

        flint_printf("%-9s %7wd ", family_name[fam], built);

        if (nn) flint_printf("%12.1f ", median(tn, nn)); else flint_printf("%12s ", "-");
        if (nb) flint_printf("%12.1f ", median(tb, nb)); else flint_printf("%12s ", "-");
        if (ns) flint_printf("%12.1f ", median(ts, ns)); else flint_printf("%12s ", "-");
        if (nc) flint_printf("%12.1f ", median(tc, nc)); else flint_printf("%12s ", "-");

        flint_printf("  %s\n", agree ? "yes" : "NO");
    }

    flint_printf("\n");
    fmpz_clear(p);
}

int main(void)
{
    flint_rand_t state;

    flint_rand_init(state);

    flint_printf("gr_ec: counting the points of E(F_q)\n");
    flint_printf("microseconds per count, median over the curves of each family\n\n");

    /*            bits  mpn  curves  naive bsgs schoof */
    run_size(state,  20,   0,      9,     1,   1,   1);
    run_size(state,  32,   0,      9,     0,   1,   1);
    run_size(state,  44,   0,      7,     0,   1,   1);
    run_size(state,  56,   0,      5,     0,   1,   1);
    run_size(state,  64,   0,      5,     0,   1,   1);
    run_size(state,  96,   1,      3,     0,   0,   1);
    run_size(state, 128,   1,      3,     0,   0,   1);

    flint_rand_clear(state);

    return 0;
}
