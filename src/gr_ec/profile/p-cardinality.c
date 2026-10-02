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
      CM D      j = j_D for one of the class number one discriminants and a
                prime that splits in that order, so the curve is ordinary
                with complex multiplication by it; a random quadratic twist
                is applied, so both signs of the trace occur
      CM -15*   j a root of the quadratic Hilbert class polynomial of
                discriminant -15. Class number two, hence outside the table
                gr_ec_ctx_cardinality_cm knows, and here to show that it
                declines rather than guesses
      ss D      the same j_D but a prime that is inert, which makes the
                reduction supersingular with #E = p + 1

    Four algorithms: the O(q) walk over the field, baby-step giant-step at
    O(q^(1/4)) group operations, Schoof, polynomial in log q, and the
    complex multiplication shortcut.

    What each needs of the base ring differs, and that decides the columns:
    naive and BSGS need square roots, to count the roots in y and to find a
    random point at all, so they are only available where gr_sqrt is; on
    mpn_mod it currently returns GR_UNABLE for a general modulus, and only
    Schoof can run. Schoof needs no square root, but it does need the short
    model and residue characteristic above 3.

    Only the families marked "full" below are timed with all four; for the
    rest, which exist to exercise the CM path, the general algorithms are
    still run wherever they are cheap, but only as a correctness check.

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

typedef struct
{
    const char * name;
    slong D;            /* the CM discriminant; 0 for the generic family */
    slong j;            /* j_D, when the class number is one */
    int split;          /* 1: p splits (ordinary), 0: inert (supersingular) */
    int h;              /* class number: 0 generic, 1 tabulated, 2 quadratic H_D */
    int full;           /* time every algorithm, not only the CM one */
}
family_t;

static const family_t fams[] = {
    { "generic",   0, WORD(0),                        1, 0, 1 },

    { "CM -3",    -3, WORD(0),                        1, 1, 1 },
    { "CM -4",    -4, WORD(1728),                     1, 1, 1 },
    { "CM -7",    -7, WORD(-3375),                    1, 1, 1 },
    { "CM -8",    -8, WORD(8000),                     1, 1, 0 },
    { "CM -11",  -11, WORD(-32768),                   1, 1, 0 },
    { "CM -12",  -12, WORD(54000),                    1, 1, 0 },
    { "CM -16",  -16, WORD(287496),                   1, 1, 0 },
    { "CM -19",  -19, WORD(-884736),                  1, 1, 0 },
    { "CM -27",  -27, WORD(-12288000),                1, 1, 0 },
    { "CM -28",  -28, WORD(16581375),                 1, 1, 0 },
    { "CM -43",  -43, WORD(-884736000),               1, 1, 0 },
    { "CM -67",  -67, WORD(-147197952000),            1, 1, 0 },
    { "CM -163",-163, WORD(-262537412640768000),      1, 1, 0 },

    { "CM -15*", -15, WORD(0),                        1, 2, 0 },

    { "ss -3",    -3, WORD(0),                        0, 1, 1 },
    { "ss -4",    -4, WORD(1728),                     0, 1, 1 },
    { "ss -7",    -7, WORD(-3375),                    0, 1, 0 },
    { "ss -8",    -8, WORD(8000),                     0, 1, 0 },
    { "ss -11",  -11, WORD(-32768),                   0, 1, 0 },
    { "ss -163",-163, WORD(-262537412640768000),      0, 1, 0 },
};

#define FAM_NUM ((slong) (sizeof(fams) / sizeof(family_t)))

/* H_{-15}(x) = x^2 + 191025 x - 121287375 */
#define H15_B WORD(191025)
#define H15_C WORD(-121287375)

/* B^2 - 4C, the discriminant of that quadratic */
static void
h15_disc(fmpz_t d)
{
    fmpz_set_si(d, H15_B);
    fmpz_mul(d, d, d);
    fmpz_sub_si(d, d, 4 * H15_C);
}

/*
    A prime of the given size in which D behaves as the family asks. The
    class number two family also needs p to be one for which H_D has a
    root, which is its discriminant being a square.
*/
static int
family_prime(fmpz_t p, flint_rand_t state, slong bits, const family_t * F)
{
    fmpz_t d;
    slong tries;
    int ok = 0;

    fmpz_init(d);

    for (tries = 0; tries < 10000 && !ok; tries++)
    {
        fmpz_randprime(p, state, bits, 0);

        if (fmpz_cmp_ui(p, 2000) <= 0)
            continue;

        if (F->h == 0)
        {
            ok = 1;
            continue;
        }

        fmpz_set_si(d, F->D);
        fmpz_mod(d, d, p);

        if (fmpz_jacobi(d, p) != (F->split ? 1 : -1))
            continue;

        if (F->h == 2)
        {
            h15_disc(d);
            fmpz_mod(d, d, p);

            if (fmpz_jacobi(d, p) != 1)
                continue;
        }

        ok = 1;
    }

    fmpz_clear(d);
    return ok;
}

/* the j-invariant of the family over F_p */
static int
family_j(fmpz_t j, const fmpz_t p, const family_t * F)
{
    if (F->h == 2)
    {
        fmpz_t d, s;
        int ok;

        fmpz_init(d);
        fmpz_init(s);

        h15_disc(d);
        fmpz_mod(d, d, p);

        ok = fmpz_sqrtmod(s, d, p);

        if (ok)
        {
            fmpz_sub_ui(j, s, (ulong) H15_B);
            fmpz_mod(j, j, p);

            if (fmpz_is_odd(j))
                fmpz_add(j, j, p);

            fmpz_fdiv_q_2exp(j, j, 1);
            fmpz_mod(j, j, p);
        }

        fmpz_clear(d);
        fmpz_clear(s);
        return ok;
    }

    fmpz_set_si(j, F->j);
    fmpz_mod(j, j, p);
    return 1;
}

/* the a-invariants of a random curve in the family */
static int
family_curve(fmpz_t a4, fmpz_t a6, flint_rand_t state, const fmpz_t p,
        const family_t * F)
{
    fmpz_t j, k, c;
    int ok = 1;

    if (F->h == 0)
    {
        fmpz_randm(a4, state, p);
        fmpz_randm(a6, state, p);
        return 1;
    }

    fmpz_init(j); fmpz_init(k); fmpz_init(c);

    if (!family_j(j, p, F))
        ok = 0;
    else if (fmpz_is_zero(j))                   /* j = 0, sextic twists */
    {
        fmpz_zero(a4);
        fmpz_randm(a6, state, p);
        if (fmpz_is_zero(a6)) fmpz_one(a6);
    }
    else
    {
        fmpz_set_ui(k, 1728);
        fmpz_sub(k, k, j);
        fmpz_mod(k, k, p);

        if (fmpz_is_zero(k))                    /* j = 1728, quartic twists */
        {
            fmpz_randm(a4, state, p);
            if (fmpz_is_zero(a4)) fmpz_one(a4);
            fmpz_zero(a6);
        }
        else
        {
            /* a4 = 3 j (1728 - j), a6 = 2 j (1728 - j)^2 has j-invariant j;
               the random quadratic twist makes both signs of the trace
               occur */
            fmpz_mul(a4, j, k);
            fmpz_mul_ui(a4, a4, 3);
            fmpz_mod(a4, a4, p);

            fmpz_mul(a6, j, k);
            fmpz_mul(a6, a6, k);
            fmpz_mul_ui(a6, a6, 2);
            fmpz_mod(a6, a6, p);

            do { fmpz_randm(c, state, p); } while (fmpz_is_zero(c));

            fmpz_mul(k, c, c);
            fmpz_mod(k, k, p);

            fmpz_mul(a4, a4, k);
            fmpz_mod(a4, a4, p);

            fmpz_mul(k, k, c);
            fmpz_mod(k, k, p);

            fmpz_mul(a6, a6, k);
            fmpz_mod(a6, a6, p);
        }
    }

    fmpz_clear(j); fmpz_clear(k); fmpz_clear(c);
    return ok;
}

static int
build_curve(gr_ec_ctx_t E, gr_ctx_t R, const fmpz_t a4, const fmpz_t a6)
{
    gr_ptr g4, g6;
    int status;

    GR_TMP_INIT2(g4, g6, R);

    status = gr_set_fmpz(g4, a4, R) | gr_set_fmpz(g6, a6, R);

    if (status == GR_SUCCESS)
        status = gr_ec_ctx_init_short_weierstrass(E, R, g4, g6);

    GR_TMP_CLEAR2(g4, g6, R);

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
    slong ifam;
    fmpz_t p, a4, a6;

    fmpz_init(p); fmpz_init(a4); fmpz_init(a6);

    flint_printf("--- %wd-bit p, %s ---\n", bits,
            use_mpn ? "mpn_mod" : (bits < 62 ? "nmod" : "fmpz_mod"));
    flint_printf("%-9s %7s %12s %12s %12s %12s   %s\n",
            "family", "curves", "naive", "bsgs", "schoof", "cm", "agree");

    for (ifam = 0; ifam < FAM_NUM; ifam++)
    {
        const family_t * F = fams + ifam;
        double tn[MAX_CURVES], tb[MAX_CURVES], ts[MAX_CURVES], tc[MAX_CURVES];
        slong nn = 0, nb = 0, ns = 0, nc = 0, built = 0, i;
        int agree = 1;

        /* the families that are here only for the CM path still get
           cross-checked against a general algorithm where that is cheap */
        int want_naive = do_naive && (F->full || bits <= 20);
        int want_bsgs = do_bsgs && (F->full || bits <= 44);
        int want_schoof = do_schoof && F->full;

        for (i = 0; i < ncurves && i < MAX_CURVES; i++)
        {
            gr_ctx_t R;
            gr_ec_ctx_t E;
            fmpz_t c_naive, c_bsgs, c_schoof, c_cm;
            int have_n = 0, have_b = 0, have_s = 0, have_c = 0;

            if (!family_prime(p, state, bits, F))
                continue;

            if (!family_curve(a4, a6, state, p, F))
                continue;

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

            if (build_curve(E, R, a4, a6) != GR_SUCCESS)
            {
                gr_ctx_clear(R);
                continue;
            }

            built++;

            fmpz_init(c_naive); fmpz_init(c_bsgs); fmpz_init(c_schoof); fmpz_init(c_cm);

            if (want_naive && gr_ec_ctx_cardinality_naive(c_naive, E) == GR_SUCCESS)
            {
                have_n = 1;

                if (F->full)
                {
                    TIME_US(tn[nn], { GR_IGNORE(gr_ec_ctx_cardinality_naive(c_naive, E)); });
                    nn++;
                }
            }

            if (want_bsgs && gr_ec_ctx_cardinality_bsgs(c_bsgs, E) == GR_SUCCESS)
            {
                have_b = 1;

                if (F->full)
                {
                    TIME_US(tb[nb], { GR_IGNORE(gr_ec_ctx_cardinality_bsgs(c_bsgs, E)); });
                    nb++;
                }
            }

            if (want_schoof && gr_ec_ctx_cardinality_schoof(c_schoof, E) == GR_SUCCESS)
            {
                have_s = 1;
                TIME_US(ts[ns], { GR_IGNORE(gr_ec_ctx_cardinality_schoof(c_schoof, E)); });
                ns++;
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
            if (!F->split && (have_b || have_s || have_n || have_c))
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

        flint_printf("%-9s %7wd ", F->name, built);

        if (nn) flint_printf("%12.1f ", median(tn, nn)); else flint_printf("%12s ", "-");
        if (nb) flint_printf("%12.1f ", median(tb, nb)); else flint_printf("%12s ", "-");
        if (ns) flint_printf("%12.1f ", median(ts, ns)); else flint_printf("%12s ", "-");
        if (nc) flint_printf("%12.1f ", median(tc, nc)); else flint_printf("%12s ", "-");

        flint_printf("  %s\n", agree ? "yes" : "NO");
    }

    flint_printf("\n");
    fmpz_clear(p); fmpz_clear(a4); fmpz_clear(a6);
}

int main(void)
{
    flint_rand_t state;

    flint_rand_init(state);

    flint_printf("gr_ec: counting the points of E(F_q)\n");
    flint_printf("microseconds per count, median over the curves of each family\n");
    flint_printf("a blank column means the algorithm was not run, not that it failed\n\n");

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
