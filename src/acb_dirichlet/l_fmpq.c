/*
    Copyright (C) 2021, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpq.h"
#include "ulong_extras.h"
#include "arb_mat.h"
#include "arb_hypgeom.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"
#include "thread_support.h"

/* number of threads available to the approximate functional equation */
static double
afe_num_threads(void)
{
    return FLINT_MAX(1, flint_get_num_threads());
}

/* Whether the approximate functional equation is estimated to be faster
   for L(s, chi), chi modulo q (the conductor if chi is imprimitive). */
int
_acb_dirichlet_l_fmpq_use_afe(ulong q, const fmpq_t s, slong prec)
{
    double cutoff;

    if ((slong) fmpz_bits(fmpq_numref(s)) - (slong) fmpz_bits(fmpq_denref(s)) > 20)
        return 0;

    if (fabs(fmpq_get_d(s)) > 10 + prec * 0.01)
        return 0;

    if (q == 1)
    {
        if (fmpz_is_one(fmpq_denref(s)))
            return 0;

        /* s = 1/2 needs a single sum; otherwise the cost grows somewhat
           with the denominator. Euler-Maclaurin summation slows down
           more than quadratically beyond 10^5 bits due to the memory
           needed for the Bernoulli numbers. The terms of the approximate
           functional equation are computed in parallel with almost
           perfect speedup, while Euler-Maclaurin summation gains less
           from threads, so the cutoffs (in wall time) decrease with the
           number of threads T. Single-threaded crossovers: 40000 bits
           for s = 1/2, about 100000 for s = 3/4; with T = 8, below
           20000 for both. */
        if (fmpz_is_one(fmpq_numref(s)) && fmpz_equal_si(fmpq_denref(s), 2))
            cutoff = 40000 * pow(afe_num_threads(), -0.5);
        else
            cutoff = 90000 * (1.0 + 0.15 * FLINT_MAX(0, (slong) fmpz_bits(fmpq_denref(s)) - 2))
                * pow(afe_num_threads(), -0.9);

        return prec > cutoff;
    }

    if (fmpq_is_zero(s))
        return 0;

    /* Compared with _acb_dirichlet_l_generic (mainly Hurwitz zeta
       functions with a precomputation). Single-threaded crossovers for
       odd characters: L(1/2) 15000-19000 bits for q = 3, 4, 12000 for
       q = 7, 8500 for q = 13, 3100 for q = 31, below 64 bits for
       q >= 101; L(1/4) 2600 bits for q = 3, 4 down to 530 for q = 31 and
       below 100 for q >= 101. With T threads the approximate functional
       equation is faster still (the generic algorithm gains little or
       even loses from threads at low precision): with T = 8, L(1/2) at
       about 500 bits for q = 3, 4 and below 64 bits for q >= 13, and
       L(1/4) at 300-800 bits. */
    if (fmpz_cmp_ui(fmpq_denref(s), 2) <= 0)
    {
        if (fmpq_is_one(s))
            return prec > 30000 || q > 1000;

        if (q > 50)
            return 1;

        cutoff = 30000.0 / (1.0 + q / 6.0);
        cutoff *= pow(afe_num_threads(), -1.5);
        return prec > cutoff;
    }

    cutoff = 8000.0 * pow(q, -0.85);
    cutoff *= pow(afe_num_threads(), -0.45);

    return prec > cutoff;
}

/* The same for the expansion of L(s+x) to length len >= 2. */
int
_acb_dirichlet_l_jet_use_afe(ulong q, const fmpq_t s, int deflate, slong len, slong prec)
{
    double cutoff;

    if ((slong) fmpz_bits(fmpq_numref(s)) - (slong) fmpz_bits(fmpq_denref(s)) > 20)
        return 0;

    if (fabs(fmpq_get_d(s)) > 10 + prec * 0.01)
        return 0;

    /* Compare with Euler-Maclaurin summation (Hurwitz zeta functions
       for q > 1), whose memory usage also grows quickly with the
       precision. Both algorithms cost roughly proportional to len. */
    if (q == 1)
    {
        double beta;

        /* Single-threaded crossovers (wall time): about 150000 bits for
           s = 1/2 and 440000 bits otherwise (zeta'(-1), and s = 1 for
           the Stieltjes constants) for len = 2, increasing slowly with
           len (660000 for len = 11). With T threads, the crossover
           decreases roughly like T^-beta since the approximate
           functional equation parallelizes better: with T = 8, below
           20000 bits for s = 1/2, 70000-100000 bits for len = 2 and
           220000 bits for len = 11. */
        if (fmpz_is_one(fmpq_numref(s)) && fmpz_equal_ui(fmpq_denref(s), 2))
        {
            cutoff = 150000;
            beta = 1.0;
        }
        else
        {
            cutoff = 440000;
            beta = 0.8;
        }

        cutoff *= pow(FLINT_MAX(len, 2) * 0.5, 0.24);
        beta = FLINT_MAX(0.4, beta - 0.03 * (FLINT_MAX(len, 2) - 2));
        cutoff *= pow(afe_num_threads(), -beta);

        /* Euler-Maclaurin summation needs O(prec^2) memory for the
           Bernoulli numbers (1.8 GB at 2.2e5 bits, 4.2 GB at 3.3e5
           bits), while the approximate functional equation needs
           O(prec) (and only some 20-80 MB there); near the crossover
           the times are close, so prefer the approximate functional
           equation beyond 350000 bits. */
        cutoff = FLINT_MIN(cutoff, 350000);
    }
    else
    {
        /* Compared with the sum of Hurwitz zeta function expansions.
           The approximate functional equation needs a single sum when
           s = 1/2 and is slower at integers, where the incomplete gamma
           functions are evaluated at integer parameters. Single-threaded
           crossovers: L(1/2+x) to O(x^2) 3600 bits for q = 3 down to
           640 for q = 31; L(1/4+x) to O(x^3) 15000 bits down to 1000 for
           q = 101; L(1+x) to O(x^2) 9600 bits down to 730 for q = 101;
           for q >= 331 the approximate functional equation is faster
           already at 64 bits. With T threads the crossovers decrease by
           roughly T^-beta with beta growing with q (the approximate
           functional equation has few terms for small q): with T = 8,
           by a factor about 1.5-2 for q = 3, 4 and 3-6 for q >= 13. */
        if (q > 200)
            return 1;

        if (fmpz_is_one(fmpq_numref(s)) && fmpz_equal_ui(fmpq_denref(s), 2))
            cutoff = 9000;
        else if (fmpz_is_one(fmpq_denref(s)))
            cutoff = 21000;
        else
            cutoff = 19000 * (1.0 + 0.2 * FLINT_MAX(0, (slong) fmpz_bits(fmpq_denref(s)) - 2));

        cutoff *= pow(FLINT_MAX(len, 2) * 0.5, 0.8);
        cutoff *= pow(q, -0.75);
        cutoff *= pow(afe_num_threads(), -FLINT_MIN(0.7, 0.17 * log2(q)));
    }

    (void) deflate;

    return prec > cutoff;
}

/* Sets x to s and returns 1 if s is an exact real number which is a
   rational number of small height (so that the approximate functional
   equation is efficient); otherwise returns 0. */
int
_acb_dirichlet_get_small_fmpq(fmpq_t x, const acb_t s)
{
    const arf_struct * m = arb_midref(acb_realref(s));

    if (!acb_is_exact(s) || !arb_is_zero(acb_imagref(s)) || !arf_is_finite(m))
        return 0;

    if (arf_is_zero(m))
    {
        fmpq_zero(x);
        return 1;
    }

    /* numerator of at most 40 bits and denominator of at most 20 bits */
    if (COEFF_IS_MPZ(ARF_EXP(m)) || arf_bits(m) > 40 || ARF_EXP(m) > 40)
        return 0;

    if ((slong) arf_bits(m) - ARF_EXP(m) > 20)
        return 0;

    arf_get_fmpq(x, m);
    return 1;
}

void
acb_dirichlet_l_fmpq(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    ulong cond;

    if (fmpz_is_one(fmpq_denref(s)) &&
        fmpz_bits(fmpq_numref(s)) < FLINT_BITS - 4)
    {
        acb_dirichlet_l_si(res, fmpz_get_si(fmpq_numref(s)), G, chi, prec);
        return;
    }

    cond = (G == NULL || G->q == 1) ? 1 : dirichlet_conductor_char(G, chi);

    if (_acb_dirichlet_l_fmpq_use_afe(cond, s, prec))
    {
        acb_dirichlet_l_fmpq_afe_jet(res, s, G, chi, 0, 1, prec);
        if (acb_is_finite(res))
            return;
    }

    /* not acb_dirichlet_l, which would come back here */
    acb_set_fmpq(res, s, prec + 10);
    _acb_dirichlet_l_generic(res, res, G, chi, prec);
}
