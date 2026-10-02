/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz.h"
#include "arb.h"
#include "thread_support.h"
#include "acb_dirichlet.h"
#include "acb_dirichlet/impl.h"
#include "dfloat.h"

/*
    The choice of the method, by a cost model measured with 2 threads
    (seconds; only the ratios matter) at the heights n = 10^d,
    d = 11, ..., 18 (extrapolated beyond), L = log2(n) and the target
    precision p (relative bits):

    * Riemann-Siegel (isolation and refinement of each zero): r(d) per
      zero up to about the default precision 64 + L.  Beyond, the
      refinement ends with four evaluations of Z at p + 12 bits (see
      _refine_hardy_z_zero_final), each about r/11 while the main sum
      stays in the dfloat range (p + 12 <= 194 - log2(t)/4) and
      1.25 r (p/192)^1.55 beyond (r + 4 of them: 6 r at p = 192 and 12 r
      at p = 320, as measured from 1e11 to 1e17).

    * The large height method (Platt) with the dfloat sums over j, at a
      working precision up to about L + 104, giving about L + 68 bits:
      C(d) per call (one call isolates thousands of zeros) plus about
      10 ms per zero; for a higher target, each zero is then finished
      by the last stage of the refinement (four evaluations of Z).

    * The large height method with arb sums at the working precision
      p + 48: Carb(d) per call (3-12 C), about C(d) (1 + (p - 194)/400),
      but its accuracy saturates at about S(d) bits (its parameters only
      depend on the height): S = 253, 234, 219, 201, 187, 170, 155 from
      1e11 to 1e17, so it is only used up to S - 6.

    At the default precision, the large height method is then used from
    158 zeros at 1e11 (it breaks even at about 140, measured directly),
    51 at 1e12 (45), 14 at 1e13, 10 at 1e14, 3 at 1e15, 2 from 1e16 and
    for one zero from 1e20 (extrapolated).  For 100 zeros at 1e15 to 180
    bits, the arb version takes 31 s against about 1200 s for
    Riemann-Siegel; to 192 bits (beyond its accuracy), the dfloat version
    and the refinement take about 5 s per zero instead of 6.4 s.  Without dfloat (x87 arithmetic), all three use
    arb throughout; then the large height method is used from 10^15, or
    from 10^11 for more than 100 zeros (the tuning that predates dfloat),
    with the arb sums when its accuracy suffices.  It is only used for
    10^11 <= n < 10^23.
*/

static const double _hzz_r[12] =      /* Riemann-Siegel per zero */
    { 0.031, 0.085, 0.26, 0.42, 1.85, 7.0, 15.6, 81, 330, 1500, 6400, 27000 };
static const double _hzz_C[12] =      /* dfloat Platt per call */
    { 3.3, 3.8, 3.3, 4.1, 4.9, 7.2, 18.1, 94, 380, 1400, 5600, 22000 };
static const double _hzz_Carb[12] =   /* arb Platt per call */
    { 12.7, 12.4, 13.7, 19.3, 34, 87, 204, 1100, 4500, 17000, 67000, 270000 };
static const slong _hzz_S[12] =       /* accuracy of the large height method */
    { 253, 234, 219, 201, 187, 170, 155, 140, 124, 109, 94, 79 };

/* the cost of the last stage of the refinement of one zero to p bits */
static double
_hzz_final_cost(double r, slong p, slong tbits)
{
    if (p + 12 <= 194 - tbits / 4)
        return 4 * r / 11;
    else
        return 4 * 1.25 * r * pow(p / 192.0, 1.55);
}

slong
_acb_dirichlet_hardy_z_zeros_use_platt(const fmpz_t n, slong len, slong prec)
{
    slong d, L, tbits, Pd, Pa, res = 0;
    fmpz_t t;

    if (len <= 0 || fmpz_sgn(n) <= 0)
        return 0;

    /* d = floor(log10(n)) */
    d = fmpz_sizeinbase(n, 10) - 1;
    fmpz_init(t);
    fmpz_ui_pow_ui(t, 10, d);
    if (fmpz_cmp(n, t) < 0)
        d--;
    fmpz_clear(t);

    if (d < 11 || d > 22)
        return 0;

    L = fmpz_clog_ui(n, 2);
    tbits = L - 2;      /* log2 of the height, about n / 5 to n / 7 */
    if (prec <= 0)
        prec = 64 + L;

    /* the working precisions: the dfloat sums (at least L + 56 bits, for
       the isolation), and the arb sums up to the accuracy */
    Pd = FLINT_MAX(prec + 40, L + 56);
    Pd = FLINT_MIN(Pd, L + 104);
    Pa = (prec > L + 68 && prec <= _hzz_S[d - 11] - 6) ? prec + 48 : 0;

    if (!dfloat_is_supported())
    {
        if (d >= 15 || len > 100)
            res = Pa ? Pa : Pd;
    }
    else
    {
        double r = _hzz_r[d - 11], fin = _hzz_final_cost(r, prec, tbits);
        double c_rs, c_pd, c_pa;

        /* Riemann-Siegel: the refinement in two stages for a high
           precision, else somewhat slower above the default precision */
        if (prec >= tbits + 96)
            c_rs = len * (r + fin);
        else
            c_rs = len * r * (1 + FLINT_MAX(0, prec - (L + 64)) / 64.0);

        c_pd = _hzz_C[d - 11] + len * (0.01 + ((prec > L + 68) ? fin : 0));
        c_pa = Pa ? _hzz_Carb[d - 11] * (1 + FLINT_MAX(0, Pa - 242) / 400.0) + len * 0.01 : HUGE_VAL;

        if (c_pd < c_rs && c_pd <= c_pa)
            res = Pd;
        else if (c_pa < c_rs)
            res = Pa;
    }

    return res;
}

typedef struct
{
    arb_ptr res;
    slong prec;
}
_hzz_work_t;

static void
_hzz_refine_worker(slong i, _hzz_work_t * work)
{
    _acb_dirichlet_refine_hardy_z_zero_ball(work->res + i, work->res + i, work->prec);
}

void
acb_dirichlet_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, slong prec)
{
    slong found = 0, wp;
    fmpz_t k;

    if (len <= 0)
        return;

    if (fmpz_sgn(n) < 1)
        flint_throw(FLINT_ERROR, "nonpositive indices of zeros are not supported\n");

    fmpz_init(k);

    wp = _acb_dirichlet_hardy_z_zeros_use_platt(n, len, prec);
    if (wp != 0)
    {
        _hzz_work_t work;

        found = acb_dirichlet_platt_hardy_z_zeros(res, n, len, wp);

        /* the zeros are rigorously isolated; those not accurate to the
           requested precision (beyond what the large height method
           delivers) are finished by Riemann-Siegel evaluations */
        work.res = res;
        work.prec = prec;
        flint_parallel_do((do_func_t) _hzz_refine_worker, &work, found, -1, FLINT_PARALLEL_STRIDED);
    }

    /* the zeros that the large height method did not find (or all) */
    if (found < len)
    {
        fmpz_add_si(k, n, found);
        _acb_dirichlet_hardy_z_zeros_rs(res + found, k, len - found, prec);
    }

    fmpz_clear(k);
}
