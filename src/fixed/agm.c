/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "fixed.h"

/* The arithmetic-geometric mean of two positive balls.

   THE ITERATION a' = (a + b)/2, b' = sqrt(a b) converges quadratically:
   with z = (a - b)/(a + b), z' = (1 - sqrt(1 - z^2))/(1 + sqrt(1 - z^2))
   ~ z^2 / 4.  Each step is a ball sum, product and square root at the
   working precision; fball carries the radii (inexact inputs lower the
   precision of every operation automatically).

   THE FINISH.  agm(a, b) = ((a + b)/2) / 2F1(1/2, 1/2; 1; z^2) exactly,
   and the Taylor coefficients of 1/2F1(1/2, 1/2; 1; x) = 1 - sum_{j>=1}
   c_j x^j are dyadic rationals with c_j > 0 (the coefficients
   ((1/2)_k / k!)^2 of 2F1 are log-convex, so by Kaluza's theorem those
   of its reciprocal beyond the constant are nonpositive) and
   sum_{j>=1} c_j = 1 (the reciprocal vanishes at x = 1).  Hence the
   series truncated after x^(m-1) has a tail below x^m for 0 <= x <= 1,
   at every order m, and once x = z^2 < 2^(-p/m) for the working
   precision p the remaining iterations -- about log2(2m) of them --
   are replaced by m - 1 terms, the j-th needed only to the relative
   precision p - j log2(1/x): powers of x by squaring at decreasing
   precision, the coefficients as integers over the common denominator
   2^e_{m-1} (a shift), one product by (a + b)/2.  arb_agm uses m = 5.

   THE FALLBACK.  For a >= b > 0, b <= agm(a, b) <= a, so
   agm in (a + b)/2 +- |a - b|/2 whatever the iteration has reached:
   used when a - b no longer determines its size (the inputs' radii
   dominate), when a + b or a b is too wide for the division and the
   square root (relative radius above 2^-40), or when the iteration
   count runs out. */

/* c_j = agm_num[j-1] / 2^agm_exp[j-1], j = 1, ..., 15 */
#define AGM_MAX_ORDER 16
static const unsigned long long agm_num[AGM_MAX_ORDER - 1] = {
    1ULL, 5ULL, 11ULL,
    469ULL, 1379ULL, 17223ULL,
    56001ULL, 11998869ULL, 41064827ULL,
    571915951ULL, 2018982161ULL, 115338112823ULL,
    415720532641ULL, 6041874952949ULL, 22103950817043ULL };
static const unsigned char agm_exp[AGM_MAX_ORDER - 1] = {
    2, 6, 8, 14, 16, 20, 22, 30, 32, 36, 38, 44, 46, 50, 52 };

/* the largest order whose scaled coefficients fit a word */
#if FLINT_BITS == 64
#define AGM_WORD_ORDER 16
#else
#define AGM_WORD_ORDER 10
#endif

/* the default order of the finish: measured flat to within a few
   percent from m = 8 to 16, 5-9% under m = 2 and 3-5% under arb's
   m = 5 at high precision */
static int
_agm_default_order(slong n)
{
    return (n < 1024) ? 8 : 12;
}

void
_fball_agm_order(fball_t res, const fball_t x, const fball_t y, slong n,
    int m)
{
    fball_t a, b, t, d;
    slong p = n + 1, iter;

    if (x->negative || y->negative)
        flint_throw(FLINT_ERROR, "fball_agm: negative argument\n");
    if (x->size == 0 || y->size == 0)
    {
        /* agm(0, y) = 0; with a radius, the fallback below is the
           enclosure 0 +- max(|x|, |y|) */
        if (fball_is_zero_exact(x) || fball_is_zero_exact(y))
        {
            fball_zero(res);
            return;
        }
    }

    if (m == 0)
        m = _agm_default_order(n);
    m = FLINT_MAX(2, FLINT_MIN(m, AGM_WORD_ORDER));

    fball_init(a);
    fball_init(b);
    fball_init(t);
    fball_init(d);
    fball_set(a, x);
    fball_set(b, y);

    for (iter = 0; ; iter++)
    {
        slong de, tl;

        fball_add(t, a, b, p);
        fball_sub(d, a, b, p);

        if (fball_is_zero_exact(d))
        {
            fball_mul_2exp_si(t, -1);
            fball_swap(res, t);
            break;
        }

        /* |z| < 2^ze with ze = de - tl: |d| < 2^de, |t| >= 2^tl */
        de = fball_mag_2exp(d);
        tl = (t->size == 0) ? WORD_MIN / 4
            : FLINT_BITS * (t->exp - 1) + FLINT_BIT_COUNT(t->d[t->size - 1]) - 1;

        if (t->size == 0 || d->size == 0 || fball_rel_2exp(d) > -2
            || fball_rel_2exp(t) > -40 || iter > 2 * FLINT_BITS + 64)
        {
            /* the fallback: agm in t/2 +- |d|/2 */
            fball_mul_2exp_si(t, -1);
            fball_add_error_2exp(t, de - 1);
            fball_swap(res, t);
            break;
        }

        if (2 * m * (de - tl) < -FLINT_BITS * p)
        {
            /* the finish: x = z^2, 1 - sum_{j<m} c_j x^j over the
               common denominator 2^E, tail below x^m */
            fball_t z, W, T;
            fball_struct pw[AGM_MAX_ORDER];
            slong xe, cp, j, E = agm_exp[m - 2];

            fball_init(z);
            fball_init(W);
            fball_init(T);

            fball_div(z, d, t, FLINT_MAX(2, p + (de - tl) / FLINT_BITS + 2));
            fball_mul(z, z, z, FLINT_MAX(2, p + 2 * (de - tl) / FLINT_BITS + 2));
            xe = fball_mag_2exp(z);

            for (j = 1; j < m; j++)
            {
                const fball_struct * xj;

                cp = FLINT_MAX(2, p + (j * xe) / FLINT_BITS + 2);
                if (j == 1)
                    xj = z;
                else
                {
                    fball_init(pw + j);
                    if (j % 2 == 0)
                        fball_mul(pw + j, (j / 2 == 1) ? z : pw + j / 2,
                            (j / 2 == 1) ? z : pw + j / 2, cp);
                    else
                        fball_mul(pw + j, pw + j - 1, z, cp);
                    xj = pw + j;
                }
                if (j == 1)
                    fball_mul_ui(W, xj, (ulong) (agm_num[0] << (E - agm_exp[0])), cp);
                else
                    fball_addmul_ui(W, W, xj, (ulong) (agm_num[j - 1] << (E - agm_exp[j - 1])),
                        FLINT_MAX(2, p + xe / FLINT_BITS + 2));
            }
            for (j = 2; j < m; j++)
                fball_clear(pw + j);

            fball_mul_2exp_si(W, -E);
            fball_set_ui(T, 1);
            fball_sub(W, T, W, p);
            fball_add_error_2exp(W, m * xe);

            fball_mul_2exp_si(t, -1);
            fball_mul(res, t, W, p);

            fball_clear(z);
            fball_clear(W);
            fball_clear(T);
            break;
        }

        /* a' = t / 2, b' = sqrt(a b); a product too wide for the
           square root (wide inputs) ends in the fallback */
        fball_mul(b, a, b, p);
        if (fball_rel_2exp(b) > -40)
        {
            fball_mul_2exp_si(t, -1);
            fball_add_error_2exp(t, de - 1);
            fball_swap(res, t);
            break;
        }
        fball_sqrt(d, b, p);
        fball_swap(b, d);
        fball_mul_2exp_si(t, -1);
        fball_swap(a, t);
    }

    fball_clear(a);
    fball_clear(b);
    fball_clear(t);
    fball_clear(d);
}

void
fball_agm(fball_t res, const fball_t x, const fball_t y, slong n)
{
    _fball_agm_order(res, x, y, n, 0);
}
