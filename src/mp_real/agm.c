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
#include "mp_real.h"
#include "impl.h"

/* The arithmetic-geometric mean of two positive balls.

   THE ITERATION a' = (a + b)/2, b' = sqrt(a b) converges quadratically:
   with z = (a - b)/(a + b), z' = (1 - sqrt(1 - z^2))/(1 + sqrt(1 - z^2))
   ~ z^2 / 4.  Each step is a ball sum, product and square root at the
   working precision; mp_real carries the radii (inexact inputs lower the
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
_mp_real_agm_order(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n,
    int m)
{
    mp_real_t a, b, t, d;
    /* at least 128 bits: the relative radii are held below 2^-40
       over the iterations (n + 1 on 64-bit machines) */
    slong p = FLINT_MAX(n + 1, 128 / FLINT_BITS), iter;

    if (x->negative || y->negative)
        flint_throw(FLINT_ERROR, "mp_real_agm: negative argument\n");
    if (x->size == 0 || y->size == 0)
    {
        /* agm(0, y) = 0; with a radius, the fallback below is the
           enclosure 0 +- max(|x|, |y|) */
        if (mp_real_is_zero(x) || mp_real_is_zero(y))
        {
            mp_real_zero(res);
            return;
        }
    }

    if (m == 0)
        m = _agm_default_order(n);
    m = FLINT_MAX(2, FLINT_MIN(m, AGM_WORD_ORDER));

    mp_real_init(a);
    mp_real_init(b);
    mp_real_init(t);
    mp_real_init(d);
    mp_real_set(a, x);
    mp_real_set(b, y);

    for (iter = 0; ; iter++)
    {
        slong de, tl;

        mp_real_add(t, a, b, p);
        mp_real_sub(d, a, b, p);

        if (mp_real_is_zero(d))
        {
            mp_real_mul_2exp_si(t, t, -1);
            mp_real_swap(res, t);
            break;
        }

        /* |z| < 2^ze with ze = de - tl: |d| < 2^de, |t| >= 2^tl */
        de = mp_real_abs_bound_lt_2exp_si(d);
        tl = (t->size == 0) ? WORD_MIN / 4
            : FLINT_BITS * (t->exp - 1) + FLINT_BIT_COUNT(t->d[t->size - 1]) - 1;

        if (t->size == 0 || d->size == 0 || mp_real_rel_radius_lt_2exp_si(d) > -2
            || mp_real_rel_radius_lt_2exp_si(t) > -40 || iter > 2 * FLINT_BITS + 64)
        {
            /* the fallback: agm in t/2 +- |d|/2 */
            mp_real_mul_2exp_si(t, t, -1);
            mp_real_add_error_2exp_si(t, de - 1);
            mp_real_swap(res, t);
            break;
        }

        if (2 * m * (de - tl) < -FLINT_BITS * p)
        {
            /* the finish: x = z^2, 1 - sum_{j<m} c_j x^j over the
               common denominator 2^E, tail below x^m */
            mp_real_t z, W, T;
            mp_real_struct pw[AGM_MAX_ORDER];
            slong xe, cp, j, E = agm_exp[m - 2];

            mp_real_init(z);
            mp_real_init(W);
            mp_real_init(T);

            mp_real_div(z, d, t, FLINT_MAX(2, p + (de - tl) / FLINT_BITS + 2));
            mp_real_mul(z, z, z, FLINT_MAX(2, p + 2 * (de - tl) / FLINT_BITS + 2));
            xe = mp_real_abs_bound_lt_2exp_si(z);

            for (j = 1; j < m; j++)
            {
                const mp_real_struct * xj;

                cp = FLINT_MAX(2, p + (j * xe) / FLINT_BITS + 2);
                if (j == 1)
                    xj = z;
                else
                {
                    mp_real_init(pw + j);
                    if (j % 2 == 0)
                        mp_real_mul(pw + j, (j / 2 == 1) ? z : pw + j / 2,
                            (j / 2 == 1) ? z : pw + j / 2, cp);
                    else
                        mp_real_mul(pw + j, pw + j - 1, z, cp);
                    xj = pw + j;
                }
                if (j == 1)
                    mp_real_mul_ui(W, xj, (ulong) (agm_num[0] << (E - agm_exp[0])), cp);
                else
                    mp_real_addmul_ui(W, W, xj, (ulong) (agm_num[j - 1] << (E - agm_exp[j - 1])),
                        FLINT_MAX(2, p + xe / FLINT_BITS + 2));
            }
            for (j = 2; j < m; j++)
                mp_real_clear(pw + j);

            mp_real_mul_2exp_si(W, W, -E);
            mp_real_set_ui(T, 1);
            mp_real_sub(W, T, W, p);
            mp_real_add_error_2exp_si(W, m * xe);

            mp_real_mul_2exp_si(t, t, -1);
            mp_real_mul(res, t, W, p);

            mp_real_clear(z);
            mp_real_clear(W);
            mp_real_clear(T);
            break;
        }

        /* a' = t / 2, b' = sqrt(a b); a product too wide for the
           square root (wide inputs) ends in the fallback */
        mp_real_mul(b, a, b, p);
        if (mp_real_rel_radius_lt_2exp_si(b) > -40)
        {
            mp_real_mul_2exp_si(t, t, -1);
            mp_real_add_error_2exp_si(t, de - 1);
            mp_real_swap(res, t);
            break;
        }
        mp_real_sqrt(d, b, p);
        mp_real_swap(b, d);
        mp_real_mul_2exp_si(t, t, -1);
        mp_real_swap(a, t);
    }

    mp_real_clear(a);
    mp_real_clear(b);
    mp_real_clear(t);
    mp_real_clear(d);
}

void
mp_real_agm(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n)
{
    _mp_real_agm_order(res, x, y, n, 0);
}
