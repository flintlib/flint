/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "impl.h"

/* Square roots and reciprocal square roots. */

void
mp_real_rsqrt_ui(mp_real_t res, ulong c, slong n)
{
    if (c <= 1)
    {
        if (c == 0)
            flint_throw(FLINT_ERROR, "mp_real_rsqrt_ui: division by zero\n");
        mp_real_set_ui(res, 1);
        return;
    }

    mp_real_fit_length(res, n);
    _mp_real_rsqrt_ui_newton(res->d, c, n);
    res->size = n;
    res->exp = 0;
    res->negative = 0;
    res->err = 1;       /* norm must know inexactness */
    _mp_real_norm(res);
    _mp_real_apply_err(res, _mp_real_dbnd(2.0, -n));
}

/* shared alignment for sqrt/rsqrt: express x = ahat * B^E with E
   even, ahat in [B^-2, 1) given as an an-limb fraction (aa, an);
   the input radius keeps its anchor: err ulps of B^(E - an) */
#define MP_REAL_SQRT_ALIGN                                        \
    slong E = x->exp, an = x->size, nd;                         \
    nn_ptr aa;                                                  \
    mp_real_dbnd_t e;                                                  \
    double ml;                                                  \
    TMP_INIT;                                                   \
    FLINT_ASSERT((const void *) res != (const void *) x);       \
    FLINT_ASSERT(x->size > 0 && !x->negative);                  \
    if (_mp_real_rel_bound(x) > 0x1p-30)                          \
        flint_throw(FLINT_ERROR,                                \
            "mp_real sqrt: relative radius above 2^-30\n");       \
    nd = FLINT_MIN(n, _mp_real_acc(x, n));                        \
    nd = FLINT_MAX(nd, 2);                                      \
    TMP_START;                                                  \
    if (E & 1)                                                  \
    {                                                           \
        aa = TMP_ALLOC((an + 1) * sizeof(ulong));               \
        flint_mpn_copyi(aa, x->d, an);                          \
        aa[an] = 0;                                             \
        an++;                                                   \
        E++;                                                    \
    }                                                           \
    else                                                        \
        aa = (nn_ptr) x->d;

/* res = 1/sqrt(x); x > 0, n-limb target */
void
mp_real_rsqrt(mp_real_t res, const mp_real_t x, slong n)
{
    MP_REAL_SQRT_ALIGN

    /* operand error through the derivative: with |x| >= ml B^(E-1-q)
       (ml the leading bits of the mantissa, q = 1 iff the alignment
       padded an odd exponent with a zero top limb),
       |d(1/sqrt x)/dx| = 1/(2 x^(3/2)) <= ml^(-3/2) B^(3(1+q-E)/2)/2
       and the radius is err B^(E - an), so
       |Delta| <= 0.51 err ml^(-3/2) B^(3(1+q)/2) B^(-an-E/2). */
    e = _mp_real_dbnd(0.0, 0);
    ml = _mp_real_mag_lo(x);
    if (x->err != 0)
    {
        int q = (int) (x->exp & 1);
        e = _mp_real_dbnd((double) x->err * 0.51 * pow(ml, -1.5)
            * (q ? MP_REAL_D_B * MP_REAL_D_B * MP_REAL_D_B : MP_REAL_D_B * MP_REAL_D_SQRTB),
            -an - E / 2);
    }

    mp_real_fit_length(res, nd + 2);
    _mp_real_rsqrt_newton(res->d, aa, an, nd);

    /* Newton: <= 4 B^-nd / sqrt(ahat) with ahat = ml B^(-1-q), so
       <= 4 ml^(-1/2) B^((1+q)/2) B^-nd, scaled B^(-E/2) */
    {
        int q = (int) (x->exp & 1);
        e = _mp_real_dbnd_add(e, _mp_real_dbnd(4.0 * sqrt(1.0 / ml)
            * (q ? MP_REAL_D_B : MP_REAL_D_SQRTB), -nd - E / 2));
    }

    res->size = nd + 2;
    res->exp = 2 - E / 2;
    res->negative = 0;
    res->err = (e.v != 0.0);    /* norm must know inexactness */
    _mp_real_norm(res);
    _mp_real_apply_err(res, e);
    TMP_END;
}

/* mp_real_sqrt takes the integer square root flint_mpn_sqrtrem (no
   remainder) of a 2m-limb input while the m-limb result is below this,
   _mp_real_sqrt_newton above: the Newton code, which needs no exact
   remainder, overtakes at about half the input length where
   flint_mpn_sqrtrem switches to its own Newton code (measured 1.6-2.5x
   faster below, crossover at m ~ 1300-2000 on x86-64) */
#define MP_REAL_BALL_SQRT_NEWTON_CUTOFF (FLINT_MPN_SQRTREM_NEWTON_CUTOFF / 4)

/* res = sqrt(x); x > 0, n-limb target */
void
mp_real_sqrt(mp_real_t res, const mp_real_t x, slong n)
{
    MP_REAL_SQRT_ALIGN

    /* operand error through the derivative: with |x| >= ml B^(E-1-q)
       as above, |d(sqrt x)/dx| = 1/(2 sqrt x) <= ml^(-1/2) B^((1+q-E)/2)/2
       and the radius err B^(E - an):
       |Delta| <= 0.51 err ml^(-1/2) B^((1+q)/2) B^(E/2 - an) */
    e = _mp_real_dbnd(0.0, 0);
    ml = _mp_real_mag_lo(x);
    if (x->err != 0)
    {
        int q = (int) (x->exp & 1);
        e = _mp_real_dbnd((double) x->err * 0.51 * sqrt(1.0 / ml)
            * (q ? MP_REAL_D_B : MP_REAL_D_SQRTB), E / 2 - an);
    }

    if (nd + 2 < MP_REAL_BALL_SQRT_NEWTON_CUTOFF)
    {
        /* integer square root of X = ahat B^(2m) (the aligned
           mantissa zero-extended, or truncated, to 2m limbs):
           S = floor(sqrt(X)) has m limbs and sqrt(x) = S B^(E/2 - m)
           to within 1 ulp, plus < 1/2 ulp when X was truncated
           (sqrt(X) - sqrt(X - t) < t / (2 sqrt(X - t)) for t < 1) */
        slong m = nd + 2, xn = 2 * m, t = xn - an;
        nn_ptr X;

        X = TMP_ALLOC(xn * sizeof(ulong));
        if (t >= 0)
        {
            flint_mpn_zero(X, t);
            flint_mpn_copyi(X + t, aa, an);
        }
        else
        {
            flint_mpn_copyi(X, aa - t, xn);
        }
        while (X[xn - 1] == 0)      /* the alignment's zero top limb */
            xn--;

        mp_real_fit_length(res, m);
        flint_mpn_zero(res->d, m);
        flint_mpn_sqrtrem(res->d, NULL, X, xn);
        e = _mp_real_dbnd_add(e, _mp_real_dbnd(1.5, E / 2 - m));

        res->size = m;
        res->exp = E / 2;
    }
    else
    {
        /* the Newton error 4 B^-nd / sqrt(ahat) is absolute while
           sqrt(ahat) can be as small as B^-1 (an odd exponent aligned
           by a zero top limb) or B^-(1/2): one or two limbs more of
           the fraction keep the result accurate to nd limbs of its
           own */
        int q = (int) (x->exp & 1);
        slong nd2 = nd + 1 + q;

        mp_real_fit_length(res, nd2 + 2);
        _mp_real_sqrt_newton(res->d, aa, an, nd2);

        /* Newton: <= 4 B^-nd2 / sqrt(ahat) with ahat = ml B^(-1-q):
           <= 4 ml^(-1/2) B^((1+q)/2) B^-nd2, scaled B^(E/2) */
        e = _mp_real_dbnd_add(e, _mp_real_dbnd(4.0 * sqrt(1.0 / ml)
            * (q ? MP_REAL_D_B : MP_REAL_D_SQRTB), -nd2 + E / 2));

        res->size = nd2 + 2;
        res->exp = 2 + E / 2;
    }
    res->negative = 0;
    res->err = (e.v != 0.0);    /* norm must know inexactness */
    _mp_real_norm(res);
    _mp_real_apply_err(res, e);
    TMP_END;
}
