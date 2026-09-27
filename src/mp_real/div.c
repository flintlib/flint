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

/* Division by a ball and by a word. */

/* res = a / c truncated to about n limbs, c >= 1: the quotient of the
   top min(n + 2, size) limbs of the mantissa over one limb of zeros
   (two for an exact dividend) by mpn_divrem_1, so that the result keeps
   the quotient's own precision: a quotient re-anchored at the
   dividend's bottom limb with its truncation error there loses
   log2(c) bits relative to its magnitude at every division, which a
   chain of divisions (a common denominator k^m) compounds to nothing.
   The radius err / c and the dropped low limbs of the dividend are
   bounded one limb below the dividend's anchor, where the truncation
   of the quotient is one unit; for an exact dividend the quotient
   runs two limbs below, its truncation one unit there. */
void
mp_real_div_ui(mp_real_t res, const mp_real_t a, ulong c, slong n)
{
    mp_real_bnd_t e;
    nn_ptr t;
    slong sa = a->size, keep, drop, anc, qxn;

    if (c == 0)
        flint_throw(FLINT_ERROR, "mp_real_div_ui: division by zero\n");

    if (sa == 0)
    {
        e = (a->err != 0) ? _mp_real_bnd_div_ui_B(a->err, c, a->exp - 1) : _mp_real_bnd_zero;
        _mp_real_zero_bnd(res, e);
        return;
    }

    keep = FLINT_MIN(n + 2, sa);
    drop = sa - keep;
    anc = a->exp - a->size + drop;  /* ulp of the kept part */

    if (a->err == 0 && drop == 0)
    {
        /* an exact dividend: the quotient to n + 2 limbs (at least two
           below the dividend), truncated one unit at its bottom */
        qxn = FLINT_MAX(2, n + 2 - keep);
        e = _mp_real_bnd(0, 1, anc - qxn);
    }
    else
    {
        /* the radius and the dropped limbs of a (below B^anc, so below
           B^anc / c after division) at B^(anc - 1), the quotient's
           truncation one unit there */
        qxn = 1;
        /* separately: a->err + 1 wraps for err = B - 1 */
        e = _mp_real_bnd_div_ui_B(a->err, c, anc - 1);
        if (drop > 0)
            e = _mp_real_bnd_add(e, _mp_real_bnd_div_ui_B(1, c, anc - 1));
    }

    mp_real_fit_length(res, keep + qxn);
    t = res->d;
    if (res != a)
        mpn_divrem_1(t, qxn, a->d + drop, keep, c);
    else
    {
        memmove(t + qxn, t + drop, keep * sizeof(ulong));
        flint_mpn_zero(t, qxn);
        mpn_divrem_1(t, 0, t, keep + qxn, c);
    }

    res->size = keep + qxn;
    res->exp = anc + keep;
    res->negative = a->negative;
    res->err = 1;       /* norm must know inexactness */
    _mp_real_norm(res);
    _mp_real_apply_bnd(res, e);
}

/* res = a / b.  From FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF limbs (of
   both the divisor and the quotient) the mantissas are divided as
   fractions by _mp_real_div_newton, which takes the divisor as it is (no
   normalization shifts; one extra guard limb absorbs a small top
   limb); below, flint_mpn_divapprox_fraction picks the register,
   schoolbook, blockwise or divide-and-conquer division and needs
   neither zero-extended nor normalized operands.  Operand errors
   enter as err_a / |b| + |a| err_b / |b|^2, with |a| and |b| bounded
   through their top limbs and the -err_b correction of the
   denominator covered by a 1.01 factor, which is rigorous because the
   relative radius of b is CHECKED to be below 2^-30 (a divisor ball
   this wide is a usage error: the mantissa would be pure noise). */
#define MP_REAL_BALL_DIV_NEWTON_CUTOFF FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF

void
mp_real_div(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
{
    slong sa = a->size, sb = b->size, nd, qn;
    int rneg;
    double mb;
    mp_real_dbnd_t e;
    TMP_INIT;

    if (sb == 0)
        flint_throw(FLINT_ERROR, "mp_real_div: division by zero ball\n");
    if (_mp_real_rel_bound(b) > 0x1p-30)
        flint_throw(FLINT_ERROR,
            "mp_real_div: divisor relative radius above 2^-30\n");

    /* magnitudes through the top limbs, as in mp_real_mul:
       |a| < ma B^(a.exp - 1), |b| >= mb B^(b.exp - 1) with
       ma = a_top + 1 rounded up and mb = b_top rounded down; the
       blanket bounds |a| < B^a.exp, |b| >= B^(b.exp - 1) overstate
       the quotient by up to a limb each, and |a| / |b|^2 by up to
       three */
    e = _mp_real_dbnd(0.0, 0);
    mb = _mp_real_mag_lo(b);
    {
        if (a->err != 0)
            e = _mp_real_dbnd_add(e, _mp_real_dbnd((double) a->err * 1.01 / mb,
                    (a->exp - sa) + 1 - b->exp));
        if (b->err != 0)
        {
            double ma = (sa == 0) ? 1.0 : _mp_real_mag_hi(a);
            slong ea = (sa == 0) ? a->exp + 1 : a->exp;
            e = _mp_real_dbnd_add(e, _mp_real_dbnd((double) b->err * 1.01 * (ma / (mb * mb)),
                    ea - 1 + (b->exp - sb) + 2
                        - 2 * b->exp));
        }
    }

    if (sa == 0)
    {
        _mp_real_zero_err(res, e);
        return;
    }

    nd = FLINT_MIN(n, FLINT_MIN(_mp_real_acc(a, n), _mp_real_acc(b, n)));
    nd = FLINT_MAX(nd, 2);
    qn = nd + 2;

    rneg = a->negative ^ b->negative;

    TMP_START;
    {
        int aliased = (res == a || res == b);
        slong rn, E;
        nn_ptr q;

        if (sb >= MP_REAL_BALL_DIV_NEWTON_CUTOFF && qn >= MP_REAL_BALL_DIV_NEWTON_CUTOFF)
        {
            /* Karp-Markstein on the mantissas read as fractions,
               (a_d B^-sa) / (b_d B^-sb): _mp_real_div_newton does not
               need a normalized divisor, only a nonzero top limb, so
               den = b_d B^-sb >= b_top / B and the error
               4 B^-(nd+1) / den is at most 4 B / b_top ulps of the
               nd + 1 fraction limbs -- one more than otherwise needed,
               which covers the unnormalized divisor */
            rn = nd + 3;
            q = aliased ? TMP_ALLOC(rn * sizeof(ulong))
                : (mp_real_fit_length(res, rn), res->d);
            _mp_real_div_newton(q, a->d, sa, b->d, sb, nd + 1);
            E = a->exp - b->exp + 2 - rn;
            e = _mp_real_dbnd_add(e, _mp_real_dbnd(4.0 / mb, E + 1));
        }
        else
        {
            /* approximate division by the mantissa of b (truncated
               to qn + 2 limbs when longer): Q = floor(a_d B^f / b_d)
               or one more, with f chosen for qn + 1 quotient limbs
               (negative f truncates the numerator, exactly as a
               floor); |Q - true| < 1 ulp of B^E */
            nn_srcptr bd = b->d;
            slong sbt = sb, f;

            if (sb > qn + 2)
            {
                /* b' = b - delta, 0 <= delta < B^(b.exp - sbt):
                   a/b' - a/b <= |a| delta / b'^2
                   <= (ma / mb^2) B^(a.exp - b.exp - sbt + 1) */
                double ma = _mp_real_mag_hi(a);
                sbt = qn + 2;
                bd = b->d + (sb - sbt);
                e = _mp_real_dbnd_add(e, _mp_real_dbnd(1.01 * ma / (mb * mb),
                    a->exp - b->exp - sbt + 1));
            }

            rn = qn + 1;
            f = qn - sa + sbt;
            q = aliased ? TMP_ALLOC(rn * sizeof(ulong))
                : (mp_real_fit_length(res, rn), res->d);
            flint_mpn_divapprox_fraction(q, a->d, sa, bd, sbt, f);
            E = (a->exp - sa) - f - (b->exp - sbt);
            e = _mp_real_dbnd_add(e, _mp_real_dbnd(1.0, E));
        }

        if (aliased)
        {
            mp_real_fit_length(res, rn);
            flint_mpn_copyi(res->d, q, rn);
        }

        res->size = rn;
        res->exp = E + rn;
    }
    res->negative = rneg;
    res->err = (e.v != 0.0);    /* norm must know inexactness */
    _mp_real_norm(res);
    if (res->size == 0)
        _mp_real_zero_err(res, e);
    else
        _mp_real_apply_err(res, e);
    TMP_END;
}
