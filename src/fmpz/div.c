/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "div_small.h"

/*
    Checked exact division: if h divides g, sets q = g / h and returns 1;
    otherwise returns 0 and sets q to 0. Division by zero is exact iff g is
    zero. The multi-limb case goes through flint_mpn_div.
*/
int
fmpz_div(fmpz_t q, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g, c2 = *h;

    if (!COEFF_IS_MPZ(c2))
        return fmpz_div_si(q, g, c2);

    if (!COEFF_IS_MPZ(c1))
    {
        /* |g| < |h|: exact iff g = 0 */
        fmpz_zero(q);
        return c1 == 0;
    }

    {
        mpz_srcptr mg = COEFF_TO_PTR(c1), mh = COEFF_TO_PTR(c2);
        mp_size_t gn = FLINT_ABS(mg->_mp_size), hn = FLINT_ABS(mh->_mp_size), qn;
        int qneg = (mg->_mp_size < 0) ^ (mh->_mp_size < 0);
        mpz_ptr mq;
        mp_ptr qd;
        int exact;

        if (gn < hn)
        {
            fmpz_zero(q);
            return 0;
        }

        qn = gn - hn + 1;

        if (q == g || q == h)
        {
            fmpz_t t;
            fmpz_init(t);
            exact = fmpz_div(t, g, h);
            fmpz_swap(q, t);
            fmpz_clear(t);
            return exact;
        }

        mq = _fmpz_promote(q);
        qd = FLINT_MPZ_REALLOC(mq, qn);
        exact = flint_mpn_div(qd, mg->_mp_d, gn, mh->_mp_d, hn);

        if (!exact)
        {
            mq->_mp_size = 0;
        }
        else
        {
            while (qn > 0 && qd[qn - 1] == 0)
                qn--;
            mq->_mp_size = qneg ? -qn : qn;
        }

        _fmpz_demote_val(q);
        return exact;
    }
}

int
fmpz_div_ui(fmpz_t q, const fmpz_t g, ulong h)
{
    fmpz c1 = *g;

    if (h == 0)
    {
        fmpz_zero(q);
        return c1 == 0;
    }

    if (!COEFF_IS_MPZ(c1))
    {
        ulong a = FLINT_ABS(c1), r;
        int exact;

        if (a < h)
        {
            fmpz_zero(q);
            return c1 == 0;
        }

        r = a % h;
        exact = (r == 0);
        if (exact)
        {
            fmpz_set_ui(q, a / h);
            if (c1 < 0)
                fmpz_neg(q, q);
        }
        else
        {
            fmpz_zero(q);
        }
        return exact;
    }
    else
    {
        mpz_srcptr mg = COEFF_TO_PTR(c1);
        mp_size_t gn = FLINT_ABS(mg->_mp_size);
        int gneg = mg->_mp_size < 0;
        mpz_ptr mq;
        mp_ptr qd;
        mp_limb_t r;

        /* one or two limbs: hardware division, no mpz output unless the
           quotient needs it */
        if (gn <= 2 && h <= (ulong) WORD_MAX)
        {
            fmpz_t rem;
            int exact;
            fmpz_init(rem);
            _fmpz_div_qr_small_divisor(q, rem, g, (slong) h, 0);
            exact = fmpz_is_zero(rem);
            fmpz_clear(rem);
            if (!exact)
                fmpz_zero(q);
            return exact;
        }

        mq = _fmpz_promote(q);   /* q may alias g: the limbs are read below
                                    from the same (unchanged) array */
        mg = COEFF_TO_PTR(*g);
        qd = FLINT_MPZ_REALLOC(mq, gn);
        r = mpn_divrem_1(qd, 0, mg->_mp_d, gn, h);

        if (r != 0)
        {
            mq->_mp_size = 0;
            _fmpz_demote_val(q);
            return 0;
        }

        while (gn > 0 && qd[gn - 1] == 0)
            gn--;
        mq->_mp_size = gneg ? -gn : gn;
        _fmpz_demote_val(q);
        return 1;
    }
}

int
fmpz_div_si(fmpz_t q, const fmpz_t g, slong h)
{
    int exact;

    if (h >= 0)
        return fmpz_div_ui(q, g, (ulong) h);

    exact = fmpz_div_ui(q, g, -(ulong) h);
    fmpz_neg(q, q);
    return exact;
}
