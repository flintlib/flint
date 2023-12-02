/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

#define TMP_ALLOC_LIMBS(__n) TMP_ALLOC((__n) * sizeof(mp_limb_t))

/*
Compute wn-limb fixed-point number w, a number of ulps error, and
an exponent q such that x = (w + xi * error) - q * log(2)
where 0 <= w < log(2) and |xi| <= 1.

Returns 0 if unsuccessful (high enough precision for log(2) to be available)
and 1 if successful.

Assumes nonspecial x; its exponent must not be an mpz.

Error analysis:

We first set

    n = |x| - e1

    d = log(2) - e2

where 0 <= e1, e2 <= 1 ulp are the errors resulting from truncating
to fixed-point numbers (note that the value of log(2) is correctly rounded).

Next, we compute q, r such that 0 <= r < d and q*d + r = n (this
is just integer division with remainder, done exactly).

The real number we want to approximate is r_exact = |x| - q * log(2).
The approximation r satisfies

    |r_exact - r| = |(|x| - q * log(2)) - (n - q * d)|
                  = ||x| - q * log(2) - (|x| - e1) + q * (log(2) - e2)|
                  = |e1 - q * e2|
                  <= (q + 1) ulp.

We select the working precision so that (q + 1) ulp is at most
1 ulp in the target precision.

It is sufficient to use tn extra limbs where (q + 1) <= 2^(FLINT_BITS * tn).

Note that q + 1 <= n / d + 1 <= |x| * 1.5 + 1 < 2^(exp+2). So it is
sufficient to choose tn = ceil((exp+2)/FLINT_BITS).

Now we round the result to the final precision: w = r - e3.
This can add 1 more ulp of error, so the error may be 2 ulp.

Finally, if x < 0, we correct the sign by setting,x = log(2) - x,
q = -(q+1). This adds 1 more ulp (from the approximation of log(2)),
for a total of 3 ulp.

*/


int
_arb_get_mpn_fixed_mod_log2(mp_ptr w, fmpz_t q, mp_limb_t * error,
                                                const arf_t x, mp_size_t wn)
{
    mp_srcptr xp;
    mp_size_t xn;
    int negative;
    slong exp;

    ARF_GET_MPN_READONLY(xp, xn, x);
    exp = ARF_EXP(x);
    negative = ARF_SGNBIT(x);

    if (exp <= -1)
    {
        /* todo: just zero top */
        flint_mpn_zero(w, wn);

        *error = _arf_get_integer_mpn(w, xp, xn, exp + wn * FLINT_BITS);

        if (!negative)
        {
            fmpz_zero(q);
        }
        else
        {
            if (wn > ARB_LOG_TAB2_LIMBS)
                return 0;

            mpn_sub_n(w, arb_log_log2_tab + ARB_LOG_TAB2_LIMBS - wn, w, wn);
            *error += 1;    /* log(2) has 1 ulp error */
            fmpz_set_si(q, -1);
        }

        return 1; /* success */
    }
    else
    {
        mp_ptr qp, rp, np;
        mp_srcptr dp;
        mp_size_t qn, rn, nn, dn, tn, alloc;
        TMP_INIT;

        tn = ((exp + 2) + FLINT_BITS - 1) / FLINT_BITS;

        dn = wn + tn;           /* denominator */
        nn = wn + 2 * tn;       /* numerator */
        qn = nn - dn + 1;       /* quotient */
        rn = dn;                /* remainder */

        if (dn > ARB_LOG_TAB2_LIMBS)
            return 0;

        TMP_START;

        alloc = qn + rn + nn;
        qp = TMP_ALLOC_LIMBS(alloc);
        rp = qp + qn;
        np = rp + rn;

        dp = arb_log_log2_tab + ARB_LOG_TAB2_LIMBS - dn;

        /* todo: prove that zeroing is unnecessary */
        flint_mpn_zero(np, nn);

        _arf_get_integer_mpn(np, xp, xn, exp + dn * FLINT_BITS);

        mpn_tdiv_qr(qp, rp, 0, np, nn, dp, dn);

        if (!negative)
        {
            flint_mpn_copyi(w, rp + tn, wn);
            *error = 2;
        }
        else
        {
            if (mpn_add_1(qp, qp, qn, 1))
            {
                /* I believe this cannot happen (should prove it) */
                flint_throw(FLINT_ERROR, "mod log(2): unexpected carry\n");
            }

            mpn_sub_n(w, dp + tn, rp + tn, wn);
            *error = 3;
        }

        /* read the exponent */
        while (qn > 1 && qp[qn-1] == 0)
            qn--;

        if (qn == 1)
        {
            if (!negative)
                fmpz_set_ui(q, qp[0]);
            else
                fmpz_neg_ui(q, qp[0]);
        }
        else
        {
            fmpz_set_mpn_large(q, qp, qn, negative);
        }

        TMP_END;

        return 1;
    }
}

