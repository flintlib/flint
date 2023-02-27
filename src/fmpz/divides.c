/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
fmpz_divides(fmpz_t q, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;
    int res, negate = 0;

    if (fmpz_is_zero(h))
    {
        res = fmpz_is_zero(g);

        fmpz_zero(q);

        return res;
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {	
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            mp_limb_t qz;

            if (c1 < 0)
            {
                c1 = -c1;
                negate ^= 1;
            }
            if (c2 < 0)
            {
                c2 = -c2;
                negate ^= 1;
            }

            res = n_divides(&qz, c1, c2);
            fmpz_set_ui(q, qz);
            if (negate)
                fmpz_neg(q, q);

            return res;
        }
        else                    /* h is large and g is small */
        {
            res = fmpz_is_zero(g);

            fmpz_zero(q);

            return res;
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mq = _fmpz_promote(q);

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            mp_limb_t r;

            if (c2 < 0)
            {
                c2 = -c2;
                negate ^= 1;
            }

            r = flint_mpz_tdiv_q_ui(mq, COEFF_TO_PTR(c1), c2);

            res = (r == 0);
        }
        else                    /* both are large */
        {
            mpz_t r;
            mq  = _fmpz_promote(q);

            mpz_init(r);

            mpz_tdiv_qr(mq, r, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));

            res = mpz_sgn(r) == 0;

            mpz_clear(r);
        }

        if (!res)
            mpz_set_ui(mq, 0);

        _fmpz_demote_val(q);

        if (negate)
            fmpz_neg(q, q);

        return res;
    }
}
