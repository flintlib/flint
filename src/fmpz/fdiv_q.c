/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "fmpz.h"

void
fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_fdiv_q). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz q = c1 / c2;       /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            if (r && (c2 ^ r) < WORD(0))
                --q;

            fmpz_set_si(f, q);
        }
        else                    /* h is large and g is small */
        {
            if ((c1 > WORD(0) && fmpz_sgn(h) < 0) || (c1 < WORD(0) && fmpz_sgn(h) > 0))  /* signs are the same */
                fmpz_set_si(f, WORD(-1));   /* quotient is negative, round down to minus one */
            else
                fmpz_zero(f);
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mf;

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            mf = _fmpz_promote(f);

            if (c2 > 0)         /* h > 0 */
            {
                flint_mpz_fdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_cdiv_q_ui(mf, COEFF_TO_PTR(c1), -c2);
                mpz_neg(mf, mf);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_fdiv_q_newton(f, g, h);
            }
            else
            {
                mf = _fmpz_promote(f);
                mpz_fdiv_q(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
                _fmpz_demote_val(f);  /* division by h may result in small value */
            }
        }
    }
}
