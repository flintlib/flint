/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
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
            fmpz q = c1 / c2;   /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            if ((c2 > WORD(0) && r < WORD(0)) || (c2 < WORD(0) && r > WORD(0)))
            {
                q--;            /* q cannot overflow as remainder implies |c2| != 1 */
                r += c2;
            }

            fmpz_set_si(f, q);
            fmpz_set_si(s, r);
        }
        else                    /* h is large and g is small */
        {
            if (c1 == WORD(0))
            {
                fmpz_set_ui(f, WORD(0)); /* g is zero */
                fmpz_set_si(s, c1);
            }
            else if ((c1 < WORD(0) && fmpz_sgn(h) < 0) || (c1 > WORD(0) && fmpz_sgn(h) > 0))  /* signs are the same */
            {
                fmpz_zero(f);   /* quotient is positive, round down to zero */
                fmpz_set_si(s, c1);
            }
            else
            {
                fmpz_add(s, g, h);
                fmpz_set_si(f, WORD(-1));    /* quotient is negative, round down to minus one */
            }
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mf, * ms;

		if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            _fmpz_promote(f); /* must not hang on to ptr whilst promoting s */
            ms = _fmpz_promote(s);
            mf  = COEFF_TO_PTR(*f);

            if (c2 > 0)         /* h > 0 */
            {
                flint_mpz_fdiv_qr_ui(mf, ms, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_cdiv_qr_ui(mf, ms, COEFF_TO_PTR(c1), -c2);
                mpz_neg(mf, mf);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
            _fmpz_demote_val(s);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_fdiv_qr_newton(f, s, g, h);
            }
            else
            {
                _fmpz_promote(f); /* must not hang on to ptr whilst promoting s */
                ms = _fmpz_promote(s);
                mf  = COEFF_TO_PTR(*f);

                mpz_fdiv_qr(mf, ms, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));

                _fmpz_demote_val(f);    /* division by h may result in small value */
                _fmpz_demote_val(s);    /* division by h may result in small value */
            }
        }
    }
}
