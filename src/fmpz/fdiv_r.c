/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

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
fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_fdiv_r). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz q = c1 / c2;   /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            if ((c2 > WORD(0) && r < WORD(0)) || (c2 < WORD(0) && r > WORD(0)))
                r += c2;

            fmpz_set_si(f, r);
        }
        else                    /* h is large and g is small */
        {
            if (c1 == WORD(0))
            {
                fmpz_set_si(f, c1);
            }
            else if ((c1 < WORD(0) && fmpz_sgn(h) < 0) || (c1 > WORD(0) && fmpz_sgn(h) > 0))  /* signs are the same */
            {
                fmpz_set_si(f, c1);
            }
            else
            {
                fmpz_add(f, g, h);
            }
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mf;

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            /* todo: should not create an mpz */
            mf = _fmpz_promote(f);

            if (c2 > 0)         /* h > 0 */
            {
                flint_mpz_fdiv_r_ui(mf, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_cdiv_r_ui(mf, COEFF_TO_PTR(c1), -c2);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_fdiv_r_newton(f, g, h);
            }
            else
            {
                mf = _fmpz_promote(f);
                mpz_fdiv_r(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
                _fmpz_demote_val(f);    /* division by h may result in small value */
            }
        }
    }
}
