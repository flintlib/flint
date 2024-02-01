/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"

void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))  /* g is small */
    {
        if (d >= 0)
            d = d >> FLINT_MIN(exp, SMALL_FMPZ_BITCOUNT_MAX);
        else
            d = -((-d) >> FLINT_MIN(exp, SMALL_FMPZ_BITCOUNT_MAX));

        fmpz_set_si(f, d);
    }
    else  /*g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);  /* g is already large */
        mpz_tdiv_q_2exp(mf, COEFF_TO_PTR(d), exp);
        _fmpz_demote_val(f);  /* division may make value small */
    }
}

void
fmpz_tdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_tdiv_q). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
            fmpz_set_si(f, c1 / c2);
        else                    /* h is large */
            fmpz_zero(f);
    }
    else                        /* g is large */
    {
        __mpz_struct * mf;

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            mf = _fmpz_promote(f);

            if (c2 > 0)         /* h > 0 */
            {
                flint_mpz_tdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_tdiv_q_ui(mf, COEFF_TO_PTR(c1), -c2);
                mpz_neg(mf, mf);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_tdiv_q_newton(f, g, h);
            }
            else
            {
                mf = _fmpz_promote(f);
                mpz_tdiv_q(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
                _fmpz_demote_val(f);    /* division by h may result in small value */
            }
        }
    }
}

void
fmpz_tdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception: division by zero in fmpz_tdiv_qr\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz q = c1 / c2;   /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            fmpz_set_si(f, q);
            fmpz_set_si(s, r);
        }
        else                    /* h is large and g is small */
        {
            fmpz_set_ui(f, WORD(0)); /* g is zero */
            fmpz_set_si(s, c1);
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
                flint_mpz_tdiv_qr_ui(mf, ms, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_tdiv_qr_ui(mf, ms, COEFF_TO_PTR(c1), -c2);
                mpz_neg(mf, mf);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
            _fmpz_demote_val(s);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_tdiv_qr_newton(f, s, g, h);
            }
            else
            {
                _fmpz_promote(f); /* must not hang on to ptr whilst promoting s */
                ms = _fmpz_promote(s);
                mf  = COEFF_TO_PTR(*f);

                mpz_tdiv_qr(mf, ms, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));

                _fmpz_demote_val(f);    /* division by h may result in small value */
                _fmpz_demote_val(s);    /* division by h may result in small value */
            }
        }
    }
}

void
fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
{
    fmpz c1 = *g;
    slong c2 = h;

    if (h == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_tdiv_q_si). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        fmpz_set_si(f, c1 / c2);
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);

        if (c2 > 0)
        {
            flint_mpz_tdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
        }
        else
        {
            flint_mpz_tdiv_q_ui(mf, COEFF_TO_PTR(c1), -(ulong) c2);
            mpz_neg(mf, mf);
        }
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}

void
fmpz_tdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
{
    fmpz c1 = *g;
    ulong c2 = h;

    if (h == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_tdiv_q_ui). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (c1 > 0)
        {
            fmpz_set_ui(f, c1 / c2);
        }
        else
        {
            ulong q = ((ulong) -c1) / c2;

            fmpz_set_si(f, - (slong) q);
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);

        flint_mpz_tdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}

void fmpz_tdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))  /* g is small */
    {
        if (d >= 0)
        {
            fmpz_set_ui(f, exp < (SMALL_FMPZ_BITCOUNT_MAX) ? d & ((WORD(1) << exp) - 1) : d);
        }
        else
        {
            d = -d;
            fmpz_neg_ui(f, exp < (SMALL_FMPZ_BITCOUNT_MAX) ? d & ((WORD(1) << exp) - 1) : d);
        }
    }
    else  /*g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);  /* g is already large */
        mpz_tdiv_r_2exp(mf, COEFF_TO_PTR(d), exp);
        _fmpz_demote_val(f);  /* division may make value small */
    }
}

ulong
fmpz_tdiv_ui(const fmpz_t g, ulong h)
{
    fmpz c1 = *g;

    if (h == UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_tdiv_ui). Division by 0.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
		/* We need the absolute value of the remainder and
		   C 90 guarantees truncation towards zero. */
		if (c1 < WORD(0))
			return -c1 % h;
		else
			return c1 % h;
    }
    else                        /* g is large */
    {
        return flint_mpz_tdiv_ui(COEFF_TO_PTR(c1), h);
    }
}
