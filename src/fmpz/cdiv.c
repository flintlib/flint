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

void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))  /* g is small */
    {
        fmpz_set_si(f, -((-d) >> FLINT_MIN(exp, SMALL_FMPZ_BITCOUNT_MAX)));
    }
    else  /*g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);  /* g is already large */
        mpz_cdiv_q_2exp(mf, COEFF_TO_PTR(d), exp);
        _fmpz_demote_val(f);  /* division may make value small */
    }
}

void
fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_cdiv_q). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))  /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz q = c1 / c2;      /* compute C quotient */
            fmpz r = c1 - c2 * q;  /* compute remainder  */

            if (r && ((c2 ^ r) > WORD(0)))  /* r != 0, c2 and r same sign */
                ++q;

            fmpz_set_si(f, q);
        }
        else  /* h is large and g is small */
        {
            if ((c1 < WORD(0) && fmpz_sgn(h) < 0) || (c1 > WORD(0) && fmpz_sgn(h) > 0))  /* signs are the same */
                fmpz_one(f);  /* quotient is positive, round up to one */
            else
                fmpz_zero(f);  /* quotient is zero, or negative, round up to zero */
        }
    }
    else  /* g is large */
    {
        __mpz_struct * mf;

        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            mf = _fmpz_promote(f);

            if (c2 > 0)  /* h > 0 */
            {
                flint_mpz_cdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_fdiv_q_ui(mf, COEFF_TO_PTR(c1), -c2);
                mpz_neg(mf, mf);
            }

            _fmpz_demote_val(f);  /* division by h may result in small value */
        }
        else  /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_cdiv_q_newton(f, g, h);
            }
            else
            {
                mf = _fmpz_promote(f);
                mpz_cdiv_q(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
                _fmpz_demote_val(f);  /* division by h may result in small value */
            }
        }
    }
}

void fmpz_cdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_cdiv_q). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz q = c1 / c2;   /* compute C quotient */
            fmpz r = c1 - c2 * q;   /* compute remainder */

            /* rounding of division is undefined, but it should satisfy: */
            FLINT_ASSERT(FLINT_ABS(r) < FLINT_ABS(c2));

            if ((c2 > WORD(0) && r > WORD(0)) || (c2 < WORD(0) && r < WORD(0)))
            {
                q += 1; /* q cannot overflow as remainder implies |c2| != 1 */
                r -= c2;
            }

            fmpz_set_si(f, q);
            fmpz_set_si(s, r);
        }
        else                    /* h is large and g is small */
        {
            int sgn_h = fmpz_sgn(h);
            if ((c1 < WORD(0) && sgn_h < 0) || (c1 > WORD(0) && sgn_h > 0))  /* signs are the same */
            {
                fmpz_sub(s, g, h);
                fmpz_one(f);   /* quotient > 0, round up to 1 */
            }
            else
            {
                fmpz_set_si(s, c1);
                fmpz_zero(f);   /* quotient <= 0, round up to 0 */
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
                flint_mpz_cdiv_qr_ui(mf, ms, COEFF_TO_PTR(c1), c2);
            }
            else
            {
                flint_mpz_fdiv_qr_ui(mf, ms, COEFF_TO_PTR(c1), -c2);
                mpz_neg(mf, mf);
            }

            _fmpz_demote_val(f);    /* division by h may result in small value */
            _fmpz_demote_val(s);    /* division by h may result in small value */
        }
        else                    /* both are large */
        {
            if (MPZ_WANT_FLINT_DIVISION(COEFF_TO_PTR(c1), COEFF_TO_PTR(c2)))
            {
                _fmpz_cdiv_qr_newton(f, s, g, h);
            }
            else
            {
                _fmpz_promote(f); /* must not hang on to ptr whilst promoting s */
                ms = _fmpz_promote(s);
                mf  = COEFF_TO_PTR(*f);

                mpz_cdiv_qr(mf, ms, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));

                _fmpz_demote_val(f);    /* division by h may result in small value */
                _fmpz_demote_val(s);    /* division by h may result in small value */
            }
        }
    }
}

void
fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
{
    fmpz c1 = *g;
    slong c2 = h;

    if (h == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_cdiv_q_si). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        fmpz q = c1 / c2;       /* compute C quotient */
        fmpz r = c1 - c2 * q;   /* compute remainder */

        if (r && ((c1 ^ c2) > WORD(0)))
            ++q;

        fmpz_set_si(f, q);
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);

        if (c2 > 0)
        {
            flint_mpz_cdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
        }
        else
        {
            flint_mpz_fdiv_q_ui(mf, COEFF_TO_PTR(c1), -(ulong) c2);
            mpz_neg(mf, mf);
        }
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}

void
fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
{
    fmpz c1 = *g;
    ulong c2 = h;

    if (h == 0)
    {
        flint_throw(FLINT_ERROR, "Exception: division by zero in fmpz_cdiv_q_ui\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (c1 > 0)
        {
            ulong q = c1 / c2;
            ulong r = c1 - c2 * q;

            if (r)
                ++q;

            fmpz_set_ui(f, q);
        }
        else
        {
            fmpz_set_si(f, - (((ulong) -c1) / c2));
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);

        flint_mpz_cdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}

void fmpz_cdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))  /* g is small */
    {
        if (d <= 0)
        {
            d = -d;
            fmpz_neg_ui(f, exp < (SMALL_FMPZ_BITCOUNT_MAX) ? d & ((WORD(1) << exp) - 1) : d);
        }
        else
        {
            if (exp <= SMALL_FMPZ_BITCOUNT_MAX)
            {
                fmpz_neg_ui(f, (-d) & ((WORD(1) << exp) - 1));
            }
            else
            {
                __mpz_struct * mf = _fmpz_promote(f);

                flint_mpz_set_ui(mf, 1);
                mpz_mul_2exp(mf, mf, exp);
                flint_mpz_sub_ui(mf, mf, d);
                mpz_neg(mf, mf);
            }
        }
    }
    else  /*g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);  /* g is already large */
        mpz_cdiv_r_2exp(mf, COEFF_TO_PTR(d), exp);
        _fmpz_demote_val(f);  /* division may make value small */
    }
}

ulong
fmpz_cdiv_ui(const fmpz_t g, ulong h)
{
    fmpz c1 = *g;
    ulong r;

    if (h == UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_cdiv_ui). Division by 0.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (c1 >= WORD(1))
            r = h - 1 - ((c1 - WORD(1)) % h);
        else
            r = (-c1) % h;

        return r;
    }
    else                        /* g is large */
    {
        return flint_mpz_cdiv_ui(COEFF_TO_PTR(c1), h);
    }
}
