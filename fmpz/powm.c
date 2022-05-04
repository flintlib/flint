/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#ifdef LONGSLONG
#define flint_mpz_init_set_si mpz_init_set_si
#define flint_mpz_init_set_ui mpz_init_set_ui
#else
#include "gmpcompat.h"
#endif

void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m)
{
    if (fmpz_sgn(m) <= 0)
    {
        flint_throw(FLINT_ERROR, "Modulus is less than 1 in fmpz_powm\n");
    }
    else if (!COEFF_IS_MPZ(*e))  /* e is small */
    {
        if (*e >= 0)
        {
            fmpz_powm_ui(f, g, *e, m);
        }
        else
        {
            fmpz_t g_inv;
            fmpz_init(g_inv);
            if (!fmpz_invmod(g_inv, g, m))
            {
                fmpz_clear(g_inv);
                flint_throw(FLINT_ERROR, "Base not invertible in fmpz_powm\n");
            }
            else
            {
                fmpz_powm_ui(f, g_inv, -*e, m);
                fmpz_clear(g_inv);
            }
        }
    }
    else  /* e is large */
    {
        if (!COEFF_IS_MPZ(*m))  /* m is small */
        {
            ulong g1 = fmpz_fdiv_ui(g, *m);
            mpz_t g2, m2;
            __mpz_struct * mf;

            flint_mpz_init_set_ui(g2, g1);
            flint_mpz_init_set_ui(m2, *m);
            mf = _fmpz_promote(f);

            mpz_powm(mf, g2, COEFF_TO_PTR(*e), m2);

            mpz_clear(g2);
            mpz_clear(m2);
            _fmpz_demote_val(f);
        }
        else  /* m is large */
        {
            if (!COEFF_IS_MPZ(*g))  /* g is small */
            {
                mpz_t g2;
                __mpz_struct * mf;

                flint_mpz_init_set_si(g2, *g);
                mf = _fmpz_promote(f);

                mpz_powm(mf, g2, COEFF_TO_PTR(*e), COEFF_TO_PTR(*m));

                mpz_clear(g2);
                _fmpz_demote_val(f);
            }
            else  /* g is large */
            {
                __mpz_struct * mf = _fmpz_promote(f);

                mpz_powm(mf, 
                    COEFF_TO_PTR(*g), COEFF_TO_PTR(*e), COEFF_TO_PTR(*m));
                _fmpz_demote_val(f);
            }
        }
    }
}

