/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m)
{
    if (fmpz_sgn(m) <= 0)
    {
        printf("Exception (fmpz_powm). Modulus is less than 1.\n");
        abort();
    }
    else if (!COEFF_IS_MPZ(*e))  /* e is small */
    {
        fmpz_powm_ui(f, g, *e, m);
    }
    else  /* e is large */
    {
        if (!COEFF_IS_MPZ(*m))  /* m is small */
        {
            ulong g1 = fmpz_fdiv_ui(g, *m);
            mpz_t g2, m2;
            __mpz_struct *mpz_ptr;

            mpz_init_set_ui(g2, g1);
            mpz_init_set_ui(m2, *m);
            mpz_ptr = _fmpz_promote(f);

            mpz_powm(mpz_ptr, g2, COEFF_TO_PTR(*e), m2);

            mpz_clear(g2);
            mpz_clear(m2);
            _fmpz_demote_val(f);
        }
        else  /* m is large */
        {
            if (!COEFF_IS_MPZ(*g))  /* g is small */
            {
                mpz_t g2;
                __mpz_struct *mpz_ptr;

                mpz_init_set_si(g2, *g);
                mpz_ptr = _fmpz_promote(f);

                mpz_powm(mpz_ptr, g2, COEFF_TO_PTR(*e), COEFF_TO_PTR(*m));

                mpz_clear(g2);
                _fmpz_demote_val(f);
            }
            else  /* g is large */
            {
                __mpz_struct *mpz_ptr = _fmpz_promote(f);

                mpz_powm(mpz_ptr, 
                    COEFF_TO_PTR(*g), COEFF_TO_PTR(*e), COEFF_TO_PTR(*m));
                _fmpz_demote_val(f);
            }
        }
    }
}

