/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_fdiv_r_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))  /* g is small */
    {
        if (d >= 0)
        {
            fmpz_set_ui(f, exp < (FLINT_BITS - 2) ? d & ((WORD(1) << exp) - 1) : d);
        }
        else
        {
            if (exp <= FLINT_BITS - 2)
            {
                fmpz_set_ui(f, d & ((WORD(1) << exp) - 1));
            }
            else
            {
                __mpz_struct * mpz_ptr = _fmpz_promote(f);

                flint_mpz_set_ui(mpz_ptr, 1);
                mpz_mul_2exp(mpz_ptr, mpz_ptr, exp);
                flint_mpz_sub_ui(mpz_ptr, mpz_ptr, -d);
            }
        }
    }
    else  /*g is large */
    {
        __mpz_struct * mpz_ptr = _fmpz_promote(f);  /* g is already large */
        mpz_fdiv_r_2exp(mpz_ptr, COEFF_TO_PTR(d), exp);
        _fmpz_demote_val(f);  /* division may make value small */
    }
}

