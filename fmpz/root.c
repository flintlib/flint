/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#undef ulong
#define ulong mp_limb_t

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_root(fmpz_t r, const fmpz_t f, slong n)
{
    fmpz c = *f;
    
    if (n == 0)
    {
        flint_printf("Exception (fmpz_root). Unable to take 0-th root.\n");
        flint_abort();
    }

    if (n == 1)
    {
        fmpz_set(r, f);
        return;
    }
    
    if (!COEFF_IS_MPZ(c)) /* f is small */
    {
        if (n == 2)
        {
            if (c < WORD(0))
            {
                flint_printf("Exception (fmpz_root). Unable to take square root of negative value.\n");
                flint_abort();
            }

            fmpz_set_ui(r, n_sqrt(c));
        } else /* n > 2 */
        {
            __mpz_struct mpz2; 
            __mpz_struct * mpz1;
            mp_limb_t cval;

            if (c == 0)
            {
                fmpz_set_ui(r, 0);
                return;
            }

            mpz1 = _fmpz_promote(r);
            
            cval = FLINT_ABS(c);
            mpz2._mp_d = &cval; /* mock up an mpz */
            mpz2._mp_size = 1;
            if (c < WORD(0))
                mpz2._mp_size = -1;
            mpz2._mp_alloc = 1;

            mpz_root(mpz1, &mpz2, n);

            _fmpz_demote_val(r); /* root may be small */
        }
    } else /* f is large */
    {
        __mpz_struct * mpz2 = COEFF_TO_PTR(c);
        __mpz_struct * mpz1 = _fmpz_promote(r);
            
        mpz_root(mpz1, mpz2, n);
        _fmpz_demote_val(r); /* root may be small */
    }
}
