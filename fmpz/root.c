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

    Copyright (C) 2011 William Hart

******************************************************************************/

#undef ulong /* avoid conflict with standard library */
#include <stdlib.h>
#define ulong unsigned long

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_root(fmpz_t r, fmpz_t f, len_t n)
{
    fmpz c = *f;
    
    if (n == 0)
    {
        printf("Exception (fmpz_root). Unable to take 0-th root.\n");
        abort();
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
            if (c < 0L)
            {
                printf("Exception (fmpz_root). Unable to take square root of negative value.\n");
                abort();
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
            if (c < 0L)
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
