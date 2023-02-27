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

int
fmpz_root(fmpz_t r, const fmpz_t f, slong n)
{
    fmpz c = *f;
    
    if (n <= 0)
    {
        flint_printf("Exception (fmpz_root). Unable to take %wd-th root.\n", n);
        flint_abort();
    }

    if (n == 1)
    {
        fmpz_set(r, f);
        return 1;
    }
    
    if (!COEFF_IS_MPZ(c)) /* f is small */
    {
        mp_limb_t rem, root;
        int sgn = c < 0;

        if (n == 2)
        {
            if (sgn)
            {
                flint_printf("Exception (fmpz_root). Unable to take square root of negative value.\n");
                flint_abort();
            }

            root = n_sqrtrem(&rem, c);
            fmpz_set_ui(r, root);
            return rem == 0;
        } else if (n == 3)
        {
            if (sgn)
                c = -c;

            root = n_cbrtrem(&rem, c);
            fmpz_set_si(r, sgn ? -root : root);
            return rem == 0;
        } else /* n > 3 */
        {
            if (sgn)
            {
                if ((n & 1) == 0) /* even root */
                {
                    flint_printf("Exception (fmpz_root). Unable to take %wd-th root of negative value.\n", n);
                    flint_abort();
                } else /* odd */
                    c = -c;
            }
            
            root = n_rootrem(&rem, c, n);
            fmpz_set_si(r, sgn ? -root : root);
            return rem == 0;
        }
    } else /* f is large */
    {
        __mpz_struct * mpz2 = COEFF_TO_PTR(c);
        __mpz_struct * mpz1 = _fmpz_promote(r);
            
        int exact = mpz_root(mpz1, mpz2, n);
        _fmpz_demote_val(r); /* root may be small */

        return exact;
    }
}
