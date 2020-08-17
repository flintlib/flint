/*
    Copyright (C) 2009 William Hart

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

void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x)
{
	fmpz c = *g;

	if (!COEFF_IS_MPZ(c))  /* g is small */
	{
        mp_limb_t sum[2];
		if (c >= WORD(0))  /* both operands non-negative */
		{
			add_ssaaaa(sum[1], sum[0], 0, c, 0, x);
            fmpz_set_uiui(f, sum[1], sum[0]);
		}
        else  /* coeff is negative, x positive */
		{
			if (-c > x)
                fmpz_set_si(f, x + c); /* can't overflow as g is small and x smaller */
			else
                fmpz_set_ui(f, x + c);  /* won't be negative and has to be less than x */
		}
	}
    else
	{	
		__mpz_struct * mpz_ptr2 = _fmpz_promote(f);  /* g is already large */
		__mpz_struct * mpz_ptr = COEFF_TO_PTR(c);
		flint_mpz_add_ui(mpz_ptr2, mpz_ptr, x);
		_fmpz_demote_val(f);  /* cancellation may have occurred */
	}
}
