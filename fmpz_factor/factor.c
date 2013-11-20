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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

void
fmpz_factor(fmpz_factor_t factor, const fmpz_t n)
{
    ulong exp;
    mp_limb_t p;
    __mpz_struct * xsrc;
    mp_ptr xd;
    mp_size_t xsize;
    slong found;
    slong trial_start, trial_stop;
    TMP_INIT;

    if (!COEFF_IS_MPZ(*n))
    {
        fmpz_factor_si(factor, *n);
        return;
    }

    _fmpz_factor_set_length(factor, 0);

    /* Get sign and size */
    xsrc = COEFF_TO_PTR(*n);
    if (xsrc->_mp_size < 0)
    {
        xsize = -(xsrc->_mp_size);
        factor->sign = -1;
    }
    else
    {
        xsize = xsrc->_mp_size;
        factor->sign = 1;
    }

    /* Just a single limb */
    if (xsize == 1)
    {
        _fmpz_factor_extend_factor_ui(factor, xsrc->_mp_d[0]);
        return;
    }

    /* Create a temporary copy to be mutated */
    TMP_START;
    xd = TMP_ALLOC(xsize * sizeof(mp_limb_t));
    flint_mpn_copyi(xd, xsrc->_mp_d, xsize);

    /* Factor out powers of two */
    xsize = flint_mpn_remove_2exp(xd, xsize, &exp);
    if (exp != 0)
        _fmpz_factor_append_ui(factor, UWORD(2), exp);

    trial_start = 1;
    trial_stop = 1000;

    while (xsize > 1)
    {
        found = flint_mpn_factor_trial(xd, xsize, trial_start, trial_stop);

        if (found)
        {
            p = n_primes_arr_readonly(found+1)[found];
            exp = 1;
            xsize = flint_mpn_divexact_1(xd, xsize, p);

            /* Check if p^2 divides n */
            if (flint_mpn_divisible_1_p(xd, xsize, p))
            {
                /* TODO: when searching for squarefree numbers
                   (Moebius function, etc), we can abort here. */
                xsize = flint_mpn_divexact_1(xd, xsize, p);
                exp = 2;
            }

            /* If we're up to cubes, then maybe there are higher powers */
            if (exp == 2 && flint_mpn_divisible_1_p(xd, xsize, p))
            {
                xsize = flint_mpn_divexact_1(xd, xsize, p);
                xsize = flint_mpn_remove_power_ascending(xd, xsize, &p, 1, &exp);
                exp += 3;
            }

            _fmpz_factor_append_ui(factor, p, exp);
            /* flint_printf("added %wu %wu\n", p, exp); */

            /* Continue using only trial division whilst it is successful.
               This allows quickly factoring huge highly composite numbers
               such as factorials, which can arise in some applications. */
            trial_start = found + 1;
            trial_stop = trial_start + 1000;
            continue;
        }
        else
        {
            /* Insert primality test, perfect power test, other factoring
               algorithms here... */
            trial_start = trial_stop;
            trial_stop = trial_start + 1000;
        }
    }

    /* Any single-limb factor left? */
    if (xd[0] != 1)
        _fmpz_factor_extend_factor_ui(factor, xd[0]);

    TMP_END;
    return;
}

