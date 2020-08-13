/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

int
fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes)
{
    ulong exp;
    mp_limb_t p;
    mpz_t x;
    mp_ptr xd;
    mp_size_t xsize;
    slong found;
    slong trial_start, trial_stop;
    int ret = 1;

    if (!COEFF_IS_MPZ(*n))
    {
        fmpz_factor_si(factor, *n);

        return ret;
    }

    _fmpz_factor_set_length(factor, 0);

    /* Make an mpz_t copy whose limbs will be mutated */
    mpz_init(x);
    fmpz_get_mpz(x, n);
    if (x->_mp_size < 0)
    {
        x->_mp_size = -(x->_mp_size);
        factor->sign = -1;
    }
    else
    {
        factor->sign = 1;
    }

    xd = x->_mp_d;
    xsize = x->_mp_size;

    /* Factor out powers of two */
    if (start == 0)
    {
       xsize = flint_mpn_remove_2exp(xd, xsize, &exp);
       if (exp != 0)
           _fmpz_factor_append_ui(factor, UWORD(2), exp);
    }

    trial_start = FLINT_MAX(1, start);
    trial_stop = FLINT_MIN(start + 1000, start + num_primes);

    do
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
            trial_stop = FLINT_MIN(trial_start + 1000, start + num_primes);
            continue;
        }
        else
        {
            /* Insert primality test, perfect power test, other factoring
               algorithms here... */
            trial_start = trial_stop;
            trial_stop = FLINT_MIN(trial_start + 1000, start + num_primes);
        }
    } while ((xsize > 1 || xd[0] != 1) && trial_start != trial_stop);

    /* Any factor left? */
    if (xsize > 1 || xd[0] != 1)
        ret = 0;

    mpz_clear(x);
    return ret;
}
