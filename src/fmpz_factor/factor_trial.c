/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_factor.h"

int
fmpz_factor_trial(fmpz_factor_t factor, const fmpz_t n, slong num_primes)
{
    ulong exp;
    mp_limb_t p;
    mpz_t x;
    mp_ptr xd;
    mp_size_t xsize;
    slong found;
    int ret = 1;
    slong * idx;
    slong bits, i;
    const mp_limb_t * primes;

    if (num_primes > 3512 || num_primes < 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_factor_trial) Number of primes must be in 0..3512\n");
    }

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
    xsize = flint_mpn_remove_2exp(xd, xsize, &exp);
    if (exp != 0)
        _fmpz_factor_append_ui(factor, UWORD(2), exp);

    bits = fmpz_sizeinbase(n, 2) - exp;
    idx = (slong *) flint_malloc((5 + bits/4)*sizeof(slong));

    found = flint_mpn_factor_trial_tree(idx, xd, xsize, num_primes);

    if (found)
    {
        primes = n_primes_arr_readonly(3512);

        for (i = 0; i < found; i++)
        {
            p = primes[idx[i]];

            if (p == 2)
                continue;

            exp = 1;
            xsize = flint_mpn_divexact_1(xd, xsize, p);

            /* Check if p^2 divides n */
            if (flint_mpn_divisible_1_odd(xd, xsize, p))
            {
                /* TODO: when searching for squarefree numbers
                   (Moebius function, etc), we can abort here. */
                xsize = flint_mpn_divexact_1(xd, xsize, p);
                exp = 2;
            }

            /* If we're up to cubes, then maybe there are higher powers */
            if (exp == 2 && flint_mpn_divisible_1_odd(xd, xsize, p))
            {
                xsize = flint_mpn_divexact_1(xd, xsize, p);
                xsize = flint_mpn_remove_power_ascending(xd, xsize, &p, 1, &exp);
                exp += 3;
            }

            _fmpz_factor_append_ui(factor, p, exp);
        }
    }

    /* Any factor left? */
    if (xsize > 1 || xd[0] != 1)
    {
        fmpz_t cofactor;
        mpz_t mockx; /* do not free */

        fmpz_init(cofactor);
        mockx->_mp_d = xd;
        mockx->_mp_size = xsize;
        mockx->_mp_alloc = x->_mp_alloc;

        fmpz_set_mpz(cofactor, mockx);
        _fmpz_factor_append(factor, cofactor, 1);

        fmpz_clear(cofactor);

        ret = 0;
    }

    mpz_clear(x);

    flint_free(idx);

    return ret;
}
