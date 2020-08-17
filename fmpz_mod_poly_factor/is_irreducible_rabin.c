/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

void
fmpz_mod_poly_powpowmod(fmpz_mod_poly_t res, const fmpz_mod_poly_t pol,
                        const fmpz_t exp, ulong exp2, const fmpz_mod_poly_t f)
{
    fmpz_mod_poly_t pow;
    ulong i;

    fmpz_mod_poly_init(pow, &f->p);
    fmpz_mod_poly_powmod_fmpz_binexp(pow, pol, exp, f);
    fmpz_mod_poly_set(res, pow);

    if (!fmpz_mod_poly_equal(pow, pol))
        for (i = 1; i < exp2; i++)
            fmpz_mod_poly_powmod_fmpz_binexp(res, res, exp, f);

    fmpz_mod_poly_clear(pow);
}

int
fmpz_mod_poly_is_irreducible_rabin(const fmpz_mod_poly_t f)
{
    int res = 1;
    
    if (fmpz_mod_poly_length(f) > 2)
    {
        const slong n = fmpz_mod_poly_degree(f);
        fmpz_mod_poly_t a, x, x_p, finv;
        fmpz_mod_poly_frobenius_powers_2exp_t pow;
        fmpz_t pk;

        fmpz_init(pk);

        fmpz_mod_poly_init(a, &f->p);
        fmpz_mod_poly_init(x, &f->p);
        fmpz_mod_poly_init(x_p, &f->p);
        fmpz_mod_poly_init(finv, &f->p);

        fmpz_mod_poly_set_coeff_ui(x, 1, 1);

        /* Compute x^q mod f */
        fmpz_mod_poly_reverse(finv, f, f->length);

        fmpz_mod_poly_inv_series_newton(finv, finv, f->length);
                   
        fmpz_mod_poly_frobenius_powers_2exp_precomp(pow, f, finv, n);
        
        fmpz_mod_poly_frobenius_power(x_p, pow, f, n);

        if (!fmpz_mod_poly_is_zero(x_p))
           fmpz_mod_poly_make_monic(x_p, x_p);

        /* Now do the irreducibility test */
        if (!fmpz_mod_poly_equal(x_p, x))
            res = 0;
        else
        {
            n_factor_t factors;
            slong i;

            n_factor_init(&factors);
            n_factor(&factors, n, 1);

            for (i = 0; i < factors.num; i++)
            {
                fmpz_mod_poly_frobenius_power(a, pow, f, n / factors.p[i]);

                fmpz_mod_poly_sub(a, a, x);

                if (!fmpz_mod_poly_is_zero(a))
                    fmpz_mod_poly_make_monic(a, a);

                fmpz_mod_poly_gcd(a, a, f);

                if (a->length != 1)
                   res = 0;
            }
        }

        fmpz_clear(pk);
        fmpz_mod_poly_frobenius_powers_2exp_clear(pow);
        fmpz_mod_poly_clear(finv);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(x);
        fmpz_mod_poly_clear(x_p);
    }

    return res;
}
