/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "fmpz.h"
#include "ulong_extras.h"

slong _fmpz_remove(fmpz_t x, const fmpz_t f, double finv)
{
    fmpz y = *x;
    fmpz q = *f;

    if (!COEFF_IS_MPZ(y))  /* x is small */
    {
        if (!COEFF_IS_MPZ(q))  /* f is small */
        {
            if (y > 0)
            {
                return n_remove2_precomp((mp_limb_t *) x, q, finv);
            }
            else
            {
                ulong z = - (ulong) y;
                slong e  = n_remove2_precomp(&z, q, finv);

                if (e > 0)
                {
                    *x = - (slong) z;
                }
                return e;
            }
        }
        else  /* f is large */
        {
            return 0;
        }
    }
    else  /* x is large */
    {
        __mpz_struct *z = COEFF_TO_PTR(y);

        if (!COEFF_IS_MPZ(q))  /* f is small */
        {
            if (!flint_mpz_divisible_ui_p(z, q))
            {
                return 0;
            }
            else
            {
                flint_mpz_divexact_ui(z, z, q);

                if (!flint_mpz_divisible_ui_p(z, q))
                {
                    _fmpz_demote_val(x);
                    return 1;
                }
                else
                {
                    mpz_t r;
                    slong e;

                    flint_mpz_divexact_ui(z, z, q);
                    flint_mpz_init_set_ui(r, q);
                    e = 2 + mpz_remove(z, z, r);
                    mpz_clear(r);
                    _fmpz_demote_val(x);

                    return e;
                }
            }
        }
        else  /* f is large */
        {
            __mpz_struct *r = COEFF_TO_PTR(q);

            if (!mpz_divisible_p(z, r))
            {
                return 0;
            }
            else
            {
                slong e;

                mpz_divexact(z, z, r);
                e = 1 + mpz_remove(z, z, r);
                _fmpz_demote_val(x);
                return e;
            }
        }
    }
}

slong fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f)
{
    double finv;

    if ((fmpz_sgn(f) <= 0) || fmpz_is_one(f))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_remove). factor f <= 1.\n");
    }

    if (rop == f)
    {
        slong ans;
        fmpz_t t;

        fmpz_init(t);
        ans = fmpz_remove(t, op, f);
        fmpz_swap(rop, t);
        fmpz_clear(t);

        return ans;
    }

    finv = (!COEFF_IS_MPZ((*f))) ? n_precompute_inverse(*f) : 0;

    fmpz_set(rop, op);
    return _fmpz_remove(rop, f, finv);
}

