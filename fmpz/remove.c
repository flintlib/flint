/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "ulong_extras.h"
#include "fmpz_mini.h"
#ifdef LONGSLONG
# define flint_mpz_divisible_ui_p mpz_divisible_ui_p
# define flint_mpz_init_set_ui mpz_init_set_ui
# define flint_mpz_divexact_ui mpz_divexact_ui
#else
# include "gmpcompat.h"
#endif

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
                return n_remove2_precomp((ulong *) x, q, finv);
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
        mpz_mock_ptr z = COEFF_TO_PTR(y);

        if (!COEFF_IS_MPZ(q))  /* f is small */
        {
            if (!flint_mpz_divisible_ui_p((mpz_ptr) z, q))
            {
                return 0;
            }
            else
            {
                flint_mpz_divexact_ui((mpz_ptr) z, (mpz_ptr) z, q);

                if (!flint_mpz_divisible_ui_p((mpz_ptr) z, q))
                {
                    _fmpz_demote_val(x);
                    return 1;
                }
                else
                {
                    mpz_t r;
                    slong e;

                    flint_mpz_divexact_ui((mpz_ptr) z, (mpz_ptr) z, q);
                    flint_mpz_init_set_ui(r, q);
                    e = 2 + mpz_remove((mpz_ptr) z, (mpz_ptr) z, r);
                    mpz_clear(r);
                    _fmpz_demote_val(x);

                    return e;
                }
            }
        }
        else  /* f is large */
        {
            mpz_mock_ptr r = COEFF_TO_PTR(q);

            if (!mpz_divisible_p((mpz_ptr) z, (mpz_ptr) r))
            {
                return 0;
            }
            else
            {
                slong e;

                mpz_divexact((mpz_ptr) z, (mpz_ptr) z, (mpz_ptr) r);
                e = 1 + mpz_remove((mpz_ptr) z, (mpz_ptr) z, (mpz_ptr) r);
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
        flint_throw(FLINT_ERROR, "factor f <= 1 in fmpz_remove\n");

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
