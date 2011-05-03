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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

long _fmpz_remove(fmpz_t x, const fmpz_t f, double finv)
{
    long e;
    fmpz y = *x;
    fmpz q = *f;

    if (!COEFF_IS_MPZ(y))  /* x is small */
    {
        if (!COEFF_IS_MPZ(q))  /* f is small */
        {
            mp_limb_t quo, rem;
            int sgn = 1;

            if (y < 0)
            {
                sgn = -1;
                y = -y;
            }

            e = 0;
            rem = n_divrem2_precomp(&quo, y, q, finv);
            while (rem == 0)
            {
                y = quo;
                e++;
                rem = n_divrem2_precomp(&quo, y, q, finv);
            }

            *x = sgn * y;

            return e;
        }
        else  /* f is large */
        {
            return 0;
        }
    }
    else  /* x is big */
    {
        if (!COEFF_IS_MPZ(q))  /* f is small */
        {
            if (!mpz_divisible_ui_p(COEFF_TO_PTR(y), q))
            {
                return 0;
            }
            else
            {
                mpz_t h;

                mpz_init_set_ui(h, q);
                e = mpz_remove(COEFF_TO_PTR(y), COEFF_TO_PTR(y), h);
                _fmpz_demote_val(x);
                mpz_clear(h);

                return e;
            }
        }
        else  /* f is large */
        {
            __mpz_struct *a = COEFF_TO_PTR(y);
            __mpz_struct *b = COEFF_TO_PTR(q);

            if (!mpz_divisible_p(a, b))
            {
                return 0;
            }
            else
            {
                e = mpz_remove(a, a, b);
                _fmpz_demote_val(x);
                return e;
            }
        }
    }
}

long fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f)
{
    double finv;
    long ans;

    if ((fmpz_sgn(f) <= 0) || fmpz_is_one(f))
    {
        printf("Exception:  factor f <= 1 in fmpz_remove\n");
        abort();
    }

    if (rop == f)
    {
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

