/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void
fmpz_ndiv_qr(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    if (fmpz_is_zero(b))
    {
        flint_printf("Exception: division by zero in fmpz_ndiv_qr\n");
        flint_abort();
    }

    if (!COEFF_IS_MPZ(*a))
    {
        if (!COEFF_IS_MPZ(*b))      /* a and b are small */
        {
            slong lquo;
            slong lrem;

            _fmpz_demote(q);
            _fmpz_demote(r);

            if (FLINT_ABS(*b) == 1) /* avoid overflow in case */
            {                       /* a = 2^(FLINT_BITS - 2) */
                fmpz_set_si(q, *a * FLINT_SGN(*b));
                fmpz_set_si(r, 0);
                return;
            }

            *q = *a / *b;
            *r = *a - *b * *q;
            lquo = *q + FLINT_SGN(*a) * FLINT_SGN(*b);
            lrem = *a - *b * lquo;

            if (FLINT_ABS(lrem) < FLINT_ABS(*r))
            {
                *q = lquo;
                *r = lrem;
            }
        }
        else                        /* a is small and b is large */
        {
            fmpz_mul_si(q, a, 2);
            if (fmpz_cmpabs(q, b) > 0) /* |a / b| > 1 / 2 */
            {
                fmpz_set_si(q, FLINT_SGN(*a) * mpz_sgn(COEFF_TO_PTR(*b)));
                fmpz_set(r, a);
                fmpz_submul(r, b, q);
            }
            else
            {
                fmpz_set_ui(q, WORD(0));
                fmpz_set(r, a);
            }
        }
    }
    else
    {
        __mpz_struct *mpz_ptr, *mpz_ptr2;
        fmpz_t lquo;
        fmpz_t lrem;

        fmpz_init(lquo);
        fmpz_init(lrem);

        _fmpz_promote(q); /* must not hang on to ptr whilst promoting s */
        mpz_ptr2 = _fmpz_promote(r);
		mpz_ptr  = COEFF_TO_PTR(*q);

		if (!COEFF_IS_MPZ(*b))  /* a is large and b is small */
        {
            if (*b > 0)
            {
                flint_mpz_tdiv_qr_ui(mpz_ptr, mpz_ptr2, COEFF_TO_PTR(*a), *b);
            }
            else
            {
                flint_mpz_tdiv_qr_ui(mpz_ptr, mpz_ptr2, COEFF_TO_PTR(*a), -*b);
                mpz_neg(mpz_ptr, mpz_ptr);
            }

            _fmpz_demote_val(q);
            _fmpz_demote_val(r);

            fmpz_add_si(lquo, q, mpz_sgn(COEFF_TO_PTR(*a)) * FLINT_SGN(*b));
        }
        else                    /* a and b are large */
        {
            mpz_tdiv_qr(mpz_ptr, mpz_ptr2, COEFF_TO_PTR(*a), COEFF_TO_PTR(*b));

            _fmpz_demote_val(q);
            _fmpz_demote_val(r);

            fmpz_add_si(lquo, q, mpz_sgn(COEFF_TO_PTR(*a))
                                * mpz_sgn(COEFF_TO_PTR(*b)));
        }

        fmpz_set(lrem, a);
        fmpz_submul(lrem, b, lquo);

        if (fmpz_cmpabs(lrem, r) < 0)
        {
            fmpz_set(q, lquo);
            fmpz_set(r, lrem);
        }

        fmpz_clear(lquo);
        fmpz_clear(lrem);
    }
}
