/*=============================================================================

    This file is part of ANTIC.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2019 Vincent Delecroix

******************************************************************************/

#include "nf.h"

void nf_init_randtest(nf_t nf, flint_rand_t state,
        slong len,
        mp_bitcnt_t bits_in)
{
    fmpq_poly_t pol;
    fmpz_poly_t q;

    if (len < 2 || bits_in < 1)
    {
        fprintf(stderr, "[nf_init_randtest] len must be >= 2 and bits_in >= 1\n");
        abort();
    }

    if (len <= 2 || n_randint(state, 10) == 0)
        len = 2; /* linear */
    else if (len <= 3 || n_randint(state, 8) == 0)
        len = 3; /* quadratic */
    else
        len = 3 + n_randint(state, len-2);

    fmpz_poly_init(q);
    fmpq_poly_init(pol);

    if (len == 3 && (n_randint(state, 8) == 0))
    {
        fmpq_poly_set_coeff_si(pol, 0, 1);
        fmpq_poly_set_coeff_si(pol, 2, 1);
    }
    else
    {
        do {
            fmpz_poly_randtest(q,
                    state,
                    len,
                    1 + n_randint(state, bits_in));
        } while (fmpz_poly_degree(q) < 1 || fmpz_is_zero(q->coeffs));

        fmpq_poly_set_fmpz_poly(pol, q);

        if (n_randint(state, 5) == 0)
            fmpz_one(pol->coeffs + pol->length - 1); /* monic */
        else
            fmpz_randtest_not_zero(fmpq_poly_denref(pol), state, bits_in);
        fmpq_poly_canonicalise(pol);
    }

    nf_init(nf, pol);
    fmpq_poly_clear(pol);
    fmpz_poly_clear(q);
}
