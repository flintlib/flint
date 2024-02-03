/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_poly_q.h"

int fmpz_poly_q_evaluate_fmpq(fmpq_t rop, const fmpz_poly_q_t f, const fmpq_t a)
{
    if (fmpz_cmp_si(fmpq_denref(a), 1))  /* a is not an integer */
    {
        fmpq_t fmpqnum, fmpqden;

        fmpq_init(fmpqden);
        fmpz_poly_evaluate_fmpq(fmpqden, f->den, a);
        if (fmpq_sgn(fmpqden) == 0)
        {
            fmpq_clear(fmpqden);
            return 1;
        }

        fmpq_init(fmpqnum);

        fmpz_poly_evaluate_fmpq(fmpqnum, f->num, a);
        fmpq_div(rop, fmpqnum, fmpqden);

        fmpq_clear(fmpqnum);
        fmpq_clear(fmpqden);
        return 0;
    }
    else  /* a is an integer */
    {
        fmpz_t num, den, a2;

        fmpz_init(num);
        fmpz_init(den);
        fmpz_init(a2);

        fmpz_set(a2, fmpq_numref(a));

        fmpz_poly_evaluate_fmpz(den, f->den, a2);
        if (fmpz_is_zero(den))
        {
            fmpz_clear(a2);
            fmpz_clear(num);
            fmpz_clear(den);
            return 1;
        }

        fmpz_poly_evaluate_fmpz(num, f->num, a2);

        fmpz_set(fmpq_numref(rop), num);
        fmpz_set(fmpq_denref(rop), den);
        _fmpq_canonicalise(fmpq_numref(rop), fmpq_denref(rop));

        fmpz_clear(a2);
        fmpz_clear(num);
        fmpz_clear(den);
        return 0;
    }
}
