/*
    Copyright (C) 2014 William Hart
    Copyright (C) 2015 Claus Fieker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "nf_elem.h"

void _nf_elem_norm_div(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf, const fmpz_t divisor, slong nbits)
{
    if (fmpz_is_zero(divisor))
    {
        fmpz_zero(rnum);
        fmpz_one(rden);
        return;
    }

    if (nf->flag & NF_LINEAR)
    {
        const fmpz * const anum = LNF_ELEM_NUMREF(a);
        const fmpz * const aden = LNF_ELEM_DENREF(a);

        fmpz_gcd(rden, anum, divisor);
        fmpz_divexact(rnum, anum, rden);
        fmpz_divexact(rden, divisor, rden);
        fmpz_mul(rden, rden, aden);
        if (fmpz_sgn(rden) == -1)
        {
            fmpz_neg(rden, rden);
            fmpz_neg(rnum, rnum);
        }
    } else if (nf->flag & NF_QUADRATIC)
    {
        const fmpz * const anum = QNF_ELEM_NUMREF(a);
        const fmpz * const aden = QNF_ELEM_DENREF(a);
        fmpz_t pow, one;

        slong alen = 2;
        while (alen > 0 && fmpz_is_zero(anum + alen - 1))
            alen--;

        if (alen == 0)
        {
            fmpz_zero(rnum);
            fmpz_one(rden);

            return;
        }

        fmpz_init_set_ui(one, 1);
        fmpz_init(pow);

        _fmpq_poly_resultant_div(rnum, rden,
            nf->pol->coeffs, one, 3, anum, aden, alen, divisor, nbits);

        if (!fmpz_is_one(nf->pol->coeffs + 2) && alen > 1)
        {
            fmpz_pow_ui(pow, nf->pol->coeffs + 2, alen - 1);
            _fmpq_mul(rnum, rden, rnum, rden, one, pow);

            if (fmpz_sgn(rden) < 0)
            {
                fmpz_neg(rnum, rnum);
                fmpz_neg(rden, rden);
            }
        }
        fmpz_clear(one);
        fmpz_clear(pow);
    } else /* generic nf_elem */
    {
        const fmpz * const anum = NF_ELEM_NUMREF(a);
        const fmpz * const aden = NF_ELEM_DENREF(a);
        fmpz_t pow, one;

        slong alen = NF_ELEM(a)->length;
        slong len = nf->pol->length;
        fmpz * coeffs = nf->pol->coeffs;

        if (alen == 0)
        {
            fmpz_zero(rnum);
            fmpz_one(rden);

            return;
        }

        fmpz_init_set_ui(one, 1);
        fmpz_init(pow);

        _fmpq_poly_resultant_div(rnum, rden,
            nf->pol->coeffs, one, len, anum, aden, alen, divisor, nbits);

        if (!fmpz_is_one(coeffs + len - 1) && alen > 1)
        {
            fmpz_pow_ui(pow, coeffs + len - 1, alen - 1);
            _fmpq_mul(rnum, rden, rnum, rden, one, pow);

            if (fmpz_sgn(rden) < 0)
            {
                fmpz_neg(rnum, rnum);
                fmpz_neg(rden, rden);
            }
        }

        fmpz_clear(one);
        fmpz_clear(pow);
    }
}

void nf_elem_norm_div(fmpq_t res, const nf_elem_t a, const nf_t nf, const fmpz_t divisor, slong nbits)
{
    _nf_elem_norm_div(fmpq_numref(res), fmpq_denref(res), a, nf, divisor, nbits);
}
