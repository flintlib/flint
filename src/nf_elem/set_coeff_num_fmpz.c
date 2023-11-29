/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nf_elem.h"

void _nf_elem_set_coeff_num_fmpz(nf_elem_t a, slong i, const fmpz_t b, const nf_t nf)
{
    if (i > 2*(fmpq_poly_degree(nf->pol) - 1))
    {
        flint_throw(FLINT_ERROR, "(%s): Degree out of range\n", __func__);
    }

    if (nf->flag & NF_LINEAR)
    {
        fmpz_set(LNF_ELEM_NUMREF(a), b);
        nf_elem_canonicalise(a, nf);
    } else if (nf->flag & NF_QUADRATIC)
    {
        fmpz_set(QNF_ELEM_NUMREF(a) + i, b);
        nf_elem_canonicalise(a, nf);
    } else
    {
        slong len = NF_ELEM(a)->length;
        const int replace = (i < len && !fmpz_is_zero(NF_ELEM(a)->coeffs + i));

        if (!replace && fmpz_is_zero(b))
            return;

        if (i + 1 > len)
        {
            fmpq_poly_fit_length(NF_ELEM(a), i + 1);
            _fmpq_poly_set_length(NF_ELEM(a), i + 1);
            flint_mpn_zero((mp_ptr) NF_ELEM(a)->coeffs + len, (i + 1) - len);
        }

        if (*NF_ELEM(a)->den == WORD(1))
        {
            fmpz_set(NF_ELEM(a)->coeffs + i, b);
            if (replace)
                _fmpq_poly_normalise(NF_ELEM(a));
        }
        else
        {
            fmpz_set(NF_ELEM(a)->coeffs + i, b);
            if (replace)
                fmpq_poly_canonicalise(NF_ELEM(a));
        }
    }
}
