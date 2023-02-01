/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 William Hart

******************************************************************************/

#include "nf_elem.h"

void nf_elem_get_fmpq_poly(fmpq_poly_t pol, const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpq_poly_set_fmpz(pol, LNF_ELEM_NUMREF(a));
        fmpz_set(fmpq_poly_denref(pol), LNF_ELEM_DENREF(a));
    } else if (nf->flag & NF_QUADRATIC)
    {
        const fmpz * const anum = QNF_ELEM_NUMREF(a);

        fmpq_poly_fit_length(pol, 2);
        _fmpq_poly_set_length(pol, 2);
        _fmpz_vec_set(pol->coeffs, anum, 2);
        _fmpq_poly_normalise(pol);
        fmpz_set(pol->den, QNF_ELEM_DENREF(a));
    } else
    {
        fmpq_poly_set(pol, NF_ELEM(a));
    }
}
