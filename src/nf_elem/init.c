/*
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013, 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nf_elem.h"

void nf_elem_init(nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_init(LNF_ELEM_NUMREF(a));
        fmpz_init(LNF_ELEM_DENREF(a));
        fmpz_one(LNF_ELEM_DENREF(a));
    } else if (nf->flag & NF_QUADRATIC)
    {
        fmpz * const anum = QNF_ELEM_NUMREF(a);
        fmpz * const aden = QNF_ELEM_DENREF(a);

        fmpz_init(anum);
        fmpz_init(anum + 1);
        fmpz_init(anum + 2);
        fmpz_init(aden);
        fmpz_one(aden);
    }
    else
    {
        fmpq_poly_init2(NF_ELEM(a), FLINT_MAX(2*nf->pol->length - 3, 0));
    }
}
