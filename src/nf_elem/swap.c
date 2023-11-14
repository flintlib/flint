/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nf_elem.h"

void
nf_elem_swap(nf_elem_t a, nf_elem_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_swap(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b));
        fmpz_swap(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *const anum = QNF_ELEM_NUMREF(a);
        fmpz *const bnum = QNF_ELEM_NUMREF(b);

        fmpz_swap(anum, bnum);
        fmpz_swap(anum + 1, bnum + 1);
        fmpz_swap(anum + 2, bnum + 2);
        fmpz_swap(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));
    }
    else
        fmpq_poly_swap(NF_ELEM(a), NF_ELEM(b));
}
