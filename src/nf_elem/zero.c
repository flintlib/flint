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
nf_elem_zero(nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_zero(LNF_ELEM_NUMREF(a));
        fmpz_one(LNF_ELEM_DENREF(a));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz * const anum = QNF_ELEM_NUMREF(a);
        fmpz_zero(anum);
        fmpz_zero(anum + 1);
        fmpz_one(QNF_ELEM_DENREF(a));
    }
    else
        fmpq_poly_zero(NF_ELEM(a));
}
