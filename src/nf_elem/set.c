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
nf_elem_set(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_set(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b));
        fmpz_set(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *const anum = QNF_ELEM_NUMREF(a);
        const fmpz *const bnum = QNF_ELEM_NUMREF(b);

        fmpz_set(anum, bnum);
        fmpz_set(anum + 1, bnum + 1);
        fmpz_set(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));
    }
    else
        fmpq_poly_set(NF_ELEM(a), NF_ELEM(b));
}

void
nf_elem_set_si(nf_elem_t a, slong c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_set_si(LNF_ELEM_NUMREF(a), c);
        fmpz_one(LNF_ELEM_DENREF(a));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *const anum = QNF_ELEM_NUMREF(a);

        fmpz_set_si(anum, c);
        fmpz_zero(anum + 1);
        fmpz_one(QNF_ELEM_DENREF(a));
    }
    else
        fmpq_poly_set_si(NF_ELEM(a), c);
}

void
nf_elem_set_ui(nf_elem_t a, ulong c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_set_ui(LNF_ELEM_NUMREF(a), c);
        fmpz_one(LNF_ELEM_DENREF(a));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *const anum = QNF_ELEM_NUMREF(a);

        fmpz_set_ui(anum, c);
        fmpz_zero(anum + 1);
        fmpz_one(QNF_ELEM_DENREF(a));
    }
    else
        fmpq_poly_set_ui(NF_ELEM(a), c);
}

void
nf_elem_set_fmpz(nf_elem_t a, const fmpz_t c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_set(LNF_ELEM_NUMREF(a), c);
        fmpz_one(LNF_ELEM_DENREF(a));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *const anum = QNF_ELEM_NUMREF(a);

        fmpz_set(anum, c);
        fmpz_zero(anum + 1);
        fmpz_one(QNF_ELEM_DENREF(a));
    }
    else
        fmpq_poly_set_fmpz(NF_ELEM(a), c);
}

void
nf_elem_set_fmpq(nf_elem_t a, const fmpq_t c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_set(LNF_ELEM_NUMREF(a), fmpq_numref(c));
        fmpz_set(LNF_ELEM_DENREF(a), fmpq_denref(c));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *const anum = QNF_ELEM_NUMREF(a);

        fmpz_set(anum, fmpq_numref(c));
        fmpz_zero(anum + 1);
        fmpz_set(QNF_ELEM_DENREF(a), fmpq_denref(c));
    }
    else
        fmpq_poly_set_fmpq(NF_ELEM(a), c);
}
