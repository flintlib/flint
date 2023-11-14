/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "nf_elem.h"

void
nf_elem_scalar_div_fmpq(nf_elem_t a, const nf_elem_t b,
                        const fmpq_t c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);
        const fmpz *const den2 = LNF_ELEM_DENREF(b);
        const fmpz *const num2 = LNF_ELEM_NUMREF(b);
        _fmpq_div(num, den, num2, den2, fmpq_numref(c), fmpq_denref(c));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);
        const fmpz *const den2 = QNF_ELEM_DENREF(b);
        const fmpz *const num2 = QNF_ELEM_NUMREF(b);

        _fmpq_poly_scalar_div_fmpq(num, den,
                                   num2, den2, 2, fmpq_numref(c),
                                   fmpq_denref(c));
    }
    else
    {
        fmpq_poly_scalar_div_fmpq(NF_ELEM(a), NF_ELEM(b), c);
    }
}

void
nf_elem_scalar_div_fmpz(nf_elem_t a, const nf_elem_t b,
                        const fmpz_t c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);
        const fmpz *const den2 = LNF_ELEM_DENREF(b);
        const fmpz *const num2 = LNF_ELEM_NUMREF(b);

        fmpz_mul(den, den2, c);
        fmpz_set(num, num2);
        _fmpq_canonicalise(num, den);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);
        const fmpz *const den2 = QNF_ELEM_DENREF(b);
        const fmpz *const num2 = QNF_ELEM_NUMREF(b);

        fmpz_mul(den, den2, c);
        _fmpz_vec_set(num, num2, 2);
        _fmpq_poly_canonicalise(num, den, 2);
    }
    else
    {
        fmpq_poly_scalar_div_fmpz(NF_ELEM(a), NF_ELEM(b), c);
    }
}

void
nf_elem_scalar_div_si(nf_elem_t a, const nf_elem_t b, slong c, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz *den = LNF_ELEM_DENREF(a);
        fmpz *num = LNF_ELEM_NUMREF(a);
        const fmpz *const den2 = LNF_ELEM_DENREF(b);
        const fmpz *const num2 = LNF_ELEM_NUMREF(b);

        fmpz_mul_si(den, den2, c);
        fmpz_set(num, num2);
        _fmpq_canonicalise(num, den);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        fmpz *den = QNF_ELEM_DENREF(a);
        fmpz *num = QNF_ELEM_NUMREF(a);
        const fmpz *const den2 = QNF_ELEM_DENREF(b);
        const fmpz *const num2 = QNF_ELEM_NUMREF(b);

        fmpz_mul_si(den, den2, c);
        _fmpz_vec_set(num, num2, 2);
        _fmpq_poly_canonicalise(num, den, 2);

    }
    else
    {
        fmpq_poly_scalar_div_si(NF_ELEM(a), NF_ELEM(b), c);
    }
}
