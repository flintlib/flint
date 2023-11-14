/*
    Copyright (C) 2018 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"
#include "nf_elem.h"

void nf_elem_rep_mat(fmpq_mat_t res, const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_set(fmpq_mat_entry_num(res, 0, 0), LNF_ELEM_NUMREF(a));
        fmpz_set(fmpq_mat_entry_den(res, 0, 0), LNF_ELEM_DENREF(a));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        nf_elem_t t;
        const fmpz * const anum = QNF_ELEM_NUMREF(a);
        const fmpz * const aden = QNF_ELEM_DENREF(a);
        fmpz * const tnum = QNF_ELEM_NUMREF(t);
        fmpz * const tden = QNF_ELEM_DENREF(t);

        nf_elem_init(t, nf);

        fmpz_set(fmpq_mat_entry_num(res, 0, 0), anum);
        fmpz_set(fmpq_mat_entry_den(res, 0, 0), aden);
        fmpq_canonicalise(fmpq_mat_entry(res, 0, 0));
        fmpz_set(fmpq_mat_entry_num(res, 0, 1), anum + 1);
        fmpz_set(fmpq_mat_entry_den(res, 0, 1), aden);
        fmpq_canonicalise(fmpq_mat_entry(res, 0, 1));

        nf_elem_mul_gen(t, a, nf);

        fmpz_set(fmpq_mat_entry_num(res, 1, 0), tnum);
        fmpz_set(fmpq_mat_entry_den(res, 1, 0), tden);
        fmpq_canonicalise(fmpq_mat_entry(res, 1, 0));
        fmpz_set(fmpq_mat_entry_num(res, 1, 1), tnum + 1);
        fmpz_set(fmpq_mat_entry_den(res, 1, 1), tden);
        fmpq_canonicalise(fmpq_mat_entry(res, 1, 1));

        nf_elem_clear(t, nf);
    }
    else
    {
        nf_elem_t t;
        slong i, j;
        slong d = fmpq_poly_degree(nf->pol);

        nf_elem_init(t, nf);
        nf_elem_set(t, a, nf);

        if (NF_ELEM(a)->length == 0)
        {
            fmpq_mat_zero(res);
            return;
        }

        for (i = 0; i <= NF_ELEM(a)->length - 1; i ++)
        {
            fmpz_set(fmpq_mat_entry_num(res, 0, i), fmpq_poly_numref(NF_ELEM(a)) + i);
            fmpz_set(fmpq_mat_entry_den(res, 0, i), fmpq_poly_denref(NF_ELEM(a)));
            fmpq_canonicalise(fmpq_mat_entry(res, 0, i));
        }

        for (i = NF_ELEM(a)->length; i <= d - 1; i++)
            fmpq_zero(fmpq_mat_entry(res, 0, i));

        for (j = 1; j <= d - NF_ELEM(a)->length; j++)
        {
            nf_elem_mul_gen(t, t, nf);
            for (i = 0; i < j; i++)
                fmpq_zero(fmpq_mat_entry(res, j, i));

            for (i = 0; i <= NF_ELEM(a)->length - 1; i++)
            {
                fmpz_set(fmpq_mat_entry_num(res, j, j + i), fmpq_poly_numref(NF_ELEM(a)) + i);
                fmpz_set(fmpq_mat_entry_den(res, j, j + i), fmpq_poly_denref(NF_ELEM(a)));
                fmpq_canonicalise(fmpq_mat_entry(res, j, j + i));
            }

            for (i = j + NF_ELEM(a)->length; i <= d - 1; i++)
                fmpq_zero(fmpq_mat_entry(res, j, i));

        }

        for (j = d - NF_ELEM(a)->length + 1; j <= d - 1; j++)
        {
            nf_elem_mul_gen(t, t, nf);
            for (i = 0; i <= d - 1; i++)
            {
                fmpz_set(fmpq_mat_entry_num(res, j, i), fmpq_poly_numref(NF_ELEM(t)) + i);
                fmpz_set(fmpq_mat_entry_den(res, j, i), fmpq_poly_denref(NF_ELEM(t)));
                fmpq_canonicalise(fmpq_mat_entry(res, j, i));
            }
        }

        nf_elem_clear(t, nf);
    }
}
