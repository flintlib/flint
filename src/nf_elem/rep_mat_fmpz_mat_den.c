/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2018 Tommy Hofmann

******************************************************************************/

#include "nf_elem.h"

void nf_elem_rep_mat_fmpz_mat_den(fmpz_mat_t res, fmpz_t den, const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_set(fmpz_mat_entry(res, 0, 0), LNF_ELEM_NUMREF(a));
        fmpz_set(den, LNF_ELEM_DENREF(a));
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        nf_elem_t t;
        const fmpz * const anum = QNF_ELEM_NUMREF(a);
        const fmpz * const aden = QNF_ELEM_DENREF(a);
        fmpz * const tnum = QNF_ELEM_NUMREF(t);
        fmpz * const tden = QNF_ELEM_DENREF(t);

        nf_elem_init(t, nf);
        nf_elem_mul_gen(t, a, nf);

        if (fmpz_equal(tden, aden))
        {
            fmpz_set(fmpz_mat_entry(res, 0, 0), anum);
            fmpz_set(fmpz_mat_entry(res, 0, 1), anum + 1);
            fmpz_set(fmpz_mat_entry(res, 1, 0), tnum);
            fmpz_set(fmpz_mat_entry(res, 1, 1), tnum + 1);

            fmpz_set(den, tden);
        }
        else
        {
            fmpz_lcm(den, tden, aden);
            fmpz_divexact(fmpz_mat_entry(res, 0, 0), den, aden);
            fmpz_mul(fmpz_mat_entry(res, 0, 1), anum + 1, fmpz_mat_entry(res, 0, 0));
            fmpz_mul(fmpz_mat_entry(res, 0, 0), anum, fmpz_mat_entry(res, 0, 0));

            fmpz_divexact(fmpz_mat_entry(res, 1, 0), den, tden);
            fmpz_mul(fmpz_mat_entry(res, 1, 1), tnum + 1, fmpz_mat_entry(res, 1, 0));
            fmpz_mul(fmpz_mat_entry(res, 1, 0), tnum, fmpz_mat_entry(res, 1, 0));
        }
        nf_elem_clear(t, nf);
    }
    else
    {
        slong i, j;
        nf_elem_t t;
        slong d = fmpq_poly_degree(nf->pol);

        nf_elem_init(t, nf);
        nf_elem_set(t, a, nf);

        if (NF_ELEM(a)->length == 0)
        {
            fmpz_mat_zero(res);
            fmpz_one(den);
        }
        else if (NF_ELEM(a)->length == 1)
        {
            fmpz_mat_zero(res);
            for (i = 0; i <= d - 1; i++)
            {
              fmpz_set(fmpz_mat_entry(res, i, i), fmpq_poly_numref(NF_ELEM(a)));
            }
            fmpz_set(den, fmpq_poly_denref(NF_ELEM(a)));
        }
        else
        {
            /* Special case if defining polynomial is monic and integral and the element also has trivial denominator */
            if (nf->flag & NF_MONIC && fmpz_is_one(fmpq_poly_denref(nf->pol)) && fmpz_is_one(fmpq_poly_denref(NF_ELEM(a))))
            {
                fmpz_one(den);

                for (i = 0; i <= NF_ELEM(a)->length - 1; i++)
                    fmpz_set(fmpz_mat_entry(res, 0, i), fmpq_poly_numref(NF_ELEM(a)) + i);

                for (i = NF_ELEM(a)->length; i <= d - 1; i++)
                    fmpz_zero(fmpz_mat_entry(res, 0, i));

                for (j = 1; j <= d - NF_ELEM(a)->length; j++)
                {
                    nf_elem_mul_gen(t, t, nf);
                    for (i = 0; i < j; i++)
                        fmpz_zero(fmpz_mat_entry(res, j, i));

                    for (i = 0; i <= NF_ELEM(a)->length - 1; i++)
                        fmpz_set(fmpz_mat_entry(res, j, j + i), fmpq_poly_numref(NF_ELEM(a)) + i);

                    for (i = j + NF_ELEM(a)->length; i <= d - 1; i++)
                        fmpz_zero(fmpz_mat_entry(res, j, i));
                }

                for (j = d - NF_ELEM(a)->length + 1; j <= d - 1; j++)
                {
                    nf_elem_mul_gen(t, t, nf);
                    for (i = 0; i <= d - 1; i++)
                        fmpz_set(fmpz_mat_entry(res, j, i), fmpq_poly_numref(NF_ELEM(t)) + i);
                }
            }
            else
            {
                /* Now the general case. For 0 <= j < d - 2 we store the
                 * denominator for row j at res[d - 1, j]. At the end we
                 * divide the lcm of all of them by the corresponding
                 * denominator of the row to get the correct multiplier for
                 * row.
                 */

                for (i = 0; i <= NF_ELEM(a)->length - 1; i++)
                    fmpz_set(fmpz_mat_entry(res, 0, i), fmpq_poly_numref(NF_ELEM(a)) + i);

                for (i = NF_ELEM(a)->length; i <= d - 1; i++)
                    fmpz_zero(fmpz_mat_entry(res, 0, i));

                fmpz_set(fmpz_mat_entry(res, d - 1, 0), fmpq_poly_denref(NF_ELEM(a)));

                for (j = 1; j <= d - NF_ELEM(a)->length; j++)
                {
                    nf_elem_mul_gen(t, t, nf);
                    for (i = 0; i < j; i++)
                        fmpz_zero(fmpz_mat_entry(res, j, i));

                    for (i = 0; i <= NF_ELEM(a)->length - 1; i++)
                        fmpz_set(fmpz_mat_entry(res, j, j + i), fmpq_poly_numref(NF_ELEM(a)) + i);

                    for (i = j + NF_ELEM(a)->length; i <= d - 1; i++)
                        fmpz_zero(fmpz_mat_entry(res, j, i));

                    fmpz_set(fmpz_mat_entry(res, d - 1, j), fmpq_poly_denref(NF_ELEM(a)));
                }

                for (j = d - NF_ELEM(a)->length + 1; j <= d - 2; j++)
                {
                    nf_elem_mul_gen(t, t, nf);
                    for (i = 0; i <= d - 1; i++)
                        fmpz_set(fmpz_mat_entry(res, j, i), fmpq_poly_numref(NF_ELEM(t)) + i);

                    fmpz_set(fmpz_mat_entry(res, d - 1, j), fmpq_poly_denref(NF_ELEM(t)));

                }

                nf_elem_mul_gen(t, t, nf);
                /* Now compute the correct denominator */

                fmpz_set(fmpz_mat_entry(res, d - 1, d - 1), fmpq_poly_denref(NF_ELEM(t)));

                fmpz_set(den, fmpq_poly_denref(NF_ELEM(t)));

                for (j = 0; j <= d - 2; j++)
                    fmpz_lcm(den, den, fmpz_mat_entry(res, d - 1, j));

                for (j = 0; j <= d - 2; j++)
                {
                    if (!fmpz_equal(den, fmpz_mat_entry(res, d - 1, j)))
                    {
                        fmpz_divexact(fmpz_mat_entry(res, d - 1, j), den, fmpz_mat_entry(res, d - 1, j));

                        for (i = 0; i <= d - 1; i++)
                            fmpz_mul(fmpz_mat_entry(res, j, i), fmpz_mat_entry(res, j, i), fmpz_mat_entry(res, d - 1, j));
                    }
                }

                if (fmpz_equal(den, fmpz_mat_entry(res, d - 1, d - 1)))
                {
                    for (i = 0; i < d; i++)
                        fmpz_set(fmpz_mat_entry(res, d - 1, i), fmpq_poly_numref(NF_ELEM(t)) + i);
                }
                else
                {
                    fmpz_divexact(fmpz_mat_entry(res, d - 1, d - 1), den, fmpq_poly_denref(NF_ELEM(t)));
                    for (i = 0; i < d; i++)
                        fmpz_mul(fmpz_mat_entry(res, d - 1, i), fmpq_poly_numref(NF_ELEM(t)) + i, fmpz_mat_entry(res, d - 1, d - 1));
                }
            }
        }
        nf_elem_clear(t, nf);
    }
}
