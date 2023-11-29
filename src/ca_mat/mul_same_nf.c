/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_mat.h"
#include "ca_mat.h"

static const fmpz * _nf_denref(const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return LNF_ELEM_DENREF(a);
    else if (nf->flag & NF_QUADRATIC)
        return QNF_ELEM_DENREF(a);
    else
        return NF_ELEM_DENREF(a);
}

static void
_nf_elem_get_fmpz_poly_lcm(fmpz_poly_t pol, fmpz_t t, const nf_elem_t a, const fmpz_t lcm, const nf_t nf)
{
    fmpz_divexact(t, lcm, _nf_denref(a, nf));

    if (nf->flag & NF_LINEAR)
    {
        fmpz_mul(t, t, LNF_ELEM_NUMREF(a));
        fmpz_poly_set_fmpz(pol, t);
    }
    else if (nf->flag & NF_QUADRATIC)
    {
        const fmpz * const anum = QNF_ELEM_NUMREF(a);

        fmpz_poly_fit_length(pol, 2);
        _fmpz_poly_set_length(pol, 2);
        _fmpz_vec_scalar_mul_fmpz(pol->coeffs, anum, 2, t);
        _fmpz_poly_normalise(pol);
    }
    else
    {
        fmpz_poly_fit_length(pol, NF_ELEM(a)->length);
        _fmpz_poly_set_length(pol, NF_ELEM(a)->length);
        _fmpz_vec_scalar_mul_fmpz(pol->coeffs, NF_ELEM(a)->coeffs, NF_ELEM(a)->length, t);
    }
}

static int
get_lcm_rowwise(fmpz * Aden, const ca_mat_t A, ca_field_t K, slong bits_limit, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(A); i++)
    {
        fmpz_one(Aden + i);

        for (j = 0; j < ca_mat_ncols(A); j++)
        {
            if (CA_IS_QQ(ca_mat_entry(A, i, j), ctx))
                fmpz_lcm(Aden + i, Aden + i, CA_FMPQ_DENREF(ca_mat_entry(A, i, j)));
            else
                fmpz_lcm(Aden + i, Aden + i, _nf_denref(CA_NF_ELEM(ca_mat_entry(A, i, j)), CA_FIELD_NF(K)));

            if (fmpz_bits(Aden + i) > bits_limit)
                return 0;
        }
    }

    return 1;
}

static int
get_lcm_colwise(fmpz * Aden, const ca_mat_t A, ca_field_t K, slong bits_limit, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_ncols(A); i++)
    {
        fmpz_one(Aden + i);

        for (j = 0; j < ca_mat_nrows(A); j++)
        {
            if (CA_IS_QQ(ca_mat_entry(A, j, i), ctx))
                fmpz_lcm(Aden + i, Aden + i, CA_FMPQ_DENREF(ca_mat_entry(A, j, i)));
            else
                fmpz_lcm(Aden + i, Aden + i, _nf_denref(CA_NF_ELEM(ca_mat_entry(A, j, i)), CA_FIELD_NF(K)));

            if (fmpz_bits(Aden + i) > bits_limit)
                return 0;
        }
    }

    return 1;
}

static void
get_mat_rowwise(fmpz_poly_mat_t Aclear, const ca_mat_t A, const fmpz * Aden, ca_field_t K, ca_ctx_t ctx)
{
    slong i, j;
    fmpz_t t;

    fmpz_init(t);

    for (i = 0; i < ca_mat_nrows(A); i++)
    {
        for (j = 0; j < ca_mat_ncols(A); j++)
        {
            if (CA_IS_QQ(ca_mat_entry(A, i, j), ctx))
            {
                fmpz_divexact(t, Aden + i, CA_FMPQ_DENREF(ca_mat_entry(A, i, j)));
                fmpz_mul(t, t, CA_FMPQ_NUMREF(ca_mat_entry(A, i, j)));
                fmpz_poly_set_fmpz(fmpz_poly_mat_entry(Aclear, i, j), t);
            }
            else
            {
                _nf_elem_get_fmpz_poly_lcm(fmpz_poly_mat_entry(Aclear, i, j), t, CA_NF_ELEM(ca_mat_entry(A, i, j)), Aden + i, CA_FIELD_NF(K));
            }
        }
    }

    fmpz_clear(t);
}

static void
get_mat_colwise(fmpz_poly_mat_t Aclear, const ca_mat_t A, const fmpz * Aden, ca_field_t K, ca_ctx_t ctx)
{
    slong i, j;
    fmpz_t t;

    fmpz_init(t);

    for (i = 0; i < ca_mat_ncols(A); i++)
    {
        for (j = 0; j < ca_mat_nrows(A); j++)
        {
            if (CA_IS_QQ(ca_mat_entry(A, j, i), ctx))
            {
                fmpz_divexact(t, Aden + i, CA_FMPQ_DENREF(ca_mat_entry(A, j, i)));
                fmpz_mul(t, t, CA_FMPQ_NUMREF(ca_mat_entry(A, j, i)));
                fmpz_poly_set_fmpz(fmpz_poly_mat_entry(Aclear, j, i), t);
            }
            else
            {
                _nf_elem_get_fmpz_poly_lcm(fmpz_poly_mat_entry(Aclear, j, i), t, CA_NF_ELEM(ca_mat_entry(A, j, i)), Aden + i, CA_FIELD_NF(K));
            }
        }
    }

    fmpz_clear(t);
}

void
_ca_set_nf_fmpz_poly_den(ca_t res, const fmpz_poly_t poly, const fmpz_t den, ca_field_t K, ca_ctx_t ctx)
{
    if (poly->length == 0)
    {
        ca_zero(res, ctx);
    }
    else if (poly->length == 1)
    {
        _ca_make_fmpq(res, ctx);

        fmpz_gcd(CA_FMPQ_DENREF(res), poly->coeffs, den);
        fmpz_divexact(CA_FMPQ_NUMREF(res), poly->coeffs, CA_FMPQ_DENREF(res));
        fmpz_divexact(CA_FMPQ_DENREF(res), den, CA_FMPQ_DENREF(res));
    }
    else
    {
        fmpq_poly_t T;

        T->coeffs = poly->coeffs;
        T->length = poly->length;
        T->alloc = poly->alloc;
        T->den[0] = den[0];
        _ca_make_field_element(res, K, ctx);

        /* work around antic bug */
        if (CA_FIELD_NF(K)->flag & NF_QUADRATIC)
        {
            fmpz_set(QNF_ELEM_NUMREF(CA_NF_ELEM(res)), T->coeffs);
            fmpz_set(QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1, T->coeffs + 1);
            if (T->length == 3)
                fmpz_set(QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 2, T->coeffs + 2);
            fmpz_set(QNF_ELEM_DENREF(CA_NF_ELEM(res)), den);
            /* todo: canonicalise before reduction? */
        }
        else
        {
            nf_elem_set_fmpq_poly(CA_NF_ELEM(res), T, CA_FIELD_NF(K));
        }

        nf_elem_reduce(CA_NF_ELEM(res), CA_FIELD_NF(K));
        /* antic is currently inconsistent about canonicalising in reduce() */
        if (CA_FIELD_NF(K)->flag & NF_LINEAR)
            nf_elem_canonicalise(CA_NF_ELEM(res), CA_FIELD_NF(K));

        /* may have reduced to a rational */
        ca_condense_field(res, ctx);
    }
}

void
ca_mat_mul_same_nf(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_field_t K, ca_ctx_t ctx)
{
    fmpz_poly_mat_t ZC, ZA, ZB;
    fmpz * Aden, * Bden;
    fmpz_t den;
    slong Ar, Ac, Br, Bc;
    slong i, j;

    Ar = ca_mat_nrows(A);
    Ac = ca_mat_ncols(A);
    Br = ca_mat_nrows(B);
    Bc = ca_mat_ncols(B);

    if (Ar == 0 || Ac == 0 || Bc == 0)
    {
        ca_mat_zero(C, ctx);
        return;
    }

    if (!CA_FIELD_IS_NF(K))
    {
        flint_throw(FLINT_ERROR, "ca_mat_mul_same_nf: expected a number field\n");
    }

    Aden = _fmpz_vec_init(Ar);
    Bden = _fmpz_vec_init(Bc);

    if (!get_lcm_rowwise(Aden, A, K, 1000, ctx) || !get_lcm_colwise(Bden, B, K, 1000, ctx))
    {
        _fmpz_vec_clear(Aden, Ar);
        _fmpz_vec_clear(Bden, Bc);
        ca_mat_mul_classical(C, A, B, ctx);
        return;
    }

    fmpz_init(den);

    fmpz_poly_mat_init(ZA, Ar, Ac);
    fmpz_poly_mat_init(ZB, Br, Bc);
    fmpz_poly_mat_init(ZC, Ar, Bc);

    get_mat_rowwise(ZA, A, Aden, K, ctx);
    get_mat_colwise(ZB, B, Bden, K, ctx);

    fmpz_poly_mat_mul(ZC, ZA, ZB);

    for (i = 0; i < ca_mat_nrows(C); i++)
    {
        for (j = 0; j < ca_mat_ncols(C); j++)
        {
            fmpz_mul(den, Aden + i, Bden + j);
            _ca_set_nf_fmpz_poly_den(ca_mat_entry(C, i, j), fmpz_poly_mat_entry(ZC, i, j), den, K, ctx);
        }
    }

    fmpz_poly_mat_clear(ZA);
    fmpz_poly_mat_clear(ZB);
    fmpz_poly_mat_clear(ZC);

    _fmpz_vec_clear(Aden, Ar);
    _fmpz_vec_clear(Bden, Bc);
    fmpz_clear(den);
}
