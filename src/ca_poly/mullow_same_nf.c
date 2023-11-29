/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

static const fmpz * _nf_denref(const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return LNF_ELEM_DENREF(a);
    else if (nf->flag & NF_QUADRATIC)
        return QNF_ELEM_DENREF(a);
    else
        return NF_ELEM_DENREF(a);
}

/* writes used coefficients; does not write padding zeros */
static void
_nf_elem_get_fmpz_poly_lcm(fmpz * pol, fmpz_t t, const nf_elem_t a, const fmpz_t lcm, const nf_t nf)
{
    fmpz_divexact(t, lcm, _nf_denref(a, nf));

    if (nf->flag & NF_LINEAR)
        fmpz_mul(pol, t, LNF_ELEM_NUMREF(a));
    else if (nf->flag & NF_QUADRATIC)
        _fmpz_vec_scalar_mul_fmpz(pol, QNF_ELEM_NUMREF(a), 2, t);
    else
        _fmpz_vec_scalar_mul_fmpz(pol, NF_ELEM(a)->coeffs, NF_ELEM(a)->length, t);
}

static int
get_lcm(fmpz_t Aden, ca_srcptr A, slong Alen, ca_field_t K, slong bits_limit, ca_ctx_t ctx)
{
    slong i;

    fmpz_one(Aden);

    for (i = 0; i < Alen; i++)
    {
        if (CA_IS_QQ(A + i, ctx))
            fmpz_lcm(Aden, Aden, CA_FMPQ_DENREF(A + i));
        else
            fmpz_lcm(Aden, Aden, _nf_denref(CA_NF_ELEM(A + i), CA_FIELD_NF(K)));

        if (fmpz_bits(Aden) > bits_limit)
            return 0;
    }

    return 1;
}

void
_ca_set_nf_fmpz_poly_den(ca_t res, const fmpz_poly_t poly, const fmpz_t den, ca_field_t K, ca_ctx_t ctx);

void
_ca_poly_mullow_same_nf(ca_ptr C, ca_srcptr A, slong Alen, ca_srcptr B, slong Blen, slong len, ca_field_t K, ca_ctx_t ctx)
{
    fmpz_t Aden, Bden;
    fmpz_t den, t;
    slong deg, m;
    fmpz *ZA, *ZB, *ZC;
    int squaring;
    slong i;

    Alen = FLINT_MIN(Alen, len);
    Blen = FLINT_MIN(Blen, len);

    if (!CA_FIELD_IS_NF(K))
    {
        flint_throw(FLINT_ERROR, "_ca_poly_mullow_same_nf: expected a number field\n");
    }

    squaring = (A == B) && (Alen == Blen);

    fmpz_init(Aden);
    fmpz_init(Bden);

    if (!get_lcm(Aden, A, Alen, K, WORD_MAX, ctx) || (!squaring && !get_lcm(Bden, B, Blen, K, WORD_MAX, ctx)))
    {
        flint_throw(FLINT_ERROR, "%s\n", __func__);
        /*
        fmpz_clear(Aden);
        fmpz_clear(Bden);
        _ca_poly_mullow_classical(C, A, Alen, B, Blen, len, ctx);
        return;
        */
    }

    fmpz_init(den);
    fmpz_init(t);

    /* Todo: could inspect elements to get degree bounds */
    deg = CA_FIELD_NF(K)->pol->length - 1;

    /* Packing length */
    m = 2 * deg - 1;

    ZA = _fmpz_vec_init(Alen * m);
    if (squaring)
        ZB = ZA;
    else
        ZB = _fmpz_vec_init(Blen * m);
    ZC = _fmpz_vec_init(len * m);

    for (i = 0; i < Alen; i++)
    {
        if (CA_IS_QQ(A + i, ctx))
        {
            fmpz_divexact(t, Aden, CA_FMPQ_DENREF(A + i));
            fmpz_mul(ZA + i * m + 0, t, CA_FMPQ_NUMREF(A + i));
        }
        else
        {
            _nf_elem_get_fmpz_poly_lcm(ZA + i * m, t, CA_NF_ELEM(A + i), Aden, CA_FIELD_NF(K));
        }
    }

    if (!squaring)
    {
        for (i = 0; i < Blen; i++)
        {
            if (CA_IS_QQ(B + i, ctx))
            {
                fmpz_divexact(t, Bden, CA_FMPQ_DENREF(B + i));
                fmpz_mul(ZB + i * m + 0, t, CA_FMPQ_NUMREF(B + i));
            }
            else
            {
                _nf_elem_get_fmpz_poly_lcm(ZB + i * m, t, CA_NF_ELEM(B + i), Bden, CA_FIELD_NF(K));
            }
        }
    }

    if (squaring)
    {
        _fmpz_poly_sqrlow(ZC, ZA, Alen * m, len * m);
        fmpz_mul(den, Aden, Aden);
    }
    else
    {
        if (Alen >= Blen)
            _fmpz_poly_mullow(ZC, ZA, Alen * m, ZB, Blen * m, len * m);
        else
            _fmpz_poly_mullow(ZC, ZB, Blen * m, ZA, Alen * m, len * m);
        fmpz_mul(den, Aden, Bden);
    }

    for (i = 0; i < len; i++)
    {
        fmpz_poly_t poly;
        poly->coeffs = ZC + i * m;
        poly->length = m;
        while (poly->length > 0 && fmpz_is_zero(poly->coeffs + poly->length - 1))
            poly->length--;
        _ca_set_nf_fmpz_poly_den(C + i, poly, den, K, ctx);
    }

    _fmpz_vec_clear(ZA, Alen * m);
    if (!squaring)
        _fmpz_vec_clear(ZB, Blen * m);
    _fmpz_vec_clear(ZC, len * m);

    fmpz_clear(Aden);
    fmpz_clear(Bden);
    fmpz_clear(den);
    fmpz_clear(t);
}
