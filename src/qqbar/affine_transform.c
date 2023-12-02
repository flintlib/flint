/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

static void
_fmpq_poly_compose2(fmpz * res, const fmpz * poly1, const fmpz_t den1,
                   slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)
{
    fmpz_t den;
    slong len;
    len = (len1 - WORD(1)) * (len2 - WORD(1)) + WORD(1);

    fmpz_init(den);

    if (fmpz_is_one(den2))
    {
        _fmpz_poly_compose(res, poly1, len1, poly2, len2);
    }
    else
    {
        fmpz_t one;
        fmpz * v = _fmpz_vec_init(len1);
        fmpz_init(one);
        fmpz_one(one);

        _fmpq_poly_rescale(v, den, poly1, den1, len1, one, den2);
        _fmpz_poly_compose(res, v, len1, poly2, len2);

        fmpz_clear(one);
        _fmpz_vec_clear(v, len1);
    }

    _fmpz_vec_content(den, res, len);
    if (fmpz_sgn(res + len - 1) < 0)
        fmpz_neg(den, den);
    _fmpz_vec_scalar_divexact_fmpz(res, res, len, den);

    fmpz_clear(den);
}

void
qqbar_scalar_op(qqbar_t res, const qqbar_t x, const fmpz_t a, const fmpz_t b, const fmpz_t c)
{
    fmpz_poly_t H;
    fmpz_t one, Gden;
    fmpz G[2];
    acb_t z, t, w;
    slong d, prec;

    if (fmpz_is_zero(c))
    {
        flint_throw(FLINT_ERROR, "qqbar_scalar_op: division by zero\n");
    }

    /* Special case: set to rational */
    if (fmpz_is_zero(a))
    {
        fmpq_t t;
        fmpq_init(t);
        fmpq_set_fmpz_frac(t, b, c);
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        return;
    }

    d = qqbar_degree(x);

    /* Special case: rational arithmetic */
    if (d == 1)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpz_neg(fmpq_numref(t), QQBAR_POLY(x)->coeffs);
        fmpz_set(fmpq_denref(t), QQBAR_POLY(x)->coeffs + 1);
        if (!fmpz_is_one(a))
            fmpq_mul_fmpz(t, t, a);
        if (!fmpz_is_zero(b))
            fmpq_add_fmpz(t, t, b);
        if (!fmpz_is_one(c))
            fmpq_div_fmpz(t, t, c);
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        return;
    }

    fmpz_poly_init2(H, d + 1);
    fmpz_init(one);
    fmpz_init(G);
    fmpz_init(G + 1);
    fmpz_init(Gden);

    /* (ax+b)/c -> inverse transformation (-b+cx)/a */
    fmpz_one(one);

    if (fmpz_sgn(a) > 0)
    {
        fmpz_neg(G, b);
        fmpz_set(G + 1, c);
        fmpz_set(Gden, a);
    }
    else
    {
        fmpz_set(G, b);
        fmpz_neg(G + 1, c);
        fmpz_neg(Gden, a);
    }

    _fmpq_poly_compose2(H->coeffs, QQBAR_POLY(x)->coeffs, one, d + 1, G, Gden, 2);
    _fmpz_poly_set_length(H, d + 1);

    acb_init(z);
    acb_init(t);
    acb_init(w);

    acb_set(z, QQBAR_ENCLOSURE(x));

    for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, prec);

        if (fmpz_is_one(a))
            acb_set(w, z);
        else if (fmpz_equal_si(a, -1))
            acb_neg(w, z);
        else
            acb_mul_fmpz(w, z, a, prec);

        if (!fmpz_is_zero(b))
            acb_add_fmpz(w, w, b, prec);

        if (!fmpz_is_one(c))
        {
            if (fmpz_equal_si(c, -1))
                acb_neg(w, w);
            else
                acb_div_fmpz(w, w, c, prec);
        }

        if (_qqbar_validate_uniqueness(t, H, w, 2 * prec)
            /* && acb_rel_accuracy_bits(t) >= QQBAR_DEFAULT_PREC / 2 */)
        {
            fmpz_poly_set(QQBAR_POLY(res), H);
            acb_set(QQBAR_ENCLOSURE(res), t);
            break;
        }
    }

    acb_clear(z);
    acb_clear(t);
    acb_clear(w);

    fmpz_poly_clear(H);
    fmpz_clear(one);
    fmpz_clear(G);
    fmpz_clear(G + 1);
    fmpz_clear(Gden);
}

