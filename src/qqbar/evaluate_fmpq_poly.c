/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_mat.h"
#include "fmpq_poly.h"
#include "nf.h"
#include "nf_elem.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

void
_qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpz * poly, const fmpz_t den, slong len, const qqbar_t x)
{
    slong d = qqbar_degree(x);

    if (len == 0)
    {
        qqbar_zero(res);
    }
    else if (len == 1)
    {
        if (fmpz_is_one(den))
        {
            qqbar_set_fmpz(res, poly);
        }
        else
        {
            fmpq_t t;
            fmpq_init(t);
            fmpq_set_fmpz_frac(t, poly, den);
            qqbar_set_fmpq(res, t);
            fmpq_clear(t);
        }
    }
    else if (qqbar_is_rational(x))
    {
        fmpq_t t, u;
        fmpq_init(t);
        fmpq_init(u);
        qqbar_get_fmpq(u, x);
        _fmpq_poly_evaluate_fmpq(fmpq_numref(t), fmpq_denref(t), poly, den, len, fmpq_numref(u), fmpq_denref(u));
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        fmpq_clear(u);
    }
    else if (len == 2)
    {
        qqbar_scalar_op(res, x, poly + 1, poly, den);
    }
    else if (len > d)
    {
        fmpz * tmp;
        fmpz_t r, one;
        slong len2;

        tmp = _fmpz_vec_init(len);
        fmpz_init(r);
        fmpz_init(one);
        fmpz_one(one);

        _fmpq_poly_rem(tmp, r, poly, den, len,
            QQBAR_COEFFS(x), one, d + 1, NULL);

        len2 = d;
        while (len2 >= 1 && fmpz_is_zero(tmp + len2 - 1))
            len2--;

        _qqbar_evaluate_fmpq_poly(res, tmp, r, len2, x);

        fmpz_clear(r);
        fmpz_clear(one);
        _fmpz_vec_clear(tmp, d);
    }
    else
    {
        fmpq_poly_t t, minpoly;
        nf_t nf;
        nf_elem_t elem;
        fmpq_mat_t mat;
        int is_power;

        /* todo: detect squaring. x^4, x^8, ...? */
        /* todo: other special cases; deflation? */
        is_power = _fmpz_vec_is_zero(poly, len - 1);

        /* nf_init wants an fmpq_poly_t, so mock up one */
        t->coeffs = QQBAR_POLY(x)->coeffs;
        t->den[0] = 1;
        t->length = QQBAR_POLY(x)->length;
        t->alloc = QQBAR_POLY(x)->alloc;

        nf_init(nf, t);
        nf_elem_init(elem, nf);

        t->coeffs = (fmpz *) poly;
        t->length = len;
        t->den[0] = *den;
        t->alloc = len;
        nf_elem_set_fmpq_poly(elem, t, nf);

        fmpq_mat_init(mat, d, d);
        nf_elem_rep_mat(mat, elem, nf);
        fmpq_poly_init(minpoly);
        fmpq_mat_minpoly(minpoly, mat);
        fmpq_mat_clear(mat);

        {
            fmpz_poly_t A;
            acb_t z, t, w;
            slong prec;
            int pure_real, pure_imag;

            A->coeffs = minpoly->coeffs;
            A->length = minpoly->length;
            A->alloc = A->length;

            acb_init(z);
            acb_init(t);
            acb_init(w);

            acb_set(z, QQBAR_ENCLOSURE(x));
            pure_real = (qqbar_sgn_im(x) == 0);
            pure_imag = (qqbar_sgn_re(x) == 0);

            for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
            {
                _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, prec);
                if (pure_real)
                    arb_zero(acb_imagref(z));
                if (pure_imag)
                    arb_zero(acb_realref(z));

                if (is_power)
                {
                    acb_pow_ui(w, z, len - 1, prec);
                    if (!fmpz_is_one(poly + len - 1))
                        acb_mul_fmpz(w, w, poly + len - 1, prec);
                    if (!fmpz_is_one(den))
                        acb_div_fmpz(w, w, den, prec);
                }
                else
                {
                    _arb_fmpz_poly_evaluate_acb(w, poly, len, z, prec);
                    if (!fmpz_is_one(den))
                        acb_div_fmpz(w, w, den, prec);
                }

                if (_qqbar_validate_uniqueness(t, A, w, 2 * prec))
                {
                    fmpz_poly_set(QQBAR_POLY(res), A);
                    acb_set(QQBAR_ENCLOSURE(res), t);
                    break;
                }
            }

            acb_clear(z);
            acb_clear(t);
            acb_clear(w);
        }

        fmpq_poly_clear(minpoly);
        nf_elem_clear(elem, nf);
        nf_clear(nf);
    }
}

void
qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpq_poly_t poly, const qqbar_t x)
{
    _qqbar_evaluate_fmpq_poly(res, poly->coeffs, poly->den, poly->length, x);
}
