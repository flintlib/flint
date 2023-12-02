/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

void
_qqbar_sqr_undeflatable(qqbar_t res, const qqbar_t x)
{
    fmpz_poly_t A, B;
    acb_t z, t, w;
    slong i, prec, d;
    int pure_real, pure_imag;

    fmpz_poly_init(A);
    fmpz_poly_init(B);
    acb_init(z);
    acb_init(t);
    acb_init(w);

    d = fmpz_poly_degree(QQBAR_POLY(x));
    for (i = 0; i <= d; i++)
    {
        if (i % 2 == 0)
            fmpz_poly_set_coeff_fmpz(A, i / 2, QQBAR_POLY(x)->coeffs + i);
        else
            fmpz_poly_set_coeff_fmpz(B, i / 2, QQBAR_POLY(x)->coeffs + i);
    }

    fmpz_poly_sqr(A, A);
    fmpz_poly_sqr(B, B);
    fmpz_poly_shift_left(B, B, 1);
    fmpz_poly_sub(A, A, B);
    if (fmpz_sgn(A->coeffs + A->length - 1) < 0)
        fmpz_poly_neg(A, A);

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

        acb_sqr(w, z, prec);

        if (_qqbar_validate_uniqueness(t, A, w, 2 * prec))
        {
            fmpz_poly_set(QQBAR_POLY(res), A);
            acb_set(QQBAR_ENCLOSURE(res), t);
            break;
        }
    }

    fmpz_poly_clear(A);
    fmpz_poly_clear(B);
    acb_clear(z);
    acb_clear(t);
    acb_clear(w);
}

/* todo: use number field arithmetic, other optimisations ... */
void
qqbar_pow_ui(qqbar_t res, const qqbar_t x, ulong n)
{
    if (n == 0)
    {
        qqbar_one(res);
    }
    else if (n == 1)
    {
        qqbar_set(res, x);
    }
    else if (qqbar_is_rational(x))
    {
        fmpq_t t;
        fmpq_init(t);
        qqbar_get_fmpq(t, x);
        fmpz_pow_ui(fmpq_numref(t), fmpq_numref(t), n);
        fmpz_pow_ui(fmpq_denref(t), fmpq_denref(t), n);
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
    }
    else
    {
        slong p;
        ulong q;
        ulong f;

        /* Fast path for roots of unity. Todo: generalize to
           rational multiples of roots of unity. */
        if (qqbar_is_root_of_unity(&p, &q, x))
        {
            if (p < 0)
                p += 2 * q;

            p = n_mulmod2(p, n, 2 * q);
            qqbar_root_of_unity(res, p, q);
            return;
        }

        /* Fast detection of perfect powers */
        f = arb_fmpz_poly_deflation(QQBAR_POLY(x));

        if (f % n == 0)
        {
            acb_t z, t, w;
            fmpz_poly_t H;
            slong prec;
            int pure_real, pure_imag;

            fmpz_poly_init(H);
            acb_init(z);
            acb_init(t);
            acb_init(w);

            arb_fmpz_poly_deflate(H, QQBAR_POLY(x), n);

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

                acb_pow_ui(w, z, n, prec);

                if (_qqbar_validate_uniqueness(t, H, w, 2 * prec))
                {
                    fmpz_poly_set(QQBAR_POLY(res), H);
                    acb_set(QQBAR_ENCLOSURE(res), t);
                    break;
                }
            }

            fmpz_poly_clear(H);
            acb_clear(z);
            acb_clear(t);
            acb_clear(w);
            return;
        }

        /* fast path for principal roots of positive rational numbers */
        if (_qqbar_fast_detect_simple_principal_surd(x))
        {
            fmpq_t t;
            fmpq_init(t);
            fmpz_neg(fmpq_numref(t), QQBAR_COEFFS(x));
            fmpz_set(fmpq_denref(t), QQBAR_COEFFS(x) + qqbar_degree(x));
            fmpq_pow_si(t, t, n);
            qqbar_fmpq_root_ui(res, t, qqbar_degree(x));
            fmpq_clear(t);
            return;
        }

        if (n == 2)
        {
            _qqbar_sqr_undeflatable(res, x);
        }
        else
        {
            /* todo: don't construct this polynomial when n is huge */
            fmpz * coeffs;
            fmpz_t den;
            coeffs = _fmpz_vec_init(n + 1);
            fmpz_one(coeffs + n);
            *den = 1;
            _qqbar_evaluate_fmpq_poly(res, coeffs, den, n + 1, x);
            _fmpz_vec_clear(coeffs, n + 1);
        }
    }
}

void
qqbar_pow_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
{
    if (fmpq_is_zero(y))
    {
        qqbar_one(res);
    }
    else if (fmpq_is_one(y))
    {
        qqbar_set(res, x);
    }
    else if (qqbar_is_one(x))
    {
        qqbar_one(res);
    }
    else if (qqbar_is_zero(x))
    {
        if (fmpq_sgn(y) <= 0)
        {
            flint_throw(FLINT_ERROR, "qqbar_pow_fmpq: division by zero\n");
        }

        qqbar_zero(res);
    }
    else
    {
        fmpq_t t;
        fmpz_t r;
        slong p;
        ulong q;

        fmpq_init(t);
        fmpz_init(r);
        fmpq_set(t, y);

        /* Fast path for roots of unity. */
        if (qqbar_is_root_of_unity(&p, &q, x))
        {
            fmpz_mul_si(fmpq_numref(t), fmpq_numref(t), p);
            fmpz_mul_ui(fmpq_denref(t), fmpq_denref(t), q);
            fmpz_mul_ui(r, fmpq_denref(t), 2);
            fmpz_fdiv_r(fmpq_numref(t), fmpq_numref(t), r);
            fmpq_canonicalise(t);

            if (COEFF_IS_MPZ(*fmpq_denref(t)))
            {
                flint_throw(FLINT_ERROR, "qqbar_pow: excessive exponent\n");
            }

            qqbar_root_of_unity(res, *fmpq_numref(t), *fmpq_denref(t));
        }
        else
        {
            if (COEFF_IS_MPZ(*fmpq_numref(t)) || COEFF_IS_MPZ(*fmpq_denref(t)))
            {
                flint_throw(FLINT_ERROR, "qqbar_pow: excessive exponent\n");
            }

            p = *fmpq_numref(t);
            q = *fmpq_denref(t);

            qqbar_root_ui(res, x, q);
            if (p >= 0)
                qqbar_pow_ui(res, res, p);
            else
            {
                qqbar_pow_ui(res, res, -p);
                qqbar_inv(res, res);
            }
        }

        fmpq_clear(t);
        fmpz_clear(r);
    }
}

void
qqbar_pow_si(qqbar_t res, const qqbar_t x, slong n)
{
    if (n >= 0)
    {
        qqbar_pow_ui(res, x, n);
    }
    else
    {
        fmpq_t t;
        fmpq_init(t);
        fmpz_set_si(fmpq_numref(t), n);
        qqbar_pow_fmpq(res, x, t);
        fmpq_clear(t);
    }
}

void
qqbar_pow_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t n)
{
    fmpq_t t;
    fmpq_init(t);
    fmpz_set(fmpq_numref(t), n);
    qqbar_pow_fmpq(res, x, t);
    fmpq_clear(t);
}

int
qqbar_pow(qqbar_t res, const qqbar_t x, const qqbar_t y)
{
    if (qqbar_is_zero(y))
    {
        qqbar_one(res);
        return 1;
    }
    else if (qqbar_is_one(y))
    {
        qqbar_set(res, x);
        return 1;
    }
    else if (qqbar_is_one(x))
    {
        qqbar_one(res);
        return 1;
    }
    else if (qqbar_is_zero(x))
    {
        if (qqbar_sgn_re(y) <= 0)
            return 0;

        qqbar_zero(res);
        return 1;
    }
    else if (qqbar_is_rational(y))
    {
        fmpq_t t;
        fmpq_init(t);
        qqbar_get_fmpq(t, y);
        qqbar_pow_fmpq(res, x, t);
        fmpq_clear(t);
        return 1;
    }
    else
    {
        return 0;
    }
}
