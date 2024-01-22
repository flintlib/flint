/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_poly.h"

void _fmpq_poly_scalar_div_fmpq(fmpz * rpoly, fmpz_t rden,
                                const fmpz * poly, const fmpz_t den, slong len,
                                const fmpz_t r, const fmpz_t s)
{
    fmpz_t gcd1;  /* GCD( poly, r ) */
    fmpz_t gcd2;  /* GCD( s, den )  */
    fmpz_init(gcd1);
    fmpz_init(gcd2);
    fmpz_one(gcd1);
    fmpz_one(gcd2);
    if (!fmpz_is_one(r))
        _fmpz_vec_content_chained(gcd1, poly, len, r);
    if (!fmpz_is_one(den) && !fmpz_is_one(s))
        fmpz_gcd(gcd2, s, den);

    if (fmpz_is_one(gcd1))
    {
        if (fmpz_is_one(gcd2))
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, s);
            fmpz_mul(rden, den, r);
        }
        else
        {
            fmpz_t s2;
            fmpz_init(s2);
            fmpz_divexact(s2, s, gcd2);
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, s2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, r);
            fmpz_clear(s2);
        }
    }
    else
    {
        fmpz_t r2;
        fmpz_init(r2);
        fmpz_divexact(r2, r, gcd1);
        if (fmpz_is_one(gcd2))
        {
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, s);
            fmpz_mul(rden, den, r2);
        }
        else
        {
            fmpz_t s2;
            fmpz_init(s2);
            fmpz_divexact(s2, s, gcd2);
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, s2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, r2);
            fmpz_clear(s2);
        }
        fmpz_clear(r2);
    }

    if (_fmpz_vec_is_zero(rpoly, len))
        fmpz_one(rden);
    if (fmpz_sgn(rden) < 0)
    {
        _fmpz_vec_neg(rpoly, rpoly, len);
        fmpz_inplace_neg(rden);
    }

    fmpz_clear(gcd1);
    fmpz_clear(gcd2);
}

void fmpq_poly_scalar_div_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
{
    if (fmpq_is_zero(c))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_scalar_div_fmpq). Division by zero.\n");
    }

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, op->length);
        _fmpq_poly_set_length(rop, op->length);

        _fmpq_poly_scalar_div_fmpq(rop->coeffs, rop->den,
                                   op->coeffs, op->den, op->length,
                                   fmpq_numref(c), fmpq_denref(c));
    }
}

void _fmpq_poly_scalar_div_fmpz(fmpz * rpoly, fmpz_t rden, const fmpz * poly,
                                const fmpz_t den, slong len, const fmpz_t c)
{
    if (fmpz_is_one(c))
    {
        if (rpoly != poly)
        {
            _fmpz_vec_set(rpoly, poly, len);
            fmpz_set(rden, den);
        }
    }
    else if (*c == WORD(-1))
    {
        _fmpz_vec_neg(rpoly, poly, len);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t d;
        fmpz_init(d);
        _fmpz_vec_content_chained(d, poly, len, c);

        if (fmpz_sgn(c) < 0)
            fmpz_inplace_neg(d);
        _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, d);
        fmpz_divexact(d, c, d);
        fmpz_mul(rden, den, d);

        fmpz_clear(d);
    }
}

void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
{
    if (*c == WORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_scalar_div_fmpz). Division by zero.\n");
    }

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);

    _fmpq_poly_scalar_div_fmpz(rop->coeffs, rop->den,
                               op->coeffs, op->den, op->length, c);
}

void _fmpq_poly_scalar_div_si(fmpz * rpoly, fmpz_t rden, const fmpz * poly,
                              const fmpz_t den, slong len, slong c)
{
    if (c == 1)
    {
        if (rpoly != poly)
        {
            _fmpz_vec_set(rpoly, poly, len);
            fmpz_set(rden, den);
        }
    }
    else if (c == -1)
    {
        _fmpz_vec_neg(rpoly, poly, len);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t d, f;

        fmpz_init(d);
        fmpz_init(f);

        fmpz_set_si(f, c);
        _fmpz_vec_content_chained(d, poly, len, f);

        if (c > 0)
        {
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, d);
            fmpz_mul_si(rden, den, c / fmpz_get_si(d));
        }
        else
        {
            ulong q = (- (ulong) c) / fmpz_get_ui(d);

            fmpz_inplace_neg(d);
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, d);
            fmpz_mul_ui(rden, den, q);
        }

        fmpz_clear(d);
        fmpz_clear(f);
    }
}

void fmpq_poly_scalar_div_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
{
    if (c == WORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_scalar_div_si). Division by zero.\n");
    }

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);

    _fmpq_poly_scalar_div_si(rop->coeffs, rop->den,
                             op->coeffs, op->den, op->length, c);
}

void _fmpq_poly_scalar_div_ui(fmpz * rpoly, fmpz_t rden, const fmpz * poly,
                              const fmpz_t den, slong len, ulong c)
{
    if (c == UWORD(1))
    {
        if (rpoly != poly)
            _fmpz_vec_set(rpoly, poly, len);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t d, fc;
        ulong ud;
        fmpz_init(d);
        fmpz_init(fc);
        fmpz_set_ui(fc, c);
        _fmpz_vec_content_chained(d, poly, len, fc);

        ud = fmpz_get_ui(d);  /* gcd of d and c fits into a ulong */

        _fmpz_vec_scalar_divexact_ui(rpoly, poly, len, ud);
        fmpz_mul_ui(rden, den, c / ud);

        fmpz_clear(d);
        fmpz_clear(fc);
    }
}

void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
{
    if (c == UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_scalar_div_ui). Division by zero.\n");
    }

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);

    _fmpq_poly_scalar_div_ui(rop->coeffs, rop->den,
                             op->coeffs, op->den, op->length, c);
}

void _fmpq_poly_scalar_mul_fmpq(fmpz * rpoly, fmpz_t rden,
                                const fmpz * poly, const fmpz_t den, slong len,
                                const fmpz_t r, const fmpz_t s)
{
    fmpz_t gcd1;  /* GCD( poly, s ) */
    fmpz_t gcd2;  /* GCD( r, den )  */

    if (fmpz_is_zero(r))
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd1);
    fmpz_init(gcd2);
    fmpz_one(gcd1);
    fmpz_one(gcd2);
    if (!fmpz_is_one(s))
        _fmpz_vec_content_chained(gcd1, poly, len, s);
    if (!fmpz_is_one(den) && !fmpz_is_one(r))
        fmpz_gcd(gcd2, r, den);

    if (fmpz_is_one(gcd1))
    {
        if (fmpz_is_one(gcd2))
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, r);
            fmpz_mul(rden, den, s);
        }
        else
        {
            fmpz_t r2;
            fmpz_init(r2);
            fmpz_divexact(r2, r, gcd2);
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, r2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, s);
            fmpz_clear(r2);
        }
    }
    else
    {
        fmpz_t s2;
        fmpz_init(s2);
        fmpz_divexact(s2, s, gcd1);
        if (fmpz_is_one(gcd2))
        {
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, r);
            fmpz_mul(rden, den, s2);
        }
        else
        {
            fmpz_t r2;
            fmpz_init(r2);
            fmpz_divexact(r2, r, gcd2);
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, r2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, s2);
            fmpz_clear(r2);
        }
        fmpz_clear(s2);
    }

    fmpz_clear(gcd1);
    fmpz_clear(gcd2);
}

void fmpq_poly_scalar_mul_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
{
    if (fmpz_is_one(fmpq_denref(c)))
    {
        fmpq_poly_scalar_mul_fmpz(rop, op, fmpq_numref(c));
    }
    else if (fmpq_is_zero(c) || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, op->length);
        _fmpq_poly_set_length(rop, op->length);

        _fmpq_poly_scalar_mul_fmpq(rop->coeffs, rop->den,
                                   op->coeffs, op->den, op->length,
                                   fmpq_numref(c), fmpq_denref(c));
    }
}

void _fmpq_poly_scalar_mul_fmpz(fmpz * rpoly, fmpz_t rden,
                                const fmpz * poly, const fmpz_t den, slong len,
                                const fmpz_t c)
{
    fmpz_t gcd;  /* GCD( den, c ) */

    if (fmpz_is_zero(c))
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd);
    fmpz_one(gcd);
    if (!fmpz_is_one(c))
        fmpz_gcd(gcd, c, den);
    if (fmpz_is_one(gcd))
    {
        _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, c);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t c2;
        fmpz_init(c2);
        fmpz_divexact(c2, c, gcd);
        _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, c2);
        fmpz_divexact(rden, den, gcd);
        fmpz_clear(c2);
    }
    fmpz_clear(gcd);
}

void fmpq_poly_scalar_mul_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
{
    if (fmpz_is_zero(c) || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);

    _fmpq_poly_scalar_mul_fmpz(rop->coeffs, rop->den,
                               op->coeffs, op->den, op->length, c);
}

void _fmpq_poly_scalar_mul_si(fmpz * rpoly, fmpz_t rden,
                              const fmpz * poly, const fmpz_t den, slong len,
                              slong c)
{
    fmpz_t gcd;  /* GCD( den, c ) */

    if (c == 0)
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd);
    fmpz_set_si(gcd, c);
    fmpz_gcd(gcd, gcd, den);
    if (fmpz_is_one(gcd))
    {
        _fmpz_vec_scalar_mul_si(rpoly, poly, len, c);
        fmpz_set(rden, den);
    }
    else
    {
        if (c > WORD_MIN || fmpz_cmp_ui(gcd, - (ulong) WORD_MIN))
        {
            slong g = fmpz_get_si(gcd);

            _fmpz_vec_scalar_mul_si(rpoly, poly, len, c / g);
            fmpz_divexact_si(rden, den, g);
        }
        else
        {
            _fmpz_vec_neg(rpoly, poly, len);
            fmpz_divexact_ui(rden, den, - (ulong) WORD_MIN);
        }
    }
    fmpz_clear(gcd);
}

void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
{
    if (c == 0 || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);

    _fmpq_poly_scalar_mul_si(rop->coeffs, rop->den,
                             op->coeffs, op->den, op->length, c);
}

void _fmpq_poly_scalar_mul_ui(fmpz * rpoly, fmpz_t rden,
                              const fmpz * poly, const fmpz_t den, slong len,
                              ulong c)
{
    fmpz_t gcd;  /* GCD( den, c ) */

    if (c == 0)
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd);
    fmpz_set_ui(gcd, c);
    fmpz_gcd(gcd, gcd, den);
    if (fmpz_is_one(gcd))
    {
        _fmpz_vec_scalar_mul_ui(rpoly, poly, len, c);
        fmpz_set(rden, den);
    }
    else
    {
        ulong gcd2 = fmpz_get_ui(gcd);
        ulong c2 = c / gcd2;
        _fmpz_vec_scalar_mul_ui(rpoly, poly, len, c2);
        fmpz_fdiv_q_ui(rden, den, gcd2);
    }
    fmpz_clear(gcd);
}

void fmpq_poly_scalar_mul_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
{
    if (c == 0 || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);

    _fmpq_poly_scalar_mul_ui(rop->coeffs, rop->den,
                             op->coeffs, op->den, op->length, c);
}
