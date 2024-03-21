/*
    Copyright (C) 2010, 2016 William Hart
    Copyright (C) 2010-2012 Sebastian Pancratz
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_scalar_abs(fmpz * vec1, const fmpz * vec2, slong len2)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_abs(vec1 + i, vec2 + i);
}

void
_fmpz_vec_scalar_addmul_fmpz(fmpz * poly1, const fmpz * poly2, slong len2,
                             const fmpz_t x)
{
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 0)
            return;
        else if (c == 1)
            _fmpz_vec_add(poly1, poly1, poly2, len2);
        else if (c == -1)
            _fmpz_vec_sub(poly1, poly1, poly2, len2);
        else
            _fmpz_vec_scalar_addmul_si(poly1, poly2, len2, c);
    }
    else
    {
        slong i;
        for (i = 0; i < len2; i++)
            fmpz_addmul(poly1 + i, poly2 + i, x);
    }
}

void
_fmpz_vec_scalar_addmul_si_2exp(fmpz * vec1, const fmpz * vec2, slong len2,
                                slong c, ulong exp)
{
    slong i;
    fmpz_t temp;

    if (c == 0)
        return;                 /* nothing to add */

    if (exp == 0)               /* just do addmul */
    {
        _fmpz_vec_scalar_addmul_si(vec1, vec2, len2, c);
        return;
    }

    fmpz_init(temp);

    if (c == 1)                 /* scalar is 1, just add c * 2^exp times c */
    {
        for (i = 0; i < len2; i++)
        {
            fmpz_mul_2exp(temp, vec2 + i, exp);
            fmpz_add(vec1 + i, vec1 + i, temp);
        }
    }
    else if (c == -1)           /* scalar is -1, subtract c * 2^exp */
    {
        for (i = 0; i < len2; i++)
        {
            fmpz_mul_2exp(temp, vec2 + i, exp);
            fmpz_sub(vec1 + i, vec1 + i, temp);
        }
    }
    else                        /* generic case */
    {
        if (c > 0)
        {
            for (i = 0; i < len2; i++)
            {
                fmpz_mul_2exp(temp, vec2 + i, exp);
                fmpz_addmul_ui(vec1 + i, temp, c);
            }
        }
        else
        {
            for (i = 0; i < len2; i++)
            {
                fmpz_mul_2exp(temp, vec2 + i, exp);
                fmpz_submul_ui(vec1 + i, temp, -c);
            }
        }
    }

    fmpz_clear(temp);
}

void
_fmpz_vec_scalar_addmul_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
{
    slong i;

    if (c >= 0)
        for (i = 0; i < len2; i++)
            fmpz_addmul_ui(vec1 + i, vec2 + i, c);
    else
        for (i = 0; i < len2; i++)
            fmpz_submul_ui(vec1 + i, vec2 + i, -c);
}

void
_fmpz_vec_scalar_addmul_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
{
    slong ix;

    for (ix = 0; ix < len2; ix++)
        fmpz_addmul_ui(vec1 + ix, vec2 + ix, c);
}

void
_fmpz_vec_scalar_divexact_fmpz(fmpz * vec1, const fmpz * vec2,
                               slong len2, const fmpz_t x)
{
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 1)
            _fmpz_vec_set(vec1, vec2, len2);
        else if (c == -1)
            _fmpz_vec_neg(vec1, vec2, len2);
        else
            _fmpz_vec_scalar_divexact_si(vec1, vec2, len2, c);
    }
    else
    {
        slong i;
        for (i = 0; i < len2; i++)
            fmpz_divexact(vec1 + i, vec2 + i, x);
    }
}

void
_fmpz_vec_scalar_divexact_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_divexact_si(vec1 + i, vec2 + i, c);
}

void
_fmpz_vec_scalar_divexact_ui(fmpz * vec1, const fmpz * vec2,
                             slong len2, ulong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_divexact_ui(vec1 + i, vec2 + i, c);
}

void
_fmpz_vec_scalar_fdiv_q_2exp(fmpz * vec1, const fmpz * vec2, slong len2,
                             ulong exp)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_fdiv_q_2exp(vec1 + i, vec2 + i, exp);
}

void
_fmpz_vec_scalar_fdiv_q_fmpz(fmpz * vec1, const fmpz * vec2, slong len2,
                             const fmpz_t c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_fdiv_q(vec1 + i, vec2 + i, c);
}

void
_fmpz_vec_scalar_fdiv_q_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_fdiv_q_si(vec1 + i, vec2 + i, c);
}

void
_fmpz_vec_scalar_fdiv_q_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_fdiv_q_ui(vec1 + i, vec2 + i, c);
}

void
_fmpz_vec_scalar_fdiv_r_2exp(fmpz * vec1, const fmpz * vec2, slong len2,
                             ulong exp)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_fdiv_r_2exp(vec1 + i, vec2 + i, exp);
}

void _fmpz_vec_scalar_mod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p)
{
    slong i;

    for (i = 0; i < len; i++)
        fmpz_mod(res + i, vec + i, p);
}

void
_fmpz_vec_scalar_mul_2exp(fmpz * vec1, const fmpz * vec2, slong len2, ulong exp)
{
    slong i;

    for (i = 0; i < len2; i++)
        fmpz_mul_2exp(vec1 + i, vec2 + i, exp);
}

void
_fmpz_vec_scalar_mul_fmpz(fmpz * poly1, const fmpz * poly2, slong len2,
                          const fmpz_t x)
{
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 0)
            _fmpz_vec_zero(poly1, len2);
        else if (c == 1)
            _fmpz_vec_set(poly1, poly2, len2);
        else if (c == -1)
            _fmpz_vec_neg(poly1, poly2, len2);
        else
            _fmpz_vec_scalar_mul_si(poly1, poly2, len2, c);
    }
    else
    {
        slong i;
        for (i = 0; i < len2; i++)
            fmpz_mul(poly1 + i, poly2 + i, x);
    }
}

void
_fmpz_vec_scalar_mul_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_mul_si(vec1 + i, vec2 + i, c);
}

void
_fmpz_vec_scalar_mul_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_mul_ui(vec1 + i, vec2 + i, c);
}

void _fmpz_vec_scalar_smod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p)
{
    slong i;
    fmpz_t pdiv2;

    fmpz_init(pdiv2);
    fmpz_fdiv_q_2exp(pdiv2, p, 1);

    for (i = 0; i < len; i++)
    {
        fmpz_mod(res + i, vec + i, p);

        if (fmpz_cmp(res + i, pdiv2) > 0)
        {
            fmpz_sub(res + i, res + i, p);
        }
    }

    fmpz_clear(pdiv2);
}

void
_fmpz_vec_scalar_submul_fmpz(fmpz * vec1, const fmpz * vec2, slong len2,
                             const fmpz_t x)
{
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 0)
            return;
        else if (c == 1)
            _fmpz_vec_sub(vec1, vec1, vec2, len2);
        else if (c == -1)
            _fmpz_vec_add(vec1, vec1, vec2, len2);
        else
            _fmpz_vec_scalar_submul_si(vec1, vec2, len2, c);
    }
    else
    {
        slong i;
        for (i = 0; i < len2; i++)
            fmpz_submul(vec1 + i, vec2 + i, x);
    }
}

void
_fmpz_vec_scalar_submul_si_2exp(fmpz * vec1, const fmpz * vec2, slong len2,
                                slong c, ulong exp)
{
    slong i;
    fmpz_t temp;

    if (c == 0)
        return;                 /* nothing to add */

    if (exp == 0)               /* just do submul */
    {
        _fmpz_vec_scalar_submul_si(vec1, vec2, len2, c);
        return;
    }

    fmpz_init(temp);

    if (c == 1)                 /* scalar is 1, just subtract c * 2^exp times c */
    {
        for (i = 0; i < len2; i++)
        {
            fmpz_mul_2exp(temp, vec2 + i, exp);
            fmpz_sub(vec1 + i, vec1 + i, temp);
        }
    }
    else if (c == -1)           /* scalar is -1, add c * 2^exp */
    {
        for (i = 0; i < len2; i++)
        {
            fmpz_mul_2exp(temp, vec2 + i, exp);
            fmpz_add(vec1 + i, vec1 + i, temp);
        }
    }
    else                        /* generic case */
    {
        if (c >= 0)
        {
            for (i = 0; i < len2; i++)
            {
                fmpz_mul_2exp(temp, vec2 + i, exp);
                fmpz_submul_ui(vec1 + i, temp, c);
            }
        }
        else
        {
            for (i = 0; i < len2; i++)
            {
                fmpz_mul_2exp(temp, vec2 + i, exp);
                fmpz_addmul_ui(vec1 + i, temp, -c);
            }
        }
    }

    fmpz_clear(temp);
}

void
_fmpz_vec_scalar_submul_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
{
    slong i;

    if (c >= 0)
        for (i = 0; i < len2; i++)
            fmpz_submul_ui(vec1 + i, vec2 + i, c);
    else
        for (i = 0; i < len2; i++)
            fmpz_addmul_ui(vec1 + i, vec2 + i, -c);
}

void
_fmpz_vec_scalar_tdiv_q_2exp(fmpz * vec1, const fmpz * vec2, slong len2,
    ulong exp)
{
    slong i;

    for (i = 0; i < len2; i++)
        fmpz_tdiv_q_2exp(vec1 + i, vec2 + i, exp);
}

void
_fmpz_vec_scalar_tdiv_q_fmpz(fmpz * vec1, const fmpz * vec2, slong len2,
                             const fmpz_t c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_tdiv_q(vec1 + i, vec2 + i, c);
}

void
_fmpz_vec_scalar_tdiv_q_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_tdiv_q_si(vec1 + i, vec2 + i, c);
}

void
_fmpz_vec_scalar_tdiv_q_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_tdiv_q_ui(vec1 + i, vec2 + i, c);
}
