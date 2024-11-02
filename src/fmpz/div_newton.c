/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"
#include "arb-impl.h"

#if 1
#define TEST_PERTURBATION
#elif 0
#define TEST_PERTURBATION fmpz_add_ui(q, q, 1);
#else
#define TEST_PERTURBATION fmpz_sub_ui(q, q, 1);
#endif

static void
fmpz_addmul_sgn(fmpz_t res, const fmpz_t a, const fmpz_t b, int sgn)
{
    if (sgn == 1)
        fmpz_add(res, a, b);
    else
        fmpz_sub(res, a, b);
}

static ulong
_fmpz_can_round(const fmpz_t x)
{
    fmpz f = *x;
    ulong c;

    if (!COEFF_IS_MPZ(f))
        c = FLINT_ABS(f);
    else
        c = COEFF_TO_PTR(f)->_mp_d[0];

    return c > 2 && c < UWORD_MAX - 2;
}

static void
_fmpz_tdiv_qr_correction(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    while (fmpz_sgn(r) < 0)
    {
        fmpz_addmul_sgn(r, r, b, fmpz_sgn(b));
        fmpz_sub_si(q, q, fmpz_sgn(b));
    }

    while (fmpz_cmpabs(r, b) >= 0)
    {
        fmpz_addmul_sgn(r, r, b, -fmpz_sgn(b));
        fmpz_add_si(q, q, fmpz_sgn(b));
    }

    if (!fmpz_is_zero(r) && fmpz_sgn(a) < 0)
    {
        fmpz_add_si(q, q, fmpz_sgn(b));
        fmpz_addmul_sgn(r, r, b, -fmpz_sgn(b));
    }
}

static void
_fmpz_fdiv_qr_correction(fmpz_t q, fmpz_t r, const fmpz_t FLINT_UNUSED(a), const fmpz_t b)
{
    while (fmpz_sgn(r) < 0)
    {
        fmpz_addmul_sgn(r, r, b, fmpz_sgn(b));
        fmpz_sub_si(q, q, fmpz_sgn(b));
    }

    while (fmpz_cmpabs(r, b) >= 0)
    {
        fmpz_addmul_sgn(r, r, b, -fmpz_sgn(b));
        fmpz_add_si(q, q, fmpz_sgn(b));
    }

    if (!fmpz_is_zero(r) && fmpz_sgn(b) != fmpz_sgn(r))
    {
        fmpz_addmul_sgn(r, r, b, -fmpz_sgn(b));
        fmpz_add_si(q, q, fmpz_sgn(b));
    }
}

static void
_fmpz_cdiv_qr_correction(fmpz_t q, fmpz_t r, const fmpz_t FLINT_UNUSED(a), const fmpz_t b)
{
    while (fmpz_sgn(r) > 0)
    {
        fmpz_addmul_sgn(r, r, b, -fmpz_sgn(b));
        fmpz_add_si(q, q, fmpz_sgn(b));
    }

    while (fmpz_sgn(r) < 0 && fmpz_cmpabs(r, b) >= 0)
    {
        fmpz_addmul_sgn(r, r, b, fmpz_sgn(b));
        fmpz_sub_si(q, q, fmpz_sgn(b));
    }

    if (!fmpz_is_zero(r) && fmpz_sgn(b) < 0)
    {
        fmpz_add_ui(q, q, 1);
        fmpz_sub(r, r, b);
    }
}

void
_fmpz_tdiv_q_newton(fmpz_t q, const fmpz_t a, const fmpz_t b)
{
    if (q == a || q == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_tdiv_q_newton(t, a, b);
        fmpz_swap(q, t);
        fmpz_clear(t);
        return;
    }

    _arb_fmpz_divapprox_newton(q, a, b, FLINT_BITS);
    TEST_PERTURBATION

    if (_fmpz_can_round(q))
    {
        fmpz_tdiv_q_2exp(q, q, FLINT_BITS);
    }
    else
    {
        fmpz_t r;
        fmpz_init(r);
        fmpz_tdiv_q_2exp(q, q, FLINT_BITS);
        fmpz_mul(r, q, b);
        fmpz_sub(r, a, r);
        _fmpz_tdiv_qr_correction(q, r, a, b);
        fmpz_clear(r);
    }
}

void
_fmpz_fdiv_q_newton(fmpz_t q, const fmpz_t a, const fmpz_t b)
{
    if (q == a || q == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_fdiv_q_newton(t, a, b);
        fmpz_swap(q, t);
        fmpz_clear(t);
        return;
    }

    _arb_fmpz_divapprox_newton(q, a, b, FLINT_BITS);
    TEST_PERTURBATION

    if (_fmpz_can_round(q))
    {
        fmpz_fdiv_q_2exp(q, q, FLINT_BITS);
    }
    else
    {
        fmpz_t r;
        fmpz_init(r);
        fmpz_fdiv_q_2exp(q, q, FLINT_BITS);
        fmpz_mul(r, q, b);
        fmpz_sub(r, a, r);
        _fmpz_fdiv_qr_correction(q, r, a, b);
        fmpz_clear(r);
    }
}

void
_fmpz_cdiv_q_newton(fmpz_t q, const fmpz_t a, const fmpz_t b)
{
    if (q == a || q == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_cdiv_q_newton(t, a, b);
        fmpz_swap(q, t);
        fmpz_clear(t);
        return;
    }

    _arb_fmpz_divapprox_newton(q, a, b, FLINT_BITS);
    TEST_PERTURBATION

    if (_fmpz_can_round(q))
    {
        fmpz_cdiv_q_2exp(q, q, FLINT_BITS);
    }
    else
    {
        fmpz_t r;
        fmpz_init(r);
        fmpz_cdiv_q_2exp(q, q, FLINT_BITS);
        fmpz_mul(r, q, b);
        fmpz_sub(r, a, r);
        _fmpz_cdiv_qr_correction(q, r, a, b);
        fmpz_clear(r);
    }
}

void
_fmpz_fdiv_qr_newton(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    if (q == NULL || q == a || q == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_fdiv_qr_newton(t, r, a, b);
        if (q != NULL)
            fmpz_swap(q, t);
        fmpz_clear(t);
        return;
    }

    if (r == a || r == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_fdiv_qr_newton(q, t, a, b);
        fmpz_swap(r, t);
        fmpz_clear(t);
        return;
    }

    _arb_fmpz_divapprox_newton(q, a, b, 0);
    TEST_PERTURBATION
    fmpz_mul(r, q, b);
    fmpz_sub(r, a, r);
    _fmpz_fdiv_qr_correction(q, r, a, b);
}

void
_fmpz_cdiv_qr_newton(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    if (q == NULL || q == a || q == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_cdiv_qr_newton(t, r, a, b);
        if (q != NULL)
            fmpz_swap(q, t);
        fmpz_clear(t);
        return;
    }

    if (r == a || r == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_cdiv_qr_newton(q, t, a, b);
        fmpz_swap(r, t);
        fmpz_clear(t);
        return;
    }

    _arb_fmpz_divapprox_newton(q, a, b, 0);
    TEST_PERTURBATION
    fmpz_mul(r, q, b);
    fmpz_sub(r, a, r);
    _fmpz_cdiv_qr_correction(q, r, a, b);
}

void
_fmpz_tdiv_qr_newton(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    if (q == NULL || q == a || q == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_tdiv_qr_newton(t, r, a, b);
        if (q != NULL)
            fmpz_swap(q, t);
        fmpz_clear(t);
        return;
    }

    if (r == a || r == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_tdiv_qr_newton(q, t, a, b);
        fmpz_swap(r, t);
        fmpz_clear(t);
        return;
    }

    _arb_fmpz_divapprox_newton(q, a, b, 0);
    TEST_PERTURBATION
    fmpz_mul(r, q, b);
    fmpz_sub(r, a, r);
    _fmpz_tdiv_qr_correction(q, r, a, b);
}

void
_fmpz_tdiv_r_newton(fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    _fmpz_tdiv_qr_newton(NULL, r, a, b);
}

void
_fmpz_fdiv_r_newton(fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    _fmpz_fdiv_qr_newton(NULL, r, a, b);
}

void
_fmpz_cdiv_r_newton(fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    _fmpz_cdiv_qr_newton(NULL, r, a, b);
}

void
_fmpz_mod_newton(fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    fmpz_t q;
    fmpz_init(q);

    if (r == a || r == b)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_mod_newton(t, a, b);
        fmpz_swap(r, t);
        fmpz_clear(t);
        return;
    }

    _arb_fmpz_divapprox_newton(q, a, b, 0);
    TEST_PERTURBATION
    fmpz_mul(r, q, b);
    fmpz_sub(r, a, r);

    if (fmpz_sgn(b) > 0)
        _fmpz_fdiv_qr_correction(q, r, a, b);
    else
        _fmpz_cdiv_qr_correction(q, r, a, b);

    fmpz_clear(q);
}

void
_fmpz_divexact_newton(fmpz_t q, const fmpz_t a, const fmpz_t b)
{
    _arb_fmpz_divapprox_newton(q, a, b, 0);
}
