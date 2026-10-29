/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"

static void
_vec_add_error(arb_ptr v, slong len, const mag_t err)
{
    slong j;
    for (j = 0; j < len; j++)
        arb_add_error_mag(v + j, err);
}

static int
_vec_overlaps(arb_srcptr a, arb_srcptr b, slong len)
{
    slong j;
    for (j = 0; j < len; j++)
        if (!arb_overlaps(a + j, b + j))
            return 0;
    return 1;
}

static void
_fail(const char * what, const fmpq_t a, const arb_t z, arb_srcptr r1, arb_srcptr r2, slong len)
{
    slong j;
    flint_printf("FAIL: %s\n\n", what);
    flint_printf("a = "); fmpq_print(a); flint_printf("\n\n");
    flint_printf("z = "); arb_printn(z, 50, 0); flint_printf("\n\n");
    for (j = 0; j < len; j++)
    {
        flint_printf("r1[%wd] = ", j); arb_printn(r1 + j, 50, 0); flint_printf("\n");
        flint_printf("r2[%wd] = ", j); arb_printn(r2 + j, 50, 0); flint_printf("\n");
    }
    flint_abort();
}

/* Gamma(a+eps, z) = Gamma(a+eps) - gamma(a+eps, z), a not a pole */
static int
_gamma_upper_via_lower(arb_ptr res, const fmpq_t a, const arb_t z, slong len, slong prec)
{
    arb_ptr G;
    mag_t tol, err;
    slong N;

    G = _arb_vec_init(len);
    mag_init(tol);
    mag_init(err);

    mag_set_ui_2exp_si(tol, 1, -prec);
    N = _arb_hypgeom_gamma_lower_fmpq_0_choose_N_jet(err, a, z, tol, len);

    if (N >= 0)
    {
        _arb_hypgeom_gamma_fmpq_jet(G, a, len, prec);
        _arb_hypgeom_gamma_lower_fmpq_0_bsplit_jet(res, a, z, NULL, N, len, prec);
        _vec_add_error(res, len, err);
        _arb_vec_sub(res, G, res, len, prec);
    }

    _arb_vec_clear(G, len);
    mag_clear(tol);
    mag_clear(err);

    return N >= 0;
}

TEST_FUNCTION_START(arb_hypgeom_gamma_upper_fmpq_jet, state)
{
    slong iter;

    /* Gamma(a+eps) against arb_poly */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t a;
        arb_ptr r1, r2, h;
        slong len, prec;
        arb_t z;

        len = 1 + n_randint(state, 10);
        prec = 2 + n_randint(state, 300);

        fmpq_init(a);
        arb_init(z);
        r1 = _arb_vec_init(len);
        r2 = _arb_vec_init(len);
        h = _arb_vec_init(2);

        do {
            fmpq_randtest(a, state, 6);
        } while (fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) <= 0);

        _arb_hypgeom_gamma_fmpq_jet(r1, a, len, prec);
        arb_set_fmpq(h, a, prec + 30);
        arb_one(h + 1);
        _arb_poly_gamma_series(r2, h, 2, len, prec + 30);

        if (!_vec_overlaps(r1, r2, len))
            _fail("gamma", a, z, r1, r2, len);

        fmpq_clear(a);
        arb_clear(z);
        _arb_vec_clear(r1, len);
        _arb_vec_clear(r2, len);
        _arb_vec_clear(h, 2);
    }

    /* asymptotic series, lower series and bit-burst step */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t a;
        arb_t z, z1, e;
        arb_ptr r1, r2, r3;
        mag_t tol, err;
        slong len, prec, N;

        len = 1 + n_randint(state, 8);
        prec = 2 + n_randint(state, 200);

        fmpq_init(a);
        arb_init(z);
        arb_init(z1);
        arb_init(e);
        mag_init(tol);
        mag_init(err);
        r1 = _arb_vec_init(len);
        r2 = _arb_vec_init(len);
        r3 = _arb_vec_init(len);

        do {
            fmpq_randtest(a, state, 5);
        } while (fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) <= 0);

        arb_set_ui(z, 1 + n_randint(state, 1000));
        arb_div_ui(z, z, 1 + n_randint(state, 20), 200);
        arb_get_mid_arb(z, z);

        if (_gamma_upper_via_lower(r1, a, z, len, prec))
        {
            mag_set_ui_2exp_si(tol, 1, -(slong) n_randint(state, 200));
            N = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N_jet(err, a, z, tol, len);

            if (N >= 0)
            {
                _arb_hypgeom_gamma_upper_fmpq_inf_bsplit_jet(r2, a, z, NULL, N, len, 2 + n_randint(state, 200));
                _vec_add_error(r2, len, err);

                if (!_vec_overlaps(r1, r2, len))
                    _fail("asymptotic", a, z, r1, r2, len);
            }

            /* step to z1 = z + small */
            arb_set_si(z1, -100 + (slong) n_randint(state, 200));
            arb_mul_2exp_si(z1, z1, -(slong) n_randint(state, 12));
            arb_mul(z1, z1, z, 200);
            arb_mul_2exp_si(z1, z1, -8);
            arb_add(z1, z1, z, 200);
            arb_get_mid_arb(z1, z1);

            if (arb_is_positive(z1) && _gamma_upper_via_lower(r2, a, z1, len, prec))
            {
                arb_neg(e, z);
                arb_exp(e, e, prec + 20);
                mag_set_ui_2exp_si(tol, 1, -(slong) n_randint(state, 200));
                if (n_randint(state, 2))
                {
                    _arb_gamma_upper_fmpq_step_bsplit_jet(r3, a, z, NULL, z1, r1, e, tol, len, prec);
                }
                else
                {
                    arb_t lz;
                    arb_init(lz);
                    arb_log(lz, z, prec + 10);
                    _arb_gamma_upper_fmpq_step_bsplit_jet(r3, a, z, lz, z1, r1, e, tol, len, prec);
                    arb_clear(lz);
                }

                if (!_vec_overlaps(r2, r3, len))
                    _fail("step", a, z1, r2, r3, len);
            }
        }

        fmpq_clear(a);
        arb_clear(z);
        arb_clear(z1);
        arb_clear(e);
        mag_clear(tol);
        mag_clear(err);
        _arb_vec_clear(r1, len);
        _arb_vec_clear(r2, len);
        _arb_vec_clear(r3, len);
    }

    /* singular case: Gamma(-n+eps, z) = (Gamma(-n+1+eps, z) - z^(-n+eps) e^(-z)) / (-n+eps) */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t a;
        arb_t z, e, lz;
        arb_ptr r1, r2, r3, d;
        mag_t tol, err;
        slong len, prec, N, n, j;

        len = 1 + n_randint(state, 8);
        prec = 2 + n_randint(state, 200);
        n = n_randint(state, 6);

        fmpq_init(a);
        arb_init(z);
        arb_init(e);
        arb_init(lz);
        mag_init(tol);
        mag_init(err);
        r1 = _arb_vec_init(len + 1);
        r2 = _arb_vec_init(len + 1);
        r3 = _arb_vec_init(len + 1);
        d = _arb_vec_init(2);

        arb_set_ui(z, 1 + n_randint(state, 100));
        arb_div_ui(z, z, 1 + n_randint(state, 20), 200);
        arb_get_mid_arb(z, z);

        /* r1 = Gamma(-n+1+eps, z), length len + 1 */
        if (n == 0)
        {
            fmpq_one(a);
            if (!_gamma_upper_via_lower(r1, a, z, len + 1, prec + 30))
                goto cleanup;
        }
        else
        {
            mag_set_ui_2exp_si(tol, 1, -prec - 30);
            N = _arb_hypgeom_gamma_upper_singular_si_choose_N_jet(err, n - 1, z, tol, len + 1);
            _arb_hypgeom_gamma_upper_singular_si_bsplit_jet(r1, n - 1, z, NULL, N, NULL, len + 1, prec + 30);
            _vec_add_error(r1, len + 1, err);
        }

        /* r2 = z^(-n+eps) e^(-z) */
        fmpq_set_si(a, -n, 1);
        arb_pow_fmpq(r2, z, a, prec + 30);
        arb_neg(e, z);
        arb_exp(e, e, prec + 30);
        arb_mul(r2, r2, e, prec + 30);
        arb_log(lz, z, prec + 30);
        for (j = 1; j <= len; j++)
        {
            arb_mul(r2 + j, r2 + j - 1, lz, prec + 30);
            arb_div_ui(r2 + j, r2 + j, j, prec + 30);
        }

        _arb_vec_sub(r1, r1, r2, len + 1, prec + 30);

        if (n == 0)
        {
            _arb_vec_set(r3, r1 + 1, len);
        }
        else
        {
            arb_set_si(d, -n);
            arb_one(d + 1);
            _arb_poly_div_series(r3, r1, len, d, 2, len, prec + 30);
        }

        mag_set_ui_2exp_si(tol, 1, -(slong) n_randint(state, 200));
        N = _arb_hypgeom_gamma_upper_singular_si_choose_N_jet(err, n, z, tol, len);
        if (n_randint(state, 2))
        {
            _arb_hypgeom_gamma_upper_singular_si_bsplit_jet(r2, n, z, NULL, N, NULL, len, 2 + n_randint(state, 200));
        }
        else
        {
            slong wp = 2 + n_randint(state, 200);
            fmpq_one(a);
            _arb_hypgeom_gamma_fmpq_jet(r1, a, len + 1, wp);
            arb_log(lz, z, wp);
            _arb_hypgeom_gamma_upper_singular_si_bsplit_jet(r2, n, z, lz, N, r1, len, wp);
        }
        _vec_add_error(r2, len, err);

        if (!_vec_overlaps(r2, r3, len))
        {
            fmpq_set_si(a, -n, 1);
            _fail("singular", a, z, r2, r3, len);
        }

cleanup:
        fmpq_clear(a);
        arb_clear(z);
        arb_clear(e);
        arb_clear(lz);
        mag_clear(tol);
        mag_clear(err);
        _arb_vec_clear(r1, len + 1);
        _arb_vec_clear(r2, len + 1);
        _arb_vec_clear(r3, len + 1);
        _arb_vec_clear(d, 2);
    }

    TEST_FUNCTION_END(state);
}
