/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "acb.h"
#include "acb_mat.h"
#include "arb_mat/impl.h"

/* The approximate solving functions switch to dfloat or nfloat with
   scaling and fall back when the exponents do not fit. The test matrices
   are A = M D 2^e and B = C F 2^f with M well-conditioned, D and F diagonal
   powers of two (column scalings, to which the algorithms are invariant
   in exact exponent arithmetic) and exponents chosen near the limits of
   the backends; the solution X of AX = B, rescaled to M^-1 C, must be
   accurate to nearly the full precision. */

static void
_rand_exp(fmpz_t e, flint_rand_t state)
{
    slong limits[] = { 0, 1000, 1022, 1074, 1100, 2000,
        ARB_MAT_APPROX_DFLOAT_RANGE(1), ARB_MAT_APPROX_DFLOAT_RANGE(2),
        ARB_MAT_APPROX_DFLOAT_RANGE(3), ARB_MAT_APPROX_DFLOAT_RANGE(4),
        ARB_MAT_APPROX_MAX_EXP / 1024, ARB_MAT_APPROX_MAX_EXP, WORD_MAX / 4 };

    switch (n_randint(state, 4))
    {
        case 0:
            fmpz_zero(e);
            break;
        case 1:
            fmpz_set_si(e, (slong) n_randint(state, 100) - 50);
            break;
        case 2:
            fmpz_set_si(e, limits[n_randint(state, sizeof(limits) / sizeof(slong))]
                + (slong) n_randint(state, 21) - 10);
            break;
        default:
            /* beyond any small exponent */
            fmpz_one(e);
            fmpz_mul_2exp(e, e, FLINT_BITS + n_randint(state, 10));
            break;
    }

    if (n_randint(state, 2))
        fmpz_neg(e, e);
}

/* column scalings: exponents in [0, spread] with both ends present */
static void
_rand_col_exps(fmpz * d, slong n, flint_rand_t state)
{
    slong j;
    fmpz_t spread;

    fmpz_init(spread);

    if (n_randint(state, 3) == 0)
        fmpz_zero(spread);
    else
        _rand_exp(spread, state);
    fmpz_abs(spread, spread);

    for (j = 0; j < n; j++)
    {
        if (j == 0 || fmpz_is_zero(spread))
            fmpz_zero(d + j);
        else if (j == 1)
            fmpz_set(d + j, spread);
        else if (fmpz_bits(spread) <= FLINT_BITS - 4)
            fmpz_set_ui(d + j, n_randint(state, fmpz_get_ui(spread) + 1));
        else
            fmpz_mul_2exp(d + j, spread, 0), fmpz_fdiv_q_2exp(d + j, d + j, n_randint(state, 8));
    }

    fmpz_clear(spread);
}

/* A = M diag(2^d) 2^e */
static void
_scale(acb_mat_t A, const acb_mat_t M, const fmpz * d, const fmpz_t e)
{
    slong i, j;
    fmpz_t t;
    fmpz_init(t);

    for (i = 0; i < acb_mat_nrows(M); i++)
        for (j = 0; j < acb_mat_ncols(M); j++)
        {
            fmpz_add(t, d + j, e);
            acb_mul_2exp_fmpz(acb_mat_entry(A, i, j), acb_mat_entry(M, i, j), t);
        }

    fmpz_clear(t);
}

/* Y = diag(2^d) X diag(2^-f) 2^(e - g), recovering M^-1 C from X */
static void
_unscale(acb_mat_t Y, const acb_mat_t X, const fmpz * d, const fmpz * f, const fmpz_t e, const fmpz_t g)
{
    slong i, j;
    fmpz_t t;
    fmpz_init(t);

    for (i = 0; i < acb_mat_nrows(X); i++)
        for (j = 0; j < acb_mat_ncols(X); j++)
        {
            fmpz_add(t, d + i, e);
            fmpz_sub(t, t, f + j);
            fmpz_sub(t, t, g);
            acb_mul_2exp_fmpz(acb_mat_entry(Y, i, j), acb_mat_entry(X, i, j), t);
        }

    fmpz_clear(t);
}

/* max |Y - Yref| <= 2^-bits max |Yref| */
static int
_close(const acb_mat_t Y, const acb_mat_t Yref, slong bits)
{
    slong i, j;
    acb_t t;
    arf_t err, nrm;
    int ok;

    acb_init(t);
    arf_init(err);
    arf_init(nrm);

    for (i = 0; i < acb_mat_nrows(Y); i++)
        for (j = 0; j < acb_mat_ncols(Y); j++)
        {
            acb_sub(t, acb_mat_entry(Y, i, j), acb_mat_entry(Yref, i, j), ARF_PREC_EXACT);
            if (arf_cmpabs(arb_midref(acb_realref(t)), err) > 0)
                arf_abs(err, arb_midref(acb_realref(t)));
            if (arf_cmpabs(arb_midref(acb_imagref(t)), err) > 0)
                arf_abs(err, arb_midref(acb_imagref(t)));
            if (arf_cmpabs(arb_midref(acb_realref(acb_mat_entry(Yref, i, j))), nrm) > 0)
                arf_abs(nrm, arb_midref(acb_realref(acb_mat_entry(Yref, i, j))));
            if (arf_cmpabs(arb_midref(acb_imagref(acb_mat_entry(Yref, i, j))), nrm) > 0)
                arf_abs(nrm, arb_midref(acb_imagref(acb_mat_entry(Yref, i, j))));
        }

    arf_mul_2exp_si(nrm, nrm, -bits);
    ok = arf_is_finite(err) && arf_cmp(err, nrm) <= 0;

    acb_clear(t);
    arf_clear(err);
    arf_clear(nrm);

    return ok;
}

static slong
_rand_prec(flint_rand_t state)
{
    slong precs[] = { 2, 30, 51, 52, 53, 54, 64, 100, 104, 105, 106, 107, 128,
                      150, 153, 154, 159, 160, 192, 200, 204, 205, 212, 213,
                      256, 512, 1000, 4224, 4225, 5000 };

    if (n_randint(state, 4) == 0)
        return 2 + n_randint(state, 300);

    return precs[n_randint(state, sizeof(precs) / sizeof(slong))];
}

TEST_FUNCTION_START(acb_mat_approx_backend, state)
{
    slong iter;

    for (iter = 0; iter < 150 * flint_test_multiplier(); iter++)
    {
        acb_mat_t M, C, A, B, X, Y, Yref, LU;
        fmpz *d, *f;
        fmpz_t e, g;
        slong n, m, prec, bits, i, j, *P;
        int which, unit, result;

        n = n_randint(state, 4) ? 4 + n_randint(state, 30) : n_randint(state, 4);
        m = 1 + n_randint(state, 8);
        which = n_randint(state, 5);   /* solve, inv, lu + lu_precomp, tril, triu */
        unit = n_randint(state, 2);
        prec = _rand_prec(state);

        if (which == 1)
            m = n;

        acb_mat_init(M, n, n);
        acb_mat_init(C, n, m);
        acb_mat_init(A, n, n);
        acb_mat_init(B, n, m);
        acb_mat_init(X, n, m);
        acb_mat_init(Y, n, m);
        acb_mat_init(Yref, n, m);
        acb_mat_init(LU, n, n);
        P = flint_malloc(sizeof(slong) * FLINT_MAX(n, 1));
        d = _fmpz_vec_init(n);
        f = _fmpz_vec_init(m);
        fmpz_init(e);
        fmpz_init(g);

        /* M diagonally dominant, with entries carrying more bits than prec */
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                arb_urandom(acb_realref(acb_mat_entry(M, i, j)), state, prec + 30);
                arb_urandom(acb_imagref(acb_mat_entry(M, i, j)), state, prec + 30);
                if (n_randint(state, 2))
                    acb_neg(acb_mat_entry(M, i, j), acb_mat_entry(M, i, j));
                if (n_randint(state, 2))
                    acb_conj(acb_mat_entry(M, i, j), acb_mat_entry(M, i, j));
                if (i == j)
                    acb_add_si(acb_mat_entry(M, i, j), acb_mat_entry(M, i, j), 2 * n, prec + 60);
            }
        acb_mat_get_mid(M, M);

        if (which == 1)
            acb_mat_one(C);
        else
        {
            acb_mat_randtest(C, state, prec + 30, 4);
            acb_mat_get_mid(C, C);
        }

        if (which >= 3)
        {
            /* triangular parts; the unit case gets a unit diagonal and
               small off-diagonal entries, and no column scaling */
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    if ((which == 3) ? (j > i) : (j < i))
                        acb_zero(acb_mat_entry(M, i, j));
            if (unit)
            {
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++)
                        if (i == j)
                            acb_one(acb_mat_entry(M, i, j));
                        else
                            acb_div_si(acb_mat_entry(M, i, j), acb_mat_entry(M, i, j), 2 * n, prec + 60);
                acb_mat_get_mid(M, M);
            }
        }

        if (!(which >= 3 && unit))
        {
            _rand_col_exps(d, n, state);
            _rand_exp(e, state);
        }
        if (which != 1)
        {
            _rand_col_exps(f, m, state);
            _rand_exp(g, state);
        }

        _scale(A, M, d, e);
        _scale(B, C, f, g);

        if (which >= 3 && unit)
        {
            /* garbage on the diagonal, which must be ignored */
            for (i = 0; i < n; i++)
                acb_set_si(acb_mat_entry(A, i, i), 3 + i);
        }

        /* the reference, with no scaling */
        if (which == 3)
            acb_mat_solve_tril(Yref, M, C, 0, prec + 100);
        else if (which == 4)
            acb_mat_solve_triu(Yref, M, C, 0, prec + 100);
        else if (!acb_mat_solve(Yref, M, C, prec + 100))
        {
            flint_printf("FAIL (reference)\n");
            fflush(stdout);
            flint_abort();
        }

        result = 1;

        switch (which)
        {
            case 0:
                if (n_randint(state, 2))
                {
                    result = acb_mat_approx_solve(X, A, B, prec);
                }
                else
                {
                    /* aliasing */
                    acb_mat_set(X, B);
                    result = acb_mat_approx_solve(X, A, X, prec);
                }
                break;
            case 1:
                if (n_randint(state, 2))
                {
                    result = acb_mat_approx_inv(X, A, prec);
                }
                else
                {
                    acb_mat_set(X, A);
                    result = acb_mat_approx_inv(X, X, prec);
                }
                break;
            case 2:
                result = acb_mat_approx_lu(P, LU, A, prec);
                if (result)
                    acb_mat_approx_solve_lu_precomp(X, P, LU, B, prec);
                break;
            case 3:
                acb_mat_approx_solve_tril(X, A, B, unit, prec);
                break;
            default:
                acb_mat_approx_solve_triu(X, A, B, unit, prec);
                break;
        }

        _unscale(Y, X, d, f, e, g);

        /* a few bits for the conditioning of M (at most 3), the size and
           the effective precision of dfloat */
        bits = prec - 12 - FLINT_BIT_COUNT(n);

        if (!result || !_close(Y, Yref, bits))
        {
            flint_printf("FAIL\n");
            flint_printf("which = %d, unit = %d, n = %wd, m = %wd, prec = %wd, result = %d\n",
                which, unit, n, m, prec, result);
            flint_printf("e = "); fmpz_print(e); flint_printf("\n");
            flint_printf("g = "); fmpz_print(g); flint_printf("\n");
            flint_printf("d = "); _fmpz_vec_print(d, n); flint_printf("\n");
            flint_printf("f = "); _fmpz_vec_print(f, m); flint_printf("\n");
            flint_printf("Y = "); acb_mat_printd(Y, 10); flint_printf("\n");
            flint_printf("Yref = "); acb_mat_printd(Yref, 10); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(M);
        acb_mat_clear(C);
        acb_mat_clear(A);
        acb_mat_clear(B);
        acb_mat_clear(X);
        acb_mat_clear(Y);
        acb_mat_clear(Yref);
        acb_mat_clear(LU);
        flint_free(P);
        _fmpz_vec_clear(d, n);
        _fmpz_vec_clear(f, m);
        fmpz_clear(e);
        fmpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}

/* max over the rows of |E_i - Eref_j| minimized over j, against 2^-bits max |Eref| */
static int
_close_eig(acb_srcptr E, acb_srcptr Eref, slong n, slong bits)
{
    slong i, j;
    acb_t t;
    arf_t err, best, nrm, u;
    int ok;

    acb_init(t);
    arf_init(err);
    arf_init(best);
    arf_init(nrm);
    arf_init(u);

    for (i = 0; i < n; i++)
    {
        arf_pos_inf(best);
        for (j = 0; j < n; j++)
        {
            acb_sub(t, E + i, Eref + j, ARF_PREC_EXACT);
            arf_abs(u, arb_midref(acb_realref(t)));
            if (arf_cmpabs(arb_midref(acb_imagref(t)), u) > 0)
                arf_abs(u, arb_midref(acb_imagref(t)));
            if (arf_cmp(u, best) < 0)
                arf_set(best, u);
        }
        if (arf_cmp(best, err) > 0)
            arf_set(err, best);
        if (arf_cmpabs(arb_midref(acb_realref(Eref + i)), nrm) > 0)
            arf_abs(nrm, arb_midref(acb_realref(Eref + i)));
    }

    arf_mul_2exp_si(nrm, nrm, -bits);
    ok = arf_is_finite(err) && arf_cmp(err, nrm) <= 0;

    acb_clear(t);
    arf_clear(err);
    arf_clear(best);
    arf_clear(nrm);
    arf_clear(u);

    return ok;
}

TEST_FUNCTION_START(acb_mat_approx_eig_qr_backend, state)
{
    slong iter;

    /* Hermitian matrices (perfectly conditioned eigenvalues) scaled by 2^e
       with e near the exponent limits of the backends */
    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        acb_mat_t H, A, L, R, Rref, T, U;
        acb_ptr E, Eref, Eapprox;
        fmpz_t e, me;
        slong n, prec, bits, i, j;
        int result, want_L, want_R;

        n = n_randint(state, 4) ? 4 + n_randint(state, 14) : n_randint(state, 4);
        prec = _rand_prec(state);
        if (prec > 1000 && n_randint(state, 4))
            prec = 2 + n_randint(state, 300);
        want_L = n_randint(state, 2);
        want_R = n_randint(state, 2);

        acb_mat_init(H, n, n);
        acb_mat_init(A, n, n);
        acb_mat_init(L, n, n);
        acb_mat_init(R, n, n);
        acb_mat_init(Rref, n, n);
        acb_mat_init(T, n, n);
        acb_mat_init(U, n, n);
        E = _acb_vec_init(n);
        Eref = _acb_vec_init(n);
        Eapprox = _acb_vec_init(n);
        fmpz_init(e);
        fmpz_init(me);

        for (i = 0; i < n; i++)
            for (j = 0; j <= i; j++)
            {
                arb_urandom(acb_realref(acb_mat_entry(H, i, j)), state, prec + 30);
                if (n_randint(state, 2))
                    arb_neg(acb_realref(acb_mat_entry(H, i, j)), acb_realref(acb_mat_entry(H, i, j)));
                if (i != j)
                {
                    arb_urandom(acb_imagref(acb_mat_entry(H, i, j)), state, prec + 30);
                    if (n_randint(state, 2))
                        arb_neg(acb_imagref(acb_mat_entry(H, i, j)), acb_imagref(acb_mat_entry(H, i, j)));
                    acb_conj(acb_mat_entry(H, j, i), acb_mat_entry(H, i, j));
                }
            }

        _rand_exp(e, state);
        fmpz_neg(me, e);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                acb_mul_2exp_fmpz(acb_mat_entry(A, i, j), acb_mat_entry(H, i, j), e);

        result = acb_mat_approx_eig_qr(E, want_L ? L : NULL, want_R ? R : NULL, A, NULL, 0, prec);

        for (i = 0; i < n; i++)
            acb_mul_2exp_fmpz(E + i, E + i, me);

        /* a rigorous reference */
        acb_mat_approx_eig_qr(Eapprox, NULL, Rref, H, NULL, 0, prec + 100);
        if (!acb_mat_eig_simple(Eref, NULL, NULL, H, Eapprox, Rref, prec + 100))
        {
            flint_printf("FAIL (reference)\n");
            fflush(stdout);
            flint_abort();
        }

        bits = prec - 12 - FLINT_BIT_COUNT(n);

        if (!result || !_close_eig(E, Eref, n, bits))
        {
            flint_printf("FAIL (eigenvalues)\n");
            flint_printf("n = %wd, prec = %wd, result = %d\n", n, prec, result);
            flint_printf("e = "); fmpz_print(e); flint_printf("\n");
            flint_printf("E = "); _acb_vec_printd(E, n, 10); flint_printf("\n");
            flint_printf("Eref = "); _acb_vec_printd(Eref, n, 10); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* residuals H R - R diag(E) and L H - diag(E) L */
        if (want_R || want_L)
        {
            for (j = 0; j < 2; j++)
            {
                if (j == 0 && !want_R) continue;
                if (j == 1 && !want_L) continue;

                if (j == 0)
                {
                    acb_mat_mul(T, H, R, prec + 100);
                    acb_mat_set(U, R);
                    for (i = 0; i < n; i++)
                    {
                        slong k;
                        for (k = 0; k < n; k++)
                            acb_mul(acb_mat_entry(U, k, i), acb_mat_entry(U, k, i), E + i, prec + 100);
                    }
                }
                else
                {
                    acb_mat_mul(T, L, H, prec + 100);
                    acb_mat_set(U, L);
                    for (i = 0; i < n; i++)
                    {
                        slong k;
                        for (k = 0; k < n; k++)
                            acb_mul(acb_mat_entry(U, i, k), acb_mat_entry(U, i, k), E + i, prec + 100);
                    }
                }

                /* compare against the scale |H| |V| */
                {
                    acb_mat_t S;
                    arb_t hn, vn;
                    acb_mat_init(S, n, n);
                    arb_init(hn);
                    arb_init(vn);
                    acb_mat_frobenius_norm(hn, H, 60);
                    acb_mat_frobenius_norm(vn, (j == 0) ? R : L, 60);
                    arb_mul(hn, hn, vn, 60);
                    acb_mat_sub(S, T, U, prec + 100);
                    acb_mat_frobenius_norm(vn, S, 60);
                    arb_mul_2exp_si(hn, hn, -bits);
                    if (arf_cmp(arb_midref(vn), arb_midref(hn)) > 0)
                    {
                        flint_printf("FAIL (residual, %s)\n", j == 0 ? "R" : "L");
                        flint_printf("n = %wd, prec = %wd\n", n, prec);
                        flint_printf("e = "); fmpz_print(e); flint_printf("\n");
                        arb_printd(vn, 10); flint_printf("   ");
                        arb_printd(hn, 10); flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }
                    acb_mat_clear(S);
                    arb_clear(hn);
                    arb_clear(vn);
                }
            }
        }

        acb_mat_clear(H);
        acb_mat_clear(A);
        acb_mat_clear(L);
        acb_mat_clear(R);
        acb_mat_clear(Rref);
        acb_mat_clear(T);
        acb_mat_clear(U);
        _acb_vec_clear(E, n);
        _acb_vec_clear(Eref, n);
        _acb_vec_clear(Eapprox, n);
        fmpz_clear(e);
        fmpz_clear(me);
    }

    TEST_FUNCTION_END(state);
}
