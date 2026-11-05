/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"
#include "arb_mat.h"
#include "acb.h"
#include "acb_mzeta.h"

/* A = L^T for G = L L^T, x = A u, y = A^-T v (solve L y = v) */
static void
gram_to_basis(arb_mat_t A, arb_ptr x, arb_ptr y, const fmpq_mat_t G,
    const fmpq * u, const fmpq * v, slong prec)
{
    slong d = fmpq_mat_nrows(G), i, j;
    arb_mat_t Ga, L;
    arb_ptr t;

    arb_mat_init(Ga, d, d);
    arb_mat_init(L, d, d);
    t = _arb_vec_init(d);
    arb_mat_set_fmpq_mat(Ga, G, prec);
    arb_mat_cho(L, Ga, prec);
    arb_mat_transpose(A, L);
    for (i = 0; i < d; i++)
        arb_set_fmpq(t + i, u + i, prec);
    for (i = 0; i < d; i++)
    {
        arb_zero(x + i);
        for (j = 0; j < d; j++)
            arb_addmul(x + i, arb_mat_entry(A, i, j), t + j, prec);
    }
    for (i = 0; i < d; i++)
    {
        arb_set_fmpq(y + i, v + i, prec);
        for (j = 0; j < i; j++)
            arb_submul(y + i, arb_mat_entry(L, i, j), y + j, prec);
        arb_div(y + i, y + i, arb_mat_entry(L, i, i), prec);
    }
    arb_mat_clear(Ga);
    arb_mat_clear(L);
    _arb_vec_clear(t, d);
}

TEST_FUNCTION_START(acb_mzeta_epstein_shell, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        slong d, i, j, k, prec;
        fmpz_mat_t U;
        fmpq_mat_t G, G2, Uq, UqT, Tm;
        fmpq * u, * v, * u2, * v2;
        arb_mat_t A;
        arb_ptr x, y;
        acb_t nu, z1, z2;
        fmpq_t c;
        int ok1, ok2;
        static const slong nus[] = { 1, 2, 3, 4, -1, -2, 6 };

        d = 1 + n_randint(state, 3);
        if (n_randint(state, 8) == 0)
            d = 4;
        prec = 64 + n_randint(state, 200);

        fmpz_mat_init(U, d, d);
        fmpq_mat_init(G, d, d);
        fmpq_mat_init(G2, d, d);
        fmpq_mat_init(Uq, d, d);
        fmpq_mat_init(UqT, d, d);
        fmpq_mat_init(Tm, d, d);
        u = _fmpq_vec_init(d);
        v = _fmpq_vec_init(d);
        u2 = _fmpq_vec_init(d);
        v2 = _fmpq_vec_init(d);
        arb_mat_init(A, d, d);
        x = _arb_vec_init(d);
        y = _arb_vec_init(d);
        acb_init(nu);
        acb_init(z1);
        acb_init(z2);
        fmpq_init(c);

        /* random positive definite G with small entries: diagonally
           dominant, sometimes scaled by 1/2 or 1/3 */
        for (i = 0; i < d; i++)
            for (j = 0; j <= i; j++)
            {
                slong e = (i == j) ? 0 : (slong) n_randint(state, 5) - 2;
                fmpq_set_si(fmpq_mat_entry(G, i, j), e, 1);
                fmpq_set_si(fmpq_mat_entry(G, j, i), e, 1);
            }
        for (i = 0; i < d; i++)
        {
            slong s = 1 + n_randint(state, 3);
            for (j = 0; j < d; j++)
                if (j != i)
                    s += FLINT_ABS(fmpz_get_si(fmpq_mat_entry_num(G, i, j)));
            fmpq_set_si(fmpq_mat_entry(G, i, i), s, 1);
        }
        if (n_randint(state, 3) == 0)
        {
            fmpq_set_si(c, 1, 2 + n_randint(state, 2));
            fmpq_mat_scalar_mul_fmpq(G, G, c);
        }

        for (i = 0; i < d; i++)
        {
            static const slong dens[] = { 1, 1, 2, 2, 3, 4, 6 };
            slong du = dens[n_randint(state, 7)], dv = dens[n_randint(state, 7)];
            fmpq_set_si(u + i, (slong) n_randint(state, 2 * du + 1) - du, du);
            fmpq_set_si(v + i, (slong) n_randint(state, 2 * dv + 1) - dv, dv);
        }

        if (n_randint(state, 4) == 0)
            acb_set_si(nu, d - 2 * (slong) n_randint(state, 2));      /* Coulomb-like */
        else if (n_randint(state, 4) == 0)
            acb_set_d(nu, (double) (1 + 2 * (slong) n_randint(state, 4)) / 2);
        else
            acb_set_si(nu, nus[n_randint(state, 7)]);

        ok1 = acb_mzeta_epstein_gram_shell(z1, nu, G, u, v, 0, prec);

        if (n_randint(state, 2))
        {
            /* against Crandall */
            gram_to_basis(A, x, y, G, u, v, prec + 64);
            ok2 = acb_mzeta_epstein_crandall(z2, nu, A, x, y, 0, NULL, prec);
        }
        else
        {
            /* basis change n = U n': G' = U^T G U, u' = U^-1 u, v' = U^T v */
            fmpz_mat_one(U);
            for (k = 0; k < 3; k++)
            {
                i = n_randint(state, d);
                j = n_randint(state, d);
                if (i != j)
                {
                    slong cc = (slong) n_randint(state, 3) - 1, l;
                    for (l = 0; l < d; l++)
                        fmpz_addmul_si(fmpz_mat_entry(U, i, l), fmpz_mat_entry(U, j, l), cc);
                }
            }
            fmpq_mat_set_fmpz_mat(Uq, U);
            fmpq_mat_transpose(UqT, Uq);
            fmpq_mat_mul(Tm, UqT, G);
            fmpq_mat_mul(G2, Tm, Uq);
            fmpq_mat_inv(Tm, Uq);
            for (i = 0; i < d; i++)
            {
                fmpq_zero(u2 + i);
                fmpq_zero(v2 + i);
                for (j = 0; j < d; j++)
                {
                    fmpq_addmul(u2 + i, fmpq_mat_entry(Tm, i, j), u + j);
                    fmpq_addmul(v2 + i, fmpq_mat_entry(UqT, i, j), v + j);
                }
            }
            ok2 = acb_mzeta_epstein_gram_shell(z2, nu, G2, u2, v2, 0, prec + 30);
        }

        if (ok1 && ok2 && !acb_overlaps(z1, z2))
        {
            flint_printf("FAIL (iter %wd)\n", iter);
            flint_printf("d = %wd, prec = %wd\n", d, prec);
            flint_printf("G = "); fmpq_mat_print(G); flint_printf("\n");
            flint_printf("u = "); _fmpq_vec_print(u, d); flint_printf("\n");
            flint_printf("v = "); _fmpq_vec_print(v, d); flint_printf("\n");
            flint_printf("nu = "); acb_printd(nu, 10); flint_printf("\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n");
            flint_printf("z2 = "); acb_printd(z2, 30); flint_printf("\n");
            flint_abort();
        }

        if (!ok1 && d <= 3)
        {
            flint_printf("FAIL (not applicable)\n");
            flint_printf("G = "); fmpq_mat_print(G); flint_printf("\n");
            flint_printf("u = "); _fmpq_vec_print(u, d); flint_printf("\n");
            flint_printf("v = "); _fmpq_vec_print(v, d); flint_printf("\n");
            flint_printf("nu = "); acb_printd(nu, 10); flint_printf("\n");
            flint_abort();
        }

        /* accuracy, away from the pole nu = d, v in Z^d */
        if (ok1 && acb_is_finite(z1) && acb_rel_accuracy_bits(z1) < prec - 8
            && !acb_contains_zero(z1))
        {
            flint_printf("FAIL (accuracy %wd, prec %wd)\n", acb_rel_accuracy_bits(z1), prec);
            flint_printf("G = "); fmpq_mat_print(G); flint_printf("\n");
            flint_printf("nu = "); acb_printd(nu, 10); flint_printf("\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(U);
        fmpq_mat_clear(G);
        fmpq_mat_clear(G2);
        fmpq_mat_clear(Uq);
        fmpq_mat_clear(UqT);
        fmpq_mat_clear(Tm);
        _fmpq_vec_clear(u, d);
        _fmpq_vec_clear(v, d);
        _fmpq_vec_clear(u2, d);
        _fmpq_vec_clear(v2, d);
        arb_mat_clear(A);
        _arb_vec_clear(x, d);
        _arb_vec_clear(y, d);
        acb_clear(nu);
        acb_clear(z1);
        acb_clear(z2);
        fmpq_clear(c);
    }

    TEST_FUNCTION_END(state);
}
