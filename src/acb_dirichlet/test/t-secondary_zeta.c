/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

void _acb_dirichlet_secondary_zeta_term_E(acb_t res, const acb_t s, const arb_t a, slong prec);
void _acb_dirichlet_secondary_zeta_term_S(acb_t res, const acb_t s, const arb_t a, slong prec);
void _acb_dirichlet_secondary_zeta_term_P(acb_t res, const acb_t s, const arb_t a, slong prec);
void _acb_dirichlet_secondary_zeta_term_A(acb_t res, const acb_t s, const arb_t a, slong prec);

TEST_FUNCTION_START(acb_dirichlet_secondary_zeta, state)
{
    acb_t s, A, P, E, S, Z, Z2;
    arb_t a;
    slong i, j, prec;
    acb_init(s); acb_init(A); acb_init(P); acb_init(E); acb_init(S);
    acb_init(Z); acb_init(Z2), arb_init(a);

    struct { double re, im; const char *zr, *zi; } C[] = {
        {2.0, 0.0, "[0.02310499311541897078893381043033901400338176039742209012 +/- 1e-50]", "0"},
        {50.0, 0.0, "[3.059367237705770705860511544941278451436737499249456979e-58 +/- 1e-108]", "0"},
        {-2.0, 10, "[188.643888367473947954809548622456913341026479844696602 +/- 1e-50]",
                     "-23.92543799796535790705975992188910266801550773379454949 +/- 1e-50"},
        {0.0, 0.0, "0.875", "0"},
        {-2.0, 0.0, "-0.28125", "0"},
        {-4.0, 0.0, "0.0234375", "0"},
        {-20.0, 0.0, "-44151686.21987307071685791015625", "0"},
    };

    double atest[] = { 0.0, 0.1, 0.01, 0.001 };

    for (i = 0; i < 7; i++)
    {
        for (j = 0; j < 4; j++)
        {
            acb_set_d_d(s, C[i].re, C[i].im);
            arb_set_d(a, atest[j]);

            for (prec = 32; prec <= (flint_test_multiplier() > 1 ? 128 : 64); prec *= 2)
            {
                if (atest[j] == 0.0)
                {
                    acb_dirichlet_secondary_zeta(Z, s, prec);
                }
                else
                {
                    _acb_dirichlet_secondary_zeta_term_A(A, s, a, prec);
                    _acb_dirichlet_secondary_zeta_term_P(P, s, a, prec);
                    _acb_dirichlet_secondary_zeta_term_E(E, s, a, prec);
                    _acb_dirichlet_secondary_zeta_term_S(S, s, a, prec);

                    /* Z = A - P + E - S */
                    acb_sub(Z, A, P, prec);
                    acb_add(Z, Z, E, prec);
                    acb_sub(Z, Z, S, prec);
                }

                arb_set_str(acb_realref(Z2), C[i].zr, prec);
                arb_set_str(acb_imagref(Z2), C[i].zi, prec);

                if (!acb_overlaps(Z, Z2))
                {
                    flint_printf("FAIL: acb_dirichlet_secondary_zeta: overlap\n");
                    flint_printf("prec = %wd\n", prec);
                    flint_printf("s = "); acb_printn(s, 50, 0); flint_printf("\n");
                    flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n");
                    flint_printf("Z = "); acb_printn(Z, 50, 0); flint_printf("\n");
                    flint_printf("Z2 = "); acb_printn(Z2, 50, 0); flint_printf("\n");
                    flint_abort();
                }
            }
        }
    }

    acb_clear(s); acb_clear(A); acb_clear(P); acb_clear(E); acb_clear(S);
    acb_clear(Z); acb_clear(Z2), arb_clear(a);

    TEST_FUNCTION_END(state);
}

