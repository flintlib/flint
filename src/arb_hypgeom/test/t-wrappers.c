/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

void
TEST(const arb_t x1, const arb_t x2, const char * s)
{
    if (!arb_overlaps(x1, x2))
    {
        flint_printf("FAIL: %s\n", s);
        arb_printn(x1, 30, 0); printf("\n\n");
        arb_printn(x2, 30, 0); printf("\n\n");
        flint_abort();
    }
}

TEST_FUNCTION_START(arb_hypgeom_wrappers, state)
{
    {
        arb_t a, b, c, d, z, r, u, v;
        slong prec;

        prec = 53;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(z);
        arb_init(r);
        arb_init(u);
        arb_init(v);

        arb_set_d(a, 0.5);
        arb_set_d(b, 0.25);
        arb_set_d(c, 0.125);
        arb_set_d(z, 0.75);

        arb_hypgeom_pfq(r, a, 1, b, 1, z, 0, prec);
        arb_set_str(v, "[3.40725820729608 +/- 6.88e-15]", prec);
        TEST(r, v, "pfq");

        arb_hypgeom_pfq(r, a, 1, b, 1, z, 1, prec);
        arb_set_str(v, "[0.93977518087904 +/- 3.26e-15]", prec);
        TEST(r, v, "pfq regularized");

        arb_hypgeom_0f1(r, a, z, 0, prec);
        arb_set_str(v, "[2.91457744017593 +/- 4.17e-15]", prec);
        TEST(r, v, "0f1");

        arb_hypgeom_0f1(r, a, z, 1, prec);
        arb_set_str(v, "[1.64437423219054 +/- 3.77e-15]", prec);
        TEST(r, v, "0f1 regularized");

        arb_hypgeom_1f1(r, a, b, z, 0, prec);
        arb_set_str(v, "[3.40725820729608 +/- 6.88e-15]", prec);
        TEST(r, v, "1f1");

        arb_hypgeom_1f1(r, a, b, z, 1, prec);
        arb_set_str(v, "[0.93977518087904 +/- 3.26e-15]", prec);
        TEST(r, v, "1f1 regularized");

        arb_hypgeom_u(r, a, b, z, prec);
        arb_set_str(v, "[0.7761320390950 +/- 3.89e-14]", prec);
        TEST(r, v, "u");

        arb_hypgeom_2f1(r, a, b, c, z, 0, prec);
        arb_set_str(v, "[3.2569394189980 +/- 8.18e-14]", prec);
        TEST(r, v, "2f1");

        arb_hypgeom_2f1(r, a, b, c, z, 1, prec);
        arb_set_str(v, "[0.4323021855542 +/- 1.67e-14]", prec);
        TEST(r, v, "2f1 regularized");

        arb_hypgeom_erf(r, z, prec);
        arb_set_str(v, "[0.711155633653515 +/- 2.35e-16]", prec);
        TEST(r, v, "erf");

        arb_hypgeom_erfc(r, z, prec);
        arb_set_str(v, "[0.288844366346485 +/- 2.34e-16]", prec);
        TEST(r, v, "erfc");

        arb_hypgeom_erfi(r, z, prec);
        arb_set_str(v, "[1.03575728441196 +/- 3.17e-15]", prec);
        TEST(r, v, "erfi");

        arb_hypgeom_fresnel(r, NULL, z, 0, prec);
        arb_set_str(v, "[0.137478632382610 +/- 1.79e-16]", prec);
        TEST(r, v, "fresnel_s");

        arb_hypgeom_fresnel(r, NULL, z, 1, prec);
        arb_set_str(v, "[0.208877111233384 +/- 4.63e-16]", prec);
        TEST(r, v, "fresnel_s normalized");

        arb_hypgeom_fresnel(NULL, r, z, 0, prec);
        arb_set_str(v, "[0.726614618304550 +/- 2.11e-16]", prec);
        TEST(r, v, "fresnel_c");

        arb_hypgeom_fresnel(NULL, r, z, 1, prec);
        arb_set_str(v, "[0.693525990787136 +/- 2.76e-16]", prec);
        TEST(r, v, "fresnel_c normalized");

        arb_hypgeom_gamma_lower(r, a, z, 0, prec);
        arb_set_str(v, "[1.38132404568612 +/- 3.66e-15]", prec);
        TEST(r, v, "gamma_lower");

        arb_hypgeom_gamma_lower(r, a, z, 1, prec);
        arb_set_str(v, "[0.77932863808015 +/- 4.17e-15]", prec);
        TEST(r, v, "gamma_lower 1");

        arb_hypgeom_gamma_lower(r, a, z, 2, prec);
        arb_set_str(v, "[0.89989119796552 +/- 2.76e-15]", prec);
        TEST(r, v, "gamma_lower 2");

        arb_hypgeom_gamma_upper(r, a, z, 0, prec);
        arb_set_str(v, "[0.39112980521940 +/- 4.21e-15]", prec);
        TEST(r, v, "gamma_upper");

        arb_hypgeom_gamma_upper(r, a, z, 1, prec);
        arb_set_str(v, "[0.22067136191985 +/- 4.08e-15]", prec);
        TEST(r, v, "gamma_upper 1");

        arb_hypgeom_beta_lower(r, a, b, z, 0, prec);
        arb_set_str(v, "[2.3363313667086 +/- 3.14e-14]", prec);
        TEST(r, v, "beta_lower");

        arb_hypgeom_beta_lower(r, a, b, z, 1, prec);
        arb_set_str(v, "[0.44551489018313 +/- 6.22e-15]", prec);
        TEST(r, v, "beta_lower 1");

        arb_hypgeom_expint(r, a, z, prec);
        arb_set_str(v, "[0.45163779666301 +/- 3.10e-15]", prec);
        TEST(r, v, "expint");

        arb_hypgeom_ei(r, z, prec);
        arb_set_str(v, "[1.20733281600122 +/- 2.72e-15]", prec);
        TEST(r, v, "ei");

        arb_hypgeom_si(r, z, prec);
        arb_set_str(v, "[0.72695424715009 +/- 3.99e-15]", prec);
        TEST(r, v, "si");

        arb_hypgeom_ci(r, z, prec);
        arb_set_str(v, "[0.152163600980330 +/- 5.49e-16]", prec);
        TEST(r, v, "ci");

        arb_hypgeom_shi(r, z, prec);
        arb_set_str(v, "[0.77383681445623 +/- 5.97e-15]", prec);
        TEST(r, v, "shi");

        arb_hypgeom_chi(r, z, prec);
        arb_set_str(v, "[0.433496001544996 +/- 7.91e-16]", prec);
        TEST(r, v, "chi");

        arb_hypgeom_li(r, z, 0, prec);
        arb_set_str(v, "[-0.93693001101265 +/- 5.67e-15]", prec);
        TEST(r, v, "li");

        arb_hypgeom_li(r, z, 1, prec);
        arb_set_str(v, "[-1.98209379113014 +/- 3.56e-15]", prec);
        TEST(r, v, "li offset");

        arb_hypgeom_bessel_j(r, a, z, prec);
        arb_set_str(v, "[0.62800587637589 +/- 4.06e-15]", prec);
        TEST(r, v, "bessel_j");

        arb_hypgeom_bessel_y(r, a, z, prec);
        arb_set_str(v, "[-0.67411792914454 +/- 5.27e-15]", prec);
        TEST(r, v, "bessel_y");

        arb_hypgeom_bessel_jy(r, u, a, z, prec);
        arb_set_str(v, "[0.62800587637589 +/- 4.06e-15]", prec);
        TEST(r, v, "bessel_jy");
        arb_set_str(v, "[-0.67411792914454 +/- 5.27e-15]", prec);
        TEST(u, v, "bessel_jy");

        arb_hypgeom_bessel_i(r, a, z, prec);
        arb_set_str(v, "[0.75761498638991 +/- 4.63e-15]", prec);
        TEST(r, v, "bessel_i");

        arb_hypgeom_bessel_i_scaled(r, a, z, prec);
        arb_set_str(v, "[0.357871979425934 +/- 8.04e-16]", prec);
        TEST(r, v, "bessel_i");

        arb_hypgeom_bessel_k(r, a, z, prec);
        arb_set_str(v, "[0.68361006034952 +/- 7.89e-15]", prec);
        TEST(r, v, "bessel_k");

        arb_hypgeom_bessel_k_scaled(r, a, z, prec);
        arb_set_str(v, "[1.4472025091165 +/- 4.23e-14]", prec);
        TEST(r, v, "bessel_k");

        arb_hypgeom_airy(r, NULL, NULL, NULL, z, prec);
        arb_set_str(v, "[0.179336305478645 +/- 2.36e-16]", prec);
        TEST(r, v, "airy ai");

        arb_hypgeom_airy(NULL, r, NULL, NULL, z, prec);
        arb_set_str(v, "[-0.193175208104376 +/- 4.67e-16]", prec);
        TEST(r, v, "airy ai'");

        arb_hypgeom_airy(NULL, NULL, r, NULL, z, prec);
        arb_set_str(v, "[1.00693090863322 +/- 3.90e-15]", prec);
        TEST(r, v, "airy bi");

        arb_hypgeom_airy(NULL, NULL, NULL, r, z, prec);
        arb_set_str(v, "[0.690299702736886 +/- 2.34e-16]", prec);
        TEST(r, v, "airy bi'");

        arb_hypgeom_coulomb(r, NULL, a, b, z, prec);
        arb_set_str(v, "[0.281871468006603 +/- 4.68e-16]", prec);
        TEST(r, v, "coulomb f");

        arb_hypgeom_coulomb(NULL, r, a, b, z, prec);
        arb_set_str(v, "[1.38897454984644 +/- 2.23e-15]", prec);
        TEST(r, v, "coulomb g");

        arb_hypgeom_chebyshev_t(r, a, z, prec);
        arb_set_str(v, "[0.935414346693485 +/- 7.02e-16]", prec);
        TEST(r, v, "chebyshev_t");

        arb_hypgeom_chebyshev_u(r, a, z, prec);
        arb_set_str(v, "[1.33630620956212 +/- 3.41e-15]", prec);
        TEST(r, v, "chebyshev_u");

        arb_hypgeom_jacobi_p(r, a, b, c, z, prec);
        arb_set_str(v, "[1.03224258095454 +/- 6.54e-15]", prec);
        TEST(r, v, "jacobi_p");

        arb_hypgeom_gegenbauer_c(r, a, b, z, prec);
        arb_set_str(v, "[0.58153237382485 +/- 3.74e-15]", prec);
        TEST(r, v, "gegenbauer_c");

        arb_hypgeom_laguerre_l(r, a, b, z, prec);
        arb_set_str(v, "[0.76858987769621 +/- 3.19e-15]", prec);
        TEST(r, v, "laguerre_l");

        arb_hypgeom_hermite_h(r, a, z, prec);
        arb_set_str(v, "[1.31600493243946 +/- 6.52e-15]", prec);
        TEST(r, v, "hermite_h");

        arb_hypgeom_legendre_p(r, a, b, z, 0, prec);
        arb_set_str(v, "[0.90435295129350 +/- 5.79e-15]", prec);
        TEST(r, v, "legendre_q");

        arb_hypgeom_legendre_q(r, a, b, z, 0, prec);
        arb_set_str(v, "[-0.3763617490859 +/- 3.83e-14]", prec);
        TEST(r, v, "legendre_q");

        arb_hypgeom_dilog(r, a, prec);
        arb_set_str(v, "[0.582240526465012 +/- 6.18e-16]", prec);
        TEST(r, v, "dilog");

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arb_clear(z);
        arb_clear(r);
        arb_clear(u);
        arb_clear(v);
    }

    TEST_FUNCTION_END(state);
}
