/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-0f1.c"
#include "t-2f1.c"
#include "t-2f1_continuation.c"
#include "t-airy_bound.c"
#include "t-airy.c"
#include "t-airy_series.c"
#include "t-bessel_i.c"
#include "t-bessel_j.c"
#include "t-bessel_k.c"
#include "t-bessel_y.c"
#include "t-beta_lower.c"
#include "t-beta_lower_series.c"
#include "t-chebyshev_t.c"
#include "t-chebyshev_u.c"
#include "t-chi.c"
#include "t-chi_series.c"
#include "t-ci.c"
#include "t-ci_series.c"
#include "t-coulomb.c"
#include "t-coulomb_series.c"
#include "t-dilog.c"
#include "t-ei.c"
#include "t-ei_series.c"
#include "t-erf.c"
#include "t-erfc.c"
#include "t-erfc_series.c"
#include "t-erfi_series.c"
#include "t-erf_series.c"
#include "t-fresnel.c"
#include "t-fresnel_series.c"
#include "t-gamma_lower.c"
#include "t-gamma_lower_series.c"
#include "t-gamma_stirling_sum.c"
#include "t-gamma_taylor.c"
#include "t-gamma_upper.c"
#include "t-gamma_upper_series.c"
#include "t-gegenbauer_c.c"
#include "t-hermite_h.c"
#include "t-jacobi_p.c"
#include "t-laguerre_l.c"
#include "t-legendre_p.c"
#include "t-legendre_q.c"
#include "t-lgamma.c"
#include "t-li_series.c"
#include "t-log_rising_ui.c"
#include "t-log_rising_ui_jet.c"
#include "t-m.c"
#include "t-pfq.c"
#include "t-pfq_series_direct.c"
#include "t-pfq_series_sum_bs.c"
#include "t-pfq_series_sum_rs.c"
#include "t-pfq_sum_bs.c"
#include "t-pfq_sum_fme.c"
#include "t-pfq_sum_invz.c"
#include "t-pfq_sum_rs.c"
#include "t-rising_ui.c"
#include "t-rising_ui_jet.c"
#include "t-shi_series.c"
#include "t-si.c"
#include "t-si_series.c"
#include "t-spherical_y.c"
#include "t-u_asymp.c"
#include "t-u.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_hypgeom_0f1),
    TEST_FUNCTION(acb_hypgeom_2f1),
    TEST_FUNCTION(acb_hypgeom_2f1_continuation),
    TEST_FUNCTION(acb_hypgeom_airy_bound),
    TEST_FUNCTION(acb_hypgeom_airy),
    TEST_FUNCTION(acb_hypgeom_airy_series),
    TEST_FUNCTION(acb_hypgeom_bessel_i),
    TEST_FUNCTION(acb_hypgeom_bessel_j),
    TEST_FUNCTION(acb_hypgeom_bessel_k),
    TEST_FUNCTION(acb_hypgeom_bessel_y),
    TEST_FUNCTION(acb_hypgeom_beta_lower),
    TEST_FUNCTION(acb_hypgeom_beta_lower_series),
    TEST_FUNCTION(acb_hypgeom_chebyshev_t),
    TEST_FUNCTION(acb_hypgeom_chebyshev_u),
    TEST_FUNCTION(acb_hypgeom_chi),
    TEST_FUNCTION(acb_hypgeom_chi_series),
    TEST_FUNCTION(acb_hypgeom_ci),
    TEST_FUNCTION(acb_hypgeom_ci_series),
    TEST_FUNCTION(acb_hypgeom_coulomb),
    TEST_FUNCTION(acb_hypgeom_coulomb_series),
    TEST_FUNCTION(acb_hypgeom_dilog),
    TEST_FUNCTION(acb_hypgeom_ei),
    TEST_FUNCTION(acb_hypgeom_ei_series),
    TEST_FUNCTION(acb_hypgeom_erf),
    TEST_FUNCTION(acb_hypgeom_erfc),
    TEST_FUNCTION(acb_hypgeom_erfc_series),
    TEST_FUNCTION(acb_hypgeom_erfi_series),
    TEST_FUNCTION(acb_hypgeom_erf_series),
    TEST_FUNCTION(acb_hypgeom_fresnel),
    TEST_FUNCTION(acb_hypgeom_fresnel_series),
    TEST_FUNCTION(acb_hypgeom_gamma_lower),
    TEST_FUNCTION(acb_hypgeom_gamma_lower_series),
    TEST_FUNCTION(acb_hypgeom_gamma_stirling_sum),
    TEST_FUNCTION(acb_hypgeom_gamma_taylor),
    TEST_FUNCTION(acb_hypgeom_gamma_upper),
    TEST_FUNCTION(acb_hypgeom_gamma_upper_series),
    TEST_FUNCTION(acb_hypgeom_gegenbauer_c),
    TEST_FUNCTION(acb_hypgeom_hermite_h),
    TEST_FUNCTION(acb_hypgeom_jacobi_p),
    TEST_FUNCTION(acb_hypgeom_laguerre_l),
    TEST_FUNCTION(acb_hypgeom_legendre_p),
    TEST_FUNCTION(acb_hypgeom_legendre_q),
    TEST_FUNCTION(acb_hypgeom_lgamma),
    TEST_FUNCTION(acb_hypgeom_li_series),
    TEST_FUNCTION(acb_hypgeom_log_rising_ui),
    TEST_FUNCTION(acb_hypgeom_log_rising_ui_jet),
    TEST_FUNCTION(acb_hypgeom_m),
    TEST_FUNCTION(acb_hypgeom_pfq),
    TEST_FUNCTION(acb_hypgeom_pfq_series_direct),
    TEST_FUNCTION(acb_hypgeom_pfq_series_sum_bs),
    TEST_FUNCTION(acb_hypgeom_pfq_series_sum_rs),
    TEST_FUNCTION(acb_hypgeom_pfq_sum_bs),
    TEST_FUNCTION(acb_hypgeom_pfq_sum_fme),
    TEST_FUNCTION(acb_hypgeom_pfq_sum_invz),
    TEST_FUNCTION(acb_hypgeom_pfq_sum_rs),
    TEST_FUNCTION(acb_hypgeom_rising_ui),
    TEST_FUNCTION(acb_hypgeom_rising_ui_jet),
    TEST_FUNCTION(acb_hypgeom_shi_series),
    TEST_FUNCTION(acb_hypgeom_si),
    TEST_FUNCTION(acb_hypgeom_si_series),
    TEST_FUNCTION(acb_hypgeom_spherical_y),
    TEST_FUNCTION(acb_hypgeom_u_asymp),
    TEST_FUNCTION(acb_hypgeom_u)
};

/* main function *************************************************************/

TEST_MAIN(tests)
