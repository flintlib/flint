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

#include "t-addmod.c"
#include "t-cbrt_binary_search.c"
#include "t-cbrt.c"
#include "t-cbrt_chebyshev_approx.c"
#include "t-cbrt_newton_iteration.c"
#include "t-cbrtrem.c"
#include "t-clog_2exp.c"
#include "t-clog.c"
#include "t-compute_primes.c"
#include "t-CRT.c"
#include "t-discrete_log_bsgs.c"
#include "t-div2_preinv.c"
#include "t-divides.c"
#include "t-divrem2_precomp.c"
#include "t-divrem2_preinv.c"
#include "t-euler_phi.c"
#include "t-factor.c"
#include "t-factor_ecm.c"
#include "t-factorial_fast_mod2_preinv.c"
#include "t-factorial_mod2_preinv.c"
#include "t-factor_lehman.c"
#include "t-factor_one_line.c"
#include "t-factor_partial.c"
#include "t-factor_pollard_brent.c"
#include "t-factor_power235.c"
#include "t-factor_pp1.c"
#include "t-factor_SQUFOF.c"
#include "t-factor_trial.c"
#include "t-factor_trial_partial.c"
#include "t-factor_trial_range.c"
#include "t-flog.c"
#include "t-gcd.c"
#include "t-gcdinv.c"
#include "t-invmod.c"
#include "t-is_oddprime_binary.c"
#include "t-is_oddprime_small.c"
#include "t-is_perfect_power235.c"
#include "t-is_perfect_power.c"
#include "t-is_prime.c"
#include "t-is_prime_pocklington.c"
#include "t-is_prime_pseudosquare.c"
#include "t-is_probabprime_BPSW.c"
#include "t-is_probabprime.c"
#include "t-is_probabprime_fermat.c"
#include "t-is_probabprime_fibonacci.c"
#include "t-is_probabprime_lucas.c"
#include "t-is_square.c"
#include "t-is_squarefree.c"
#include "t-is_strong_probabprime2_preinv.c"
#include "t-is_strong_probabprime_precomp.c"
#include "t-jacobi.c"
#include "t-lll_mod_preinv.c"
#include "t-ll_mod_preinv.c"
#include "t-mod2_precomp.c"
#include "t-mod2_preinv.c"
#include "t-mod_precomp.c"
#include "t-moebius_mu.c"
#include "t-mulmod2.c"
#include "t-mulmod2_preinv.c"
#include "t-mulmod_precomp.c"
#include "t-mulmod_preinv.c"
#include "t-mulmod_shoup.c"
#include "t-nextprime.c"
#include "t-nth_prime_bounds.c"
#include "t-urandint.c"
#include "t-pow.c"
#include "t-powmod2.c"
#include "t-powmod2_preinv.c"
#include "t-powmod2_ui_preinv.c"
#include "t-powmod.c"
#include "t-powmod_precomp.c"
#include "t-powmod_ui_precomp.c"
#include "t-powmod_ui_preinv.c"
#include "t-preinvert_limb_prenorm.c"
#include "t-prime_pi_bounds.c"
#include "t-prime_pi.c"
#include "t-primes.c"
#include "t-primes_jump_after.c"
#include "t-primitive_root_prime.c"
#include "t-remove2_precomp.c"
#include "t-remove.c"
#include "t-revbin.c"
#include "t-root.c"
#include "t-rootrem.c"
#include "t-sizeinbase.c"
#include "t-sqrt.c"
#include "t-sqrtmod.c"
#include "t-sqrtmodn.c"
#include "t-sqrtmod_primepow.c"
#include "t-sqrtrem.c"
#include "t-submod.c"
#include "t-xgcd.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(n_addmod),
    TEST_FUNCTION(n_cbrt_binary_search),
    TEST_FUNCTION(n_cbrt),
    TEST_FUNCTION(n_cbrt_chebyshev_approx),
    TEST_FUNCTION(n_cbrt_newton_iteration),
    TEST_FUNCTION(n_cbrtrem),
    TEST_FUNCTION(n_clog_2exp),
    TEST_FUNCTION(n_clog),
    TEST_FUNCTION(compute_primes),
    TEST_FUNCTION(n_CRT),
    TEST_FUNCTION(n_discrete_log_bsgs),
    TEST_FUNCTION(n_div2_preinv),
    TEST_FUNCTION(n_divides),
    TEST_FUNCTION(n_divrem2_precomp),
    TEST_FUNCTION(n_divrem2_preinv),
    TEST_FUNCTION(n_euler_phi),
    TEST_FUNCTION(n_factor),
    TEST_FUNCTION(n_factor_ecm),
    TEST_FUNCTION(n_factorial_fast_mod2_preinv),
    TEST_FUNCTION(n_factorial_mod2_preinv),
    TEST_FUNCTION(n_factor_lehman),
    TEST_FUNCTION(n_factor_one_line),
    TEST_FUNCTION(n_factor_partial),
    TEST_FUNCTION(n_factor_pollard_brent),
    TEST_FUNCTION(n_factor_power235),
    TEST_FUNCTION(n_factor_pp1),
    TEST_FUNCTION(n_factor_SQUFOF),
    TEST_FUNCTION(n_factor_trial),
    TEST_FUNCTION(n_factor_trial_partial),
    TEST_FUNCTION(n_factor_trial_range),
    TEST_FUNCTION(n_flog),
    TEST_FUNCTION(n_gcd),
    TEST_FUNCTION(n_gcdinv),
    TEST_FUNCTION(n_invmod),
    TEST_FUNCTION(n_is_oddprime_binary),
    TEST_FUNCTION(n_is_oddprime_small),
    TEST_FUNCTION(n_is_perfect_power235),
    TEST_FUNCTION(n_is_perfect_power),
    TEST_FUNCTION(n_is_prime),
    TEST_FUNCTION(n_is_prime_pocklington),
    TEST_FUNCTION(n_is_prime_pseudosquare),
    TEST_FUNCTION(n_is_probabprime_BPSW),
    TEST_FUNCTION(n_is_probabprime),
    TEST_FUNCTION(n_is_probabprime_fermat),
    TEST_FUNCTION(n_is_probabprime_fibonacci),
    TEST_FUNCTION(n_is_probabprime_lucas),
    TEST_FUNCTION(n_is_square),
    TEST_FUNCTION(n_is_squarefree),
    TEST_FUNCTION(n_is_strong_probabprime2_preinv),
    TEST_FUNCTION(n_is_strong_probabprime_precomp),
    TEST_FUNCTION(n_jacobi),
    TEST_FUNCTION(n_lll_mod_preinv),
    TEST_FUNCTION(n_ll_mod_preinv),
    TEST_FUNCTION(n_mod2_precomp),
    TEST_FUNCTION(n_mod2_preinv),
    TEST_FUNCTION(n_mod_precomp),
    TEST_FUNCTION(n_moebius_mu),
    TEST_FUNCTION(n_mulmod2),
    TEST_FUNCTION(n_mulmod2_preinv),
    TEST_FUNCTION(n_mulmod_precomp),
    TEST_FUNCTION(n_mulmod_preinv),
    TEST_FUNCTION(n_mulmod_shoup),
    TEST_FUNCTION(n_nextprime),
    TEST_FUNCTION(n_nth_prime_bounds),
    TEST_FUNCTION(n_urandint),
    TEST_FUNCTION(n_pow),
    TEST_FUNCTION(n_powmod2),
    TEST_FUNCTION(n_powmod2_preinv),
    TEST_FUNCTION(n_powmod2_ui_preinv),
    TEST_FUNCTION(n_powmod),
    TEST_FUNCTION(n_powmod_precomp),
    TEST_FUNCTION(n_powmod_ui_precomp),
    TEST_FUNCTION(n_powmod_ui_preinv),
    TEST_FUNCTION(n_preinvert_limb_prenorm),
    TEST_FUNCTION(n_prime_pi_bounds),
    TEST_FUNCTION(n_prime_pi),
    TEST_FUNCTION(n_primes),
    TEST_FUNCTION(n_primes_jump_after),
    TEST_FUNCTION(n_primitive_root_prime),
    TEST_FUNCTION(n_remove2_precomp),
    TEST_FUNCTION(n_remove),
    TEST_FUNCTION(n_revbin),
    TEST_FUNCTION(n_root),
    TEST_FUNCTION(n_rootrem),
    TEST_FUNCTION(n_sizeinbase),
    TEST_FUNCTION(n_sqrt),
    TEST_FUNCTION(n_sqrtmod),
    TEST_FUNCTION(n_sqrtmodn),
    TEST_FUNCTION(n_sqrtmod_primepow),
    TEST_FUNCTION(n_sqrtrem),
    TEST_FUNCTION(n_submod),
    TEST_FUNCTION(n_xgcd)
};

/* main function *************************************************************/

TEST_MAIN(tests)
