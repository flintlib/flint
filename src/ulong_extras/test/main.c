/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>

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

int (*test_functions[])(void) =
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

char n_addmod_name[] = "n_addmod";
char n_cbrt_binary_search_name[] = "n_cbrt_binary_search";
char n_cbrt_name[] = "n_cbrt";
char n_cbrt_chebyshev_approx_name[] = "n_cbrt_chebyshev_approx";
char n_cbrt_newton_iteration_name[] = "n_cbrt_newton_iteration";
char n_cbrtrem_name[] = "n_cbrtrem";
char n_clog_2exp_name[] = "n_clog_2exp";
char n_clog_name[] = "n_clog";
char compute_primes_name[] = "compute_primes";
char n_CRT_name[] = "n_CRT";
char n_discrete_log_bsgs_name[] = "n_discrete_log_bsgs";
char n_div2_preinv_name[] = "n_div2_preinv";
char n_divides_name[] = "n_divides";
char n_divrem2_precomp_name[] = "n_divrem2_precomp";
char n_divrem2_preinv_name[] = "n_divrem2_preinv";
char n_euler_phi_name[] = "n_euler_phi";
char n_factor_name[] = "n_factor";
char n_factor_ecm_name[] = "n_factor_ecm";
char n_factorial_fast_mod2_preinv_name[] = "n_factorial_fast_mod2_preinv";
char n_factorial_mod2_preinv_name[] = "n_factorial_mod2_preinv";
char n_factor_lehman_name[] = "n_factor_lehman";
char n_factor_one_line_name[] = "n_factor_one_line";
char n_factor_partial_name[] = "n_factor_partial";
char n_factor_pollard_brent_name[] = "n_factor_pollard_brent";
char n_factor_power235_name[] = "n_factor_power235";
char n_factor_pp1_name[] = "n_factor_pp1";
char n_factor_SQUFOF_name[] = "n_factor_SQUFOF";
char n_factor_trial_name[] = "n_factor_trial";
char n_factor_trial_partial_name[] = "n_factor_trial_partial";
char n_factor_trial_range_name[] = "n_factor_trial_range";
char n_flog_name[] = "n_flog";
char n_gcd_name[] = "n_gcd";
char n_gcdinv_name[] = "n_gcdinv";
char n_invmod_name[] = "n_invmod";
char n_is_oddprime_binary_name[] = "n_is_oddprime_binary";
char n_is_oddprime_small_name[] = "n_is_oddprime_small";
char n_is_perfect_power235_name[] = "n_is_perfect_power235";
char n_is_perfect_power_name[] = "n_is_perfect_power";
char n_is_prime_name[] = "n_is_prime";
char n_is_prime_pocklington_name[] = "n_is_prime_pocklington";
char n_is_prime_pseudosquare_name[] = "n_is_prime_pseudosquare";
char n_is_probabprime_BPSW_name[] = "n_is_probabprime_BPSW";
char n_is_probabprime_name[] = "n_is_probabprime";
char n_is_probabprime_fermat_name[] = "n_is_probabprime_fermat";
char n_is_probabprime_fibonacci_name[] = "n_is_probabprime_fibonacci";
char n_is_probabprime_lucas_name[] = "n_is_probabprime_lucas";
char n_is_square_name[] = "n_is_square";
char n_is_squarefree_name[] = "n_is_squarefree";
char n_is_strong_probabprime2_preinv_name[] = "n_is_strong_probabprime2_preinv";
char n_is_strong_probabprime_precomp_name[] = "n_is_strong_probabprime_precomp";
char n_jacobi_name[] = "n_jacobi";
char n_lll_mod_preinv_name[] = "n_lll_mod_preinv";
char n_ll_mod_preinv_name[] = "n_ll_mod_preinv";
char n_mod2_precomp_name[] = "n_mod2_precomp";
char n_mod2_preinv_name[] = "n_mod2_preinv";
char n_mod_precomp_name[] = "n_mod_precomp";
char n_moebius_mu_name[] = "n_moebius_mu";
char n_mulmod2_name[] = "n_mulmod2";
char n_mulmod2_preinv_name[] = "n_mulmod2_preinv";
char n_mulmod_precomp_name[] = "n_mulmod_precomp";
char n_mulmod_preinv_name[] = "n_mulmod_preinv";
char n_mulmod_shoup_name[] = "n_mulmod_shoup";
char n_nextprime_name[] = "n_nextprime";
char n_nth_prime_bounds_name[] = "n_nth_prime_bounds";
char n_urandint_name[] = "n_urandint";
char n_pow_name[] = "n_pow";
char n_powmod2_name[] = "n_powmod2";
char n_powmod2_preinv_name[] = "n_powmod2_preinv";
char n_powmod2_ui_preinv_name[] = "n_powmod2_ui_preinv";
char n_powmod_name[] = "n_powmod";
char n_powmod_precomp_name[] = "n_powmod_precomp";
char n_powmod_ui_precomp_name[] = "n_powmod_ui_precomp";
char n_powmod_ui_preinv_name[] = "n_powmod_ui_preinv";
char n_prime_pi_bounds_name[] = "n_prime_pi_bounds";
char n_prime_pi_name[] = "n_prime_pi";
char n_primes_name[] = "n_primes";
char n_primes_jump_after_name[] = "n_primes_jump_after";
char n_primitive_root_prime_name[] = "n_primitive_root_prime";
char n_remove2_precomp_name[] = "n_remove2_precomp";
char n_remove_name[] = "n_remove";
char n_revbin_name[] = "n_revbin";
char n_root_name[] = "n_root";
char n_rootrem_name[] = "n_rootrem";
char n_sizeinbase_name[] = "n_sizeinbase";
char n_sqrt_name[] = "n_sqrt";
char n_sqrtmod_name[] = "n_sqrtmod";
char n_sqrtmodn_name[] = "n_sqrtmodn";
char n_sqrtmod_primepow_name[] = "n_sqrtmod_primepow";
char n_sqrtrem_name[] = "n_sqrtrem";
char n_submod_name[] = "n_submod";
char n_xgcd_name[] = "n_xgcd";

char * test_names[] =
{
    n_addmod_name,
    n_cbrt_binary_search_name,
    n_cbrt_name,
    n_cbrt_chebyshev_approx_name,
    n_cbrt_newton_iteration_name,
    n_cbrtrem_name,
    n_clog_2exp_name,
    n_clog_name,
    compute_primes_name,
    n_CRT_name,
    n_discrete_log_bsgs_name,
    n_div2_preinv_name,
    n_divides_name,
    n_divrem2_precomp_name,
    n_divrem2_preinv_name,
    n_euler_phi_name,
    n_factor_name,
    n_factor_ecm_name,
    n_factorial_fast_mod2_preinv_name,
    n_factorial_mod2_preinv_name,
    n_factor_lehman_name,
    n_factor_one_line_name,
    n_factor_partial_name,
    n_factor_pollard_brent_name,
    n_factor_power235_name,
    n_factor_pp1_name,
    n_factor_SQUFOF_name,
    n_factor_trial_name,
    n_factor_trial_partial_name,
    n_factor_trial_range_name,
    n_flog_name,
    n_gcd_name,
    n_gcdinv_name,
    n_invmod_name,
    n_is_oddprime_binary_name,
    n_is_oddprime_small_name,
    n_is_perfect_power235_name,
    n_is_perfect_power_name,
    n_is_prime_name,
    n_is_prime_pocklington_name,
    n_is_prime_pseudosquare_name,
    n_is_probabprime_BPSW_name,
    n_is_probabprime_name,
    n_is_probabprime_fermat_name,
    n_is_probabprime_fibonacci_name,
    n_is_probabprime_lucas_name,
    n_is_square_name,
    n_is_squarefree_name,
    n_is_strong_probabprime2_preinv_name,
    n_is_strong_probabprime_precomp_name,
    n_jacobi_name,
    n_lll_mod_preinv_name,
    n_ll_mod_preinv_name,
    n_mod2_precomp_name,
    n_mod2_preinv_name,
    n_mod_precomp_name,
    n_moebius_mu_name,
    n_mulmod2_name,
    n_mulmod2_preinv_name,
    n_mulmod_precomp_name,
    n_mulmod_preinv_name,
    n_mulmod_shoup_name,
    n_nextprime_name,
    n_nth_prime_bounds_name,
    n_urandint_name,
    n_pow_name,
    n_powmod2_name,
    n_powmod2_preinv_name,
    n_powmod2_ui_preinv_name,
    n_powmod_name,
    n_powmod_precomp_name,
    n_powmod_ui_precomp_name,
    n_powmod_ui_preinv_name,
    n_prime_pi_bounds_name,
    n_prime_pi_name,
    n_primes_name,
    n_primes_jump_after_name,
    n_primitive_root_prime_name,
    n_remove2_precomp_name,
    n_remove_name,
    n_revbin_name,
    n_root_name,
    n_rootrem_name,
    n_sizeinbase_name,
    n_sqrt_name,
    n_sqrtmod_name,
    n_sqrtmodn_name,
    n_sqrtmod_primepow_name,
    n_sqrtrem_name,
    n_submod_name,
    n_xgcd_name
};

/* main function *************************************************************/

#define NUMBER_OF_TESTS (sizeof(test_functions) / sizeof(int (*)(void)))

int
main(int argc, char ** argv)
{
    int ix, jx;

    if (argc < 2)
    {
        for (ix = 0; ix < NUMBER_OF_TESTS; ix++)
            if (test_functions[ix]())
                flint_abort();
    }
    else
    {
        for (ix = 1; ix < argc; ix++)
        {
            for (jx = 0; jx < NUMBER_OF_TESTS; jx++)
            {
                /* If argument equals to test name, run it */
                if (strcmp(argv[ix], test_names[jx]) == 0)
                {
                    if (test_functions[jx]())
                        flint_abort();
                    break;
                }
            }

            if (jx == NUMBER_OF_TESTS)
            {
                fprintf(stderr, "Error: Could not find test function for %s\n", argv[ix]);
                flint_abort();
            }
        }
    }

    return 0;
}
