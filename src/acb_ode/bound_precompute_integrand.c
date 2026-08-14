#include "arf.h"
#include "acb_ode.h"
#include "arb_poly.h"


/* [M19], Algorithm 6.1, step 4. */

void
acb_ode_bound_precompute_integrand(arb_poly_t itg_pol, arb_poly_t itg_num,
                                   const acb_ode_group_t group,
                                   const acb_ode_bound_t bound,
                                   const acb_ode_group_bound_t gbound,
                                   slong n0, slong nlogs, slong prec)
{
    mag_ptr maj_all_nums = _mag_vec_init(bound->all_nums_len);

    /* In view of [M19], Eq. (5.9), one can take here for nlogs the actual
       length wrt log(x) of the solutions of interest truncated after their term
       in x^n, instead of τ(n0) as Algo 6.1, step 3 would prescribe. */
    acb_ode_bound_rat_vec(maj_all_nums, bound->all_nums + 1, bound->all_nums_len - 1,
                           gbound->ind, group, gbound->ind_lbound, gbound->stairs,
                           n0, nlogs, prec);

    // flint_printf("n0=%ld maj_all_nums = %{mag*}\n", n0, maj_all_nums, bound->all_nums_len - 1);

    arb_poly_fit_length(itg_pol, bound->pol_part_len - 1);
    for (slong i = 0; i < bound->pol_part_len - 1; i++)
        arf_set_mag(arb_midref(itg_pol->coeffs + i), maj_all_nums + i);
    _arb_poly_set_length(itg_pol, bound->pol_part_len - 1);
    _arb_poly_normalise(itg_pol);

    slong num_len = bound->all_nums_len - bound->pol_part_len;
    arb_poly_fit_length(itg_num, num_len);
    _arb_poly_set_length(itg_num, num_len);
    for (slong i = 0; i < num_len; i++)
        arf_set_mag(arb_midref(itg_num->coeffs + i),
                    maj_all_nums + bound->pol_part_len - 1 + i);
    _arb_poly_normalise(itg_num);

    // flint_printf("n0=%ld itg_pol=%{arb_poly} itg_num=%{arb_poly}\n", n0, itg_pol, itg_num);

    _mag_vec_clear(maj_all_nums, bound->all_nums_len);
}
