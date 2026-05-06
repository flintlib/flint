#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"

/* todo: incremental version for increased pol_part_len (and maybe prec) but
 * same group */

void
acb_ode_bound_precompute_group(acb_ode_bound_t bound,
                               const acb_poly_struct * dop, slong dop_len,
                               const acb_ode_exponents_t expos,
                               slong grp,  // group index in expos
                               slong prec)
{
    /* flint_printf("== acb_ode_bound_precompute_group %{acb} ==\n", expos->groups[grp].leader); */

    /* Monic indicial polynomial, centered at the current group leader */
    acb_poly_struct * ind = bound->all_nums;
    acb_poly_taylor_shift(bound->ind, ind, expos->groups[grp].leader, prec);
    FLINT_ASSERT(expos->groups[grp].shifts[0].n == 0);
    _acb_vec_zero(bound->ind->coeffs, expos->groups[grp].shifts[0].mult);
    FLINT_ASSERT(bound->ind->length == dop_len);
    _acb_vec_scalar_div(bound->ind->coeffs, bound->ind->coeffs, dop_len - 1,
                        bound->ind->coeffs + dop_len - 1, prec);
    acb_one(bound->ind->coeffs + dop_len - 1);

    /* flint_printf("ind=%{acb_poly}\n", bound->ind); */

    /* Rat(n) denominator data */

    acb_ode_ind_lbound_precompute(bound->ind_lbound, expos, grp, prec);

    /* Ignore the first entry of all_nums, anticipating the factor w⁻¹ under the
     * integral in Algo 6.11, step 3 of [M19] */
    acb_ode_stairs_precompute(bound->stairs,
            bound->all_nums + 1, bound->all_nums_len - 1, bound->ind,
            expos->groups + grp, bound->ind_lbound, prec);

    /* flint_printf("stairs = %{mag*}\n", bound->stairs->h, bound->stairs->length); */
}

