#include "acb_ode.h"
#include "acb_poly.h"

void
acb_ode_group_bound_init(acb_ode_group_bound_t gbound)
{
    acb_poly_init(gbound->ind);
    acb_ode_ind_lbound_init(gbound->ind_lbound);
    acb_ode_stairs_init(gbound->stairs);
}

void
acb_ode_group_bound_clear(acb_ode_group_bound_t gbound)
{
    acb_poly_clear(gbound->ind);
    acb_ode_ind_lbound_clear(gbound->ind_lbound);
    acb_ode_stairs_clear(gbound->stairs);
}

/* todo: incremental version for increased pol_part_len (and maybe prec) but
 * same group */
void
acb_ode_group_bound_precompute(acb_ode_group_bound_t gbound,
                               const acb_poly_struct * dop, slong dop_len,
                               const acb_ode_exponents_t expos,
                               slong grp,  // group index in expos
                               const acb_ode_bound_t bound,
                               slong prec)
{
    /* Monic indicial polynomial, centered at the current group leader */
    acb_poly_struct * ind = bound->all_nums;
    acb_poly_taylor_shift(gbound->ind, ind, expos->groups[grp].leader, prec);
    FLINT_ASSERT(expos->groups[grp].shifts[0].n == 0);
    _acb_vec_zero(gbound->ind->coeffs, expos->groups[grp].shifts[0].mult);
    FLINT_ASSERT(gbound->ind->length == dop_len);
    _acb_vec_scalar_div(gbound->ind->coeffs, gbound->ind->coeffs, dop_len - 1,
                        gbound->ind->coeffs + dop_len - 1, prec);
    acb_one(gbound->ind->coeffs + dop_len - 1);

    /* flint_printf("ind=%{acb_poly}\n", gbound->ind); */

    /* Rat(n) denominator data */

    acb_ode_ind_lbound_precompute(gbound->ind_lbound, expos, grp, prec);

    /* Ignore the first entry of all_nums, anticipating the factor w⁻¹ under the
     * integral in Algo 6.11, step 3 of [M19] */
    acb_ode_stairs_precompute(gbound->stairs,
            bound->all_nums + 1, bound->all_nums_len - 1, gbound->ind,
            expos->groups + grp, gbound->ind_lbound, prec);

    /* flint_printf("stairs = %{mag*}\n", gbound->stairs->h, gbound->stairs->length); */
}

