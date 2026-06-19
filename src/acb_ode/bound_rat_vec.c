#include "acb_ode.h"
#include "mag.h"


static void
acb_ode_bound_rat_exceptional_vec(mag_ptr res, slong len,
                                   const acb_ode_group_t group,
                                   const acb_ode_stairs_t stairs, slong n)
{
    slong s = 0;
    while (s < group->nshifts && group->shifts[s].n < n)
        s++;
    if (s < group->nshifts)
        for (slong i = 0; i < len; i++)
            mag_set(res + i, stairs->h + i * group->nshifts + s);
    else
        for (slong i = 0; i < len; i++)
            mag_zero(res + i);
}

/* [M19], Algorithm 7.4, steps 2 to 4, with the slight refinement that
 * τ(n₀) = ord can differ from sum(m(k), k=0..n₀). */

void
acb_ode_bound_rat_vec(mag_ptr res,
                      const acb_poly_struct * num, slong len,
                      const acb_poly_t ind, const acb_ode_group_t group,
                      const acb_ode_ind_lbound_t ind_lbound,
                      const acb_ode_stairs_t stairs,
                      slong n0, slong ord, slong prec)
{
    /* Bound at indices >= (smallest exceptional index > n0, if any). This term
     * does not depend on ord (possible values are precomputed in stairs using
     * worst-case bounds for the degrees in log(x) of the solutions). */
    acb_ode_bound_rat_exceptional_vec(res, len, group, stairs, n0);

    for (slong s = 0; s < group->nshifts; s++)
        if (group->shifts[s].n == n0)
            return;

    mag_ptr ordinary = _mag_vec_init(len);

    /* Bound at indices up < (smallest exceptional index > n0, if any).
     * Past the last exceptional index, this is the only contribution. */
    acb_ode_bound_rat_ordinary_vec(ordinary, num, len, ind, ind_lbound, n0, ord,
                                   prec);

    for (slong i = 0; i < len; i++)
        mag_max(res + i, res + i, ordinary + i);

    _mag_vec_clear(ordinary, len);
}

