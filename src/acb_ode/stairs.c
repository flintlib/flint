#include "acb_ode.h"
#include "mag.h"


void
acb_ode_stairs_init(acb_ode_stairs_t stairs)
{
    stairs->length = 0;
    stairs->h = NULL;
}

void
acb_ode_stairs_clear(acb_ode_stairs_t stairs)
{
    _mag_vec_clear(stairs->h, stairs->length);
}

/* For each p in num and for each shift s in group, computes an upper bound of
   F_{\tau(n)}(n) valid for n >= s, where

   F_T(n) = n \cdot \sum_{t=0}^{T-1} \left| [X^t] \frac{p(n+X)}{X^{-m(n)} q(n+X)} \right|.

   [M19], Algorithm 7.4, step 1 */

void
acb_ode_stairs_precompute(acb_ode_stairs_t stairs,
                          const acb_poly_struct * num, slong len,
                          const acb_poly_t ind, const acb_ode_group_t group,
                          const acb_ode_ind_lbound_t ind_lbound,
                          slong prec)
{
    /* flint_printf("_acb_ode_stairs_precompute %{acb}\n", group->leader); */

    if (stairs->h != NULL)
        acb_ode_stairs_clear(stairs);
    stairs->length = len /* +1 ? */ * group->nshifts;
    stairs->h = _mag_vec_init(stairs->length);

    mag_ptr bound = _mag_vec_init(2*len);
    mag_ptr bound1 = bound + len;

    /* This could in principle be refined when it is known that the degree in
     * log(x) of the solutions will not reach the bound. */
    slong ord = 0;
    for (slong s = 0; s < group->nshifts; s++)
        ord += group->shifts[s].mult;

    for (slong s = group->nshifts - 1; s >= 0; s--) {
        slong n = group->shifts[s].n;
        slong mult = acb_ode_group_multiplicity(group, n);
        /* flint_printf("s=%ld n=%ld mult=%ld ord=%ld\n", s, n, mult, ord); */
        acb_ode_bound_rat_ref_vec(bound, num, len, ind, n, mult, ord, prec);
        /* flint_printf("  n=%ld ref = %{mag*}\n", n, bound + i, len); */
        if (s == group->nshifts - 1 || group->shifts[s + 1].n != n + 1) {
            acb_ode_bound_rat_ordinary_vec(bound1, num, len, ind, ind_lbound,
                                           n + 1, ord, prec);
            /* flint_printf("  n=%ld bound1 = %{mag*}\n", n, bound1 + i, len); */
            for (slong i = 0; i < len; i++)
                mag_max(bound + i, bound + i, bound1 + i);
        }
        for (slong i = 0; i < len; i++) {
            slong j = i * group->nshifts + s;
            mag_swap(stairs->h + j, bound + i);
            if (s + 1 < group->nshifts)
                mag_max(stairs->h + j, stairs->h + j, stairs->h + j + 1);
        }
        ord -= group->shifts[s].mult;
    }

    _mag_vec_clear(bound, 2*len);
}
