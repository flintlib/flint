#include "acb_ode.h"
#include "acb_poly.h"
#include "arb_poly.h"


void
acb_ode_bound_init(acb_ode_bound_t bound)
{
    bound->prec = -1;
    arb_init(bound->cst);
    arb_poly_init(bound->den_lbound);
    bound->den_rt = NULL;
    bound->den_rt_len = 0;
    bound->all_nums = NULL;
    bound->all_nums_len = 0;
    bound->pol_part_len = 0;
    acb_poly_init(bound->ind);
    acb_ode_ind_lbound_init(bound->ind_lbound);
    acb_ode_stairs_init(bound->stairs);
}

void
acb_ode_bound_clear(acb_ode_bound_t bound)
{
    arb_clear(bound->cst);
    arb_poly_clear(bound->den_lbound);
    _arb_vec_clear(bound->den_rt, bound->den_rt_len);
    _acb_poly_vec_clear(bound->all_nums, bound->all_nums_len);
    acb_poly_clear(bound->ind);
    acb_ode_ind_lbound_clear(bound->ind_lbound);
    acb_ode_stairs_clear(bound->stairs);
}
