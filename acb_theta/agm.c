
#include "acb_theta.h"

void
acb_theta_agm(acb_t r, acb_srcptr a, acb_srcptr roots, slong nb_bad,
        slong nb_good, slong g, slong prec)
{
    acb_ptr v;
    acb_t scal;
    arf_t err;
    slong n = 1<<g;
    slong k;

    v = _acb_vec_init(n);
    acb_init(scal);
    arf_init(err);

    arf_one(err);
    arf_mul_2exp_si(err, err, -prec);
    _acb_vec_set(v, a, n);
  
    for (k = 0; k < nb_bad; k++)
    {
        acb_theta_agm_step_bad(v, v, roots + k*n, g, prec);
    }
    
    acb_set(scal, &v[0]);
    _acb_vec_scalar_div(v, v, n, scal, prec);

    flint_printf("(agm_with_err) Starting %wd good steps\n", nb_good);
    for (k = 0; k < n; k++)
    {
        acb_printd(&v[k], 10); flint_printf("\n");
    }
        
    for (k = 0; k < nb_good-1; k++)
    {
        acb_theta_agm_step_good(v, v, g, prec);
    }
    if (nb_good > 0) acb_theta_agm_step_last(r, v, g, prec);
    else acb_set(r, &v[0]);

    acb_add_error_arf(r, err);
    acb_mul(r, r, scal, prec);
    
    flint_printf("(agm_with_err) Reached agm\n");
    acb_printd(r, 10); flint_printf("\n");

    _acb_vec_clear(v, n);
    acb_clear(scal);
    arf_clear(err);
}
