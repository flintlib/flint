
#include "acb_theta.h"

void
acb_theta_agm(acb_t r, acb_srcptr a, acb_srcptr all_roots, const arf_t rel_err,
        slong nb_bad, slong nb_good, slong g, slong prec)
{
    acb_ptr v;
    acb_t scal;
    arb_t abs;
    arf_t err;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong n = 1<<g;
    slong k;

    v = _acb_vec_init(n);
    acb_init(scal);
    arb_init(abs);
    arf_init(err);
  
    _acb_vec_set(v, a, n);
  
    for (k = 0; k < nb_bad; k++)
    {
        acb_theta_agm_step_bad(v, v, all_roots + k*n, g, prec);
    }
    
    acb_set(scal, &v[0]);
    _acb_vec_scalar_div(v, v, n, scal, prec);
        
    for (k = 0; k < nb_good; k++)
    {
        acb_theta_agm_step_good(v, v, g, prec);
    }
    
    acb_mul(r, &v[0], scal, prec);
    
    acb_abs(abs, r, lowprec);
    arb_get_ubound_arf(err, abs, lowprec);
    arf_mul(err, err, rel_err, lowprec, ARF_RND_CEIL);  
    acb_add_error_arf(r, err);
    
    _acb_vec_clear(v, n);
    acb_clear(scal);
    arb_clear(abs);
    arf_clear(err);
}
