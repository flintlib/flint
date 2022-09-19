
#include "acb_theta.h"

void
acb_theta_agm_ext(acb_t r, acb_t s, acb_srcptr a, acb_srcptr all_roots,
	const arf_t rel_err, slong nb_bad, slong nb_good, slong g, slong prec)
{  
    acb_ptr v;
    acb_t scal;
    fmpz_t exp;
    arb_t abs;
    arf_t err;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong n = 1<<g;
    slong k;

    v = _acb_vec_init(2*n);
    acb_init(scal);
    fmpz_init(exp);
    arb_init(abs);
    arf_init(err);
  
    _acb_vec_set(v, a, 2*n);
  
    for (k = 0; k < nb_bad; k++)
    {
	acb_theta_agm_ext_step_bad(v, v, all_roots + k*2*n, g, prec);
    }

    acb_set(scal, &v[0]);
    _acb_vec_scalar_div(v, v, 2*n, scal, prec);
    
    for (k = 0; k < nb_good; k++)
    {
	acb_theta_agm_ext_step_good(v, v, g, prec);
    }

    acb_div(r, &v[0], &v[n], prec);
    fmpz_one(exp);
    fmpz_mul_2exp(exp, exp, nb_good);
    acb_pow_fmpz(r, r, exp, prec);

    acb_abs(abs, r, lowprec);
    arb_get_ubound_arf(err, abs, lowprec);
    arf_mul(err, err, rel_err, lowprec, ARF_RND_CEIL);
    acb_add_error_arf(r, err);

    fmpz_one(exp);
    fmpz_mul_2exp(exp, exp, nb_bad);
    acb_pow_fmpz(r, r, exp, prec);

    _acb_vec_clear(v, 2*n);
    acb_clear(scal);
    fmpz_clear(exp);
    arb_clear(abs);
    arf_clear(err);
}
