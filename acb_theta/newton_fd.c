
#include "acb_theta.h"

void acb_theta_newton_fd(acb_ptr v, acb_mat_t fd, acb_srcptr th, const arb_t eta,
			 const acb_theta_agm_ctx_t ctx, slong prec)
{
  slong g = acb_theta_agm_ctx_g(ctx);
  slong n = acb_theta_agm_ctx_nb(ctx);
  acb_ptr thproj, thmod;
  acb_ptr r0, r;
  slong k, j;

  thproj = _acb_vec_init(1<<g);
  thmod = _acb_vec_init(1<<g);
  r0 = _acb_vec_init(n-1);
  r = _acb_vec_init(n-1);
  
  if (!acb_is_one(&th[0]))
    {
      for (k = 1; k < (1<<g); k++) acb_div(&thproj[k], &th[k], &th[0], prec);
      acb_one(&thproj[0]);
    }
  else _acb_vec_set(thproj, th, 1<<g);    
  
  acb_theta_newton_eval(r0, thproj, ctx, prec);
  
  for (k = 1; k < (1<<g); k++)
    {
      _acb_vec_set(thmod, thproj, 1<<g);
      acb_add_arb(&thmod[k], &thmod[k], eta, prec);
      acb_theta_newton_eval(r, thmod, ctx, prec);
      _acb_vec_sub(r, r, r0, n-1, prec);
      _acb_vec_scalar_div_arb(r, r, n-1, eta, prec);
      for (j = 0; j < n-1; j++)
	{
	  acb_set(acb_mat_entry(fd, j, k-1), &r[j]);
	}
    }
  _acb_vec_set(v, r0, n-1);
  
  _acb_vec_clear(thproj, 1<<g);
  _acb_vec_clear(thmod, 1<<g);
  _acb_vec_clear(r0, n-1);
  _acb_vec_clear(r, n-1);
}
