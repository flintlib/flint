
#include "acb_theta.h"

void acb_theta_newton_eval(acb_ptr r, acb_scrptr th, const acb_theta_newton_t ctx, slong prec)
{
  slong g = acb_theta_newton_g(ctx);
  slong n = acb_theta_newton_nb(ctx);
  acb_ptr dupl;
  acb_ptr transf;
  acb_ptr agm;
  slong k, j;

  dupl = _acb_vec_init(1<<(2*g));
  transf = _acb_vec_init(1<<g);
  agm = _acb_vec_init(n);
  
  acb_theta_duplication_all(dupl, th, prec);
  for (k = 0; k < n; k++)
    {      
      acb_theta_transform_sqr_proj(transf, dupl, acb_theta_newton_matrix(ctx, k), prec);
      for (j = 1; j < n; j++) acb_div(&transf[j], &transf[j], &transf[0], prec);
      acb_one(&transf[0]);
      acb_theta_agm(&agm[k], transf, acb_theta_newton_roots(ctx, k),
		    acb_theta_newton_nb_bad_steps(ctx, k),
		    acb_theta_newton_nb_bad_steps(ctx, k) + acb_theta_agm_nb_good_steps(g, k),
		    g, prec);
    }
  for (k = 0; k < n-1; k++)
    {
      acb_div(&r[k], &agm[0], &agm[k+1], prec);
    }

  _acb_vec_clear(dupl, 1<<(2*g));
  _acb_vec_clear(transf, 1<<g);
  _acb_vec_clear(agm, n);
}
