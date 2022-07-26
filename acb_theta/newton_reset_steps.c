
#include "acb_theta.h"

void acb_theta_newton_reset_steps(acb_theta_newton_t ctx, slong k, slong m)
{
  slong n = acb_theta_newton_nb(ctx);
  slong prev = acb_theta_newton_nb_bad_steps(ctx, k);
  slong j;

  acb_theta_newton_nb_bad_steps(ctx, k) = m;
  if (prev == 0 && m > 0)
    {
      acb_theta_newton_roots(ctx, k) = _acb_vec_init(m * (n+1));
      acb_theta_newton_mi(ctx, k) = flint_malloc(m * sizeof(arf_struct));
      for (j = 0; j < m; j++) arf_init(&acb_theta_newton_mi(ctx, k)[j]);
    }
  else if (prev > 0 && m == 0)
    {
      _acb_vec_clear(acb_theta_newton_roots(ctx, k), prev * (n+1));
      for (j = 0; j < prev; j++) arf_clear(&acb_theta_newton_mi(ctx, k)[j]);
      flint_free(acb_theta_newton_mi(ctx, k));      
    }
  else if (prev > 0 && m > 0)
    {
      acb_theta_newton_reset_steps(ctx, k, 0);
      acb_theta_newton_reset_steps(ctx, k, m);
    }
}
