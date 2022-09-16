
#include "acb_theta.h"

void
acb_theta_agm_ctx_reset_steps(acb_theta_agm_ctx_t ctx, slong k, slong m)
{
  slong n = acb_theta_agm_ctx_nb(ctx);
  slong prev = acb_theta_agm_ctx_nb_bad_steps(ctx, k);
  slong nb_th;
  slong j;

  nb_th = n;
  if (acb_theta_agm_ctx_is_ext(ctx)) nb_th *= 2;

  acb_theta_agm_ctx_nb_bad_steps(ctx, k) = m;
  if (prev == 0 && m > 0)
    {
      acb_theta_agm_ctx_roots(ctx, k) = _acb_vec_init(m * nb_th);
      acb_theta_agm_ctx_mi(ctx, k) = flint_malloc(m * sizeof(arf_struct));
      for (j = 0; j < m; j++) arf_init(&acb_theta_agm_ctx_mi(ctx, k)[j]);
    }
  else if (prev > 0 && m == 0)
    {
      _acb_vec_clear(acb_theta_agm_ctx_roots(ctx, k), prev * nb_th);
      for (j = 0; j < prev; j++) arf_clear(&acb_theta_agm_ctx_mi(ctx, k)[j]);
      flint_free(acb_theta_agm_ctx_mi(ctx, k));      
    }
  else if (prev > 0 && m > 0)
    {
      acb_theta_agm_ctx_reset_steps(ctx, k, 0);
      acb_theta_agm_ctx_reset_steps(ctx, k, m);
    }
}
