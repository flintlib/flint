
#include "acb_theta.h"

void acb_theta_newton_logs(slong* log_max, slong* log_rho, slong* log_B1, slong* log_B2,
			   slong* log_B3, const acb_theta_agm_ctx_t ctx)
{
  arf_t c;
  fmpz_t e;
  slong lowprec = ACB_THETA_AGM_LOWPREC;

  arf_init(c);
  fmpz_init(e);
  
  arf_frexp(c, e, acb_theta_agm_ctx_max(ctx));
  *log_max = fmpz_get_si(e);
  arf_frexp(c, e, acb_theta_agm_ctx_rho(ctx));
  *log_rho = fmpz_get_si(e) - 1;  
  arf_mul_si(c, acb_theta_agm_ctx_max(ctx), 2*n, lowprec, ARF_RND_CEIL);
  arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
  arf_frexp(c, e, c);
  *log_B1 = fmpz_get_si(e);
  arf_mul_si(c, acb_theta_agm_ctx_max(ctx), 2*n*(n+1), lowprec, ARF_RND_CEIL);
  arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
  arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
  arf_frexp(c, e, c);
  *log_B2 = fmpz_get_si(e);  
  arf_frexp(c, e, acb_theta_agm_ctx_inv_der(ctx));
  *log_B3 = fmpz_get_si(e);

  arf_clear(c);
  arf_clear(e);
}
