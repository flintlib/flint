
#include "acb_theta.h"

int
acb_theta_agm_ctx_set_inv_der(acb_theta_agm_ctx_t ctx, acb_srcptr th,
        slong prec)
{
    slong dim = acb_theta_agm_ctx_nb(ctx) - 1;
    arf_t B2;
    fmpz_t e;
    slong exp;
    arb_t eta;
    acb_ptr r;
    acb_mat_t fd, fdinv;
    arb_t norm, bound, test;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    int res;
    
    if (acb_theta_agm_ctx_is_ext(ctx)) dim *= 2;
    
    arf_init(B2);
    fmpz_init(e);
    arb_init(eta);
    r = _acb_vec_init(dim);
    acb_mat_init(fd, dim, dim);
    acb_mat_init(fdinv, dim, dim);
    arb_init(norm);
    arb_init(bound);
    arb_init(test);
        
    /* Evaluate finite difference */
    acb_theta_cauchy(B2, acb_theta_agm_ctx_rho(ctx),
            acb_theta_agm_ctx_max(ctx), 2, dim, lowprec);
    arf_frexp(B2, e, B2);
    exp = fmpz_get_si(e);
    arb_one(eta);
    arb_mul_2exp_si(eta, eta, FLINT_MIN(- exp - n_clog(dim, 2), - prec/2));
    acb_theta_newton_fd(r, fd, th, eta, ctx, prec);
    res = acb_mat_inv(fdinv, fd, prec);

    if (!res) arb_pos_inf(norm);
    else acb_mat_ninf(norm, fdinv, lowprec);
      
    /* Is ||FD^-1||*n*B2*eta less than 1? */        
    arb_mul_arf(bound, norm, B2, lowprec);
    arb_mul_si(bound, bound, dim, lowprec);
    arb_mul(bound, bound, eta, lowprec);
    arb_sub_si(test, bound, 1, lowprec);
    if (!arb_is_negative(test)) res = 0;

    arb_mul(bound, bound, norm, lowprec);
    arb_add(bound, bound, norm, lowprec);
    arb_get_ubound_arf(acb_theta_agm_ctx_inv_der(ctx), bound, lowprec);

    arf_clear(B2);
    fmpz_clear(e);
    arb_clear(eta);
    _acb_vec_clear(r, dim);
    acb_mat_clear(fd);
    acb_mat_clear(fdinv);
    arb_clear(norm);
    arb_clear(bound);
    arb_clear(test);
    return res;
}
