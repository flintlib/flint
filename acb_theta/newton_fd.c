
#include "acb_theta.h"

void
acb_theta_newton_fd(acb_ptr v, acb_mat_t fd, acb_srcptr th, const arb_t eta,
                    const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong dim = acb_theta_agm_ctx_dim(ctx);
    slong nb_th = 1 << g;
    acb_ptr thmod;
    acb_ptr r0, r;
    slong k, j, c;

    if (acb_theta_agm_ctx_is_ext(ctx))
        nb_th *= 2;
    thmod = _acb_vec_init(nb_th);
    r0 = _acb_vec_init(dim);
    r = _acb_vec_init(dim);

    acb_theta_newton_eval(r0, th, ctx, prec);

    c = 0;
    for (k = 1; k < nb_th; k++)
    {
        if (k == (1 << g))
            continue;

        _acb_vec_set(thmod, th, nb_th);
        acb_add_arb(&thmod[k], &thmod[k], eta, prec);
        acb_theta_newton_eval(r, thmod, ctx, prec);
        _acb_vec_sub(r, r, r0, dim, prec);
        _acb_vec_scalar_div_arb(r, r, dim, eta, prec);
        for (j = 0; j < dim; j++)
        {
            acb_set(acb_mat_entry(fd, j, c), &r[j]);
        }
        c++;
    }
    _acb_vec_set(v, r0, dim);

    _acb_vec_clear(thmod, nb_th);
    _acb_vec_clear(r0, dim);
    _acb_vec_clear(r, dim);
}
