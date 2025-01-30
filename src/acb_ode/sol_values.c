#include "acb_types.h"
#include "acb_poly.h"
#include "acb_ode.h"


void
_acb_ode_sol_value(acb_poly_struct * val, acb_srcptr expo,
                   const acb_poly_struct * sums, slong nlogs,
                   acb_srcptr pt, slong nder,
                   slong nfrobshifts, slong prec)
{
    acb_poly_t expt, exptpow, log, logpow, term;

    acb_poly_init(expt);
    acb_poly_init(exptpow);
    acb_poly_init(log);
    acb_poly_init(logpow);
    acb_poly_init(term);

    acb_poly_set_coeff_si(expt, 1, 1);
    acb_poly_set_coeff_acb(expt, 0, pt);

    /* XXX does this work when pt contains 0 but expo is a nonneg int? */
    acb_poly_pow_acb_series(exptpow, expt, expo, nder, prec);
    acb_poly_log_series(log, expt, nder, prec);

    acb_poly_one(logpow);

    if (nlogs >= 1)
    {
        slong len = FLINT_MIN(nlogs, nfrobshifts);
        _acb_poly_vec_set(val, sums, len);
        _acb_poly_vec_zero(val + len, nlogs - len);
    }

    for (slong p = 1; p < nlogs; p++)
    {
        acb_poly_mullow(logpow, logpow, log, nder, prec);
        _acb_vec_scalar_div_ui(logpow->coeffs, logpow->coeffs, logpow->length,
                               p, prec);

        for (slong d = 0; d < FLINT_MIN(nlogs - p, nfrobshifts); d++)
        {
            acb_poly_mullow(term, sums + d + p, logpow, nder, prec);
            acb_poly_add(val + d, val + d, term, prec);
        }
    }

    for (slong d = 0; d < nfrobshifts; d++)
        acb_poly_mullow(val + d, val + d, exptpow, nder, prec);

    acb_poly_clear(expt);
    acb_poly_clear(exptpow);
    acb_poly_clear(log);
    acb_poly_clear(logpow);
    acb_poly_clear(term);
}


/* XXX rename? move? */

void
acb_ode_sum_value(acb_poly_struct * val, slong nfrobshifts,
                  const acb_ode_sum_context_struct * ctx,
                  slong i, slong j, slong prec)  /* sol #i, pt #j */
{
    /* flint_printf("i=%wd j=%wd nlogs=%wd\n", i, j, ctx->sol[i].nlogs); */
    /* flint_printf("%{acb_poly}\n", acb_ode_sol_sum_ptr(ctx->sol + i, j, 0)); */
    _acb_ode_sol_value(val, ctx->expo,
                       acb_ode_sol_sum_ptr(ctx->sol + i, j, 0),
                       ctx->sol[i].nlogs, ctx->pts + j, ctx->nder,
                       nfrobshifts, prec);
}
