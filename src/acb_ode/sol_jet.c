#include "acb_types.h"
#include "acb_poly.h"
#include "acb_ode.h"

/* XXX should val be a plain acb array? acb functions that compute jets usually
 * return them as acb_ptrs, not acb_polys, but here we are computing several
 * jets at once, so the user cannot have these functions write into an acb_poly
 * for further processing */

void
_acb_ode_sol_jet(acb_poly_struct * val, const acb_t expo,
                 acb_srcptr sums, slong stride, mag_srcptr tb, const acb_t pt,
                 slong nlogs, slong nder, slong nfrobshifts,
                 slong prec)
{
    acb_poly_t expt, exptpow, log, logpow, term;

    if (nlogs <= 0)
    {
        _acb_poly_vec_zero(val, nfrobshifts);
        return;
    }

    acb_poly_init(expt);
    acb_poly_init(exptpow);
    acb_poly_init(log);
    acb_poly_init(logpow);
    acb_poly_init(term);

    acb_poly_fit_length(term, nder);

    acb_poly_set_coeff_si(expt, 1, 1);
    acb_poly_set_coeff_acb(expt, 0, pt);

    /* XXX does this work when pt contains 0 but expo is a nonneg int? */
    acb_poly_pow_acb_series(exptpow, expt, expo, nder, prec);
    acb_poly_log_series(log, expt, nder, prec);

    acb_poly_one(logpow);

    slong len = FLINT_MIN(nlogs, nfrobshifts);
    for (slong d = 0; d < len; d++)
    {
        acb_poly_fit_length(val + d, nder);
        _acb_vec_set((val + d)->coeffs, sums + d * stride, nder);
        for (slong i = 0; i < nder; i++)
            acb_add_error_mag((val + d)->coeffs + i, tb + i);
        _acb_poly_set_length(val + d, nder);
        _acb_poly_normalise(val + d);
    }
    _acb_poly_vec_zero(val + len, nfrobshifts - len);

    for (slong p = 1; p < nlogs; p++)
    {
        acb_poly_mullow(logpow, logpow, log, nder, prec);
        _acb_vec_scalar_div_ui(logpow->coeffs, logpow->coeffs, logpow->length,
                               p, prec);

        for (slong d = 0; d < FLINT_MIN(nlogs - p, nfrobshifts); d++)
        {
            _acb_poly_mullow(term->coeffs,
                             sums + (d + p) * stride, nder,
                             logpow->coeffs, nder,
                             nder, prec);
            _acb_poly_set_length(term, nder);
            _acb_poly_normalise(term);
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


void
acb_ode_sol_jet(acb_poly_struct * val,
                const acb_t expo, const acb_ode_sol_t sol, slong j,
                const acb_t pt, slong ord, slong nfrobshifts, slong prec)
{
    FLINT_ASSERT(ord <= sol->nder);

    if (ord == 0)
        return;

    _acb_ode_sol_jet(val, expo,
                     acb_ode_sol_sum_ptr(sol, j, 0, 0), sol->nder, sol->tb, pt,
                     sol->nlogs, ord, nfrobshifts,
                     prec);
}
