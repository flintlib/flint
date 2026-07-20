#include "acb_types.h"
#include "acb_poly.h"
#include "acb_ode.h"

/* XXX should val be a plain acb array? acb functions that compute jets usually
 * return them as acb_ptrs, not acb_polys, but here we are computing several
 * jets at once, so the user cannot have these functions write into an acb_poly
 * for further processing */

/* Here nlogs is (a bound on) the number of logs in the **mathematical**
   solution, not in the current partial sum, because we need the tail bound to
   be added to the series in front of each log(x)^k/k!, including those not
   visible in the current partial sum. */
void
_acb_ode_sol_jet(acb_poly_struct * val, const acb_t expo,
                 acb_srcptr sums, slong stride, mag_srcptr tb, const acb_t pt,
                 slong nlogs, slong ord, slong nfrobshifts, slong prec)
{
    acb_poly_t expt, exptpow, log, logpow, term;

    FLINT_ASSERT(ord > 0);

    if (nlogs == 0)
    {
        _acb_poly_vec_zero(val, nfrobshifts);
        return;
    }

    acb_poly_init(expt);
    acb_poly_init(exptpow);
    acb_poly_init(log);
    acb_poly_init(logpow);
    acb_poly_init(term);

    acb_poly_struct * sum_with_err = _acb_poly_vec_init(nlogs);

    acb_poly_fit_length(term, ord);

    acb_poly_set_coeff_si(expt, 1, 1);
    acb_poly_set_coeff_acb(expt, 0, pt);

    /* XXX does this work when pt contains 0 but expo is a nonneg int? */
    acb_poly_pow_acb_series(exptpow, expt, expo, ord, prec);
    acb_poly_log_series(log, expt, ord, prec);

    acb_poly_one(logpow);

    for (slong k = 0; k < nlogs; k++)
    {
        acb_poly_fit_length(sum_with_err + k, ord);
        _acb_vec_set((sum_with_err + k)->coeffs, sums + k * stride, ord);
        for (slong i = 0; i < ord; i++)
            acb_add_error_mag((sum_with_err + k)->coeffs + i, tb + i);
        _acb_poly_set_length(sum_with_err + k, ord);
        _acb_poly_normalise(sum_with_err + k);
    }

    slong len = FLINT_MIN(nlogs, nfrobshifts);
    for (slong fs = 0; fs < len; fs++)
        acb_poly_set(val + fs, sum_with_err + fs);

    for (slong dk = 1; dk < nlogs; dk++)
    {
        acb_poly_mullow(logpow, logpow, log, ord, prec);
        _acb_vec_scalar_div_ui(logpow->coeffs, logpow->coeffs, logpow->length,
                               dk, prec);

        for (slong fs = 0; fs < FLINT_MIN(nlogs - dk, nfrobshifts); fs++)
        {
            acb_poly_mullow(term, sum_with_err + fs + dk, logpow, ord, prec);
            acb_poly_add(val + fs, val + fs, term, prec);
        }
    }

    for (slong fs = 0; fs < len; fs++)
        acb_poly_mullow(val + fs, val + fs, exptpow, ord, prec);

    _acb_poly_vec_zero(val + len, nfrobshifts - len);

    _acb_poly_vec_clear(sum_with_err, nlogs);
    acb_poly_clear(expt);
    acb_poly_clear(exptpow);
    acb_poly_clear(log);
    acb_poly_clear(logpow);
    acb_poly_clear(term);
}


void
acb_ode_sol_jet(acb_poly_struct * val,
                const acb_t expo, const acb_ode_sol_t sol, slong p,
                const acb_t pt, slong ord, slong nfrobshifts, slong prec)
{
    FLINT_ASSERT(ord <= sol->nder);

    if (ord == 0)
        return;

    _acb_ode_sol_jet(val, expo,
                     acb_ode_sol_sum_ptr(sol, p, 0, 0), sol->nder, sol->tb, pt,
                     sol->nlogs + sol->future_logs, ord, nfrobshifts,
                     prec);
}
