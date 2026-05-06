#include "acb.h"
#include "acb_ode.h"
#include "arf.h"
#include "mag.h"

slong
acb_ode_sol_estimate_terms(mag_t est,
                           const acb_ode_sol_t sol, slong off, slong len,
                           const mag_t radpow)
{
    mag_t m;

    mag_init(m);

    slong accuracy = ARF_PREC_EXACT;
    mag_zero(est);

    for (slong k = 0; k < sol->nlogs; k++)
    {
        for (slong i = 0; i < len; i++)
        {
            acb_ptr c = (sol->series + k)->coeffs + off + i;

            acb_get_mag(m, c);
            mag_add(est, est, m);

            accuracy = FLINT_MIN(accuracy, acb_rel_accuracy_bits(c));
        }
    }

    mag_mul(est, est, radpow);

    mag_clear(m);
    return accuracy;
}

slong
acb_ode_sol_estimate_sums(mag_t mag, mag_t rad, const acb_ode_sol_t sol)
{
    slong len = sol->npts * sol->nlogs * sol->nder;

    _acb_vec_get_mag(mag, sol->sums, len);  /* XXX good enough? */

    mag_zero(rad);
    for (slong i = 0; i < len; i++)
    {
        mag_max(rad, rad, &sol->sums[i].real.rad);
        mag_max(rad, rad, &sol->sums[i].imag.rad);
    }

    slong accuracy = ARF_PREC_EXACT;
    for (slong i = 0; i < len; i++)
        accuracy = FLINT_MIN(accuracy, acb_rel_accuracy_bits(sol->sums + i));

    return accuracy;
}
