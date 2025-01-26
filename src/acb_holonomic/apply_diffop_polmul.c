#include "acb_types.h"
#include "acb_poly.h"


static void
xD_logcoeff_inplace(acb_poly_struct * f, acb_srcptr expo, slong offset, slong k,
                    slong nlogs, slong prec)
{
    acb_poly_t tmp;
    acb_poly_init(tmp);

    /* XXX handle rational expo more efficiently (-> gr?) */
    if (expo != NULL)
        acb_poly_scalar_mul(tmp, f + k, expo, prec);
    for (slong j = 0; j < (f + k)->length; j++)
        acb_mul_ui((f + k)->coeffs + j, (f + k)->coeffs + j, offset + j, prec);
    acb_poly_add(f + k, f + k, tmp, prec);
    if (k + 1 < nlogs)
        acb_poly_add(f + k, f + k, f + k + 1, prec);

    acb_poly_clear(tmp);
}


/* Computes the coefficients of x^start to x^{start+len-1} (inclusive) of the
 * log-polynomial g such that dop(x^{expo+offset}·f) = x^{expo+offset}·g.
 *
 * expo may be NULL (treated as zero)
 *
 * The underscore version operates on the chunk of coefficients of length
 * flen starting at offset foff in each of the components of the input vector,
 * and _adds_ each component of g at offset goff in the corresponding component
 * of the output vector.
 */

/* XXX Take acb_struct**s instead of acb_poly_struct*s? More generally, how can
 * we make the interface more natural / consistent with FLINT? */

void
_acb_holonomic_apply_diffop_polmul(
        acb_poly_struct * g, slong goff,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f, slong foff, slong flen,
        slong nlogs,
        slong start, slong len,
        slong prec)
{
    acb_poly_t tmp;
    acb_poly_init(tmp);
    acb_poly_fit_length(tmp, start + len);

    /* XXX leave it to the caller to do this, and destroy f? */
    acb_poly_struct * curder = _acb_poly_vec_init(nlogs);
    _acb_poly_vec_set_block(curder, f, nlogs, foff, flen);

    for (slong i = 0; i < dop_len; i++)
    {
        for (slong k = 0; k < nlogs; k++)
        {
            /* should be a mulmid in the typical case (and should accept a
             * pretransformed first operand) */
            acb_poly_mullow(tmp, dop + i, curder + k, start + len, prec);
            acb_poly_shift_right(tmp, tmp, start);

            FLINT_ASSERT (tmp->length <= len);
            _acb_poly_add((g + k)->coeffs + goff,
                          (g + k)->coeffs + goff, len,
                          tmp->coeffs, tmp->length,
                          prec);

            /* curder[k] ← (x·d/dx(prev curder))[k] */
            xD_logcoeff_inplace(curder, expo, offset, k, nlogs, prec);
        }
    }

    _acb_poly_vec_clear(curder, nlogs);
    acb_poly_clear(tmp);
}


/* XXX support aliasing */
void
acb_holonomic_apply_diffop_polmul(
        acb_poly_struct * g,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f,
        slong nlogs,
        slong start, slong len,
        slong prec)
{
    _acb_poly_vec_fit_length(g, nlogs, len);
    _acb_holonomic_apply_diffop_polmul(
            g, 0,
            dop, dop_len,
            expo, offset,
            f, 0, start + len,
            nlogs,
            start, len,
            prec);
    _acb_poly_vec_set_length(g, nlogs, len);
    _acb_poly_vec_normalise(g, nlogs);
}
