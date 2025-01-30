#include "acb_types.h"
#include "acb.h"
#include "acb_poly.h"
#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"


/* Notation:
 *
 * We are considering an operator
 *
 *   \sum_{i,j} a[i,j] x^j (x d/dx)^i
 *
 * applied to a series
 *
 *   f = x^{expo + offset} \sum_{p,k} f[k,p] x^p log(x)^k/k!,
 *
 * resulting in
 *
 *   g = x^{expo + offset} \sum{m1,k1} g[k1,m1] x^m1 log(x)^k1/k1!.
 *
 * We want to compute the terms of index m1 = start + p1 of g for 0 <= p1 < len.
 */


/* TODO specialized version for dop with fmpz coefficients (and expo == 0) */


static void
acb_addmul_binom(acb_ptr c, acb_srcptr b, fmpz_mat_t binom, slong i, slong t,
        slong prec)
{
    if (t == 0)
        acb_add(c, c, b, prec);
    else if (t == 1)
        acb_addmul_si(c, b, i, prec);
    else
        acb_addmul_fmpz(c, b, fmpz_mat_entry(binom, i, t), prec);
}


void
_acb_holonomic_apply_diffop_basecase_weights(
        acb_ptr weights,  /* len * nlogs * flen */
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        slong flen, slong nlogs, slong start, slong len,
        slong prec)
{

    fmpz_mat_t binom;
    if (nlogs >= 2)  /* XXX: worth caching? */
    {
        gr_ctx_t fmpz_ctx;
        gr_ctx_init_fmpz(fmpz_ctx);
        fmpz_mat_init(binom, dop_len, dop_len);
        GR_MUST_SUCCEED(gr_mat_pascal((gr_mat_struct *) binom, -1, fmpz_ctx));
        gr_ctx_clear(fmpz_ctx);
    }

    acb_t expo1;
    acb_init(expo1);

    if (expo != NULL && acb_is_zero(expo))
        expo = NULL;

    for (slong p1 = 0; p1 < len; p1++)  /* n1 = offset + start + p1 */
    {
        for (slong t = 0; t < nlogs; t++)  /* t = k - k1 */
        {
            for (slong p = 0; p < flen; p++)
            {
                slong j = start + p1 - p;
                slong n = offset + p;

                if (j < 0)
                    break;

                if (expo != NULL)
                    acb_add_si(expo1, expo, n, prec);

                acb_ptr c = weights + p1 * nlogs * flen + t * flen + p;
                acb_zero(c);
                for (slong i = dop_len - 1; i >= t; i--)  /* Horner */
                {
                    if (expo == NULL)
                        acb_mul_si(c, c, n, prec);
                    else
                        acb_mul(c, c, expo1, prec);
                    if (j >= (dop + i)->length)
                        continue;
                    acb_ptr aij = (dop + i)->coeffs + j;
                    acb_addmul_binom(c, aij, binom, i, t, prec);
                }
            }
        }
    }

    if (nlogs >= 2)
        fmpz_mat_clear(binom);
}


void
_acb_holonomic_apply_diffop_basecase_precomp(
        acb_poly_struct * g, slong goff,
        acb_srcptr weights, slong weights_nlogs,
        const acb_poly_struct * f, slong foff, slong flen,
        slong nlogs,
        slong start, slong len,
        slong prec)
{
    for (slong p1 = 0; p1 < len; p1++)  /* n1 = offset + start + p1 */
    {
        for (slong k1 = 0; k1 < nlogs; k1++)
        {
            acb_ptr dest = (g + k1)->coeffs + goff + p1;

            /* XXX integrate the loop on k in acb_dot? */
            for (slong k = k1; k < nlogs; k++)
            {
                acb_srcptr src = (f + k)->coeffs + foff;
                acb_srcptr cofac = weights + p1 * weights_nlogs * flen + (k - k1) * flen;
                slong fklen = FLINT_MIN(flen, (f + k)->length - foff);
                /* loop on p */
                acb_dot(dest, dest, 0, cofac, 1, src, 1, fklen, prec);
            }
        }
    }
}


void
acb_holonomic_apply_diffop_basecase(
        acb_poly_struct * g,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f,
        slong nlogs,
        slong start, slong len,
        slong prec)
{
    slong flen = 0;
    for (int k = 0; k < nlogs; k++)
        flen = FLINT_MAX(flen, (f + k)->length);
    /* XXX could we take flen = FLINT_MIN(flen, len) or similar? */

    acb_ptr weights = _acb_vec_init(len * nlogs * flen);

    _acb_holonomic_apply_diffop_basecase_weights(
            weights, dop, dop_len, expo, offset,
            flen, nlogs, start, len, prec);

    _acb_poly_vec_fit_length(g, nlogs, len);
    _acb_holonomic_apply_diffop_basecase_precomp(
            g, 0,
            weights, nlogs, f, 0, flen,
            nlogs, start, len, prec);
    _acb_poly_vec_set_length(g, nlogs, len);
    _acb_poly_vec_normalise(g, nlogs);

    _acb_vec_clear(weights, len * nlogs * flen);
}
