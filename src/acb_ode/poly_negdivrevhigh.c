#include "acb.h"
#include "acb_ode.h"

typedef enum
{
    LC_UNKNOWN,
    LC_ACB,
    LC_FMPZ
}
lc_status_t;

// XXX maybe accept separate lengths
void
_acb_ode_poly_negdivrevhigh(acb_ptr res, acb_srcptr a, acb_srcptr cst, acb_srcptr b,
                        slong len, slong prec)
{
    fmpz_t lc_fmpz;
    acb_t invlc;
    lc_status_t lc_status = LC_UNKNOWN;

    fmpz_init(lc_fmpz);
    acb_init(invlc);

    for (slong k = len - 1; k >= 0; k--)
    {
        /* cst*b[k] + sum(a[1+j] * res[k+1+j], j≥0) */
        if (cst == NULL)
            acb_dot(res + k, b + k, 0, a + 1, 1, res + k + 1, 1, len - k - 1, prec);
        else
        {
            acb_mul(res + k, cst, b + k, prec);
            acb_dot(res + k, res + k, 0, a + 1, 1, res + k + 1, 1, len - k - 1, prec);
        }

        if (acb_is_zero(res + k))
            continue;

        /* only compute invlc if needed */
        // XXX share between calls?
        if (lc_status == LC_UNKNOWN)
        {
            if (acb_is_exact(a) && acb_get_unique_fmpz(lc_fmpz, a))
            {
                lc_status = LC_FMPZ;
                acb_indeterminate(invlc);
            }
            else
            {
                lc_status = LC_ACB;
                acb_inv(invlc, a, prec);
            }
        }

        if (lc_status == LC_FMPZ)
            acb_div_fmpz(res + k, res + k, lc_fmpz, prec);
        else
            acb_mul(res + k, res + k, invlc, prec);

        acb_neg(res + k, res + k);
    }

    acb_clear(invlc);
    fmpz_clear(lc_fmpz);
}
