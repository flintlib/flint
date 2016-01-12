#include<gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_mul_linear(fmpz * res, const fmpz * poly1,
                      slong len1, const fmpz * poly2, slong len2)
{
    slong i, j, c, rlen = len1 + len2 - 1;
    fmpz_one(res);
    fmpz_one(res + rlen - 1);

    /* lowest term */
    fmpz_mul(res, poly1, poly2);

    /* highest term */
    fmpz_mul(res + rlen - 1, res + rlen - 1, poly1 + len1 - 1);

    for (i = 1; i < rlen - 1; i++)
    {
        res[i] = poly1[i - 1] + poly2[0] * poly1[i];
    }
}

void
fmpz_poly_mul_linear(fmpz_poly_t res,
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong rlen;
    
    /* if none of the polynomials are linear*/
    if (len1 != 1 && len2 != 1)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    rlen = len1 + len2 - 1;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t = t;
        fmpz_poly_init2(t, rlen);
        if (len1 >= len2)
            _fmpz_poly_mul_linear(t->coeffs, poly1->coeffs, len1,                                      poly2->coeffs, len2);
        else
            _fmpz_poly_mul_linear(t->coeffs, poly2->coeffs, len2,                                      poly1->coeffs, len1);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, rlen);
        if(len1 >= len2)
            _fmpz_poly_mul_linear(res->coeffs, poly1->coeffs, len1                                   poly2->coeffs, len2);
        else
            _fmpz_poly_mul_linear(res->coeffs, poly2->coeffs, len2                                   poly1->coeffs, len1);
    }

    _fmpz_poly_set_length(res, rlen);
}


