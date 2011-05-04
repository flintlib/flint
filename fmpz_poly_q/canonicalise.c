#include "fmpz_poly_q.h"

void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop)
{
    fmpz_poly_t gcd;

    if (fmpz_poly_is_zero(rop->den))
    {
        printf("ERROR (fmpz_poly_q_canonicalise).  Denominator is zero.\n");
        abort();
    }

    if (fmpz_poly_is_one(rop->den))
        return;

    fmpz_poly_init(gcd);
    fmpz_poly_gcd(gcd, rop->num, rop->den);
    if (!fmpz_poly_is_unit(gcd))
    {
        fmpz_poly_div(rop->num, rop->num, gcd);
        fmpz_poly_div(rop->den, rop->den, gcd);
    }
    fmpz_poly_clear(gcd);

    if (fmpz_sgn(fmpz_poly_lead(rop->den)) < 0)
    {
        fmpz_poly_neg(rop->num, rop->num);
        fmpz_poly_neg(rop->den, rop->den);
    }
}

