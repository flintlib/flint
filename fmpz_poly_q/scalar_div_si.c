#include "fmpz_poly_q.h"

void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, long x)
{
    fmpz_t cont, fx, gcd;
    
    if (x == 0)
    {
        printf("ERROR (fmpz_poly_q_scalar_div_si).  Division by zero.\n");
        abort();
    }
    if (x == 1)
    {
        fmpz_poly_q_set(rop, op);
        return;
    }
    if (fmpz_poly_q_is_zero(op))
    {
        fmpz_poly_q_zero(rop);
        return;
    }

    fmpz_init(cont);
    fmpz_poly_content(cont, op->num);

    if (fmpz_is_one(cont))
    {
        if (x > 0)
        {
            fmpz_poly_set(rop->num, op->num);
            fmpz_poly_scalar_mul_si(rop->den, op->den, x);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->num);
            fmpz_poly_scalar_mul_si(rop->den, op->den, -x);
        }
        fmpz_clear(cont);
        return;
    }

    fmpz_init(fx);
    fmpz_init(gcd);

    fmpz_set_si(fx, x);
    fmpz_gcd(gcd, cont, fx);

    if (fmpz_is_one(gcd))
    {
        if (x > 0)
        {
            fmpz_poly_set(rop->num, op->num);
            fmpz_poly_scalar_mul_si(rop->den, op->den, x);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->num);
            fmpz_poly_scalar_mul_si(rop->den, op->den, -x);
        }
    }
    else
    {
        fmpz_poly_scalar_divexact_fmpz(rop->num, op->num, gcd);
        fmpz_divexact(fx, fx, gcd);
        fmpz_poly_scalar_mul_fmpz(rop->den, op->den, fx);
        if (x < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_poly_neg(rop->den, rop->den);
        }
    }

    fmpz_clear(cont);
    fmpz_clear(fx);
    fmpz_clear(gcd);
}
