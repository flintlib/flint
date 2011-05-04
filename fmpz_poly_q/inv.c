#include "fmpz_poly_q.h"

void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    if (fmpz_poly_is_zero(op->num))
    {
        printf("ERROR (fmpz_poly_q_inv).  Denominator is zero.\n");
        abort();
    }
    
    if (rop == op)
    {
        fmpz_poly_swap(rop->num, rop->den);
        if (fmpz_sgn(fmpz_poly_lead(rop->den)) < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_poly_neg(rop->den, rop->den);
        }
    }
    else
    {
        if (fmpz_sgn(fmpz_poly_lead(op->num)) > 0)
        {
            fmpz_poly_set(rop->num, op->den);
            fmpz_poly_set(rop->den, op->num);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->den);
            fmpz_poly_neg(rop->den, op->num);
        }
    }
}

