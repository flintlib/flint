#include "fmpz_poly_q.h"

void fmpz_poly_q_derivative(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    fmpz_poly_t d, lhs, rhs;
    
    if (fmpz_poly_q_is_zero(op))
    {
        fmpz_poly_q_zero(rop);
        return;
    }
    
    if (fmpz_poly_length(op->den) == 1)
    {
        fmpz_poly_derivative(rop->num, op->num);
        fmpz_poly_set(rop->den, op->den);
        fmpz_poly_q_canonicalise(rop);
        return;
    }

    fmpz_poly_init(d);
    fmpz_poly_init(rhs);

    fmpz_poly_derivative(rhs, op->den);
    fmpz_poly_gcd(d, rhs, op->den);
    if (!fmpz_poly_is_one(d))
        fmpz_poly_div(rhs, rhs, d);
    fmpz_poly_mul(rhs, op->num, rhs);

    fmpz_poly_derivative(rop->num, op->num);
    if (fmpz_poly_is_one(d))
    {
        fmpz_poly_mul(rop->num, rop->num, op->den);
        fmpz_poly_pow(rop->den, op->den, 2);
    }
    else
    {
        fmpz_poly_init(lhs);
        fmpz_poly_div(lhs, op->den, d);
        fmpz_poly_mul(rop->num, rop->num, lhs);
        fmpz_poly_mul(rop->den, op->den, lhs);
        fmpz_poly_clear(lhs);
    }
    fmpz_poly_sub(rop->num, rop->num, rhs);

    /* Canonicalise:  there can be at most a constant factor */
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_poly_content(a, rop->num);
        fmpz_poly_content(b, rop->den);
        fmpz_gcd(c, a, b);

        if (!fmpz_is_one(c))
        {
            fmpz_poly_scalar_divexact_fmpz(rop->num, rop->num, c);
            fmpz_poly_scalar_divexact_fmpz(rop->den, rop->den, c);
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    
    fmpz_poly_clear(d);
    fmpz_poly_clear(rhs);
}
