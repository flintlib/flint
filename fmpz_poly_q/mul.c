#include "fmpz_poly_q.h"

void fmpz_poly_q_mul_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    fmpz_poly_t d, e;

    if (fmpz_poly_q_is_zero(rop) || fmpz_poly_q_is_zero(op))
    {
        fmpz_poly_q_zero(rop);
        return;
    }
    
    if (rop == op)
    {
        fmpz_poly_pow(rop->num, op->num, 2);
        fmpz_poly_pow(rop->den, op->den, 2);
        return;
    }

    /*
        From here on, rop and op point to two different object in memory, 
        and these are both non-zero rational functions
     */

    fmpz_poly_init(d);
    fmpz_poly_gcd(d, rop->num, op->den);

    if (fmpz_poly_is_one(d))
    {
        fmpz_poly_gcd(d, rop->den, op->num);
        if (fmpz_sgn(fmpz_poly_lead(d)) < 0)
            fmpz_poly_neg(d, d);
        
        if (fmpz_poly_is_one(d))
        {
            fmpz_poly_mul(rop->num, rop->num, op->num);
            fmpz_poly_mul(rop->den, rop->den, op->den);
        }
        else
        {
            fmpz_poly_div(rop->den, rop->den, d);
            fmpz_poly_mul(rop->den, rop->den, op->den);
            fmpz_poly_div(d, op->num, d);
            fmpz_poly_mul(rop->num, rop->num, d);
        }
    }
    else
    {
        fmpz_poly_init(e);
        fmpz_poly_gcd(e, rop->den, op->num);
        if (fmpz_sgn(fmpz_poly_lead(e)) < 0)
            fmpz_poly_neg(e, e);
        
        if (fmpz_poly_is_one(e))
        {
            fmpz_poly_div(rop->num, rop->num, d);
            fmpz_poly_mul(rop->num, rop->num, op->num);
            fmpz_poly_div(d, op->den, d);
            fmpz_poly_mul(rop->den, rop->den, d);
        }
        else
        {
            fmpz_poly_div(rop->num, rop->num, d);
            fmpz_poly_div(rop->den, rop->den, e);
            fmpz_poly_div(e, op->num, e);
            fmpz_poly_div(d, op->den, d);
            fmpz_poly_mul(rop->num, rop->num, e);
            fmpz_poly_mul(rop->den, rop->den, d);
        }
        fmpz_poly_clear(e);
    }
    fmpz_poly_clear(d);
}

void fmpz_poly_q_mul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    fmpz_poly_t t, u;
    
    if (op1 == op2)
    {
        fmpz_poly_pow(rop->num, op1->num, 2);
        fmpz_poly_pow(rop->den, op1->den, 2);
        return;
    }
    
    if (rop == op1)
    {
        fmpz_poly_q_mul_in_place(rop, op2);
        return;
    }
    if (rop == op2)
    {
        fmpz_poly_q_mul_in_place(rop, op1);
        return;
    }

    if (fmpz_poly_q_is_zero(op1) || fmpz_poly_q_is_zero(op2))
    {
        fmpz_poly_q_zero(rop);
        return;
    }

    /*
        From here on, we may assume that rop, op1 and op2 refer to distinct 
        objects in memory, and that op1 and op2 are non-zero
     */
    
    /* Are both denominators equal to one? */
    if (fmpz_poly_is_one(op1->den) && fmpz_poly_is_one(op2->den))
    {
        fmpz_poly_mul(rop->num, op1->num, op2->num);
        fmpz_poly_zero(rop->den);
        fmpz_poly_set_coeff_si(rop->den, 0, 1);
        return;
    } 
    
    fmpz_poly_gcd(rop->num, op1->num, op2->den);
    if (fmpz_sgn(fmpz_poly_lead(rop->num)) < 0)
        fmpz_poly_neg(rop->num, rop->num);
    
    if (fmpz_poly_is_one(rop->num))
    {
        fmpz_poly_gcd(rop->den, op2->num, op1->den);
        if (fmpz_sgn(fmpz_poly_lead(rop->den)) < 0)
            fmpz_poly_neg(rop->den, rop->den);
        
        if (fmpz_poly_is_one(rop->den))
        {
            fmpz_poly_mul(rop->num, op1->num, op2->num);
            fmpz_poly_mul(rop->den, op1->den, op2->den);
        }
        else
        {
            fmpz_poly_div(rop->num, op2->num, rop->den);
            fmpz_poly_mul(rop->num, op1->num, rop->num);
            fmpz_poly_div(rop->den, op1->den, rop->den);
            fmpz_poly_mul(rop->den, rop->den, op2->den);
        }
    }
    else
    {
        fmpz_poly_gcd(rop->den, op2->num, op1->den);
        if (fmpz_sgn(fmpz_poly_lead(rop->den)) < 0)
            fmpz_poly_neg(rop->den, rop->den);
        
        if (fmpz_poly_is_one(rop->den))
        {
            fmpz_poly_div(rop->den, op2->den, rop->num);
            fmpz_poly_mul(rop->den, op1->den, rop->den);
            fmpz_poly_div(rop->num, op1->num, rop->num);
            fmpz_poly_mul(rop->num, rop->num, op2->num);
        }
        else
        {
            fmpz_poly_init(t);
            fmpz_poly_init(u);
            fmpz_poly_div(t, op1->num, rop->num);
            fmpz_poly_div(u, op2->den, rop->num);
            fmpz_poly_div(rop->num, op2->num, rop->den);
            fmpz_poly_mul(rop->num, t, rop->num);
            fmpz_poly_div(rop->den, op1->den, rop->den);
            fmpz_poly_mul(rop->den, rop->den, u);
            fmpz_poly_clear(t);
            fmpz_poly_clear(u);
        }
    }
}
