#include "fmpz_poly_q.h"

void fmpz_poly_q_div_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    if (fmpz_poly_q_is_zero(op))
    {
        printf("ERROR (fmpz_poly_q_div_in_place).  Division by zero.\n");
        abort();
    }
    
    if (fmpz_poly_q_is_zero(rop))
        return;
    
    fmpz_poly_q_inv(rop, rop);
    fmpz_poly_q_mul_in_place(rop, op);
    fmpz_poly_q_inv(rop, rop);
}

void fmpz_poly_q_div(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    fmpz_poly_t t, u;
    
    if (fmpz_poly_q_is_zero(op2))
    {
        printf("ERROR (fmpz_poly_q_div).  Division by zero.\n");
        abort();
    }
    
    if (op1 == op2)
    {
        fmpz_poly_q_one(rop);
        return;
    }
    
    if (fmpz_poly_q_is_zero(op1))
    {
        fmpz_poly_q_zero(rop);
        return;
    }
    
    if (rop == op1)
    {
        fmpz_poly_q_div_in_place(rop, op2);
        return;
    }
    if (rop == op2)
    {
        fmpz_poly_q_div_in_place(rop, op1);
        fmpz_poly_q_inv(rop, rop);
        return;
    }
    
    /*
        From here on, we know that rop, op1 and op2 refer to distinct objects 
        in memory, and that op1 and op2 are non-zero rational functions
     */
    
    /*
        XXX:  Do not maintain the remaining part of the function separately!!!
              Instead, note that this is the same as the corresponding part of 
              the multiplication code, with op2->num and op2->den swapped.

              The only caveat to this is that we cannot assume the leading 
              coefficient of op2->num to be positive, and thus check for this 
              in the very end.
     */
    
    /* Are both denominators equal to one? */
    if (fmpz_poly_is_one(op1->den) && fmpz_poly_is_one(op2->num))
    {
        fmpz_poly_mul(rop->num, op1->num, op2->den);
        fmpz_poly_set_si(rop->den, 1);
        return;
    } 
    
    fmpz_poly_gcd(rop->num, op1->num, op2->num);
    
    if (fmpz_poly_is_one(rop->num))
    {
        fmpz_poly_gcd(rop->den, op2->den, op1->den);
        
        if (fmpz_poly_is_one(rop->den))
        {
            fmpz_poly_mul(rop->num, op1->num, op2->den);
            fmpz_poly_mul(rop->den, op1->den, op2->num);
        }
        else
        {
            fmpz_poly_div(rop->num, op2->den, rop->den);
            fmpz_poly_mul(rop->num, op1->num, rop->num);
            fmpz_poly_div(rop->den, op1->den, rop->den);
            fmpz_poly_mul(rop->den, rop->den, op2->num);
        }
    }
    else
    {
        fmpz_poly_gcd(rop->den, op2->den, op1->den);
        
        if (fmpz_poly_is_one(rop->den))
        {
            fmpz_poly_div(rop->den, op2->num, rop->num);
            fmpz_poly_mul(rop->den, op1->den, rop->den);
            fmpz_poly_div(rop->num, op1->num, rop->num);
            fmpz_poly_mul(rop->num, rop->num, op2->den);
        }
        else
        {
            fmpz_poly_init(t);
            fmpz_poly_init(u);
            fmpz_poly_div(t, op1->num, rop->num);
            fmpz_poly_div(u, op2->num, rop->num);
            fmpz_poly_div(rop->num, op2->den, rop->den);
            fmpz_poly_mul(rop->num, t, rop->num);
            fmpz_poly_div(rop->den, op1->den, rop->den);
            fmpz_poly_mul(rop->den, rop->den, u);
            fmpz_poly_clear(t);
            fmpz_poly_clear(u);
        }
    }
    
    /* XXX:  Check that the numerator has the appropriate sign. */
    if (fmpz_sgn(fmpz_poly_lead(rop->den)) < 0)
    {
        fmpz_poly_neg(rop->num, rop->num);
        fmpz_poly_neg(rop->den, rop->den);
    }
}

