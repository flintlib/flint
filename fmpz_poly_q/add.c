#include "fmpz_poly_q.h"

void fmpz_poly_q_add_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    fmpz_poly_t d, poly, r2, s2;

    if (rop == op)
    {
        fmpz_poly_q_scalar_mul_si(rop, rop, 2);
        return;
    }

    /* Is one of the numerators zero? */
    if (fmpz_poly_is_zero(rop->num))
    {
        fmpz_poly_q_set(rop, op);
        return;
    }
    if (fmpz_poly_is_zero(op->num))
    {
        return;
    }
    
    /* Is one or are both denominators equal to one? */
    if (fmpz_poly_is_one(rop->den))
    {
        if (fmpz_poly_is_one(op->den))
        {
            fmpz_poly_add(rop->num, rop->num, op->num);
            return;
        }
        else
        {
            fmpz_poly_mul(rop->num, rop->num, op->den);
            fmpz_poly_add(rop->num, rop->num, op->num);
            fmpz_poly_set(rop->den, op->den);
            return;
        }
    }
    if (fmpz_poly_is_one(op->den))
    {
        fmpz_poly_init(poly);
        fmpz_poly_mul(poly, rop->den, op->num);
        fmpz_poly_add(rop->num, rop->num, poly);
        fmpz_poly_clear(poly);
        return;
    }
    
    /* Henrici's algorithm for summation in quotient fields */
    fmpz_poly_init(d);
    fmpz_poly_gcd(d, rop->den, op->den);
    if (fmpz_sgn(fmpz_poly_lead(d)) < 0)
        fmpz_poly_neg(d, d);
    
    if (fmpz_poly_is_one(d))
    {
        fmpz_poly_mul(rop->num, rop->num, op->den);
        fmpz_poly_mul(d, rop->den, op->num);  /* Using d as temp */
        fmpz_poly_add(rop->num, rop->num, d);
        fmpz_poly_mul(rop->den, rop->den, op->den);
    }
    else
    {
        fmpz_poly_init(r2);
        fmpz_poly_init(s2);
        
        fmpz_poly_div(r2, rop->den, d);
        fmpz_poly_div(s2, op->den, d);
        
        fmpz_poly_mul(rop->num, rop->num, s2);
        fmpz_poly_mul(s2, op->num, r2);  /* Using s2 as temp */
        fmpz_poly_add(rop->num, rop->num, s2);
        
        if (fmpz_poly_is_zero(rop->num))
        {
            fmpz_poly_zero(rop->den);
            fmpz_poly_set_coeff_si(rop->den, 0, 1);
        }
        else
        {
            fmpz_poly_mul(rop->den, r2, op->den);
            
            fmpz_poly_gcd(r2, rop->num, d);  /* Using r2 as temp */
            if (fmpz_sgn(fmpz_poly_lead(r2)) < 0)
                fmpz_poly_neg(r2, r2);
            
            if (!fmpz_poly_is_one(r2))
            {
                fmpz_poly_div(rop->num, rop->num, r2);
                fmpz_poly_div(rop->den, rop->den, r2);
            }
        }
        fmpz_poly_clear(r2);
        fmpz_poly_clear(s2);
    }
    fmpz_poly_clear(d);
}

/**
 * \ingroup  Addition
 * 
 * Sets \c rop to the sum of \c op1 and \c op2.
 */
void 
fmpz_poly_q_add(fmpz_poly_q_t rop, 
                const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    fmpz_poly_t d, r2, s2;
    
    /* If op1 == op2, return 2 * op1 */
    if (op1 == op2)
    {
        fmpz_poly_q_scalar_mul_si(rop, op1, 2);
        return;
    }
    
    /* Is rop == op1 or rop == op2? */
    if (rop == op1)
    {
        fmpz_poly_q_add_in_place(rop, op2);
        return;
    }
    if (rop == op2)
    {
        fmpz_poly_q_add_in_place(rop, op1);
        return;
    }
    
    /*
        From here on, we may assume that rop, op1 and op2 all refer to 
        distinct objects in memory, although they may still be equal
     */
    
    if (fmpz_poly_is_zero(op1->num))
    {
        fmpz_poly_q_set(rop, op2);
        return;
    }
    if (fmpz_poly_is_zero(op2->num))
    {
        fmpz_poly_q_set(rop, op1);
        return;
    }
    
    /* Is one or are even both denominators equal to one? */
    if (fmpz_poly_is_one(op1->den))
    {
        if (fmpz_poly_is_one(op2->den))
        {
            fmpz_poly_add(rop->num, op1->num, op2->num);
            fmpz_poly_zero(rop->den);
            fmpz_poly_set_coeff_si(rop->den, 0, 1);
            return;
        }
        else
        {
            fmpz_poly_mul(rop->num, op1->num, op2->den);
            fmpz_poly_add(rop->num, rop->num, op2->num);
            fmpz_poly_set(rop->den, op2->den);
            return;
        }
    }
    if (fmpz_poly_is_one(op2->den))
    {
        fmpz_poly_mul(rop->num, op2->num, op1->den);
        fmpz_poly_add(rop->num, op1->num, rop->num);
        fmpz_poly_set(rop->den, op1->den);
        return;
    }
    
    /* Henrici's algorithm for summation in quotient fields */
    
    /*
        We begin by using rop->num as a temporary variable for 
        the gcd of the two denominators' greatest common divisor
     */
    fmpz_poly_gcd(rop->num, op1->den, op2->den);
    if (fmpz_sgn(fmpz_poly_lead(rop->num)) < 0)
        fmpz_poly_neg(rop->num, rop->num);
    
    if (fmpz_poly_is_one(rop->num))
    {
        fmpz_poly_mul(rop->num, op1->num, op2->den);
        fmpz_poly_mul(rop->den, op1->den, op2->num);  /* Using rop->den as temp */
        fmpz_poly_add(rop->num, rop->num, rop->den);
        fmpz_poly_mul(rop->den, op1->den, op2->den);
    }
    else
    {
        /*
            We now copy rop->num into a new variable d, so we 
            no longer need rop->num as a temporary variable
         */
        fmpz_poly_init(d);
        fmpz_poly_swap(d, rop->num);
        
        fmpz_poly_init(r2);
        fmpz_poly_init(s2);
        
        fmpz_poly_div(r2, op1->den, d);  /* +ve leading coeff */
        fmpz_poly_div(s2, op2->den, d);  /* +ve leading coeff */
        
        fmpz_poly_mul(rop->num, op1->num, s2);
        fmpz_poly_mul(rop->den, op2->num, r2);  /* Using rop->den as temp */
        fmpz_poly_add(rop->num, rop->num, rop->den);
        
        if (fmpz_poly_is_zero(rop->num))
        {
            fmpz_poly_zero(rop->den);
            fmpz_poly_set_coeff_si(rop->den, 0, 1);
        }
        else
        {
            fmpz_poly_mul(rop->den, op1->den, s2);
            
            fmpz_poly_gcd(r2, rop->num, d);
            if (fmpz_sgn(fmpz_poly_lead(r2)) < 0)
                fmpz_poly_neg(r2, r2);
            
            if (!fmpz_poly_is_one(r2))
            {
                fmpz_poly_div(rop->num, rop->num, r2);
                fmpz_poly_div(rop->den, rop->den, r2);
            }
        }
        
        fmpz_poly_clear(d);
        fmpz_poly_clear(r2);
        fmpz_poly_clear(s2);
    }
}
