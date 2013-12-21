#ifndef FMPZ_CONVERSIONS_H
#define FMPZ_CONVERSIONS_H

/* turn a pointer to an __mpz_struct into a fmpz_t */
#define PTR_TO_COEFF(x) (((ulong) (x) >> 2) | (WORD(1) << (FLINT_BITS - 2)))

/* turns an fmpz into a pointer to an mpz */
#define COEFF_TO_PTR(x) ((__mpz_struct *) ((x) << 2))

#endif /* FMPZ_CONVERSIONS_H */
