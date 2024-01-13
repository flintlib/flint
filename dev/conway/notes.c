/* prime < 260, degrees are decaying but always less than 410. Coefficients can
   always be stored in uint8_t. Coefficient sequences are always on the form

       [a0, a1, a2, ..., an, 0, 0, ..., 0, 1].

   Hence, store prime, degree, number of "non-trivial" coefficients, as well as
   the "non-trivial" coefficients. */

/* 260 < prime < 300, degrees = 1, ..., 12. Polynomials are on the form

       [b0,   1],
       [a0,  b1,   1],
       [b0,  a1,   0,   1],
       [a0,  b2,  a2,   0,   1],
       [b0,  a3,   0,   0,   0,   1],
       [a0,  b3,  b4,  a4,  a5,   0,   1],
       [b0,  a6,   0,   0,   0,   0,   0,   1],
       [a0,  b6,  b7,  b8,  a7,   0,   0,   0, 1],
       [b0,  b9, b10,  a8,   0,   0,   0,   0, 0, 1],
       [a0, b11, b12, b13, b14, b15,   0,   0, 0, 0, 1],
       [b0,  a9,   0,   0,   0,   0,   0,   0, 0, 0, 0, 1],
       [a0, a10, b16, a11, b17, a12, a13, a14, 0, 0, 0, 0, 1],

   where a0, ..., a14 < 2^8 and b0, ..., b10 < 2^16. Store prime as uint8_t with
   offset of 2^8, and store coefficients as 15 uint8_t with 11 uint16_t. */

/* 300 < prime < 1000, degrees = 1, ..., 9.
   Polynomials are on the form

       [b0,  1],
       [a0, b1,   1],
       [b0, a1,   0,  1],
       [a0, b2,  a2,  0,  1],
       [b0, a3,   0,  0,  0, 1],
       [a0, b3,  b4, b5, a4, 0, 1],
       [b0, a5,   0,  0,  0, 0, 0, 1],
       [a0, b6,  b7, b8, a6, 0, 0, 0, 1],
       [b0, b9, b10, a7,  0, 0, 0, 0, 0, 1],

   where a0, ..., a7 < 2^8 and b0, ..., b10 < 2^16. Store prime as uint16_t, and
   store coefficients as 8 uint8_t with 11 uint16_t. */

/* 1000 < prime < 3371, degrees = 1, ..., 7, 9. Some exceptions exists where
   degree 7 and 9 does not exist. These are p = 2689, 2797, 2833, 3019, 3163,
   3209 and 3331. Polynomials are on the form

       [b0,  1],
       [a0, b1,  1],
       [b0, a1,  0,  1],
       [a0, b2, a2,  0,  1],
       [b0, a3,  0,  0,  0, 1],
       [a0, b3, b4, b5, a4, 0, 1],
       [b0, a5,  0,  0,  0, 0, 0, 1],
       [b0, b6, b7, a6,  0, 0, 0, 0, 0, 1],

   where a0, ..., a6 < 2^8 and b0, ..., b7 < 2^16. Store prime as uint16_t, and
   store coefficients as 7 uint8_t with 8 uint16_t. */

/* 3371 <= prime < 11000, degrees = 1, ..., 6.
   Polynomials are on the form

       [b0,  1],
       [a0, b1,  1],
       [b0, a1,  0,  1],
       [a0, b2, a2,  0,  1],
       [b0, a3,  0,  0,  0, 1],
       [a0, b3, b4, b5, a4, 0, 1],

   where a0, ..., a4 < 2^8 and b0, ..., b5 < 2^16. Store prime as uint16_t, and
   store coefficients as 5 uint8_t with 6 uint16_t. */

/* 11000 < prime < 2^16, degrees = 1, ..., 4.
   Polynomials are on the form

       [b0, 1],
       [a0, b1, 1],
       [b0, a1, 0, 1],
       [a0, b2, a2, 0, 1],

   where a0, ..., a2 < 2^8 and b0, ..., b2 < 2^16. Store prime as uint16_t, and
   store each prime as 3 uint8_t with 3 uint16_t. */

/* 2^16 < prime, then degrees are all four.
   Polynomials are on the form

       [a0, b0, a1, 0, 1],

   where a0, a1 < 2^8 and b0 < 2^32. Store prime as uint16_t where everything
   primes are offset by 2^16, and store coefficients as 2 uint8_t with 1
   uint32_t. */
