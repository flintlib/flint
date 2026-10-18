#include <math.h>
#include "double_extras.h"
#include "fft_small.h"

/*
    For arithmetic mod n, the product a*b is represented exactly as the usual
    double-double format h + l where h and l are prec-bit floating point
    numbers. a*b will be roughly the size of n^2, but since one or both of a or
    b may not be reduced modulo n, it is common for the product a*b to be a
    reasonable multiple of n^2. Since this is calculated as

      h = mul(a, b)
      l = fmsub(a, b, h),

    if we have |a*b| < 2^e <= 2^(2*prec), then

      |l| <= 2^(e - prec - 1).

    First, ninv is an floating approximation of 1/n with some error epsilon

      ninv := 1/n + epsilon.

    Next, we would like to approximate the true quotient a*b/n as

      q = round_int(h*ninv),

    so that h + l can be reduced mod n. Since

      h*ninv = (a*b - l)*(1/n + epsilon)
               = a*b/n + a*b*epsilon - l*ninv,

    we have (with |h*ninv - q| <= q_error)

(*)   |a*b/n - q| <= |a*b*epsilon| + 2^(e - prec - 1)*ninv + q_error.

    When the quantity on the rhs is <= k, we can guarantee that

      |h + l - q*n| <= k*n,

    while h + l - q*n can be calculated exactly as add(l, fnmadd(q, n, h)).

    *** q_error ***
    when mul(h, ninv) does rounding B bits after the units place, the additional
    rounding done by round(mul(h, ninv)) can introduce an error as large as

      q_error = 1/2 + 1/2^(B+1)

    For example, suppose

             h*ninv = XXXXXX1.0100001....
       mul(h, ninv) = XXXXXX1.1        rounded to 1 bit to the right of the dot

    Now mul(h,ninv) rounds up although h*ninv rounds down. The difference between
    the final rounded integer and the original h*ninv is at most 1/2 + 1/4 = 3/4

    We can get the optimal q_error = 1/2 by calculating round_int(h*ninv) as

      sub(fmadd(h, ninv, r), r)     with r = 3*2^(prec - 2)

    However, this requires |h*ninv| < 2^(prec - 2), which is not always true.
    Furthermore, round(mul(.,.)) is fast on AMD :)

    Illustration of the truth of 1/2 + 1/2^(B+1): suppose mul(h, ninv) rounds
    to B=2 bits after the dot (XXXX means some non-zero bits):
      h*ninv      mul(h*ninv) round(mul(h*ninv))   bound on |h*ninv - round(mul(h*ninv))|
      0.000       0.00        0.                   0
      0.000XXXX   0.00        0.                   1/8
      0.001       0.00        0.                   1/8
      0.001XXXX   0.01        0.                   2/8
      0.010       0.01        0.                   2/8
      0.010XXXX   0.01        0.                   3/8
      0.011       0.10        0.  or  1.           3/8 or 5/8
      0.011XXXX   0.10        0.  or  1.           4/8 or 5/8
      0.100       0.10        0.  or  1.           4/8
      0.100XXXX   0.10        0.  or  1.           5/8 or 4/8
      0.101       0.10        0.  or  1.           5/8 or 3/8
      0.101XXXX   0.11        1.                   3/8
      0.110       0.11        1.                   2/8
      0.110XXXX   0.11        1.                   2/8
      0.111       1.00        1.                   1/8
      0.111XXXX   1.00        1.                   1/8

    We see q_error = 5/8 = 1/2 + 1/2^(B+1)
*/
/*
Analysis for 2^(prec - 4) < n < 2^(prec - 3) (aka 50 bit primes)
----------------------------------------------------------------

    In this case we have |epsilon| < 2^(3 - 2*prec) and ninv < 2^(4 - prec).

    For |a*b| < 2*n^2, we can take e = 2*prec - 5.

    Then, |h*ninv| < 2^(prec-2) and mul(h, ninv) does rounding at least
    two bits to the right of the dot. Then, we may take q_error = 1/2 + 2^-3:

      |a*b*epsilon| + 2^(e - prec - 1)*ninv + 1/2+2^-3 < 2^-2 + 2^-2 + 1/2 + 2^-3 = 1+1/8

    Similarly for |a*b| < 4*n^2, we can take e = 2*prec - 4 and q_error = 1/2 + 2^-2:

      |a*b*epsilon| + 2^(e - prec - 1)*ninv + 1/2+2^-2 < 2^-1 + 2^-1 + 1/2 + 2^-2 = 3/2+1/4

    Therefore,

     Products in the range (-2*n^2, 2*n^2) are reduced to the range (-9/8*n, 9/8*n).
     Products in the range (-4*n^2, 4*n^2) are reduced to the range (-7/4*n, 7/4*n).

    These pessimistic bounds are not good enough. For specific values of n, epsilon can be
    much smaller, and it is not very restrictive to assume what mulmod needs
    (checked by fft_small_mulmod_satisfies_bounds below, and required by
    sd_fft_ctx_init_prime):

(1)  Products in the range (-2*n^2, 2*n^2) are reduced to the range (-1*n, 1*n).
(2)  Products in the range (-4*n^2, 4*n^2) are reduced to the range (-3/2*n, 3/2*n).
*/
/*
    define e as |a*b| < 2^e
    bound = |a*b*(ninv-1/n)| + 2^(e - prec - 1)*ninv + q_error
    want bound < 1   for |a*b| < 2*n^2  where  q_error = 1/2 + 2^-3
    want bound < 1.5 for |a*b| < 4*n^2  where  q_error = 1/2 + 2^-2
*/

int fft_small_mulmod_satisfies_bounds(ulong nn)
{
    double n = nn;
    double ninv = 1.0/n;
    double t1 = fabs(fma(n, ninv, -1.0));  /* epsilon ~= t1/n  good enough */
    double limit2, limit4;
    int B, ok, n1bits, n2bits;
    ulong n2hi, n2lo;

    n1bits = n_nbits(nn);
    umul_ppmm(n2hi, n2lo, nn, nn);
    if (n2hi != 0)
       n2bits = FLINT_BITS + n_nbits(n2hi);
    else
       n2bits = n_nbits(n2lo);

    /* for |a*b| < 2*n^2*/
    /* |h*n_inv| < 2*n, so rounding in mul(h, ninv) at least B bits after the .*/
    B = D_BITS - n1bits - 1;
    if (B < 2)
       return 0;
    limit2 = 2*n*t1 + ldexp(ninv, 1+n2bits-D_BITS-1) + 0.5 + ldexp(1.0, -(B+1));

    /* for |a*b| < 4*n^2*/
    B -= 1;
    limit4 = 4*n*t1 + ldexp(ninv, 2+n2bits-D_BITS-1) + 0.5 + ldexp(1.0, -(B+1));

    /* fudge the limits 1 and 3/2 because the above is double arithmetic */
    ok = (limit2 < 0.99) && (limit4 < 1.49);
    return ok;
}

/*

assuming satisfies_bounds(n):

(1)  Products in the range (-2*n^2, 2*n^2) are reduced to the range (-1*n, 1*n).
(2)  Products in the range (-4*n^2, 4*n^2) are reduced to the range (-3/2*n, 3/2*n).

(1) also holds for mulmod_fast, whose quotient is at least as accurate,
provided |a*b| <= 2^51*n (implied by |a*b| < 2*n^2 for n < 2^50).
(2) does not: 4*n^2 exceeds 2^51*n for n > 2^49, see machine_vectors.h.


Analysis of the forward butterflies for 50 bit primes
-----------------------------------------------------

        b0 = 1*a0 + w^2*a2 +   w*(a1 + w^2*a3)
        b1 = 1*a0 + w^2*a2 -   w*(a1 + w^2*a3)
        b2 = 1*a0 - w^2*a2 + i*w*(a1 - w^2*a3)
        b3 = 1*a0 - w^2*a2 - i*w*(a1 - w^2*a3)

    Each of w^2, w, and i*w is in (-n/2, n/2). As part of this operation, the 
    multiplication 1*a0 reduces a0 to the range (-n, n) and the other
    multiplications use the above. In practice, a0 does not get large even
    without this reduction, but this is needed to prove things stay bounded.
    The following claim is trivial from (1):

(*)     If the ai in (-3*n, 3*n), then the bi are also in (-3*n, 3*n).

    One exception: LENGTH8_ZERO_J in sd_fft.c, which is the whole transform
    for depth 3 (sd_fft_basecase_3_1), returns in x0 and x1 sums of four
    outputs of reduce_to_pm1n, so in (-4*n, 4*n) by its documented range
    [-n, n]. With the actual range of reduce_to_pm1n (|result| <= n/2 + 2,
    see machine_vectors.h) they are below 2*n + 8, so sd_fft_trunc does
    map (-3*n, 3*n) into itself at every depth.


Analysis of the truncated forward butterflies for 50 bit primes
---------------------------------------------------------------

    Trivial because some of the ai are just 0.


Analysis of the pointwise muls
------------------------------

    The output of one fft in the range (-3*n, 3*n) is multiplied with a
    multiplier (2^-depth) in the range (-n/2, n/2). This produces a product x
    in the range (-n, n). Then, x is multiplied by the output of another fft
    in the range (-3*n, 3*n). Thus,

(*)     The outputs of the pointwise muls are in the range (-3/2*n, 3/2*n).


Analysis of the reverse butterflies for 50 bit primes
-----------------------------------------------------

        4a0 =    1*(     (b0 + b1) +        (b2 + b3))
        4a2 = w^-2*(     (b0 + b1) -        (b2 + b3))
        4a1 =       w^-1*(b0 - b1) - i*w^-1*(b2 - b3)
        4a3 = w^-2*(w^-1*(b0 - b1) + i*w^-1*(b2 - b3))

    As before, each of w^-2, w^-1, and i*w^-1 is in (-n/2, n/2). As part of
    this operation, 1*(b0 + b1 + b2 + b3) reduces 4a0 to the range (-n, n),
    and the other multiplications are as before.
    The following claim is trivial from (1) and (2):

(*)     If the bi are in (-2*n, 2*n), then the 4ai are also in (-2*n, 2*n).

    The bound n < 2^50 guarantees that |(b0 + b1) +- (b2 + b3)| < 2^53 so that
    no bits are lost.


Analysis of the truncated reverse butterflies for 50 bit primes
---------------------------------------------------------------

    The analogue of the claim above, outputs in (-2*n, 2*n) if the inputs
    are, does *not* follow from (1), (2) and the documented range [-n, n]
    of reduce_to_pm1n. With these bounds alone, five truncated butterflies
    of sd_ifft.c leave (-2*n, 2*n):

        radix_4_moth_inv_trunc_block_3_4_1      X0 up to 5n, X3 up to 3n
        radix_4_moth_inv_trunc_block_3_4_0      X0 up to 5n
        radix_4_moth_inv_trunc_block_3_3_1      X0 up to 4n
        radix_4_moth_inv_trunc_block_3_3_0      X0 up to 4n
        radix_4_moth_inv_trunc_block_1_4_0      X0 up to 3n

    (e.g. X0 = 2*c - v + u + p in 3_4_x, with c and u outputs of
    reduce_to_pm1n, v and p outputs of mulmod_fast), and through the
    recursion these values reach mulmod_fast with products up to 3*n^2
    and sums up to 11*n, beyond what mulmod_fast and exact addition
    (below 2^53, that is 8*n for n close to 2^50) allow.
    So the documented ranges do not prove sd_ifft_trunc correct.

    It is correct thanks to the actual ranges of the primitives, see
    machine_vectors.h: reduce_to_pm1n returns at most n/2 + 2^-52*|a|, and
    mulmod and mulmod_fast at most n/2 + 3*2^-53*|a*b| (n/2 + 2^-52*|a*b|
    for mulmod_fast with the magic constant).

    This was checked by running sd_ifft_trunc on bounds instead of values:
    each entry holds a bound on the magnitude of the integer it stands for,
    additions and multiplications act on the bounds, and each primitive
    checks its precondition and returns the smaller of its documented
    range ((1), (2), [-n, n]) and its actual range above. The check does
    not depend on the prime: every twiddle factor other than +-1 is
    replaced by the bound (n - 1)/2, the error |n*ninv - 1| is set to its
    maximum 2^-53, and n is taken close to 2^50, which maximizes the
    rounding errors relative to n; the conditions |a*b| < 2*n^2 and sums
    below 8*n are then sufficient for every n < 2^50. So the result
    holds for every prime that fft_small accepts. For inputs of magnitude
    at most 2*n - 1, every truncation, depths up to 25 and both variants
    of round_mul:

     - every product given to mulmod_fast is below 2*n^2 (hence below
       2^51*n), and every product given to mulmod below 4*n^2;
     - every floating-point addition is below 8*n, hence below 2^53; the
       only larger results are multiplications by 4, which are exact and
       immediately followed by reduce_to_pm1n;
     - the values passed between butterflies are below 3.5*n, and the
       outputs at most 2.5*n (about 2.2*n for the default primes, whose
       twiddle factors and error |n*ninv - 1| are smaller).

    So the results are correct, but outputs outside (-2*n, 2*n) do occur:
    for instance 2.04*n in radix_4_moth_inv_trunc_block_3_3_0 at depth 10
    with trunc = 768, n = 0x0003f00000000001 and suitable inputs in [0, n),
    when mulmod_fast returns slightly more than n/2 and the two reduced
    values are close to n/2. The documented output range of sd_ifft_trunc
    is therefore (-3*n, 3*n), which is also the input range of
    sd_fft_trunc; the callers in FLINT multiply the outputs by a reduced
    scaling factor or reduce them before anything else.

    Reducing the five outputs listed above with reduce_to_pm1n would give
    (-2*n, 2*n) from the documented ranges alone, at a cost of about 1%
    of sd_ifft_trunc (median over truncations, up to 4%).
*/

int fft_small_prime_supports_depth(ulong n, ulong depth)
{
    if (n <= 2 || n > (UWORD(1) << 50))
        return 0;

    if (!fft_small_mulmod_satisfies_bounds(n))
        /* not implied by the 50-bit bound: on more than 100000 tests with
           random 50-bit primes, about 1 in 900 fails, e.g. 1104645587927041
           (however, none of 49 bits or less failed) */
        return 0;

    if (depth > SD_FFT_CTX_W2TAB_SIZE)
        /* unlikely, the transform length would have to be massive */
        return 0;

    /* sd_fft_ctx_init_prime also requires 2^4 | n - 1 */
    return n_trailing_zeros(n - 1) >= FLINT_MAX(depth, 4);
}
