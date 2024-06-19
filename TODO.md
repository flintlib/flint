# TODO

TODO: Update this TODO

## `fmpz`

* [maybe] figure out how to write robust test code for `fmpz_read` (which reads
  from `stdin`), perhaps using a pipe

* [maybe] Avoid the double allocation of both an `__mpz_struct` and limb data,
  having an `fmpz` point directly to a combined structure. This would require
  writing replacements for most `mpz` functions.


## `ulong_extras`

* in `is_prime_pocklington` allow the cofactor to be a perfect power not just
  prime

* factor out some common code between `n_is_perfect_power235` and
  `n_factor_power235`

* `n_mod2_preinv` may be slower than the chip on Core2 due to the fact that it
  can pipeline 2 divisions. Check all occurrences of this function and replace
  with divisions where it will speed things up on Core2. Beware, this will slow
  things down on AMD, so it is necessary to do this per architecture. The macros
  in `nmod_vec` will also be faster than the functions in `ulong_extras`, thus
  they should be tried first.

* add profile code for `factor_trial`, `factor_one_line`, `factor_SQUFOF`

* [maybe] make `n_factor_t` an array of length 1 so it can be passed by reference
  automatically, as per `mpz_t`'s, etc

* [enhancement] Implement a primality test which only requires factoring of
  `n - 1` up to `n^1/3` or `n^1/4`

* [enhancement] Implement a combined `p - 1` and `p + 1` primality test as per
  http://primes.utm.edu/prove/prove3_3.html

* [enhancement] Implement a quadratic sieve and use it in `n_factor` once things
  get too large for SQUFOF

* Claus Fieker suggested that in BPSW we should do Miller-Rabin rather than
  Fermat test. The advantage is it will throw out more composites before the
  Fibonacci or Lucas tests are run, but still ensures you have a base-2 probable
  prime. This should speed the test up.

## `long_extras`

* write and use `z_invert` in `fmpz_invert`


## `fmpz_factor`

* Tune Brent-Pollard Rho for optimal values

* Find optimal values for B1, B2 for ECM.

* Use multipoint polynomial evaluation in ECM stage II.

## `fmpz_mpoly`

* Implement array version of `divrem_ideal`.

* Don't restart division if upgrading to multiprecision.

* Check quotient and remainder coefficients are normalised the same way for
  small vs multiprecision division functions. Beware one normalisation may
  require an extra bit for either q or r in signed integer division in small
  case. This is probably why C division is normalised the way it is.

* Implement `fmpz_mpoly_ui/si/fmpz_sub` functions for Nemo.

* Add references to Monagan and Pearce papers to Bibtex and reference in heap
  multiplication, powering and division functions.

* Move functions in `mpoly/MISSING_FXNS.c` and decarations in mpoly.h to their
  appropriate places.

## `nmod_poly`

* Make some assembly optimisations to `nmod_poly` module. (FIXME: Which ones?)

* Add basecase versions of `log`, `sqrt`, `invsqrt` series.

* Add $O(M(n))$ powering $\mod x^n$ based on `exp` and `log`.

* Implement fast `mulmid` and use to improve Newton iteration

* Determine cutoffs for `compose_series_divconquer` for default use in
  `compose_series` (only when one polynomial is small).

* Optimise, write an underscore version of, and test `nmod_poly_remove`

* Improve `powmod` and `powpowmod` using precomputed Newton inverse and
  $2^k$-ary/sliding window powering.

* Maybe restructure the code in `factor.c`

* Add a (fast) function to convert an `nmod_poly_factor_t` to an expanded
  polynomial


## `fmpz_poly`

* add test code for `fmpz_poly_max_limbs`

* Improve the implementations of `fmpz_poly_[divrem/div/rem]`, check that the
  documentations still apply, and write test code for this --- all of this makes
  more sense once there is a choice of algorithms

* Include test code for `fmpz_poly_inv_series`, once this method does anything 
  better than aliasing `fmpz_poly_inv_newton`

* Sort out the `fmpz_poly_pseudo_[div/rem]` code.  Currently this is just a hack
  to call `fmpz_poly_pseudo_divrem`

* Fix the inefficient cases in `CRT_ui`, and move the relevant parts of this
  function to the `fmpz_vec` module

* Avoid redundant work and memory allocation in `fmpz_poly_bit_unpack`
  and `fmpz_poly_bit_unpack_unsigned`.

* Add functions for composition, multiplication and division by a monic linear
  factor, i.e. $P(x \pm c)$, $P \cdot (x \pm c)$, $P / (x \pm c)$.

* `xgcd_modular` is really slow. But it is not clear why. 1/3 of the time is
  spent in resultant, but the CRT code or the `nmod_poly_xgcd` code may also be
  to blame.

* Make resultants use fast GCD?

* In `fmpz_poly_pseudo_divrem_divconquer`, fix the repeated memory allocation of
  size $O(\mathtt{lenA})$ in the case when $\mathtt{lenA} \gg \mathtt{lenB}$.

## `fmpq_poly`

* add `fmpq_poly_fprint_pretty`

* Rewrite `_fmpq_poly_interpolate_fmpz_vec` to use the Newton form as done in
  the `fmpz_poly` interpolation function. In general this should be much faster.

* Add versions of `fmpq_poly_interpolate_fmpz_vec` for `fmpq` `y` values, and
  `fmpq` values of both `x` and `y`.

* Add `mulhigh`

* Add `gcdinv`

* Add discriminant

## `fmpz_mod_poly`

* Replace `fmpz_mod_poly_rem` by a proper implementation which actually saves
  work over `divrem`.  Then, also add test code.

* Implement a faster GCD function then the euclidean one, then make the wrapping
  GCD function choose the appropriate algorithm, and add some test code

## `fmpz_poly_mat`

* Tune multiplication cutoffs.

* Take sparseness into account when selecting between algorithms.

* Investigate more clever pivoting strategies in row reduction.


## `arith`

* Write a faster `arith_divisors` using the merge sort algorithm (see Sage's
  implementation). Alternatively (or as a complement) write a supernaturally
  fast `_fmpz_vec_sort`.

* Write tests for the `arith_hrr_expsum_factored` functions.

## `fmpz_mat`

* Add `fmpz_mat/randajtai2.c` based on Stehle's `matrix.cpp` in `fpLLL`
  (requires `mpfr_t`'s).

* Add element getter and setter methods.

* Write multiplication functions optimised for sparse matrices by changing the
  loop order and discarding zero multipliers.

* The Dixon $p$-adic solver currently spends most of the time computing integer
  matrix-vector products. Instead of using a single prime, it is likely to be
  faster to use a product of primes to increase the proportion of time spent on
  modular linear algebra. The code should also use fast radix conversion instead
  of building up the result incrementally to improve performance when the output
  gets large.

* Maybe optimise multimodular multiplication by pre-transposing so that
  transposed `nmod_mat` multiplication can be used directly instead of creating
  a transposed copy in `nmod_mat_mul`. However, this doesn't help in the
  Strassen range unless there also is a transpose version of
  `nmod_mat_mul_strassen`.

* Take sparseness into account when selecting between algorithms.

* Maybe simplify the interface for row reduction by guaranteeing
  that the denominator is the determinant.

## `fmpz_lll`

* Improve the wrapper strategy, if possible.

* Add an `mpf` version of `is_reduced` functions so that the dependency on MPFR
  can be removed.

* Componentize the `is_reduced` code so that the `R` computed during LLL can be
  reused.

## `nmod_mat`

* Write a Strassen for packed entries that does additions faster.

* Investigate why the constant of solving/rref/inverse compared to
  multiplication appears to be worse than in theory (recent paper by Jeannerod,
  Pernet and Storjohann).

* See if Strassen can be improved using combined addmul operations.

* Consider getting rid of the row pointer array, using offsets instead of window
  matrices. The row pointer is only useful for Gaussian elimination, but there
  we end up working with a separate permutation array anyway.

* Add functions for computing $A B^\mathrm{T}$ and $A^\mathrm{T} B$, using
  transpose multiplications directly to avoid creating a temporary copy.

* Maybe: add asserts to check that the modulus is a prime where this is assumed.

* The current `addmul`/`submul` functions are misnamed since they implement a
  more general operation.

* Improve `rref` and inverse to perform everything in-place.


## `fmpq`

* Add more functions for generating random numbers.


## `fmpq_mat`

* Add more random functions.

* Add a user-friendly function for LUP decomposition.

* Add a nullspace function.

## `padic`

* Add test code for the various output formats; perhaps in the form of examples?

* Implement `padic_val_fac` for generic inputs

## `mpf_[vec/mat]`

* `_mpf_vec_approx_equal` uses the bizarre `mpf` notion of equal; it should be
  either renamed `equal_bits` or the absolute value of the difference should be
  compared with zero

* The conditions used in `mpf_mat_qr/gso` (similarly for `d_mat` module) work,
  but don't match those used in the algorithm in the paper referenced in the
  docs. This is possibly because `mpf` doesn't do exact rounding. The tests could
  probably be improved.

## `d_mat`

* `d_mat_transpose` tries to be clever with the cache, but it could use L1 sized
  blocks and optimise for sequential reading/writing. It could also handle
  in-place much more efficiently.

* `d_mat_is_reduced` doesn't seem to make any guarantees about reducedness, just
  approximately checks the Lovasz conditions in double arithmetic. I think this
  is superseded now and can be removed. It is not used anywhere.

