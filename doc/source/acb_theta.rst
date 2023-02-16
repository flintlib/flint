.. _acb-theta:

**acb_theta.h** -- Theta functions in any dimension
===============================================================================

This module provides methods for numerical evaluation of theta functions in any
dimension `g`. All methods also accept `g=1`, duplicating functionality from
:ref:`acb_modular.h <acb-modular>` (without the specific speed-ups).

In the context of this module, *tau* or `\tau` always denotes an element of the
Siegel complex upper half-space `\mathbb{H}_g = \{\tau \in
\operatorname{Mat}_{g\times g}(\mathbb{C}) : \tau^t = \tau, \quad
\operatorname{Im}(\tau) \text{ is positive definite}\}`.

For each `a,b\in \{0,1\}^g`, the Riemann theta function is the following
analytic function in two variables `\tau\in \mathbb{H}_g` and `z\in
\mathbb{C}^g`:

    .. math ::
    
        \theta_{a,b}(z,\tau) = \sum_{n\in a/2 + \mathbb{Z}^{g}} \exp(\pi i n^T\tau n + 2\pi i n^T (z + b/2))

considering `a, b, z` as column vectors. The pair `(a,b)` is called a theta
characteristic.

When handling vectors of theta values, the value of `\theta_{a,b}` always
appear at index *ab* (concatenation). Note that this convention is *not* the
same as the one chosen in :ref:`acb_modular.h <acb-modular>`: indeed we order
the vector of genus 1 theta values as `\theta_3,\theta_4,\theta_2,\theta_1` in
this order. We encode *ab* as an :type:`ulong` of length *2g*, allowing us to
work with theta functions up to genus at least 32 on 64-bit machines.

The main focus of this module is the efficient evaluation in different
situations, indicated by combinations of suffixes from the following
categories:

1. Choice of algorithm:
    * Naive algorithm: suffix :func:`naive`.
    * Newton's method and the AGM (quasi-linear in the required precision):
      suffix :func:`newton`.
    * Uniform algorithm (when available): suffix :func:`unif`.
2. Number of theta values:
    * All values `\theta_{0,b}` for `b\in \{0,1\}^g`: default (no suffix).
    * All values `\theta_{a,b}` for all *a,b*: suffix :func:`all`.
    * Individual value `\theta_{a,b}` for specified *a,b*: suffix :func:`ind`.
3. Value of *z*:
    * `z=0` (theta constants): suffix :func:`const`. The result is zero
      whenever `a^T b` is odd.
    * Specified *z*: default (no suffix). Some functions accept several vectors
      *z* simultaneously: in this case an extra argument *nb_z* is provided.
4. Theta values taken at `\tau/2` instead of `tau`: suffix :func:`half`.
5. Projective theta values (i.e., the result is defined up to simultaneous
   multiplication by a nonzero complex number): suffix :func:`proj`.
6. Squared theta values: suffix :func:`sqr`.
7. Also compute derivatives of theta functions up to some order: suffix
   :func:`jet`.

As usual, the numerical functions in this module compute strict error bounds:
if *tau* is represented by an :type:`acb_mat_t` which is not certainly positive
definite, the output will have an infinite radius.

The Siegel modular group
-------------------------------------------------------------------------------

We use the type `fmpz_mat_t` directly for matrices in `\operatorname{Sp}_{2g}(\mathbb{Z})`
or `\operatorname{GSp}_{2g}(\mathbb{Z})`. We always assume that the input
matrix *mat* is square of even size `2g`.

.. function:: void sp2gz_dim(const fmpz_mat_t mat)

    Returns the dimension `g`, which is half the number of rows (or columns)
    of *mat*.

.. function:: void sp2gz_get_a(fmpz_mat_t res, const fmpz_mat_t mat)

.. function:: void sp2gz_get_b(fmpz_mat_t res, const fmpz_mat_t mat)

.. function:: void sp2gz_get_c(fmpz_mat_t res, const fmpz_mat_t mat)

.. function:: void sp2gz_get_d(fmpz_mat_t res, const fmpz_mat_t mat)

    Sets *res* to the corresponding block of *mat*, written as `\left(\begin{textmatrix} a&b\\c&d \end{textmatrix}\right)`.

.. function:: void sp2gz_set_abcd(fmpz_mat_t m, const fmpz_mat_t a, const
              fmpz_mat_t b, const fmpz_mat_t c, const fmpz_mat_t d)

    Sets *mat* to `\left(\begin{textmatrix} a&b\\c&d \end{textmatrix}\right)`,
    where `a,b,c,d` are `g\times g` blocks.

.. function:: int sp2gz_is_correct(const fmpz_mat_t mat)

.. function:: int sp2gz_is_gsp(const fmpz_mat_t mat)

    Returns whether *mat* is an element of `\operatorname{Sp}_{2g}(\mathbb{Z})`
    or `\operatorname{GSp}_{2g}(\mathbb{Z})`, respectively.

.. function:: int sp2gz_is_scalar(const fmpz_mat_t mat)

    Returns whether *mat* is a scalar matrix, i.e. diagonal with equal entries.
    
.. function:: void sp2gz_j(fmpz_mat_t mat)

    Sets *mat* to the symplectic matrix
    `J = \left(\begin{textmatrix} 0&I_g\\-I_g&0 \end{textmatrix}\right)`.

.. function:: void sp2gz_block_diag(fmpz_mat_t mat, const fmpz_mat_t U)

    Sets *mat* to the symplectic matrix
    `\left(\begin{textmatrix} U&0\\0&U^{-T} \end{textmatrix}\right)`.
    Requires that `U\in \operatorname{GL}_g(\mathbb{Z})`.

.. function:: void sp2gz_trig(fmpz_mat_t mat, const fmpz_mat_t S)

    Sets *mat* to `\left(\begin{textmatrix} I_g&S\\0&I_g \end{textmatrix}\right)`,
    which is symplectic if and only if *S* is symmetric.

.. function:: void sp2gz_randtest(fmpz_mat_t mat, flint_rand_t state, slong bits)

    Sets *mat* to a random symplectic matrix whose coefficients have length
    approximately *bits*.

.. function:: sp2gz_nb_fundamental(slong g)

    Returns the number of fundamental symplectic matrices used in the reduction
    algorithm on `\mathbb{H}_g`. This number is currently `19` when `g=2` and
    `1` otherwise.

.. function:: void sp2gz_fundamental(fmpz_mat_t mat, slong j)

    Sets *mat* to the `j^{\text{th}}` fundamental symplectic matrix as defined
    above.

The Siegel upper half space
-------------------------------------------------------------------------------

The Siegel upper half space `\mathbb{H}_g` contains the standard fundamental
domain `\mathcal{F}_g`, defined in..., as a closed subset.

.. function:: void acb_siegel_cocycle(acb_mat_t res, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Sets *res* to `c\tau+d` where *c,d* are the lower `g\times g` blocks of
    *mat*.

.. function:: void acb_siegel_transform(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t tau, slong prec)

    Sets *res* to `(a\tau + b)(c\tau + d)^{-1}` where *a,b,c,d* are the
    `g\times g` blocks of *mat*.

.. function:: void acb_siegel_transform_ext(acb_ptr r, acb_mat_t w, const fmpz_mat_t mat, acb_srcptr z, const acb_mat_t tau, slong prec)

    Sets *r* and *w* to `(c\tau + d)^{-T} z` and `(a\tau + b)(c\tau + d)^{-1}`
    respectively, where *a,b,c,d* are the `g\times g` blocks of *mat*.

.. function:: void acb_siegel_reduce_imag(fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Reduces the imaginary part of *tau* by calling :func:`arb_mat_spd_lll_reduce`
    and sets *mat* to the corresponding unimodular transformation.

.. function:: void acb_siegel_reduce_real(fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Computes a symmetric, integral matrix *mat* such that *tau+mat* has a small
    real part, ideally at most `1/2` in absolute value for each coefficient.

.. function:: void acb_siegel_reduce(acb_mat_t res, fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Computes a symplectic matrix *mat* such that the result *res* of *mat*
    acting on *tau* is closer to `\mathcal{F}_g`, by repeatedly reducing its
    real and imaginary parts and applying fundamental symplectic matrices.
    
.. function:: void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random matrix *tau* in `\mathbb{H}_g`, possibly far from the
    fundamental domain.

.. function:: void acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random matrix *tau* in `\mathbb{H}_g` that is close to the
    fundamental domain by calling :func:`acb_siegel_reduce` on a random matrix.

.. function:: void acb_siegel_randtest_nice(acb_mat_t tau, flint_rand_t state, slong prec)

    Generates a random matrix that is well in the interior of `\mathcal{F}_g`.

AGM sequences
-------------------------------------------------------------------------------

The classical arithmetic-geometric mean (AGM) of two positive real numbers
admits a generalization to tuples of `2^g` complex numbers: see for
instance... We look at sequences in which each step takes the form

    .. math::

        (x_b)_{b\in (\mathbb{Z}/2\mathbb{Z})^g \mapsto (y_b)_{b\in (\mathbb{Z}/2\mathbb{Z})^g}

where

    .. math::
        
        y_b = \sum_{b'\in (\mathbb{Z}/2\mathbb{Z})^g} r_{b'} r_{b+b'}

for some choice of square roots `(r_b)` of the tuple `(x_b)`. In this
generality, AGM sequences converge quadratically if and only if the chosen
square roots `r_b` are eventually always in *good position*, i.e. they all
belong to a common quarter plane seen from the origin.

Following..., we also compute *extended Borchardt sequences*, defined by
similar formulas for a tuple of `2^{g+1}` complex numbers.

.. function:: void acb_theta_agm_hadamard(acb_ptr r, acb_srcptr a, slong g, slong prec)

    Sets *r* to the image of *a* under multiplication by *H*, the `2^g\times
    2^g` Hadamard matrix (see ...). Requires that `g\geq 0` and *r* and *a* are
    initialized with at least `2^g` elements.

.. function:: void acb_theta_agm_sqrt_lowprec(acb_t r, const acb_t a, const acb_t root, slong prec)

    Sets *r* to a square root of *a*. Unlike :func:`acb_sqrt`, no special
    precision losses happen when *a* touches the negative real axis. The sign
    of the output is determined: it must overlap *root*, which is a
    (low-precision) complex ball containing either `\sqrt{a}` or `-\sqrt{a}`.
    Returns indeterminate if the correct sign cannot be determined.

.. function:: void acb_theta_agm_step_sqrt(acb_ptr r, acb_srcptr a, slong g,
              slong prec)

.. function:: void acb_theta_agm_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr
              roots, slong g, slong prec)

.. function:: void acb_theta_agm_step_good(acb_ptr r, acb_srcptr a, slong g,
              slong prec)

    Sets *r* to the result of an AGM step starting from *a*. In the
    :func:`sqrt` version, *a* is the vector of square roots. In the :func:`bad`
    version, a low-precision approximation of the roots is given. In the
    :func:`good` version, we assume that all entries of *a* have positive real
    parts, and a good choice of square roots is made. We require that `g\geq 0`
    and all vectors are initialized with at least `2^g` elements.

.. function:: void acb_theta_agm_ext_step_sqrt(acb_ptr r, acb_srcptr a, slong
              g, slong prec)

.. function:: void acb_theta_agm_ext_step_bad(acb_ptr r, acb_srcptr a,
              acb_srcptr roots, slong g, slong prec)

.. function:: void acb_theta_agm_ext_step_good(acb_ptr r, acb_srcptr a, slong
              g, slong prec)
    
    Analogous functions for extended Borchardt sequences. All vectors must be
    initialized with at least `2^{g+1}` elements.
    
.. function:: void acb_theta_agm_step_last(acb_t r, acb_srcptr a, slong g, slong prec)

    Sets *r* to the average of the first `2^g` entries of *a*.

.. function:: void acb_theta_agm_ext_step_last(acb_t r, const acb_t s, acb_srcptr a, slong g, slong prec)

    Computes an extended Borchardt mean *r* given the last term of the
    associated AGM sequence and the associated (regular) Borchardt mean *s*.
    
.. function:: void acb_theta_agm_max_abs(arb_t max, acb_srcptr a, slong nb, slong prec)
              
.. function:: void acb_theta_agm_min_abs(arb_t min, acb_srcptr a, slong nb, slong prec)

    Sets *max* (resp. *min*) to the maximum (resp. minimum) absolute value of
    the first *nb* entries of *a*.    
              
.. function:: void acb_theta_agm_abs_dist(arb_t eps, acb_srcptr a, slong nb, slong lowprec, slong prec)
    
    Computes `\varepsilon = \max_{0< i< nb} |a_i - a_0|`. Differences are
    computed at precision *prec* and absolute values at precision *lowprec*.

.. function:: void acb_theta_agm_rel_dist(arb_t eps, acb_srcptr a, slong nb, slong lowprec, slong prec)

    Computes `1|a_0|` times the output of :func:`acb_theta_agm_abs_dist`.

.. function:: void acb_theta_agm_conv_rate(arf_t c, arf_t r, const arf_t eps, slong prec)

    Computes the convergence rate of an AGM sequence consisting of good steps
    only, i.e. *c* and *r<1* such that the `i\text{th}` term of the sequence
    satisfies `|a_0 - m|\leq c r^i |a_0|` for all `i\geq 0`. The input *eps* is
    an upper bound on the relative distance for the term `i=0`, as computed by
    :func:`acb_theta_agm_rel_dist`, and must be less than *1/4*. Otherwise
    *c,r* are set to infinite values.

.. function:: slong acb_theta_agm_nb_good_steps(const arf_t c, const arf_t r, slong prec)

    Given the convergence rate *c,r* of an AGM sequence with good steps as,
    above, returns the (nonnegative) number of steps to compute before the
    equality `|a_0-m|\leq 2^{-\mathrm{prec}}|a_0|` holds. Returns negative if
    this number is infinite or cannot be computed from the given
    *c,r*. Computations are performed at a low precision specified by...

.. function:: void acb_theta_agm(acb_t r, acb_srcptr a, acb_srcptr roots, slong nb_bad, slong g, slong prec)

    Computes the limit of an AGM sequence starting from `2^g` complex
    numbers. The input data is as follows: *a* is the first term; *nb_bad* is
    the number of (possibly) bad steps; and *roots* consists of low-precision
    approximations of the correct roots for the first *nb_bad* steps, as in
    :func:`acb_theta_agm_sqrt_lowprec`. Returns an indeterminate result if a
    suitable convergence rate cannot be determined after *nb_bad* steps.

.. function:: void acb_theta_agm_ext_conv_rate(arf_t c1, arf_t c2, arf_t r, const arf_t eps, const arf_t m, const arf_t M, slong prec)

    Computes the convergence rate of an extended AGM sequence consisting of
    good steps only, i.e. *c1, c2* and *r<1* such that *c1,r* is the
    convergence rate of the regular AGM and for all `n\geq 1`, the inequality
    `|q_{n+1}-1|\leq c_2 r^{2^{n-1}}` holds. The input is as follows: *eps* is
    an upper bound on the relative distance for the term `i=0`, as computed by
    :func:`acb_theta_agm_rel_dist`, and must be less than *1/4*; and *m*
    (resp. *M*) is a lower (resp. upper) bound on the modulus of all entries of
    the initial extended AGM vector, which must be finite, with *m>0*. If these
    conditions are not satisfied then *c1, c2, r* are set to infinite values.

.. function:: void acb_theta_agm_ext_rel_err(arf_t err, const arf_t c2, const arf_t r, slong nb_good, slong prec)
    
    Computes the relative error for an extended AGM computation with
    convergence rate given by *c1, c2, r* and *nb_good* steps: the extended AGM
    is equal to `(u_n^{(0)}/\mu\cdot (1+\delta))^{2^n}` where *n* is given by
    *nb_good* and *err* is an upper bound on `|\delta|`. Requires that
    *nb_good* is at least `1` and `r < 1/2`, otherwise sets *err* to an
    infinite value.
              
.. function:: void acb_theta_agm_ext(acb_t r, acb_t s, acb_srcptr a, acb_srcptr roots, slong nb_bad, slong g, slong prec)

    Computes the extended Borchardt mean starting from `2^(g+1)` complex
    numbers. The input data is as follows: *a* is the first term; *nb_bad* is
    the number of (possibly) bad steps; and *roots* consists of low-precision
    approximations of the correct roots for the first *nb_bad* steps, as in
    :func:`acb_theta_agm_sqrt_lowprec`. Returns an indeterminate result if a
    suitable convergence rate cannot be determined after *nb_bad* steps.

.. function:: slong acb_theta_agm_nb_bad_steps(const acb_mat_t tau, slong prec)

.. function:: slong acb_theta_agm_ext_nb_bad_steps(acb_srcptr z, const acb_mat_t tau, slong prec)

    Given `\tau\in \mathcal{H}_g` and `z\in \mathbb{C}^g`, computes a
    nonnegative upper bound on the number of bad steps for the (extended) AGM
    sequence formed by theta values at `(z, 2^n\tau)` as *n* grows. A return
    value of -1 indicates that this bound cannot be computed.

.. function:: void acb_theta_agm_roots(acb_ptr roots, const acb_mat_t tau, slong nb_bad, slong prec)

.. function:: void acb_theta_agm_ext_roots(acb_ptr roots, acb_srcptr z, const acb_mat_t tau, slong nb_bad, slong prec)
    
    Given `\tau\in \mathcal{H}_g`, `z\in \mathbb{C}^g` and a number of bad
    steps *nb_bad*, computes an approximation of the required square root as
    required by :func:`acb_theta_agm` and :func:`acb_theta_agm_ext`
    respectively, using the naive algorithm for theta functions.

.. function:: void acb_theta_agm_radius(arf_t rad, const arf_struct* mi, const arf_struct* Mi, const arf_t abs_dist, slong nb, slong prec)

    Sets *rad* to the radius of a polydisk where a certain Borchardt mean
    function is surely analytic. The input data is as follows: *nb* is the
    number of (possibly) bad steps; *abs_dist* is the output of
    :func:`acb_theta_agm_abs_dist` for the vector obtained after *nb* steps;
    and *mi* (resp. *Mi*) contains a lower (resp. upper) bound for the absolute
    values of all entries in the `i\text{th}` term of the sequence for each *i*
    between *0* and *nb-1*.

Transformation formulas
-------------------------------------------------------------------------------

.. function:: slong acb_theta_char_dot(ulong a, ulong b, slong g)

    Returns *a^T b* mod *2*.

.. function:: slong acb_theta_dot(ulong a, slong* n, slong g)

    Returns *a^T n* mod *8*.

.. function:: void acb_theta_dupl_const(acb_ptr th2, acb_srcptr th, slong g, slong prec)

    Applies the duplication formula to compute `(\theta_{0,b}^2(0,2\tau))_{b\in
    \{0,1\}^g}` from `(\theta_{0,b}(0,\tau))_{b\in \{0,1\}^g}`. If the input is
    projective (i.e. given up to a common scalar factor), then so is the
    output.

    This function simply calls :func:`acb_theta_agm_step_sqrt`.

.. function:: void acb_theta_dupl_all_const(acb_ptr th2, acb_srcptr th, slong g, slong prec)

    Applies the duplication formula to compute to
    `(\theta_{a,b}^2(0,2\tau))_{a,b\in \{0,1\}^g}` from
    `(\theta_{0,b}(0,\tau))_{b\in \{0,1\}^g}`. If the input is projective, then
    so is the output.

.. function:: void acb_theta_dupl(acb_ptr th2, acb_srcptr th, slong g, slong prec)

.. function:: void acb_theta_dupl_all(acb_ptr th2, acb_srcptr th, slong g, slong prec)

    Analogues of the above to compute `(theta^2(z,2\tau), \theta^2(0,2\tau))`
    from `(theta(z,\tau),\theta(0,\tau))`. The first function simply calls
    :func:`acb_theta_agm_ext_step_sqrt`.

.. function:: void acb_theta_dupl_z(acb_ptr r, acb_srcptr th, slong g, slong prec)

    Computes `(\theta_{a,b}(2z,\tau))` from `(\theta_{a,b}(z,\tau))`.
        
.. function:: ulong acb_theta_transform_image_char(fmpz_t eps, ulong ab, const fmpz_mat_t mat)

    Computes the theta characteristic *a',b'* and an integer `\varepsilon` such
    that `\theta_{a,b}(0,N\tau) = \exp(i\pi \varepsilon/4) \theta_{a',b'}(0,\tau)`
    up to a scalar factor depending only on *N* and `\tau`. The matrix *N* must
    be symplectic. See also :func:`acb_modular_theta_transform`.

.. function:: void acb_theta_transform_proj(acb_ptr res, acb_srcptr th, const fmpz_mat_t mat, slong prec)

.. function:: void acb_theta_transform_sqr_proj(acb_ptr res, acb_srcptr th2, const fmpz_mat_t mat, slong prec)

.. function:: void acb_theta_transform_all_sqr_proj(acb_ptr res, acb_srcptr th2, const fmpz_mat_t mat, slong prec)
    
    Computes projective vectors of theta values at `(Nz,N\tau)` starting from
    the projective vector `(\theta_{a,b}(0,\tau))_{a,b\in \{0,1\}^g}`. Exactly
    what is computed depends on the suffix, as explained above.

.. function:: void acb_theta_transform_scal_const(acb_t scal, const acb_mat_t tau, const fmpz_mat_t mat, slong k2, slong prec)

.. function:: void acb_theta_transform_scal(acb_t scal_z, acb_t scal_0, acb_srcptr z, const acb_mat_t tau, const fmpz_mat_t mat, slong k2, slong prec)

    Computes the scalar factor appearing in the transformation formula for
    theta values at `(z,\tau)`. The input `k2` can be computed by
    :func:`sp2gz_k2`.

.. function:: void acb_theta_dupl_radius(arf_t rho, const arf_t r, acb_srcptr th, slong nb, slong prec)
              
.. function:: void acb_theta_transform_radius(arf_t rho, const arf_t r, acb_srcptr th, const fmpz_mat_t mat, slong prec)

.. function:: void acb_theta_dupl_transform_const_radius(arf_t rho, const arf_t r, acb_srcptr th, const fmpz_mat_t mat, slong prec)

.. function:: void acb_theta_dupl_transform_radius(arf_t rho, const arf_t r, acb_srcptr th, const fmpz_mat_t mat, slong prec)

    Computes a radius *rho* such that adding a deformation of entrywise modulus
    at most *rho* to the input vector leads to a deformation of radius at most
    *r* for the output. The operation is: either duplication, transformation,
    duplication+transformation for either theta constants or all theta values.

Ellipsoids
-------------------------------------------------------------------------------

The principle in naive algorithms to compute theta constants is to compute
partial sums of the theta series, with a strict error bound on the tail of the
series. Following..., we consider partial sums over points `n` in the lattice
`2\mathbb{Z}^g + a` contained in certain ellipsoids.

In the :func:`acb_theta_naive` functions, we first compute the relevant
ellipsoid using low-precision computations; our representation uses
`O(R^{g-1})` space for an ellipsoid of radius `R`, containing approximately
`R^g` points, gathered in one-dimensional lines.

.. type:: acb_theta_eld_struct

.. type:: acb_theta_eld_t

    Represents a *d*-dimensional sheet in an ellipsoid of ambient dimension
    *g*, i.e. a set of points of the form `n = (n_0,\ldots,n_{g-1})\in
    2\mathbb{Z}^g + a` such that `v + Yn` has `L^2` norm bounded by `R`, for
    some (upper-triangular) Cholesky matrix `Y`, some radius `R>0`, and some
    offset `v\in \mathbb{R}^g`, and finally `(n_{d},\ldots,n_{g-1})` have fixed
    values. This is a recursive type: we store
    * the interval of values for `n_{d-1}`,
    * the midpoint of that interval,
    * in the case `d\geq 2`, a number of *d-1* dimensional children of *E*,
    split between left and right children depending on the position of `n_{d-1}`
    relative to the center of the interval.

    Full ellipsoids correspond to the special case `d=g`. We always require
    `1\leq d \leq g`. Coordinates of lattice points are integers of type
    :type:`slong`.

.. function::  void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)

    Initializes *E* as a *d*-dimensional ellipsoid in ambient dimension *g*.

.. function:: void acb_theta_eld_clear(acb_theta_eld_t E)

    Clears *E* as well as any recursive data contained in it.

.. function:: void acb_theta_eld_interval(slong* min, slong* mid, slong* max, const arb_t ctr, const arf_t rad, int a, slong prec)

    Computes the minimum, middle point, and maximum of a subinterval of
    `2\mathbb{Z} + a` that is guaranteed to contain all points within a
    distance *rad* of the real number *ctr*. Both *ctr* and *rad* must be
    finite values, otherwise an error is thrown.

.. function:: void acb_theta_eld_round(slong* r, const arb_mat_t v)

    Given a `g\times 1` matrix *v*, computes a vector *r* of length *g* with
    integer entries that is close to *v*. The entries of *v* must be finite,
    otherwise an error is thrown.

.. function:: void acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t Y, const arf_t R2, arb_srcptr offset, slong* last_coords, ulong a, slong prec)

    Sets *E* to represent an ellipsoid as defined above, where *R2* indicates
    `R^2` and *offset* contains the vector `v`. The matrix *Y* must be a valid
    Cholesky matrix, i.e. an upper triangular matrix with positive diagonal
    entries, and *R2* must be finite, otherwise an error is thrown.

The following macros return meaningful values after the function
:func:`arb_eld_fill` has been called, with no computational cost.

.. macro:: acb_theta_eld_dim(E)

    Returns *d*.    

.. macro:: acb_theta_eld_ambient_dim(E)

    Returns *g*.

.. macro:: acb_theta_eld_coord(E, k)

    For `d <= k < g`, returns the common coordinate `n_k` of all lattice
    points in the ellipsoid sheet *E*.

.. macro:: acb_theta_eld_min(E)
.. macro:: acb_theta_eld_mid(E)
.. macro:: acb_theta_eld_max(E)
    
    Returns the minimum, midpoint, and maximum of `n_{d-1}` in the ellipsoid
    sheet `E`.

.. macro:: acb_theta_eld_nr(E) ((E)->nr)
.. macro:: acb_theta_eld_nl(E) ((E)->nl)

    Returns the number of right and left children of *E*, respectively.

.. macro:: acb_theta_eld_rchild(E, k)
.. macro:: acb_theta_eld_lchild(E, k)
    
    Macro giving a pointer to the `k^{\text{th}}` right (resp. left) child of
    *E*.

.. macro:: acb_theta_eld_nb_pts(E) ((E)->nb_pts)

    Returns the number of lattice points contained in *E*.

.. macro:: acb_theta_eld_box(E, k)

    Returns an integer `M_k` such that all lattice points `n` inside the
    ellipsoid sheet *E* satisfy `|n_k|\leq M_k`.

Finally, the following functions are available for convenience.

.. function:: void acb_theta_eld_points(slong* pts, const acb_theta_eld_t E)

    Sets *pts* to the list of lattice points contained in *E* (as a
    concatenation of vectors of length *g*).

.. function:: int acb_theta_eld_contains(const acb_theta_eld_t E, slong* pt)

    Returns nonzero iff *pt* is contained in the ellipsoid sheet *E*.

.. function:: void acb_theta_eld_print(const acb_theta_eld_t E)

    Prints a compact representation of *E* to :type:`stdout`.

Naive algorithms
-------------------------------------------------------------------------------

After computing a suitable ellipsoid, we can evaluate partial sums of the
series defining theta functions at high precisions. Some precomputation occurs
for each line in the ellipsoid, so that, on average as `R\to\infty`, the code
uses only two multiplications per exponential term. Further, many of these
multiplications are performed only at a fraction of the full precision,
resulting in considerable speedups. Note that using short addition sequences as
in :func:`acb_modular_addseq_theta` does not seem to further accelerate the
computations in genus `g\geq 2`.

The different :func:`theta_naive` functions only differ by their way of
handling individual lattice points. Using function pointers thus allows us to
factor out significant amounts of code.

.. function:: void acb_theta_naive_tail(arf_t bound, const arf_t R2, const
              arb_mat_t Y, slong ord, slong prec)

    Computes an upper bound for the following sum, where `p` stands for *ord*:

    .. math::

        \sum_{n\in Y\Z^g + v, \lVert n\rVert^2 \geq R^2} \lVert n\rVert^{2p} e^{-\lVert n\rVert^2)}

    using the following upper bound, valid after replacing `R^2` by
    `{\operatorname{max}\{R^2, 4, 2p\}}`

    .. math::

        2^{2g+2} R^{g-1+2p} e^{-R^2} \prod_{i=1}^g (1 + \gamma_i^{-1})

    where the `gamma_i` are the entries on the diagonal of `Y`.

.. function:: void acb_theta_naive_radius(arf_t R2, const arb_mat_t Y, slong ord,
              const arf_t eps, slong prec)

    Returns `R^2` such that the above upper bound is at most `\varepsilon`.

.. function:: void acb_theta_naive_ellipsoid(acb_theta_eld_t E, arf_struct*
              eps, acb_ptr c, acb_ptr new_z, ulong ab, int all, slong ord,
              acb_srcptr z, slong nb_z, const acb_mat_t tau, slong prec)

    Sets the ellipsoid *E* and `\varepsilon` *c*, *new_z*, `\varepsilon` such
    that summing exponential terms involving *new_z* over points of *E* and
    multiplying by *c* will yield an approximation of theta values at *z* up to
    an error at most `\varepsilon`, resulting in theta values at relative
    precision roughly *prec*.

    A value *nb_z > 1* indicates that several vectors *z* are provided. In this
    case, a unique ellipsoid is chosen for all of them, but *new_z*, *c* and
    *epsilon* will vary (hence vectors as return values).

    If *all=0*, the ellipsoid consists of lattice points in `2\mathbb{Z}^g+a`
    only, where *a* is specified by the theta characteristic *ab*. If *all* is
    nonzero, the ellipsoid consists of lattice points in `2\mathbb{Z}^g` and
    the radius is doubled, making *E* suitable for evaluating
    `\theta_{a,b}(z,\tau)` for all *a*.

.. function:: slong acb_theta_naive_newprec(slong prec, slong coord, slong
              dist, slong max_dist, slong ord)

    Returns a good choice of precision to process the next ellipsoid
    sheet. Here *coord* should be `n_{d-1}`, *dist* should be the distance to the
    midpoint of the interval, *max_dist* the half-length of the interval, and
    *ord* is the order of derivation.

.. function:: slong acb_theta_naive_fullprec(const acb_theta_eld_t E, slong
              prec)

    Returns a good choice of full precision for the summation phase.

.. type:: acb_theta_precomp_struct

.. type:: acb_theta_precomp_t

    Data structure containing precomputed data in the context of naive
    algorithms.

.. function:: void acb_theta_precomp_init(acb_theta_precomp_t D, slong nb_z,
              slong g)

    Initializes *D* to contain precomputations about *nb_z* vectors `z\in
    \mathbb{C}^g`.

.. function:: void acb_theta_precomp_clear(acb_theta_precomp_t D)

    Clears *D*.

.. function:: void acb_theta_precomp_set(acb_theta_precomp_t D, acb_srcptr z,
              const acb_mat_t tau, const acb_theta_eld_t E, slong prec)

    Precomputes the necessary data to evaluate theta functions at `(z,tau)` for
    all the vectors *z* in the provided list, using naive algorithms with
    lattice points contained in the ellipsoid *E*.

After :func:`acb_theta_precomp_set` has been called, the following macros are
available.

.. macro:: acb_theta_precomp_exp_mat(D)

    Macro giving a pointer to the matrix whose entry `(j,k)` contains
    `\exp(i\pi/4 \tau_{j,j})` if `j=k`, and `\exp(i\pi/2 \tau_{j,k})`
    otherwise.

.. macro:: acb_theta_precomp_sqr_pow(D, k, j)

    Macro giving a pointer to the complex number `\exp(i\pi/4 (2j + t)^2
    \tau_{k,k})`, where `t=1` if the lattice points in *E* has odd coordinates
    `n_k`, and `t=0` if these coordinates are even.

.. macro:: acb_theta_precomp_nb_z(D)

    Macro giving the number of vectors *z* stored in *D*.

.. macro:: acb_theta_precomp_exp_z(D, k, j)

    Macro giving a pointer to the complex number `exp(\pi i z_j)`, where *z* is
    the `k^\text{th}` vector stored in *D*.

.. type:: acb_theta_naive_worker_t

    Represents a function pointer to the "dimension 0" worker in different
    kinds of naive algorithm. A function :func:`worker_dim0` of this type has
    the following signature:

    .. function:: void worker_dim0(acb_ptr th, const acb_t term, slong* coords,
                  slong g, ulong ab, slong ord, slong prec, slong fullprec)

    where
    * *th* denotes the output vector of theta values,
    * *term* denotes the exponential term that has been computed for the
      current lattice point,
    * *coods* denotes the coordinates of that lattice point,
    * *g* is the genus,
    * *ab* is the theta characteristic, if applicable,
    * *ord* is the order of derivation, if applicable,
    * *prec* is the (relative) precision at which *term* was computed,
    * *fullprec* is the desired full precision in the summation phase.

.. function:: acb_theta_naive_worker(acb_ptr th, slong nb, const acb_t c, const
              arf_t eps, const acb_theta_eld_t E, const acb_theta_precomp_t D,
              slong k, ulong ab, slong ord, slong prec,
              acb_theta_naive_worker_t worker_dim0)

    Run the naive algorithm on the ellipsoid *E* to evaluate `\theta(z,\tau)`
    using precomputed data stored in *D*, where *z* is the `k^\text{th}` vector
    in the data structure.

.. function:: void acb_theta_naive(acb_ptr th, acb_srcptr z, slong nb_z, const
              acb_mat_t tau, slong prec)

.. function:: void acb_theta_naive_const(acb_ptr th, const acb_mat_t tau, slong
              prec)

.. function:: void acb_theta_naive_const_proj(acb_ptr th, const acb_mat_t tau,
              slong prec)

.. function:: void acb_theta_naive_all(acb_ptr th, acb_srcptr z, slong nb_z,
              const acb_mat_t tau, slong prec)

.. function:: void acb_theta_naive_all_const(acb_ptr th, const acb_mat_t tau,
              slong prec)

.. function:: void acb_theta_naive_ind(acb_t th, ulong ab, acb_srcptr z, const
              acb_mat_t tau, slong prec)

.. function:: void acb_theta_naive_ind_const(acb_t th, ulong ab, const
              acb_mat_t tau, slong prec)

    Evaluates theta functions using the naive algorithm. See above for the
    meaning of different suffixes.

Conversions
-------------------------------------------------------------------------------

.. function:: void acb_theta_renormalize_const_sqr(acb_t scal, acb_srcptr th2,
              const acb_mat_t tau, slong prec)

    Renormalizes the projective vector of squared theta constants at `tau`,
    computing *scal* such that multiplication by *scal* yields the actual theta
    values.

.. function:: void acb_theta_renormalize_sqr(acb_t scal_z, acb_t scal_0,
              acb_srcptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)

    Renormalizes the projective vectors `(\theta_{0,b}^2(z,\tau))` and
    `(\theta_{0,b}^2(0,\tau))` (concatenated in *th2*), computing the
    multiplicative factors *scal_z* and *scal_0* necessary to reach the actual
    theta values.


Newton/AGM algorithms
-------------------------------------------------------------------------------

We implement certified Newton iterations for the computation of theta functions
as detailed in...

The code first attempts to collect the necessary data to perform Newton
iterations in a dedicated data structure. If such data cannot be collected (due
to insufficient precision, or singular points in the algorithm), we fall back
to naive methods.

In the specific case of genus *1* theta functions and genus *2* theta
constants, Newton's method results in a uniform, quasi-linear time algorithm
for all inputs in the Siegel fundamental domain.

.. function:: void acb_theta_bound(arf_t rad, arf_t bound, acb_srcptr z, const
              acb_mat_t tau, slong prec)

    Computes *rad* and *bound* such that for any point `(z',\tau')` at a
    distance of at most *rad* from `(z,\tau)` entrywise, the absolute value
    `|\theta_{a,b}(z',\tau')|` is at most *bound*.

.. function:: void acb_theta_bound_const(arf_t rad, arf_t bound, const
              acb_mat_t tau, slong prec)

    Computes *rad* and *bound* such that for any point `\tau'` at a distance of
    at most *rad* from `\tau` entrywise, the absolute value
    `|\theta_{a,b}(0,\tau')|` is at most *bound*.

.. function:: void acb_theta_cauchy(arf_t bound_der, const arf_t rad, const
              arf_t bound, slong ord, slong dim, slong prec)

    Applies Cauchy's formula to compute *bound_der* with the following
    property: if *f* is an analytic function defined on a disk of radius *rad*
    around *x* and bounded in absolute value by *bound* on that disk, then the
    derivative of order *ord* of *f* at *x* is bounded by *bound_der* (in the
    sense of the infinity-operator norm for multilinear maps).

.. type:: acb_theta_agm_ctx_struct

.. type:: acb_theta_agm_ctx_t

    Data structure used to set up certified Newton iterations for theta
    functions. The following macros are available:

.. macro:: acb_theta_agm_ctx_g(ctx)

    Macro giving access to the genus *g*.

.. macro:: acb_theta_agm_ctx_nb(ctx)
    
    Macro giving access to the number of symplectic matrices used in the AGM
    method.

.. macro:: acb_theta_agm_ctx_matrix(ctx, k)

    Macro giving access to the `k^\text{th}` symplectic matrix stored in *ctx*.

.. macro:: acb_theta_agm_ctx_nb_bad_steps(ctx, k)    
.. macro:: acb_theta_agm_ctx_roots(ctx, k)
.. macro:: acb_theta_agm_ctx_mi(ctx, k)
.. macro:: acb_theta_agm_ctx_M0(ctx, k)
.. macro:: acb_theta_agm_ctx_minf(ctx, k)

    Macros giving access to the number of bad steps, precomputed choices of
    square roots, the vector of lower bounds `m_i` (as an :type:`arf_struct*`),
    the upper bound `M_0`, and the lower bound `m_\infty` (of type
    :type:`arf_t`) for the Borchardt sequence attached to the `k^\text{th}`
    symplectic matrix in *ctx*.

.. macro:: acb_theta_agm_ctx_rho(ctx)
.. macro:: acb_theta_agm_ctx_max(ctx)
.. macro:: acb_theta_agm_ctx_inv_der(ctx)

    Macros giving access to the quantities `rho`, `M`, `B_3` (in the notation
    of...) for the Newton scheme encoded by *ctx*.

.. function:: void acb_theta_agm_ctx_init(acb_theta_agm_ctx_t ctx, slong g, slong nb)
    
    Initializes *ctx* to contain data for *nb* symplectic matrices in genus *g*.

.. function:: void acb_theta_agm_ctx_clear(acb_theta_agm_ctx_t ctx)

    Clears *ctx*.

.. function:: void acb_theta_agm_ctx_set_all(acb_theta_agm_ctx_t ctx, const
              acb_mat_t tau, slong prec)

    Attempts to set *ctx* to a valid Newton scheme for the computation of theta
    constants at *tau*.

.. function:: int acb_theta_agm_ctx_is_valid(const acb_theta_agm_ctx_t ctx)

    Returns nonzero iff *ctx* encodes a valid Newton scheme, characterized by
    having nonzero `\rho` and finite `M, B_3`.

.. function:: void acb_theta_newton_eval(acb_ptr r, acb_srcptr th, const
              acb_theta_agm_ctx_t ctx, slong prec)

    Evaluates *F(th)*, where *F* is the analytic function encoded by the Newton
    scheme *ctx*.

.. function:: void acb_theta_newton_fd(acb_ptr r, acb_mat_t fd, acb_srcptr th,
              const arb_t eta, const acb_theta_agm_ctx_t ctx, slong prec)

    Evaluates *F(th)* as above and stores the result in *r*. Additionally stores
    the directional finite differences of *F* at *th* with radius *eta* in the
    columns of the matrix *fd*.

.. function:: void acb_theta_newton_run(acb_ptr r, const acb_mat_t tau, const
              acb_theta_agm_ctx_t ctx, slong prec)

    Run the Newton scheme encoded in *ctx* to compute theta values to a high
    precision *prec*. The context *ctx* must be valid.

.. function:: void acb_theta_newton_const_half_proj(acb_ptr th, const acb_mat_t
              tau, slong prec)

.. function:: void acb_theta_newton_all_sqr(acb_ptr th, const acb_mat_t tau,
              acb_srcptr z, slong prec)

.. function:: void acb_theta_newton_const_sqr(acb_ptr th2, const acb_mat_t tau,
              slong prec)

.. function:: void acb_theta_newton_all_const_sqr(acb_ptr th, const acb_mat_t
              tau, slong prec)

    Compute theta values using Newton iterations. Suffixes follow the same
    conventions as for naive algorithms above.
