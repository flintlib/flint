.. _acb-modular:

**acb_theta.h** -- theta functions and modular forms in genus 2 and above
===============================================================================

This module provides methods for numerical evaluation of theta functions and
modular forms in genus `g=2` and above. All methods also accept `g=1`,
duplicating functionality from :ref:`acb_modular.h <acb-modular>` (without the
specific speed-ups).

In the context of this module, *tau* or `\tau` always denotes an element of the
Siegel complex upper half-space `\mathbb{H}_g = \{\tau \in
\operatorname{Mat}_{g\times g}(\mathbb{C}) : \tau^t = \tau, \quad
\operatorname{Im}(\tau) \text{ is positive definite}\}`.

As usual, the numerical functions in this module compute strict error bounds:
if *tau* is represented by an :type:`acb_mat_t` which is not certainly positive
definite, the output will have an infinite radius.

Helper functions for real/complex scalars and matrices
-------------------------------------------------------------------------------

.. function:: void arb_randtest_pos(arb_t x, flint_rand_t state, slong prec,
              slong mag_bits)

    Generates a random positive number with radius around `2^{-\text{prec}}`
    the magnitude of the midpoint.
    
.. function:: void acb_randtest_disk(acb_t x, const acb_t ctr, const arf_t rad,
              flint_rand_t state, slong prec)

    Generates a random complex number with radius around `2^{-\text{prec}}` the
    magnitude of the midpoint, that is guaranteed to lie in a disk of radius
    *rad* centered at the midpoint of *ctr*.

.. function:: void acb_mat_get_real(arb_mat_t re, const acb_mat_t mat)

    Sets *re* to the real part of *mat*.

.. function:: void acb_mat_get_imag(arb_mat_t im, const acb_mat_t mat)

    Sets *im* to the imaginary part of *mat*.

.. function:: void acb_mat_set_arb_arb(acb_mat_t mat, const arb_mat_t re, const
              arb_mat_t im)

    Sets *mat* to the complex matrix with real and imaginary parts *re*, *im*.

.. function:: void arb_mat_randtest_cho(arb_mat_t mat, flint_rand_t state,
              slong prec, slong mag_bits)
    
    Sets the square matrix *mat* to a random upper triangular real matrix with
    positive diagonal entries, calling :func:`arb_randtest_precise` or
    :func:`arb_randtest_pos` on each relevant entry.

.. function:: void arb_mat_randtest_sym_pos(arb_mat_t mat, flint_rand_t state,
              slong prec, slong mag_bits)

    Sets *mat* to a random symmetric, positive definite real matrix with
    precise entries.

.. function:: int arb_mat_is_nonsymmetric(const arb_mat_t mat)

    Returns nonzero iff *mat* is certainly not symmetric.

.. function:: void arb_mat_pos_lambda(arb_t lambda, const arb_mat_t mat, slong
              prec)

    Given a symmetric, positive definite real matrix *mat*, sets *lambda* to a
    lower bound for the smallest eigenvalue of *mat*.

.. function:: void arb_mat_pos_radius(arf_t rad, const arb_mat_t mat, slong prec)

    Given a symmetric, positive definite real matrix *m*, computes a
    nonnegative *rad* such that any symmetric matrix obtained from *m* by
    adding an error of at most *rad* to each coefficient will still be positive
    definite.

.. function:: void arb_mat_reduce(arb_mat_t R, fmpz_mat_t U, const arb_mat_t M,
              slong prec)

    Given a symmetric, positive definite `g\times g` real matrix *M*, look for
    `U \in \operatorname{GL}_g(\mathbb{Z})` such that `R = U^T M U` is "more
    reduced" than *M*.

    If `g=2`, uses the Minkowski reduction algorithm; otherwise, relies on
    Flint's implementation of the LLL algorithm.

.. function:: void acb_mat_ninf(arb_t norm, const acb_mat_t mat, slong prec)

    Returns the infinity-operator norm of *mat*, defined as the maximum sum of
    absolute values of all entries on any line of *mat*.

Helper functions for integral matrices
-------------------------------------------------------------------------------

We implement matrices in `\operatorname{GSp}_{2g}(\mathbb{Z})` acting on the Siegel
upper half space as elements of type :type:`fmpz_mat_t`. As is usual in that
context, we allow single lowercase letters as matrix names when convenient.

.. function:: void fmpz_mat_get_a(fmpz_mat_t res, const fmpz_mat_t mat)

.. function:: void fmpz_mat_get_b(fmpz_mat_t res, const fmpz_mat_t mat)

.. function:: void fmpz_mat_get_c(fmpz_mat_t res, const fmpz_mat_t mat)

.. function:: void fmpz_mat_get_d(fmpz_mat_t res, const fmpz_mat_t mat)

    Sets *res* to the corresponding block of the `2g\times 2g` square matrix `m
    = \left(\begin{textmatrix} a&b\\c&d \end{textmatrix}\right)`.

.. function:: void fmpz_mat_set_abcd(fmpz_mat_t m, const fmpz_mat_t a, const
              fmpz_mat_t b, const fmpz_mat_t c, const fmpz_mat_t d)

    Sets the `2g\times 2g` matrix *mat* to `\left(\begin{textmatrix} a&b\\c&d
    \end{textmatrix}\right)`, where `a,b,c,d` are `g\times g` blocks.

.. function:: void fmpz_mat_J(fmpz_mat_t mat)

    Sets the `2g\times 2g` matrix *mat* to the symplectic matrix
    `\left(\begin{textmatrix} 0&I_g\\-I_g&0 \end{textmatrix}\right)`.

.. function:: int fmpz_mat_is_scalar(const fmpz_mat_t mat)

    Returns nonzero iff *m* is a square scalar matrix.

.. function:: int fmpz_mat_is_sp(const fmpz_mat_t mat)

.. function:: int fmpz_mat_is_gsp(const fmpz_mat_t mat)

    Returns nonzero iff the `2g\times 2g` matrix *m* is symplectic,
    resp. general symplectic.

.. function:: void fmpz_mat_diag_sp(fmpz_mat_t mat, const fmpz_mat_t U)

    Sets the `2g\times 2g` matrix *mat* to the symplectic matrix
    `\left(\begin{textmatrix} U&0\\0&U^{-T} \end{textmatrix}\right)`. We
    require `U\in \operatorname{GL}_g(\mathbb{Z})`.

.. function:: void fmpz_mat_trig_sp(fmpz_mat_t mat, const fmpz_mat_t S)

    Sets the `2g\times 2g` matrix *mat* to `\left(\begin{textmatrix}
    I_g&S\\0&I_g \end{textmatrix}\right)`, which is symplectic iff *S* is
    symmetric.

.. function:: void fmpz_mat_randtest_sp(fmpz_mat_t mat, flint_rand_t state,
              slong bits)

    Sets *mat* to a random symplectic matrix whose coefficients have length
    around *bits*.

.. function:: void fmpz_mat_siegel_fund(fmpz_mat_t mat, slong j)

    Sets the `2g\times 2g` matrix *mat* to the `j^{\text{th}}` matrix defining
    the boundary of the Siegel fundamental domain (in an arbitrary
    numbering). For `g=1`, we require `j=0`; for `g=2`, we require `0\leq j\leq
    18`; results in an error for `g\geq 3` where such a set of matrices is not
    explicitly known.

The Siegel upper half space
-------------------------------------------------------------------------------

We denote the Siegel upper half space by `\mathbb{H}_g`. It contains the
standard fundamental domain `\mathbb{F}_g` as a closed subset, defined
in... For `\varepsilon\geq 0`, closed neighborhoods `\mathcal{F}_g^\varepsilon`
can be defined following...

.. function:: void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong
              prec, slong mag_bits)

.. function:: void acb_siegel_randtest_fund(acb_mat_t tau, flint_rand_t state,
              slong prec)

    Sets the `g\times g` matrix *tau* to a random element of *\mathbb{H}_g*. In
    the second version, *tau* is guaranteed to belong to *\mathcal{F}_g*.

.. function:: void acb_siegel_cocycle(acb_mat_t res, const fmpz_mat_t mat,
              const acb_mat_t tau, slong prec)

    Sets *res* to `c\tau+d` where *c,d* are the lower `g\times g` blocks of
    *mat*.

.. function:: void acb_siegel_transform(acb_mat_t w, const fmpz_mat_t m, const
              acb_mat_t tau, slong prec)

    Sets *res* to `(a\tau + b)(c\tau + d)^{-1}` where *a,b,c,d* are the
    `g\times g` blocks of *mat*.

.. function:: int acb_siegel_is_real_reduced(const acb_mat_t tau, const arf_t
              eps, slong prec)

    Returns nonzero if each entry *z* of the square matrix *tau* satisfies
    `|\operatorname{Re}(z)|\leq 1/2+\varepsilon`. Returns 0 if this is false or
    cannot be determined.

.. function:: int acb_siegel_not_real_reduced(const acb_mat_t tau, slong prec)

    Returns nonzero if some entry *z* of the square matrix *tau* satisfies
    `|\operatorname{Re}(z)|> 1/2`. Returns 0 if this is false or cannot be
    determined.

.. function:: void acb_siegel_reduce_real(acb_mat_t res, fmpz_mat_t mat, const
              acb_mat_t tau, slong prec)

    Given a `g\times g` square matrix *tau*, computes a symmetric integer
    matrix *M* approximating `\operatorname{Re}(tau)`, sets *mat* to
    `\left(\begin{textmatrix} U_g&-M\\0&I_g \end{textmatrix}\right)`, and sets
    *res* to the image of *tau* under the action of *mat*, which should have a
    more reduced real part.

.. function:: void acb_siegel_reduce(acb_mat_t res, fmpz_mat_t mat, const
              acb_mat_t tau, slong prec)

    Given `\tau\in \mathbb{H}_g`, attempts to compute a symplectic matrix *mat*
    such that the image *res* of *tau* under this matrix is closer to the
    fundamental domain `\mathcal{F}_g`. We require `g\leq 2`.

    As in :func:`acb_modular_fundamental_domain_approx`, the output *mat* is
    always a valid symplectic matrix, but it us up to the user to check that
    the output *res* is close enough to the fundamental domain.

.. function:: int acb_siegel_is_reduced(const acb_mat_t tau, const arf_t eps,
              slong prec)

    Returns nonzero if the `g\times g` matrix *tau* belongs to
    `\mathcal{F}_g^\varepsilon`. We require `g\leq 2`. Returns 0 if this is
    false or cannot be determined.


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

Following..., we will also be interested in *extended Borchardt sequences*,
defined by similar formulas for a tuple of `2^{g+1}` complex numbers.

The formulas for steps in (extended) AGM sequences replicate the duplication
formulas for theta functions (see below). This remark is at the heart of
quasi-linear algorithms to evaluate theta functions; see below.

.. function:: void acb_theta_agm_hadamard(acb_ptr r, acb_srcptr s, slong g,
              slong prec)

    Sets *r* to the image of *s* under multiplication by *H*, the `2^g\times
    2^g` Hadamard matrix. We require `g\geq 0`; moreover *r* and *s* must be
    initialized with at least `2^g` elements.

.. function:: void acb_theta_agm_sqrt_lowprec(acb_t r, const acb_t x, const
              acb_t root, slong prec)

    Sets *r* to a square root of *x* to high precision that is contained in the
    (low-precision) approximation *root*. Unlike :func:`acb_sqrt`, no special
    precision losses happen when *x* touches the negative real axis.

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
    parts, and a good choice of square roots is made. We require `g\geq 0`; all
    vectors must be initialized with at least `2^g` elements.

.. function:: void acb_theta_agm_ext_step_sqrt(acb_ptr r, acb_srcptr a, slong
              g, slong prec)

.. function:: void acb_theta_agm_ext_step_bad(acb_ptr r, acb_srcptr a,
              acb_srcptr roots, slong g, slong prec)

.. function:: void acb_theta_agm_ext_step_good(acb_ptr r, acb_srcptr a, slong
              g, slong prec)
    
    Analogous functions for extended Borchardt sequences. All vectors must be
    initialized with at least `2^{g+1}` elements.

.. function:: void acb_theta_agm(acb_t r, acb_srcptr a, acb_srcptr all_roots,
              const arf_t rel_err, slong nb_bad, slong nb_good, slong g,
              slong prec)

.. function:: void acb_theta_agm_ext(acb_t r, acb_srcptr a, acb_srcptr
              all_roots, const arf_t rel_err, slong nb_bad, slong nb_good,
              slong g, slong prec)

    Evaluates the limit of an AGM sequence starting from *a*. First takes
    *nb_bad* bad steps using low-precision square roots stored in *all_roots*
    of length *nb_bad* `\times 2^g`; then, renormalizes and takes *nb_good*
    good steps.

    The first entry of the resulting vector is an approximation of the
    limit. We finally add some relative error specified by *rel_err* to account
    for the mathematical convergence error; this error must be computed by the
    user in terms of the starting data: while general formulas predict suitable
    values of *nb_bad*, *nb_good* and *rel_err* in terms of *a*, they are
    overly pessimistic for our applications.

.. function:: slong acb_theta_agm_nb_bad_steps(const acb_mat_t tau, slong prec)

    Given `\tau\in \mathcal{H}_g`, computes *n\geq 0* such that theta constants
    at `2^n\tau` lie in a disk centered at `1` with radius `1/20`. The result
    is intended for use as *nb_bad* in :func:`acb_theta_agm`.

.. function:: slong acb_theta_agm_nb_good_steps(arf_t rel_err, slong g, slong prec)

    Computes the number of good AGM steps, starting from a configuration of
    complex numbers within the disk centered at `1` with radius `1/20`, to
    approximate the limit value up to a relative error of
    `2^{-\text{prec}}`. Also sets *rel_err* to this value. The result is
    intended for use as *nb_good* and *rel_err* in :func:`acb_theta_agm`.


Conventions on theta functions
-------------------------------------------------------------------------------

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
2. Number of theta values:
    * All values `\theta_{0,b}` for `b\in \{0,1\}^g`: default (no suffix).
    * All values `\theta_{a,b}` for all *a,b*: suffix :func:`all`.
    * Individual value `\theta_{a,b}` for specified *a,b*: suffix :func:`ind`.
3. Value of *z*:
    * `z=0` (theta constants): suffix :func:`const`. The result is zero
      whenever `a^T b` is odd.
    * Specified *z*: default (no suffix). In this case, also compute theta
      constants.
    * Vector of values of *z*: not yet implemented.
4. Theta values taken at `\tau/2` instead of `tau`: suffix :func:`half`.
5. Projective theta values (i.e., the result is defined up to simultaneous
   multiplication by a nonzero complex number): suffix :func:`proj`. (This
   scalar factor may not be the same for theta constants and theta values at
   `z\neq 0`).
6. Squared theta values: suffix :func:`sqr`.
7. Also compute derivatives of theta functions up to some order: suffix
   :func:`jet`.

In the following,
* *th* is a vector of theta values, projective or not,
* *th2* is a vector of squared theta values,

Transformation formulas
-------------------------------------------------------------------------------

.. function:: void acb_theta_duplication(acb_ptr th2, acb_srcptr th, slong g,
              slong prec)

    Applies the duplication formula to compute `(\theta_{0,b}^2(0,2\tau))_{b\in
    \{0,1\}^g}` from `(\theta_{0,b}(0,\tau))_{b\in \{0,1\}^g}`. If the input is
    given up to a common scalar factor, so is the output.

    This function simply calls :func:`acb_theta_agm_step_sqrt`.

.. function:: void acb_theta_duplication_all(acb_ptr th2, acb_srcptr th, slong
              g, slong prec)

    Applies the duplication formula to compute to
    `(\theta_{a,b}^2(0,2\tau))_{a,b\in \{0,1\}^g}` from
    `(\theta_{0,b}(0,\tau))_{b\in \{0,1\}^g}`. If the input is given up to a
    common scalar factor, so is the output.
    
    
.. function:: ulong acb_theta_transform_image_char(fmpz_t epsilon, ulong ab,
              const fmpz_mat_t N)

    Computes the theta characteristic *a',b'* and an integer `\varepsilon` such
    that `\theta_{a,b}(0,N\tau) = \exp(i\pi \varepsilon/4) \theta_{a',b'}(0,\tau)`
    up to a scalar factor depending only on *N* and `\tau`. The matrix *N* must
    be symplectic. See also :func:`acb_modular_theta_transform`.

.. function:: void acb_theta_transform_sqr_proj(acb_ptr res, acb_srcptr th,
              const fmpz_mat_t N, slong prec)

    Applies the transformation formula from to compute the projective vector
    `(\theta_{0,b}^2(0,N\tau)_{b\in \{0,1\}^g}` from the projective vector
    `(\theta_{a,b}(0,\tau))_{a,b\in \{0,1\}^g}`.

Ellipsoids for naive algorithms
-------------------------------------------------------------------------------

The principle in naive algorithms to compute theta constants is to compute
partial sums of the theta series, with a strict error bound on the tail of the
series. Following..., we consider partial sums over points `n` in the lattice
`2\mathbb{Z}^g + a` contained in certain ellipsoids. The inequalities used to
bound the tail of the series are detailed here. (todo: add link).

In the :func:`acb_theta_naive` functions, we first compute the relevant
ellipsoid using low-precision computations; our representation uses
`O(R^{g-1})` space for an ellipsoid of radius `R`, containing approximately
`R^g` points, gathered in one-dimensional lines. The partial sum of exponential
terms is then computed at high precision. Some precomputation occurs for each
line so that, on average as `R\to\infty`, the code uses only two
multiplications per exponential term. Further, many of these multiplications
are performed only at a fraction of the full precision, resulting in
considerable speedups. Note that using short addition sequences as in
:func:`acb_modular_addseq_theta` does not seem to further accelerate the
computations in genus `g\geq 2`.

Many similar :func:`theta_naive` functions are provided; they essentially
differ by their way of handling individual lattice points. Using function
pointers for this last step allows us to factor out significant amounts of
code.

.. type:: acb_theta_eld_struct

.. type:: acb_theta_eld_t

    Represents a *d*-dimensional sheet in an ellipsoid of ambient dimension
    *g*, i.e. a set of points of the form `n = (n_1,\ldots,n_g)\in
    2\mathbb{Z}^g + a` such that `v + Yn` has `L^2` norm bounded by `R`, for
    some Cholesky matrix `Y`, some radius `R>0`, and some offset `v\in
    \mathbb{R}^g`, and finally `(n_{d+1},\ldots,n_g)` have fixed values. This is
    a recursive type: we store
    * the interval of values for `n_d`,
    * the midpoint of that interval,
    * in the case `d\geq 2`, a number of *d-1* dimensional children of *E*,
    split between left and right children depending on the position of `n_d`
    relative to the center of the interval.

    Full ellipsoids correspond to the special case `d=g`. We always require
    `1\leq d \leq g`. Coordinates of lattice points are encoded as
    :type:`slong` 's.

.. function::  void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)

    Initializes *E* as a *d*-dimensional ellipsoid sheet in ambient dimension
    *g*.

.. function:: void acb_theta_eld_clear(acb_theta_eld_t E)

    Clears *E* as well as any recursive data contained in it.

.. function:: void acb_theta_eld_interval(slong* min, slong* mid, slong* max,
              const arb_t ctr, const arb_t rad, int a, slong prec)

    Computes the minimum, middle point, and maximum of a subinterval of
    `2\mathbb{Z} + a` that is guaranteed to contain all points within a
    distance *rad* of the real number *ctr*. Both *ctr* and *rad* must be
    finite values.

.. function:: void acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t Y,
              const arb_t normsqr, arb_srcptr offset, slong* last_coords, ulong
              a, slong prec)

    Sets *E* to represent lattice points in an ellipsoid as defined above,
    where *normsqr* indicates `R^2` and *offset* contains the vector `v`. The
    matrix *Y* must be a valid Cholesky matrix, i.e. an upper triangular matrix
    with positive diagonal entries, and *normsqr* must be finite.

.. function:: void acb_theta_eld_points(slong* pts, const acb_theta_eld_t E)

    Sets *pts* to the list of lattice points contained in *E*.

.. function:: int acb_theta_eld_contains(const acb_theta_eld_t E, slong* pt)

    Returns nonzero iff *pt* is contained in the ellipsoid sheet *E*.

.. function:: void acb_theta_eld_print(const acb_theta_eld_t E)

    Prints a compact representation of *E* to :type:`stdout`.

In addition, the following macros are available after the function
:func:`arb_eld_fill` has been called, with no computational cost.

.. macro:: acb_theta_eld_dim(E)

    Returns *d*.    

.. macro:: acb_theta_eld_ambient_dim(E)

    Returns *g*.

.. macro:: acb_theta_eld_coord(E, k)

    For `d < k \leq g`, returns the common coordinate `n_k` of all lattice
    points in the ellipsoid sheet *E*.

.. macro:: acb_theta_eld_min(E)
.. macro:: acb_theta_eld_mid(E)
.. macro:: acb_theta_eld_max(E)
    
    Returns the minimum, midpoint, and maximum of `n_d` in the ellipsoid sheet `E`.

.. macro:: acb_theta_eld_nr(E) ((E)->nr)
.. macro:: acb_theta_eld_nl(E) ((E)->nl)

    Returns the number of right and left children of *E*, respectively

.. macro:: acb_theta_eld_rchild(E, k)
.. macro:: acb_theta_eld_lchild(E, k)
    
    Macro giving a pointer to the `k^{\text{th}}` right (resp. left) child of
    *E*.

.. macro:: acb_theta_eld_nb_pts(E) ((E)->nb_pts)

    Returns the number of lattice points contained in *E*.

.. macro:: acb_theta_eld_box(E, k)

    Returns an integer `M_k` such that all lattice points `n` inside the
    ellipsoid sheet *E* satisfy `|n_k|\leq M_k`.

Precomputations for naive algorithms
-------------------------------------------------------------------------------

.. function:: void acb_theta_naive_tail(arf_t B, const arf_t R, const arb_mat_t
              Y, slong p, slong prec)

    Computes an upper bound for the sum

    .. math::

        \sum_{n\in Y\Z^g + v, \lVert n\rVert^2 \geq R^2} \lVert n\rVert^{2p} e^{-\lVert n\rVert^2)}

    using the following upper bound, valid after replacing `R^2` by
    `{\operatorname{max}\{R^2, 4, 2p\}}`

    .. math::

        2^{2g+2} R^{g-1+2p} e^{-R^2} \prod_{i=1}^g (1 + \gamma_i^{-1})

    where the `gamma_i` are the entries on the diagonal of `Y`.

.. function:: void acb_theta_naive_ellipsoid(acb_theta_eld_t E, arf_t epsilon,
              ulong ab, int all, int unif, slong ord, acb_srcptr z, const
              acb_mat_t tau, slong prec)

    Sets the ellipsoid *E* and `\varepsilon` such that summing over points of
    *E* will yield an approximation of theta values up to an error of at most
    `\varepsilon` in naive algorithms for theta functions, resulting in theta
    values at relative precision roughly *prec*. (Todo: explain all, unif).

.. function:: void acb_theta_naive_radius(arf_t R, const arb_mat_t Y, slong p,
              const arf_t epsilon, slong prec)

    Returns `R^2` such that the above upper bound is at most `\varepsilon`.

.. function:: slong acb_theta_naive_newprec(slong prec, slong coord, slong
              dist, slong max_dist, slong ord)

    Returns a good choice of precision to process the next ellipsoid
    sheet. Here *coord* should be `n_d`, *dist* should be the distance to the
    midpoint of the interval, *max_dist* the half-length of the interval, and
    *ord* is the order of derivation.

.. function:: slong acb_theta_naive_fullprec(const acb_theta_eld_t E, slong
              prec)

    Returns a good choice of full precision for the summation phase.

.. type:: acb_theta_precomp_struct

.. type:: acb_theta_precomp_t

    Data structure containing precomputed data in the context of naive
    algorithms.

.. function:: void acb_theta_precomp_init(acb_theta_precomp_t D, slong g)

    Initializes *D*.

.. function:: void acb_theta_precomp_clear(acb_theta_precomp_t D)

    Clears *D*.

.. function:: void acb_theta_precomp_set(acb_theta_precomp_t D, const acb_mat_t
              tau, const acb_theta_eld_t E, slong prec)

    Precomputes the necessary data to evaluate theta functions at *tau* using
    naive algorithms with lattice points contained in the ellipsoid *E*.

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

Naive algorithms
-------------------------------------------------------------------------------

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

.. function:: void acb_theta_naive_worker_dim1(acb_ptr th, const
              acb_theta_eld_t E, const acb_theta_precomp_t D, const acb_t lin,
              const acb_t cofactor, ulong ab, slong ord, slong prec, slong
              fullprec, acb_theta_naive_worker_t worker_dim0)

    Process a dimension 1 ellipsoid sheet *E* (i.e. a line), calling
    :func:`worker_dim0` on each point of the line. In such a line, the
    exponential terms take the form `C \cdot L^{n_0} \cdot \exp(i\pi
    n_0^2\tau_{0,0}/4)` for some complex numbers `C,L` given by *cofactor* and
    *lin* respectively.

    This function uses only two multiplications per exponential term.

.. function:: acb_theta_naive_worker_rec(acb_ptr th, acb_mat_t lin_powers,
              const acb_theta_eld_t E, const acb_theta_precomp_t D, acb_srcptr
              exp_z, const acb_t cofactor, ulong ab, slong ord, slong prec,
              slong fullprec, acb_theta_naive_worker_t worker_dim0)

    Recursively process a general ellipsoid sheet *E* of dimension *d*, calling
    :func:`acb_theta_naive_worker_dim1` when dimension 1 is reached. Other
    arguments are as follows:
    * The matrix *lin_powers* contains `\exp(\pi i n_j \tau_{k,j}/4)` as its
    `(k,j)^{\text{th}}` entry for `k\leq d` and `j > d`; these entries are
    considered :type:`const`, but others may me modified during the recursion
    process.
    * The vector *exp_z* is assumed to contain `\exp(\pi i z_k)` as its
    `k^{\text{th}}` entry, where *z* is the argument of the theta function.
    * The scalar *cofactor* is assumed to contain the common cofactor for all
    exponential terms in *E*.

.. function:: void acb_theta_naive_term(acb_t exp, const acb_mat_t tau,
              acb_srcptr z, ulong ab, slong* coords, slong prec)

    Computes the exponential term in the series `\theta_{a,b}(z,\tau)`
    corresponting to the lattice point *coods* individually.

.. function:: void acb_theta_naive(acb_ptr th, acb_srcptr z, const acb_mat_t
              tau, slong prec);

.. function:: void acb_theta_naive_const(acb_ptr th, const acb_mat_t tau, slong
              prec);

.. function:: void acb_theta_naive_const_proj(acb_ptr th, const acb_mat_t tau,
              slong prec);

.. function:: void acb_theta_naive_all(acb_ptr th, acb_srcptr z, const
              acb_mat_t tau, slong prec);

.. function:: void acb_theta_naive_all_const(acb_ptr th, const acb_mat_t tau,
              slong prec);

.. function:: void acb_theta_naive_ind(acb_t th, ulong ab, acb_srcptr z, const
              acb_mat_t tau, slong prec);

.. function:: void acb_theta_naive_ind_const(acb_t th, ulong ab, const
              acb_mat_t tau, slong prec);

    Evaluates theta functions using the naive algorithm. See above for the
    meaning of different suffixes.

Context for Newton iterations
-------------------------------------------------------------------------------

.. function:: void acb_theta_bound(arf_t rad, arf_t bound, acb_srcptr z, const
              acb_mat_t tau, slong prec)

    Computes *rad* and *bound* such that for any point `(z',\tau')` at a
    distance of at most *rad* from `(z,\tau)` entrywise, the absolute value
    `|\theta_{a,b}(z',\tau')|` is at most *bound*. Todo: document the
    inequalities used.

.. function:: void acb_theta_bound_const(arf_t rad, arf_t bound, const
              acb_mat_t tau, slong prec)

    Computes *rad* and *bound* such that for any point `\tau'` at a distance of
    at most *rad* from `\tau` entrywise, the absolute value
    `|\theta_{a,b}(0,\tau')|` is at most *bound*. Todo: document the
    inequalities used.

.. function:: void acb_theta_cauchy(arf_t bound_der, const arf_t rad, const
              arf_t bound, slong ord, slong dim, slong prec);

    Applies Cauchy's formula to compute *bound_der* with the following
    property: if *f* is an analytic function defined on a disk of radius *rad*
    around *x* and bounded in absolute value by *bound* on that disk, then the
    derivative of order *ord* of *f* at *x* is bounded by *bound_der* (in the
    sense of the infinity-operator norm for multilinear maps).

Issues related to compatibility with existing code
-------------------------------------------------------------------------------

* Cholesky is upper rather than lower triangular.
* Some functions are not used? nonsymmetric,
* in agm.c: missing renormalization, since agm_step_good assumes positive real part.
* need a nb_bad_steps for extended B means.

/* General comments:
   - In case a computation fails, output values are set to NaNs if possible, otherwise abort
   - A suffix sqr indicates squares of theta values
   - A suffix proj indicates theta values up to common scaling, and derivatives of those
   - A suffix half indicates theta values taken at tau/2
   - A suffix all indicates theta values for all characteristics (a,b), not only a=0
   - A suffix ind indicates a single theta value
   - A suffix const indicates theta constants (z=0). If not present, we compute both theta constants and regular theta values; "proj" is understood for each half independently.
   - A suffix jet indicates successive derivatives with respect to z. Return a vector of matrices as follows: one matrix per derivation order; in each of these, a row of the matrix contains partial derivatives of a fixed theta value
   - Order: naive/newton, all/ind, const, half,  proj, sqr, jet
   - Characteristics (a,b) are encoded as ulongs; first half is a, second half is b
*/
