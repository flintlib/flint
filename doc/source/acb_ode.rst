.. _acb-ode:

**acb_ode.h** -- linear differential equations with polynomial coefficients
==================================================================================

.. note::

    This module is experimental; the API is subject to change.

This module provides functions for solving linear homogeneous differential
equations

.. math::

   a_r(x) y^{(r)}(x) + \dots + a_1(x) y'(x) + a_0(x) y(x) = 0

with polynomial coefficients `a_i(x) \in \mathbb C[x]`.
When working locally at 0, is is often convenient to rewrite the equation as

.. math::

   \tilde a_r(x) (\theta^r y)(x) + \dots + \tilde a_1(x) (\theta y)(x) + \tilde a_0(x) y(x) = 0

where `\theta = x \frac{\mathrm d}{\mathrm d x}` is the *Euler derivation*.
The integer `r` is called the *order* of the equation.

Some high-level functions take as input differential operators represented as
``gr_ore_poly_t`` objects and may accept both of the above forms, as well as
various coefficient types.
Lower-level functions typically expect pointers to ``acb_poly_struct`` objects
instead. Unless stated otherwise, these correspond to coefficient vectors with
respect to the Euler derivation.

This module currently supports evaluating solutions given by initial conditions
at an *ordinary* or a *regular singular* point.
Initial conditions at a point are coordinates in a given *local basis* of the
solution space.
At an ordinary point, the equation admits an `r`-dimensional vector space of
analytic solutions, all of valuation less than `r`.
At a regular singular point, it admits an `r`-dimensional space of solutions
of the form

.. math::

   \sum_\lambda x^\lambda \sum_{k=0}^{K-1} y_{\lambda,k}(x) \frac{\log(x)^k}{k!}

where the `y_{\lambda,k}` are analytic and `\lambda` ranges over the local
*exponents* of the differential equation, which are the roots of the local
*indicial polynomial*.
A maximal set of exponents that differ by integers is called a *group*
(this has nothing to do with the algebraic structure of the same name).
We call the exponent of minimal real part in a group the *leader* of the group,
and the difference of an exponent with its group leader its *shift*.
Solutions that can be written with a single common `\lambda` are also said to
belong to the same group.

There is currently no support for inhomogeneous equations.

Fundamental matrices
-------------------------------------------------------------------------------

See also the :ref:`usage example <ode-example-basic>` below.

.. type:: acb_ode_basis_t

    Constants used to identify families of local bases of solutions.
    To specify a basis, one needs a constant of this type, an expansion
    point, and possibly some additional information.

    Currently the following two families are supported.
    In both cases, the basis is the concatenation of blocks corresponding to
    groups of local exponents.
    Within each block, the columns are ordered by increasing shift relative to
    the group leader (that is, by increasing real part of the exponent), and for
    a given shift by decreasing degree in `\log(x)` of the coefficient of `x` of
    minimal valuation.
    The order of the blocks is specified separately.

    Let `E` denote the set of pairs `(\alpha, k)` where `\alpha` is a local
    exponent of multiplicity `m` and `k` ranges from `0` to `m - 1`.
    At a regular singular point, the cardinality of `E` is equal to the order of
    the differential equation.

    .. macro:: ACB_ODE_BASIS_ECHELON

        The basis is made of a solution of leading term
        `x^{\alpha} \log(x)^k/k!` for each `(\alpha, k) \in E`.
        The series expansion of the basis element attached to
        `(\alpha, k)` involves no term in `x^{\beta} \log(x)^\ell`
        with `(\beta, \ell) \in E \setminus \{ (\alpha, k) \}`.

        This is the default when no basis choice can be specified.

    .. macro:: ACB_ODE_BASIS_CASCADE

        The basis contains a solution of leading term `x^{\alpha} \log(x)^m/m!`
        for each exponent `\alpha`, with `m` the multiplicity of `\alpha`.
        The series expansions of these basis elements involve no term in
        `x^{\beta} \log(x)^\ell`
        with `(\beta, \ell) \in E \setminus \{ (\alpha, m) \}`.
        The remaining elements of the basis are obtained by repeatedly changing
        each `\log(x)^k/k!` to `\log(x)^{k-1}/(k-1)!` and dropping the terms not
        involving `\log(x)` in previously constructed elements.

    More constants of this type may be added in the future, including some not
    sharing the block structure of the previous two.

.. function:: int acb_ode_fundamental_matrix(acb_mat_t mat, const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop, const acb_ode_exponents_t expos, acb_srcptr lcrt, const acb_t pt, acb_ode_basis_t basis, slong prec)

    Sets *mat* to the value at *pt* of the fundamental matrix in the
    basis *basis* of the differential operator *dop*.

    This function is currently limited to operators *dop* of algebra type
    :macro:`GR_ORE_POLY_EULER_DERIVATIVE`.

    The origin may be an ordinary point or a regular singular point of *dop*.
    Operators whose coefficients with respect to the Euler derivative have a
    common factor with a root at the origin are currently not supported even if
    they fall in the previous class.

    The evaluation point *pt* must lie within an open disk centered at the
    origin and containing no singular point of *dop* except, possibly, for the
    origin itself.

    The *expos* parameter corresponds to the exponent structure of *dop* at the
    origin in the format returned by :func:`acb_ode_exponents`.
    In general, it cannot be automatically determined from *dop* and must be
    provided by the caller (see ``examples/ode_exponents``).
    When the exponent structure can be computed by :func:`acb_ode_exponents`
    (typically, over an exact field of constants with a complete zero test),
    though, *expos* may be set to ``NULL``.
    Note however that, in this case, the order of the exponent groups is
    unspecified.

    Similarly, the user must provide in *lcrt* enclosures of the roots of the
    leading coefficient of *dop* (with multiple roots repeated), except when the
    ground field supports computing complex ball roots using
    :func:`gr_poly_roots_other`, in which case *lcrt* may be set to ``NULL``.

    The columns of the output matrix correspond to the basis of type *basis* of
    the solution space, with solution groups ordered as in *expos*.
    (Note that when *expos* is omitted on input, the order is unspecified.)
    The row dimension of the output matrix determines the number of derivatives
    of the solutions to be evaluated.
    (The current implementation is not designed to evaluate a number of
    derivatives exceeding the order of *dop* efficiently.)

    The return value is *gr*-style status code. A return value other that
    ``GR_SUCCESS`` indicates a failure in the algebraic part of the computation.
    It may happen, however, that the return value is ``GR_SUCCESS`` but the
    output matrices are indeterminate.

.. function:: void _acb_ode_fundamental_matrix_vec(acb_mat_struct * mat, const acb_poly_struct * dop, slong dop_len, const acb_ode_exponents_t expos, acb_srcptr lcrt, acb_ode_bound_t bound, acb_srcptr pts, slong npts, acb_ode_basis_t basis, acb_ode_sum_worker_t sum_worker, void * worker_data, slong prec)
              int acb_ode_fundamental_matrix_vec(acb_mat_struct * mat, const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop, const acb_ode_exponents_t expos, acb_srcptr lcrt, acb_srcptr pts, slong npts, acb_ode_basis_t basis, slong prec)

     Like :func:`acb_ode_fundamental_matrix` but evaluates the fundamental
     matrix at multiple points.

     The underscore method assumes that *dop* has order at least one and does
     not allow *expos* and *lcrt* parameters to be omitted.
     Additionally, the caller must provide an object *bound* pre-filled using
     :func:`acb_ode_bound_precompute` and a summation callback *sum_worker*
     with an arbitrary user data pointer *worker_data*.
     The user data is accessible to *sum_worker* in the *data* field of its
     *sum* argument.

Solution groups and exponents
-------------------------------------------------------------------------------

See ``examples/ode_exponents`` example of how to use these functions in
combination with :func:`acb_ode_fundamental_matrix`.

.. type:: acb_ode_shift_struct
          acb_ode_shift_t
          acb_ode_group_struct
          acb_ode_group_t
          acb_ode_exponents_struct
          acb_ode_exponents_t

    Types used to represent the structure of a local basis of solutions, that
    is, the roots of the indicial polynomial and the integer differences between
    these roots.

    An *acb_ode_exponents_struct* contains an array *groups* of length *ngroups*
    of *acb_ode_group_struct* objects.
    Each *acb_ode_group_struct* contains an exponent *leader* of type
    :type:`acb_t` and an array of length *nshifts* of *acb_ode_shift_struct*
    objects.
    An *acb_ode_shift_struct* is a pair of an integer shift *n* and a
    multiplicity *mult*.
    The shifts are sorted in increasing order.

.. function:: void acb_ode_exponents_init(acb_ode_exponents_t expos)
              void acb_ode_exponents_clear(acb_ode_exponents_t expos)

.. function:: void acb_ode_exponents_ordinary(acb_ode_exponents_t expos, slong order)

    Sets *expos* to the exponent structure of an operator of order *order* at
    an ordinary point.

.. function:: void acb_ode_exponents_randtest(acb_ode_exponents_t expos, flint_rand_t state, slong len, slong disp, slong prec, slong mag_bits)

    Sets *expos* to a random exponent structure with *len* exponents in total
    (counting multiplicities) and a maximum dispersion of *disp*.

.. function:: void acb_ode_exponents_println(const acb_ode_exponents_t expos)

.. function:: int acb_ode_exponents(acb_ode_exponents_t expos, const gr_ore_poly_t dop, gr_ctx_t Dop, gr_ctx_t CC)

    Sets *expos* to the exponent structure of the operator *dop* at the origin,
    computed as elements of *CC*. The only currently supported *CC* are complex
    ball fields.

.. function:: slong acb_ode_exponents_length(const acb_ode_exponents_t expos)

    Returns the total number of solutions across all groups in *expos*, counted
    with multiplicities.

.. function:: void acb_ode_group_init(acb_ode_group_t group, slong len)

    Initializes *group* with space for up to *len* shifts.

.. function:: void acb_ode_group_clear(acb_ode_group_t group)

    Clears *group*.

.. function:: void acb_ode_group_set(acb_ode_group_t dest, const acb_ode_group_t src)

    Sets *dest* (which must have space for *src->nshifts* shifts) to a copy of
    *src*.

.. function:: slong acb_ode_group_length(const acb_ode_group_t group)

    Returns the number of solutions of *group*, i.e., the sum of
    multiplicities of its shifts.

.. function:: slong acb_ode_group_multiplicity (const acb_ode_group_t group, slong n)

    Returns the multiplicity of *n* as a shift in group *group*.

.. function:: slong acb_ode_group_nlogs (const acb_ode_group_t group, slong n)

    Returns a bound on the length of the term of index *n* in any solution of
    the group *group* when that term is viewed as a polynomial in `\log(x)`.

.. function:: void acb_ode_indicial_polynomial_from_exponents(acb_poly_t ind, const acb_ode_exponents_t expos, slong prec)

    Sets *ind* to the indicial polynomial of an equation with exponents *expos*.

Partial sums of logarithmic series solutions
-------------------------------------------------------------------------------

This section describes low-level functions for computing series expansions
and/or values of partial sums or derivatives of partial sums of solutions.
These functions can handle several solutions from the same group (or possibly,
in some cases, from conjugate groups) in parallel, and can evaluate them at
zero, one, or more points.

The general usage pattern is as follows:

1. initialize a context object of type *acb_ode_sum_t*;
2. configure it by setting a differential operator
   (*acb_ode_sum_set_diffop*), solution group (*acb_ode_sum_set_group* or
   similar) and some vectors of initial values (*acb_ode_sum_set_ini_\**), in
   this order;
3. optionally set some evaluation points and other options;
4. call one or more summation functions such as *acb_ode_sum_divconquer* to
   compute partial sums of solutions;
5. extract the value and/or coefficients of the solution (XXX this step
   currently requires using *sol_* functions).

.. type:: acb_ode_sum_struct
          acb_ode_sum_t

    Context object.

    The following options can be set by writing to the corresponding fields
    after initializing the object:

    .. member:: ulong flags

        .. macro:: ACB_ODE_WANT_SERIES

            Keep all series coefficients.

        .. macro:: ACB_ODE_APPROX

            Skip the computation of rigorous error bounds.

.. function:: void acb_ode_sum_init(acb_ode_sum_t sum, slong dop_len, slong npts, slong nsols, slong nder)

    Initializes a context object suitable for computing the derivatives of order
    0 to *nder - 1* of *nsols* solutions of an equation *dop_len - 1* at
    *npts* points.
    When only the coefficients of the series are of interest, *npts* can be set
    to zero.

    XXX Clarify the case *nder = 0*. It makes sense to have *npts = 0* with
    *nder = 1* if one wants to compute the coefficients of a rigorous
    approximation valid in a certain disk. For computing a fixed number of
    coefficients without any kind of tail bound, though, we may want to allow
    *nder = 0*.

.. function:: void acb_ode_sum_clear(acb_ode_sum_t sum)

    Clears the context object.

.. function:: void acb_ode_sum_set_diffop(acb_ode_sum_t sum, const acb_poly_struct * dop, slong dop_len, const mag_t cvrad)

    Sets the differential operator whose solutions are to be computed.
    The operator must be nonzero.
    The *cvrad* parameter is a lower bound on the radius of convergence of the
    solutions and will limit the magnitude of the evaluation points.

.. function:: void acb_ode_sum_set_group(acb_ode_sum_t sum, const acb_ode_group_t grp)
              void acb_ode_sum_set_ordinary(acb_ode_sum_t  sum)

    Sets the solution group.

.. function:: void acb_ode_sum_set_ini_echelon(acb_ode_sum_t sum)
              void acb_ode_sum_set_ini_highest(acb_ode_sum_t sum)

    Sets initial values defining the solutions to be computed.
    The solution group must have been set before.

.. function:: void acb_ode_sum_set_points(acb_ode_sum_t sum, acb_srcptr pts, slong npts)

    Sets the evaluation points.

.. function:: slong acb_ode_sum_max_nlogs(acb_ode_sum_struct * sum)

.. function:: slong acb_ode_sum_max_nlogs_xn(acb_ode_sum_struct * sum, slong off)

.. function:: void acb_ode_sum_precompute(acb_ode_sum_t sum)

    Performs various precomputations.

.. type:: acb_ode_sum_worker_t

    Pointer to a worker function
    ``void worker(acb_ode_sum_t sum, slong nterms, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)``.

.. function:: void acb_ode_sum_divconquer(acb_ode_sum_t sum, slong nterms, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)

    Sums the series expansions of the solutions contained in *sum*, updating the
    *sum* object with their coefficients and/or partial sums.
    The *nterms* parameter can be used to limit the number of terms (-1 means
    unlimited).
    The called must supply two objects *bound* and *gbound* precomputed using
    :func:`acb_ode_bound_precompute` and :func:`acb_ode_group_bound_precompute`.

    This function computes the sum from scratch even if *sum* already contains a
    partial sum from a previous run.

    XXX Clarify what happens with no evaluation points and how to evaluate the
    *coefficients* to a given precision.

All *_acb_ode_sum_forward_\** functions assume that
:func:`acb_ode_sum_precompute` has been applied and that the *sum* object's
component *acb_ode_sol_t* all have space for the new coefficients.
Additionally, these functions may leave some of these components in a state
that requires :func:`_acb_ode_sum_fix` to be called before using other
functions on them.

.. function:: void _acb_ode_sum_forward_1(acb_ode_sum_struct * sum)

    Adds one term to *sum*, whose components are assumed to have space for the
    update.

.. function:: void _acb_ode_sum_forward_divconquer(acb_ode_sum_struct * sum, slong high, slong block_len, slong stride, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)

    Extends the sum until the partial sum of order *high* is reached, some error
    tolerance based on *prec* is met, or the working precision and/or tail bound
    quality are determined to be insufficient (whichever comes first).

    The algorithm computes strides of *stride* terms (checking for convergence
    and making some internal adjustments after each stride) by chaining blocks
    of *block_len* terms that are computed using a divide-and-conquer method.
    It is assumed that *block_len* is at least the degree of the coefficients of
    the differential operator written in terms the Euler derivation.
    Currently it is also assumed that *stride* is a multiple of *block_len* and
    that the starting index (value of ``sum->n``) of a top-level recursive call
    is a multiple of *stride*.

    The called must supply an object *gbound* precomputed using
    :func:`acb_ode_group_bound_precompute`.

.. function:: void _acb_ode_sum_fix(acb_ode_sum_struct * sum)

    Restores some invariants of *sum* broken by the *sum_foward_\** functions.
    After a sequence of calls to those functions, *_acb_ode_sum_fix* must be
    called before accessing the computed sum.

.. function:: int acb_ode_sum_done(acb_ode_sum_struct * sum, slong stride, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)

    Returns nonzero when either the sum has converged (within a tolerance based
    on *prec*) or the current working precision is insufficient to significantly
    improve the current estimate.
    Updates some heuristic convergence estimates stored in *sum*, assuming that
    *stride* terms have been added since the last update.
    When returning nonzero, also updates the rigorous tail bounds (unless the
    :macro:`ACB_ODE_APPROX` flag is set in *sum*).
    The called must supply two objects *bound* and *gbound* precomputed using
    :func:`acb_ode_bound_precompute` and :func:`acb_ode_group_bound_precompute`.
    This function may be called on a *sum* object modified by *sum_foward_\**
    without first calling *sum_fix*.

Partial sums of individual solutions
-------------------------------------------------------------------------------

.. type:: acb_ode_sol_struct
          acb_ode_sol_t

    Represents a logarithmic series solution in the process of being computed.
    Mainly makes sense as a component of an :type:`acb_ode_sum_t` object or
    given such an object as context.

.. function:: void acb_ode_sol_init(acb_ode_sol_t sol, slong nshifts, slong nlogs, slong npts, slong nder)

    Initializes an :type:`acb_ode_sol_t` object suitable for a solution
    belonging to a solution group with *nshifts* shifts, involving `\log(x)` to
    powers less than *nlogs*, the derivatives of order up to *nder - 1* of
    which are to be evaluated at *npts* points.

.. function:: void acb_ode_sol_clear(acb_ode_sol_t sol)

    Clears *sol*.

.. function:: void acb_ode_sol_fit_length(acb_ode_sol_t sol, slong len)

    Ensures that the *series* field of *sol* has space for *len* coefficients.

.. function:: void acb_ode_sol_zero(acb_ode_sol_t sol)

    Zeroes the series coefficients and partial sums (without touching the
    initial conditions), resets bounds and termination flag.

.. function:: acb_ptr acb_ode_sol_sum_ptr(const acb_ode_sol_t sol, slong j, slong k, slong i)

    Returns a pointer to the coefficient of index *i* in the Taylor expansion
    at the evaluation point of index *j* of the partial sum of the coefficient
    of `\log(x)^k/k!`.

.. function:: void _acb_ode_sol_jet(acb_poly_struct * val, const acb_t expo, acb_srcptr sums, slong stride, mag_srcptr tb, const acb_t pt, slong nlogs, slong nder, slong nfrobshifts, slong prec)
              void acb_ode_sol_jet(acb_poly_struct * val, const acb_t expo, const acb_ode_sol_t sol, slong j, const acb_t pt, slong ord, slong nfrobshifts, slong prec)

    Computes the jets of order *ord* at the evaluation point *pt* of
    the Frobenius shifts of order 0 to *nfrobshifts* - 1 of the solution *sol*.
    The integer *j* is the index of the partial sum at *pt* in *sol*.
    The Frobenius shift of order *p* of `y(x)` is defined as the series obtained
    by replacing `\log(x)^k/k!` for `k \geq p` by `\log(x)^{k-p}/(k-p)!` in the
    series expansion and dropping the terms with `k < p`.

    The underscore version assumes *ord* > 0.

.. function:: slong acb_ode_sol_estimate_terms(mag_t est, const acb_ode_sol_t sol, slong off, slong len, const mag_t radpow)

    Sets *est* to heuristic estimate of the magnitude of the coefficients at
    positions *off* to *off + len* in the coefficient buffer of *sol*,
    multiplied by *radpow*.
    Returns an estimate of the relative accuracy in bits of these coefficients.

.. function:: slong acb_ode_sol_estimate_sums(mag_t mag, mag_t rad, const acb_ode_sol_t sol)

    Sets *mag* to a heuristic estimate of the magnitude of the partial sums
    encoded by *sol* and *rad* to an estimate of their radii.
    Returns an estimate of the minimum relative accuracy in bits of the sums.

Functions for working with differential operators
-------------------------------------------------------------------------------

.. function:: void _acb_ode_apply_diffop_basecase_precomp(acb_poly_struct * g, slong goff, acb_srcptr weights, slong weights_nlogs, const acb_poly_struct * f, slong foff, slong flen, slong nlogs, slong start, slong len, slong prec)
              void acb_ode_apply_diffop_basecase(acb_poly_struct * g, const acb_poly_struct * dop, slong dop_len, acb_srcptr expo, slong offset, const acb_poly_struct * f, slong nlogs, slong start, slong len, slong prec)
              void _acb_ode_apply_diffop_polmul(acb_poly_struct * g, slong goff, const acb_poly_struct * dop, slong dop_len, acb_srcptr expo, slong offset, const acb_poly_struct * f, slong foff, slong flen, slong nlogs, slong start, slong len, slong prec)
              void acb_ode_apply_diffop_polmul(acb_poly_struct * g, const acb_poly_struct * dop, slong dop_len, acb_srcptr expo, slong offset, const acb_poly_struct * f, slong nlogs, slong start, slong len, slong prec)

    Computes a slice of the coefficients of the image of a series by a
    differential operator.

    More precisely, the operator *dop* is applied to the logarithmic series

    .. math::

        f = x^{\mathit{expo} + \mathit{offset}} \sum_{p,k} f[k,p] x^p \log(x)^k/k!,

    resulting in a logarithmic series

    .. math::

        x^{\mathit{expo} + \mathit{offset}}
            \sum_{m',k'} g[k',m'] x^{m'} \log(x)^{k'}/k'!

    whose terms of index `m' = start + p'` for `0 \leq p' < \mathit{len}` are
    written to *g*.

    The non-underscore versions represent *f* and *g* as arrays of polynomials
    corresponding to the coefficients of `\log(x)^k/k!`.
    The underscore versions read the coefficients of `f[k]` (resp. write those
    of `g[k]`) starting at offset *foff* (resp. *goff*) in the coefficient array
    of the element of index `k` of the array *f* (resp. *g*).

    The *polmul* version uses fast polynomial multiplication. The *basecase*
    version uses a quadratic-time (w.r.t. *len*) algorithm. The *precomp*
    version is similar to the *basecase* version but takes as input, instead of
    the differential operator, an array of weights that can be precomputed using
    :func:`_acb_ode_apply_diffop_basecase_weights`.

.. function:: void _acb_ode_apply_diffop_basecase_weights(acb_ptr weights, const acb_poly_struct * dop, slong dop_len, acb_srcptr expo, slong offset, slong flen, slong nlogs, slong start, slong len, slong prec)

    Sets *weights* to an array of *len* × *nlogs* × *flen* weights that can be
    used to apply the differential operator *dop* to multiple series using
    :func:`_acb_ode_apply_diffop_basecase_precomp`.

Other summation utilities
-------------------------------------------------------------------------------

.. function:: void _acb_ode_poly_negdivrevhigh(acb_ptr res, acb_srcptr a, acb_srcptr cst, acb_srcptr b,  slong len, slong prec)

    Solves the equation

    .. math::

        a(1/X) \cdot \mathit{res}(X) + \mathit{cst} \cdot b(X) = (\text{terms of negative degree}),

    where *a*, *res*, *b* are polynomials of length *len*. The *cst* pointer may
    be null; this is understood as *cst = 1*.

.. function:: void acb_ode_poly_taylor_shift_aps_trunc(acb_poly_t g, const acb_poly_t f, acb_srcptr a, slong n, slong len, slong prec);

    Shifts the variable by *a + n*, truncates to length *len*.

Random generation
-------------------------------------------------------------------------------

The following function generate differential operators with a high probability
of being regular singular with a nontrivial exponent structure.

.. function:: void acb_ode_randtest_acb(acb_poly_struct * dop, acb_struct * lcroots, acb_ode_exponents_t expos, flint_rand_t state, slong disp, slong lcdeg, slong clen, slong len, slong prec)

    Writes to *dop* a random differential operator, to *lcroots* the roots of
    its leading coefficient, and to *expos* its exponent structure.
    The generated operator has length *len*, a leading coefficient of degree
    *lcdeg*, remaining coefficients of length at most *clen*, and exponent dispersion at most *disp*.

Error bounds
-------------------------------------------------------------------------------

The functions in this section operate simultaneously on one or more logarithmic
series solutions of the differential equation `L(y) = 0` belonging to the same
*group*, that is, of the form

.. math::

   y(x) = x^\lambda \sum_{k=0}^{K-1} y_k(x) \frac{\log(x)^k}{k!}
        = x^\lambda \sum_{k=0}^{K-1} \sum_{j=0}^{\infty} c_{k,j} x^j \frac{\log(x)^k}{k!}.

with the same `\lambda`.
Given one such solution, the main goal is to bound the tails

.. math::

   y_{k,n:}(a) = \sum_{j=n}^{\infty} c_{k,j} a^j, \qquad 0 \leq k < K,

of the power series `y_k` truncated to some order `n` and evaluated at
`a \in \mathbb C`.

The tail bound is computed from the *normalized residual* `q(x)` defined as
follows. Let `Q_0` be the monic indicial polynomial of `L` at the origin.
The equation

.. math::

    (Q_0(x\cdot\mathrm d/\mathrm dx))(q) = L(y_{n:})

has a unique solution of the form

.. math::

    q(x) = x^{\lambda+n} \sum_{k} q_{k}(x) \frac{\log(x)^k}{k!}
         = x^{\lambda+n} \sum_{k,j} q_{k,j} x^j \frac{\log(x)^k}{k!}

where the sums are finite and `q_{k,j} = 0` whenever `\lambda + n + j` is a root
of multiplicity strictly more than `k` of `Q_0`.

All functions in this section assume that the leading coefficient with respect
to the Euler derivative of the differential equation of interest does not
vanish at the origin. (For an equation with coprime coefficients, this is
equivalent to the origin being an ordinary point or a regular singular point.)

.. type:: acb_ode_bound_struct
          acb_ode_bound_t

    Stores precomputed bounds independent of the particular solution and the
    truncation order.

.. type:: acb_ode_group_bound_struct
          acb_ode_group_bound_t

    Stores additional precomputed bounds common to solutions belonging to the
    same group.

.. function:: void acb_ode_bound_init(acb_ode_bound_t bound)
              void acb_ode_bound_clear(acb_ode_bound_t bound)

.. function:: void acb_ode_group_bound_init(acb_ode_group_bound_t gbound)
              void acb_ode_group_bound_clear(acb_ode_group_bound_t gbound)

.. function:: void acb_ode_bound_precompute(acb_ode_bound_t bound, const acb_poly_struct * dop, slong dop_len, acb_srcptr lcrt, slong pol_part_len, slong prec)

    Fills *bound* with data used when computing bounds on the tails of series
    solutions of a given equation.

    The input consists of the operator *dop* of length *dop_len*,
    a vector *lcrt* containing enclosures of the roots of the leading
    coefficient (with multiple roots repeated),
    a minimum number *pol_part_len* of initial terms of the local series
    expansion of *dop* to be considered individually,
    and a working precision *prec*.

    As a general rule, increasing *pol_part_len* leads to sharper bounds but
    makes both the precomputation and subsequent evaluation of the bound slower.
    Values of *prec* beyond 30–50 bits are rarely useful except for equations
    with extremely large coefficients and/or clusters of singular points.

    When this function is called repeatedly with the same destination object,
    *dop* must correspond to the same operator (possibly represented at
    different precision) for all calls.
    In addition, changing *pol_part_len* invalidates any
    :type:`acb_ode_group_bound_t` objects derived from *bound*.

.. function:: void acb_ode_group_bound_precompute(acb_ode_group_bound_t gbound, const acb_poly_struct * dop, slong dop_len, const acb_ode_exponents_t expos, slong grp, const acb_ode_bound_t bound, slong prec)

    Fills *gbound* with data used for computing bounds on solutions from a given
    group.

    The input consists of the operator *dop* of length *dop_len*,
    its local exponent structure *expos*,
    the index *grp* of the group in the exponent structure,
    a group-independent bound data object computed with
    :func:`acb_ode_bound_precompute`,
    and a numeric working precision *prec*.

.. function:: void acb_ode_tail_bound_jet(arb_poly_t res, const acb_ode_group_t group, const acb_ode_bound_t bound, const acb_ode_group_bound_t gbound, slong n, slong nlogs, const arb_poly_t nres_maj, const arb_t rad, slong ord, slong prec)
              void acb_ode_tail_bound_jet_precomp(arb_poly_t res, const acb_ode_bound_t bound, slong n, const arb_poly_t itg_pol, const arb_poly_t itg_num, const arb_poly_t nres_maj, const arb_t rad, slong ord, slong prec)

    Computes upper bounds on the first *ord* coefficients of the Taylor
    expansion with respect to `z` of `y_{k,n:}(a + z)` where `|a|` is bounded
    by *rad*.
    The returned bounds are valid uniformly for all *k*, for all evaluation
    points `a` in the disk of radius *rad* around zero, and for all
    solutions `y` that belong to the given *group* and satisfy the following
    conditions:

    * The polynomial *nres_maj* is a common majorant of the
      coefficients `q_k(x)` of the normalized residual, i.e., for all `j`
      and `k`, the absolute value `|q_{k,j}|` is bounded by the coefficient
      of `x^j` in *nres_maj*.

    * The coefficients `c_{k,j}` with `j \geq n` such that `\lambda + j` is a root of the
      indicial polynomial `Q_0` of multiplicity *k* or more are zero.

    * The degree in `\log(x)` of the truncation at order `x^n` of `y(x)` is
      strictly smaller than *nlogs*. In other words, `c_{k,j} = 0` when `j < n`
      and `k \geq nlogs`.

    The caller must supply global and group-specific bound objects *bound*,
    *gbound* computed as described above.

    The *precomp* variant is suitable for repeated calls with the same *group*,
    the same value of *nlogs*, and values of *n* bounded from below.
    It takes two polynomials *itg_pol* and *itg_num* obtained by calling
    :func:`acb_ode_bound_precompute_integrand` with `n_0 \leq n`.
    Typically, values of *n0* closer to *n* yield sharper bounds, but the effect
    becomes less pronounced as *n* grows.

.. function:: void acb_ode_bound_precompute_integrand(arb_poly_t itg_pol, arb_poly_t itg_num, const acb_ode_group_t group, const acb_ode_bound_t bound, const acb_ode_group_bound_t gbound, slong n0, slong nlogs, slong prec)

    Computes polynomials *itg_pol* and *itg_num* that can be passed to
    :func:`acb_ode_tail_bound_jet_precomp` to obtain bounds on the tails of
    solutions from the group *group* at truncation orders >= *n0*.

    The remaining input arguments have the same meaning as the corresponding
    arguments of :func:`acb_ode_tail_bound_jet`.

Bounds on rational sequences
-------------------------------------------------------------------------------

The functions in this section compute upper bounds of the expression

.. math::

   F_T(n) = n \cdot \sum_{t=0}^{T-1} \left| [X^t] \frac{p(n+X)}{X^{-m(n)} q(n+X)} \right|

where `q` is a monic polynomial, `m(n)` denotes the multiplicity of `n` as a
root of `q`, and `p` is a polynomial of strictly smaller degree than `q`.
All functions work with a vector of numerators `p` and a single denominator `q`.

The function and parameter names reflect that, in the intended application,
`q` is the monic indicial polynomial of the differential equation being solved.

Some functions accept a working precision *prec* for intermediate computations
even though their outputs are represented in fixed precision.

.. type:: acb_ode_ind_lbound_struct
          acb_ode_ind_lbound_t

    A precomputed lower bound on `q(n)/n^r` where `q` is a monic polynomial
    of degree `r`.
    The data structure contains an integer *length* and a vector *r* of
    length *length* of structures each containing
    a root *root* of `q`,
    its multiplicity *mult*,
    an integer *n_min* such that *|1-root/n|* is nondecreasing for *n >= n_min*,
    and a lower bound *global_lbound* of *|1-root/n|^mult* for *n >= n_min*.

.. function:: void acb_ode_ind_lbound_init(acb_ode_ind_lbound_t ind_lbound)
              void acb_ode_ind_lbound_clear(acb_ode_ind_lbound_t ind_lbound)

.. function:: void acb_ode_ind_lbound_precompute(acb_ode_ind_lbound_t ind_lbound, const acb_ode_exponents_t expos, slong grp, slong prec)

    Fills *ind_lbound* with data corresponding to `q(X) = q_0(\lambda + X)` where
    `q_0` is the monic polynomial whose roots are given by *expos* and `\lambda` is
    the leader of the group of index *grp*.

In the following descriptions `\tau(n) = \sum_{k=0}^{n} m(k)`.

.. type:: acb_ode_stairs_struct
          acb_ode_stairs_t

    A vector *h* of length *length* of *mag_struct*.

.. function:: void acb_ode_stairs_init(acb_ode_stairs_t stairs)
              void acb_ode_stairs_clear(acb_ode_stairs_t stairs)

.. function:: void acb_ode_stairs_precompute(acb_ode_stairs_t stairs, const acb_poly_struct * num, slong len, const acb_poly_t ind, const acb_ode_group_t group, const acb_ode_ind_lbound_t ind_lbound, slong prec)

    For each `p` in the vector *num* of length *len* and for each shift `s` in
    *group*, computes an upper bound of `F_{\tau(n)}(n)` valid for *n >= s*.
    The parameter *ind* corresponds to `q`.
    The parameter *group* must describe the integer roots of `q` and their
    multiplicities (its *leader* field is ignored).
    Additionally, the caller must supply an *ind_lbound* structure precomputed
    as described above.

.. function:: void acb_ode_bound_rat_vec(mag_ptr res, const acb_poly_struct * num, slong len, const acb_poly_t ind, const acb_ode_group_t group, const acb_ode_ind_lbound_t ind_lbound, const acb_ode_stairs_t stairs, slong n0, slong ord, slong prec)

    For each `p` in the vector *num* of length *len*, computes a bound on
    `F_{\rho(n)}(n)` valid for all `n \geq n_0`, where `\rho(n) = \tau(n)` when
    there is no integer root of `q` between `n_0` (exclusive) and `n`
    (inclusive), and `\rho(n)` is equal to *ord* otherwise.
    The parameter *ind* specifies the common denominator `q`.

    The caller must supply the integer roots of `q` in the parameter *group*
    (the *leader* field is ignored), as well as two data structures *ind_lbound*
    and *stairs* precomputed as described above.

.. function:: void acb_ode_bound_rat_ordinary_vec(mag_ptr res, const acb_poly_struct * num, slong len, const acb_poly_t ind, const acb_ode_ind_lbound_t ind_lbound, slong n0, slong ord, slong prec)

    Similar but assumes `m(n_0) = 0` and computes a bound on `F_{T}(n)` valid
    for all `n \geq n_0` with `m(n) = 0`.
    The value of `T` is given by the parameter *ord*.

.. function:: void acb_ode_bound_rat_ref_vec(mag_ptr res, const acb_poly_struct * num, slong len, const acb_poly_t ind, slong n, slong mult, slong ord, slong prec)

    Similar but assumes that `m(n)` is equal to *mult* and computes a bound on
    `F_{T}(n)` at the given *n*.

Miscellaneous bounds and estimates
-------------------------------------------------------------------------------

.. function:: void _acb_ode_solution_growth(mag_t order, mag_t base, const acb_poly_struct * dop, slong dop_len)

   Estimates parameters *order* and *base* such that the coefficient sequences
   of series solutions of *dop* are
   `\mathrm O(\mathit{base}^n/n!^{1/\mathit{order}})`.
   Assumes that the leading coefficient of *dop* is a constant.
   The return values are currently not rigorous bounds.

.. function:: slong acb_ode_choose_prec(slong * rec_prec, const acb_poly_struct * dop, slong dop_len, mag_srcptr rad, mag_srcptr cvrad, slong tgt_prec)

.. type:: acb_ode_cvest_struct
          acb_ode_cvest_t

.. function:: void acb_ode_cvest_init(acb_ode_cvest_t cvest)
              void acb_ode_cvest_clear(acb_ode_cvest_t cvest)

.. function:: void acb_ode_cvest_update(acb_ode_cvest_t cvest, const acb_ode_cvest_t old, mag_t est, slong accuracy, slong stride, slong prec, slong work_prec)

Examples
-------------------------------------------------------------------------------

.. _ode-example-basic:

Basic usage for computing a fundamental matrix:

.. literalinclude:: ../../examples/ode_fundamental_matrix.c
   :language: c
