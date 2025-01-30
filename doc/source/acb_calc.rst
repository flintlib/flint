.. _acb-calc:

**acb_calc.h** -- calculus with complex-valued functions
===============================================================================

This module provides functions for operations of calculus
over the complex numbers (intended to include root-finding,
integration, and so on).
The numerical integration code is described in [Joh2018a]_.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: acb_calc_func_t

    Typedef for a pointer to a function with signature::

        int func(acb_ptr out, const acb_t inp, void * param, slong order, slong prec)

    implementing a univariate complex function `f(z)`.
    The *param* argument may be used to pass through
    additional parameters to the function.
    The return value is reserved for future use as an
    error code. It can be assumed that *out* and *inp* are not aliased.

    When called with *order* = 0, *func* should write to *out* the value
    of `f(z)` at the point *inp*, evaluated at a precision of *prec* bits.
    In this case, *f* can be an arbitrary complex function, which
    may have branch cuts or even be non-holomorphic.

    When called with *order* = *n* for `n \ge 1`, *func* should write to
    *out* the first *n* coefficients in the Taylor series expansion of `f(z)`
    at the point *inp*, evaluated at a precision of *prec* bits.
    In this case, the implementation of *func* must verify that *f*
    is holomorphic on the complex interval defined by *z*, and set the
    coefficients in *out* to non-finite values otherwise.

    For algorithms that do not rely on derivatives, *func* will always
    get called with *order* = 0 or *order* = 1, in which case the user
    only needs to implement evaluation of the direct function value `f(z)`
    (without derivatives). With *order* = 1, *func* must verify
    holomorphicity (unlike the *order* = 0 case).

    If *f* is built from field operations and meromorphic functions, then
    no special action is necessary when *order* is positive
    since division by zero or evaluation
    of builtin functions at poles automatically produces infinite enclosures.
    However, manual action is needed for bounded functions with branch cuts.
    For example, when evaluating `\sqrt{z}`, the output must be set to
    an non-finite value if `z` overlaps with the branch cut `[-\infty,0]`.
    The easiest way to accomplish this is to use versions of basic
    functions (sqrt, log, pow, etc.) that test holomorphicity of their
    arguments individually.

    Some functions with branch cut detection are available as builtins:
    see :func:`acb_sqrt_analytic`,
    :func:`acb_rsqrt_analytic`, :func:`acb_log_analytic`,
    :func:`acb_pow_analytic`. It is not difficult to write custom functions
    of this type, using the following pattern::

        /* Square root function on C with detection of the branch cut. */
        void sqrt_analytic(acb_t res, const acb_t z, int analytic, slong prec)
        {
            if (analytic &&
                arb_contains_zero(acb_imagref(z)) &&
                arb_contains_nonpositive(acb_realref(z))))
            {
                acb_indeterminate(res);
            }
            else
            {
                acb_sqrt(res, z, prec);
            }
        }

    The built-in methods :func:`acb_real_abs`, :func:`acb_real_sgn`,
    :func:`acb_real_heaviside`, :func:`acb_real_floor`, :func:`acb_real_ceil`,
    :func:`acb_real_max`, :func:`acb_real_min` provide piecewise holomorphic
    functions that are useful for integrating piecewise-defined real functions.

    For example, here we define a piecewise holomorphic extension
    of the function
    `f(z) = \sqrt{\lfloor z \rfloor}` (for simplicity, without implementing
    derivatives)::

        int func(acb_ptr out, const acb_t inp, void * param, slong order, slong prec)
        {
            if (order > 1) flint_abort();  /* derivatives not implemented */

            acb_real_floor(out, inp, order != 0, prec);
            acb_sqrt_analytic(out, out, order != 0, prec);
            return 0;
        }

    (Here, :func:`acb_real_sqrtpos` may be slightly better if it is
    known that *z* will be nonnegative on the path.)

    See the demo program ``examples/integrals.c`` for more examples.

Integration
-------------------------------------------------------------------------------

.. function:: int acb_calc_integrate(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, slong rel_goal, const mag_t abs_tol, const acb_calc_integrate_opt_t options, slong prec)

    Computes a rigorous enclosure of the integral

    .. math::

        I = \int_a^b f(t) dt

    where *f* is specified by (*func*, *param*), following a straight-line
    path between the complex numbers *a* and *b*.
    For finite results, *a*, *b* must be finite and *f* must be bounded
    on the path of integration.
    To compute improper integrals, the user should therefore truncate the path
    of integration manually (or make a regularizing change of variables,
    if possible).
    Returns *ARB_CALC_SUCCESS* if the integration converged to the
    target accuracy on all subintervals, and returns
    *ARB_CALC_NO_CONVERGENCE* otherwise.

    By default, the integrand *func* will only be called with *order* = 0
    or *order* = 1; that is, derivatives are not required.

    - The integrand will be called with *order* = 0 to evaluate *f*
      normally on the integration path (either at a single point
      or on a subinterval). In this case, *f* is treated as a pointwise defined
      function and can have arbitrary discontinuities.

    - The integrand will be called with *order* = 1 to evaluate *f*
      on a domain surrounding a segment of the integration path for the purpose
      of bounding the error of a quadrature formula. In this case, *func* must
      verify that *f* is holomorphic on this domain (and output a non-finite
      value if it is not).

    The integration algorithm combines direct interval enclosures,
    Gauss-Legendre quadrature where *f* is holomorphic,
    and adaptive subdivision. This strategy supports integrands with
    discontinuities while providing exponential convergence for typical
    piecewise holomorphic integrands.

    The following parameters control accuracy:

    - *rel_goal* - relative accuracy goal as a number of bits, i.e.
      target a relative error less than `\varepsilon_{rel} = 2^{-r}`
      where *r* = *rel_goal*
      (note the sign: *rel_goal* should be nonnegative).

    - *abs_tol* - absolute accuracy goal as a :type:`mag_t` describing
      the error tolerance, i.e.
      target an absolute error less than `\varepsilon_{abs}` = *abs_tol*.

    - *prec* - working precision. This is the working precision used to
      evaluate the integrand and manipulate interval endpoints.
      As currently implemented, the algorithm does not attempt to adjust the
      working precision by itself, and adaptive
      control of the working precision must be handled by the user.

    For typical usage, set *rel_goal* = *prec* and *abs_tol* = `2^{-prec}`.
    It usually only makes sense to have *rel_goal* between 0 and *prec*.

    The algorithm attempts to achieve an error of
    `\max(\varepsilon_{abs}, M \varepsilon_{rel})` on each subinterval,
    where *M* is the magnitude of the integral.
    These parameters are only guidelines; the cumulative error may be larger
    than both the prescribed
    absolute and relative error goals, depending on the number of
    subdivisions, cancellation between segments of the integral, and numerical
    errors in the evaluation of the integrand.

    To compute tiny integrals with high relative accuracy, one should set
    `\varepsilon_{abs} \approx M \varepsilon_{rel}` where *M* is a known
    estimate of the magnitude. Setting `\varepsilon_{abs}` to 0 is also
    allowed, forcing use of a relative instead of an absolute tolerance goal.
    This can be handy for exponentially small or
    large functions of unknown magnitude. It is recommended to avoid
    setting `\varepsilon_{abs}` very small
    if possible since the algorithm might need many extra
    subdivisions to estimate *M* automatically; if the approximate
    magnitude can be estimated by some external means (for example if
    a midpoint-width or endpoint-width estimate is known to be accurate),
    providing an appropriate `\varepsilon_{abs} \approx M \varepsilon_{rel}`
    will be more efficient.

    If the integral has very large magnitude, setting the absolute
    tolerance to a corresponding large value is recommended for best
    performance, but it is not necessary for convergence since the absolute
    tolerance is increased automatically during the execution of the
    algorithm if the partial integrals are found to have larger error.

    Additional options for the integration can be provided via the *options*
    parameter (documented below). To use all defaults, *NULL* can be passed
    for *options*.

Options for integration
...............................................................................

.. type:: acb_calc_integrate_opt_struct

.. type:: acb_calc_integrate_opt_t

    This structure contains several fields, explained below.
    An *acb_calc_integrate_opt_t* is defined as an array of
    *acb_calc_integrate_opt_struct*
    of length 1, permitting it to be passed by reference.
    An *acb_calc_integrate_opt_t* must be initialized before use, which sets
    all fields to 0 or *NULL*. For fields that have not been set to other
    values, the integration algorithm will choose defaults automatically
    (based on the precision and accuracy goals).
    This structure will most likely be extended in the future to
    accommodate more options.

    .. member:: slong deg_limit

        Maximum quadrature degree for each subinterval.
        If a zero or negative value is provided, the limit is set to a default
        value which currently equals `0.5 \cdot \min(prec, rel\_goal) + 60` for
        Gauss-Legendre quadrature.
        A higher quadrature degree can be beneficial for functions that
        are holomorphic on a large domain around the integration path
        and yet behave irregularly, such as oscillatory entire functions.
        The drawback of increasing the degree is that
        the precomputation time for quadrature nodes increases.

    .. member:: slong eval_limit

        Maximum number of function evaluations.
        If a zero or negative value is provided, the limit is set to a default
        value which currently equals `1000 \cdot prec + prec^2`.
        This is the main parameter used to limit the amount of work before
        aborting due to possible slow convergence or non-convergence.
        A lower limit allows aborting faster. A higher limit may be needed
        for integrands with many discontinuities or many singularities
        close to the integration path.
        This limit is only taken as a rough guideline, and the actual number of
        function evaluations may be slightly higher depending on the
        actual subdivisions.

    .. member:: slong depth_limit

        Maximum search depth for adaptive subdivision. Technically, this is not
        the limit on the local bisection depth but the limit on the number
        of simultaneously queued subintervals.
        If a zero or negative value is provided, the limit is set to the
        default value `2 \cdot \text{prec}`.
        Warning: memory usage may increase in proportion to this limit.

    .. member:: int use_heap

        By default (if set to 0), new subintervals generated by adaptive
        bisection will be appended to the top of a stack.
        If set to 1, a binary heap will be used to maintain a priority queue
        where the subintervals with larger error have higher priority.
        This sometimes gives better results
        in case of convergence failure, but can
        lead to a much larger array of subintervals (requiring a higher
        *depth_limit*) when many global bisections are needed.

    .. member:: int verbose

        If set to 1, some information about the overall integration process
        is printed to standard output. If set to 2, information about each
        subinterval is printed.

.. function:: void acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options)

    Initializes *options* for use, setting all fields to 0 indicating
    default values.

Local integration algorithms
-------------------------------------------------------------------------------

.. function:: int acb_calc_integrate_gl_auto_deg(acb_t res, slong * num_eval, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const mag_t tol, slong deg_limit, int flags, slong prec)

    Attempts to compute `I = \int_a^b f(t) dt` using a single application
    of Gauss-Legendre quadrature with automatic determination of the
    quadrature degree so that the error is smaller than *tol*.
    Returns *ARB_CALC_SUCCESS* if the integral has been evaluated successfully
    or *ARB_CALC_NO_CONVERGENCE* if the tolerance could not be met.
    The total number of function evaluations is written to *num_eval*.

    For the interval `[-1,1]`, the error of the *n*-point Gauss-Legendre
    rule is bounded by

    .. math::

        \left| I - \sum_{k=0}^{n-1} w_k f(x_k) \right| \le \frac{64 M}{15 (\rho-1) \rho^{2n-1}}

    if `f` is holomorphic with `|f(z)| \le M` inside the ellipse *E*
    with foci `\pm 1` and semiaxes
    `X` and `Y = \sqrt{X^2 - 1}` such that `\rho = X + Y`
    with `\rho > 1` [Tre2008]_.

    For an arbitrary interval, we use `\int_a^b f(t) dt = \int_{-1}^1 g(t) dt`
    where `g(t) = \Delta f(\Delta t + m)`,
    `\Delta = \tfrac{1}{2}(b-a)`, `m = \tfrac{1}{2}(a+b)`.
    With `I = [\pm X] + [\pm Y]i`, this means that we evaluate
    `\Delta f(\Delta I + m)` to get the bound `M`.
    (An improvement would be to reduce the wrapping effect of rotating the
    ellipse when the path is not rectilinear).

    We search for an `X` that makes the error small by trying steps `2^{2^k}`.
    Larger `X` will give smaller `1 / \rho^{2n-1}` but larger `M`. If we try
    successive larger values of `k`, we can abort when `M = \infty`
    since this either means that we have hit a singularity or a branch cut or
    that overestimation in the evaluation of `f` is becoming too severe.

Integration (old)
-------------------------------------------------------------------------------

.. function:: void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func, void * param, const acb_t x, const arb_t radius, slong maxdepth, slong prec)

    Sets *bound* to a ball containing the value of the integral

    .. math::

        C(x,r) = \frac{1}{2 \pi r} \oint_{|z-x| = r} |f(z)| dz
               = \int_0^1 |f(x+re^{2\pi i t})| dt

    where *f* is specified by (*func*, *param*) and *r* is given by *radius*.
    The integral is computed using a simple step sum.
    The integration range is subdivided until the order of magnitude of *b*
    can be determined (i.e. its error bound is smaller than its midpoint),
    or until the step length has been cut in half *maxdepth* times.
    This function is currently implemented completely naively, and
    repeatedly subdivides the whole integration range instead of
    performing adaptive subdivisions.

.. function:: int acb_calc_integrate_taylor(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const arf_t inner_radius, const arf_t outer_radius, slong accuracy_goal, slong prec)

    Computes the integral

    .. math::

        I = \int_a^b f(t) dt

    where *f* is specified by (*func*, *param*), following a straight-line
    path between the complex numbers *a* and *b* which both must be finite.

    The integral is approximated by piecewise centered Taylor polynomials.
    Rigorous truncation error bounds are calculated using the Cauchy integral
    formula. More precisely, if the Taylor series of *f* centered at the point
    *m* is `f(m+x) = \sum_{n=0}^{\infty} a_n x^n`, then

    .. math::

        \int f(m+x) = \left( \sum_{n=0}^{N-1} a_n \frac{x^{n+1}}{n+1} \right)
                  + \left( \sum_{n=N}^{\infty} a_n \frac{x^{n+1}}{n+1} \right).

    For sufficiently small *x*, the second series converges and its
    absolute value is bounded by

    .. math::

        \sum_{n=N}^{\infty} \frac{C(m,R)}{R^n} \frac{|x|^{n+1}}{N+1}
            = \frac{C(m,R) R x}{(R-x)(N+1)} \left( \frac{x}{R} \right)^N.

    It is required that any singularities of *f* are
    isolated from the path of integration by a distance strictly
    greater than the positive value *outer_radius* (which is the integration
    radius used for the Cauchy bound). Taylor series step lengths are
    chosen so as not to
    exceed *inner_radius*, which must be strictly smaller than *outer_radius*
    for convergence. A smaller *inner_radius* gives more rapid convergence
    of each Taylor series but means that more series might have to be used.
    A reasonable choice might be to set *inner_radius* to half the value of
    *outer_radius*, giving roughly one accurate bit per term.

    The truncation point of each Taylor series is chosen so that the absolute
    truncation error is roughly `2^{-p}` where *p* is given by *accuracy_goal*
    (in the future, this might change to a relative accuracy).
    Arithmetic operations and function
    evaluations are performed at a precision of *prec* bits. Note that due
    to accumulation of numerical errors, both values may have to be set
    higher (and the endpoints may have to be computed more accurately)
    to achieve a desired accuracy.

    This function chooses the evaluation points uniformly rather
    than implementing adaptive subdivision.

