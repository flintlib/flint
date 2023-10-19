.. _acb-theta:

**acb_theta.h** -- Riemann theta functions in any dimension
===============================================================================

This module provides methods for the numerical evaluation of theta functions in
any dimension `g`. The algorithms will be detailed in the forthcoming paper
[EK]_. In the case `g=1`, we rely, but also improve on functionality from
:ref:`acb_modular.h <acb-modular>`.

In the context of this module, *tau* or `\tau` always denotes an element of the
Siegel complex upper half-space `\mathbb{H}_g`, which consists of all symmetric
`g\times g` complex matrices with positive definite imaginary part. The letter
`z` denotes an element of `\mathbb{C}^g`. For each `a,b\in \{0,1\}^g`, the
Riemann theta function of characteristic `(a,b)` is the following analytic
function in two variables `\tau\in \mathbb{H}_g` and `z\in \mathbb{C}^g`:

    .. math ::

        \theta_{a,b}(z,\tau) = \sum_{n\in \mathbb{Z}^{g} + \tfrac a2} \exp(\pi i n^T\tau n + 2\pi i n^T (z + \tfrac b2)),

considering `a`, `b` and `z` as column vectors.

We encode a theta characteristic `a\in \{0,1\}^g` as the :type:`ulong` between
`0` and `2^g-1` that has the corresponding expansion in base 2: thus `a =
(1,0,0)` for `g = 3` will be numbered `8`. We also use this encoding to order
vectors of theta values throughout. Similarly, a pair of characteristics
`(a,b)` is encoded as an :type:`ulong` between `0` and `2^{2g}-1`, where `a`
corresponds to the `g` more significant bits. Note that with these conventions,
the output of :func:`acb_modular_theta` is
`(-\theta_3,\theta_2,\theta_0,\theta_1)`.

The main user-facing function to evaluate theta functions is
:func:`acb_theta_all`. This function first reduces the input `(z,\tau)` using
the action of the Siegel modular group `\mathrm{Sp}_{2g}(\mathbb{Z})` on
`\mathbb{C}^g\times \mathbb{H}_g`, then uses a quasi-linear algorithm to
compute theta values on the reduced domain. At low precisions and when `\tau`
is reasonably reduced, one may also consider using "naive algorithms" directly,
which consist in evaluating a partial sum of the theta series. The main
functions to do so are `acb_theta_naive_00` and `acb_theta_naive_all`. We also
provide functionality to evaluate Siegel modular forms in terms of theta
functions when `g=2`.

As usual, the numerical functions in this module compute certified error
bounds: for instance, if `\tau` is represented by an :type:`acb_mat_t` which is
not certainly positive definite at the chosen working precision, the output
will have an infinite radius. Some internal functions may abort on
ill-conditioned input as indicated in the documentation below.

Main user functions
-------------------------------------------------------------------------------



The Siegel modular group
-------------------------------------------------------------------------------

We use the type `fmpz_mat_t` to handle matrices in
`\operatorname{Sp}_{2g}(\mathbb{Z})`. A large number of methods from
:ref:`fmpz_mat.h <fmpz_mat>`, e.g. :func:`fmpz_mat_equal`, can thus be used
directly on symplectic matrices.

In the following functions (with the exception of :func:`sp2gz_is_correct`) we
always assume that the input matrix *mat* is square of even size `2g`, and
write it as `m = \left(\begin{smallmatrix} \alpha&\beta\\ \gamma&\delta
\end{smallmatrix}\right)` where `\alpha,\beta,\gamma,\delta` are `g\times g`
blocks.

.. function:: slong sp2gz_dim(const fmpz_mat_t mat)

    Returns the dimension `g`, which is half the number of rows (or columns)
    of *mat*.

.. function:: void sp2gz_set_blocks(fmpz_mat_t mat, const fmpz_mat_t alpha, const fmpz_mat_t beta, const fmpz_mat_t gamma, const fmpz_mat_t delta)

    Sets *mat* to `\left(\begin{smallmatrix} \alpha&\beta\\ \gamma&\delta
    \end{smallmatrix}\right)`.

.. function:: void sp2gz_j(fmpz_mat_t mat)

    Sets *mat* to the symplectic matrix `J = \left(\begin{smallmatrix}
    0&I_g\\-I_g&0 \end{smallmatrix}\right)`.

.. function:: void sp2gz_block_diag(fmpz_mat_t mat, const fmpz_mat_t U)

    Sets *mat* to the symplectic matrix `\left(\begin{smallmatrix}
    U&0\\0&U^{-T} \end{smallmatrix}\right)`. We require that `U\in
    \operatorname{GL}_g(\mathbb{Z})`.

.. function:: void sp2gz_trig(fmpz_mat_t mat, const fmpz_mat_t S)

    Sets *mat* to `\left(\begin{smallmatrix} I_g&S\\0&I_g
    \end{smallmatrix}\right)`, where *S* is a symmetric `g\times g` matrix.

.. function:: void sp2gz_embed(fmpz_mat_t res, const fmpz_mat_t mat)

    Assuming that *mat* is a symplectic matrix of size `2r\times 2r` and *res*
    is square of size `2g\times 2g` for some `g\geq r`, sets *res* to the symplectic matrix

    .. math ::

        \begin{pmatrix} \alpha && \beta & \\ & I_{g-r} && 0_{g-r} \\ \gamma &&\delta &\\ & 0_{g-r} && I_{g-r} \end{pmatrix}

    where `\alpha,\beta,\gamma,\delta` are the `r\times r` blocks of *mat*.

.. function:: void sp2gz_restrict(fmpz_mat_t res, const fmpz_mat_t mat)

    Assuming that *mat* is a symplectic matrix of size `2g\times 2g` and *res*
    is square of size `2r\times 2r` for some `r\leq g`, sets *res* to the
    matrix whose `r\times r` blocks are the upper left corners of the
    corresponding `g\times g` block of *mat*. The result may not be a
    symplectic matrix.

.. function:: slong sp2gz_nb_fundamental(slong g)

    Returns the number of fundamental symplectic matrices used in the reduction
    algorithm on `\mathbb{H}_g`. This number is 1 when `g=1` (the `J` matrix)
    and 19 when `g=2` [Got1959]_. When `g>2`, a complete set of matrices
    defining the boundary of a fundamental domain for the action of
    `\mathrm{Sp}_{2g}(\mathbb{Z})` is not currently known. As a substitute, we
    consider two types of matrices: the `19 g(g-1)/2` matrices obtained by
    mimicking the `g=2` matrices on any pair of indices between 0 and `g-1`,
    and the `2^g` matrices obtained by embedding a copy of a lower-dimensional
    `J` matrix on any subset of indices.

.. function:: void sp2gz_fundamental(fmpz_mat_t mat, slong j)

    Sets *mat* to the `j^{\text{th}}` fundamental symplectic matrix as defined
    above.

.. function:: int sp2gz_is_correct(const fmpz_mat_t mat)

    Returns true (nonzero) iff *mat* is a symplectic matrix.

.. function:: int sp2gz_is_j(const fmpz_mat_t mat)

    Returns true (nonzero) iff the symplectic matrix *mat* is the `J` matrix.

.. function:: int sp2gz_is_block_diag(const fmpz_mat_t mat)

    Returns true (nonzero) iff the symplectic matrix *mat* is of block-diagonal
    form as in :func:`sp2gz_block_diag`.

.. function:: int sp2gz_is_trig(const fmpz_mat_t mat)

    Returns true (nonzero) iff the sympletic matrix *mat* is of trigonal form
    as in :func:`sp2gz_trig`.

.. function:: int sp2gz_is_embedded(fmpz_mat_t res, const fmpz_mat_t mat)

    Assuming that *mat* is a `2g\times 2g` symplectic matrix and *res* is
    square of size `2r` for some `r\leq g`, returns true (nonzero) iff *mat*
    can be obtained as the result of :func:`sp2gz_embed` from a `2r\times 2r`
    symplectic matrix, and store this matrix in *res*. Otherwise, returns false
    and leaves *res* undefined.

.. function:: void sp2gz_inv(fmpz_mat_t inv, const fmpz_mat_t mat)

    Sets *inv* to the inverse of the symplectic matrix *mat*.

.. function:: fmpz_mat_struct* sp2gz_decompose(slong* nb, const fmpz_mat_t mat)

    Returns a vector *res* of symplectic matrices and store its length in *nb*
    such that the following holds: *mat* is the product of the elements of
    *res* from left to right, and each element of *res* is block-diagonal,
    trigonal, the `J` matrix, an embedded `J` matrix from a lower dimension, or
    an embedded matrix from dimension 1 (i.e. `\mathrm{SL}_2(\mathbb{Z})`). The
    output vector *res* will need to be freed by the user as follows:

    .. code-block:: c

        slong k;
        for (k = 0; k < *nb; k++)
        {
            fmpz_mat_clear(&res[k]);
        }
        flint_free(res);

.. function:: void sp2gz_randtest(fmpz_mat_t mat, flint_rand_t state, slong bits)

    Sets *mat* to a random symplectic matrix whose coefficients have length
    approximately *bits*, obtained as a product of block-diagonal and trigonal
    symplectic matrices and the `J` matrix.

The Siegel half space
-------------------------------------------------------------------------------

We continue to denote by `\alpha,\beta,\gamma,\delta` the `g\times g` blocks of
*mat*, which is always assumed to be symplectic.

.. function:: void acb_siegel_cocycle(acb_mat_t c, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Sets *c* to `\gamma\tau + \delta`.

.. function:: void acb_siegel_transform_cocycle_inv(acb_mat_t w, acb_mat_t c, acb_mat_t cinv, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Sets *w*, *c* and *cinv* to `(\alpha\tau + \beta)(\gamma\tau +
    \delta)^{-1}`, `\gamma\tau + \delta` and `(\gamma\tau + \delta)^{-1}`
    respectively.

.. function:: void acb_siegel_transform(acb_mat_t w, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Sets *w* to `(\alpha\tau + \beta)(\gamma\tau + \delta)^{-1}`.

.. function:: void acb_siegel_transform_z(acb_ptr r, acb_mat_t w, const fmpz_mat_t mat, acb_srcptr z, const acb_mat_t tau, slong prec)

    Sets *w* to `(\alpha\tau + \beta)(\gamma\tau + \delta)^{-1}` and *r* to
    `(\gamma\tau + \delta)^{-T}z`.

.. function:: void acb_siegel_reduce(fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Sets *mat* to a symplectic matrix such that `\mathit{mat}\cdot\tau` is as
    reduced as possible, repeatedly reducing the imaginary and real parts of
    `\tau` and applying fundamental symplectic matrices. If the coefficients of
    `\tau` do not have a reasonable size or if `\det \mathrm{Im}(\tau)` is
    vanishingly small, we simply set *mat* to the identity.

.. function:: int acb_siegel_is_reduced(const acb_mat_t tau, slong tol_exp, slong prec)

    Returns true (nonzero) iff it is certainly true that *tau* belongs to the
    reduced domain defined by the tolerance parameter `\varepsilon =
    2^{-\mathit{tol_exp}}`. This means the following:
    `|\mathrm{Re}(\tau_{j,k})| < \frac12 + \varepsilon` for all `0\leq j,k <
    g`, the imaginary part of *tau* passes \func{arb_mat_spd_is_lll_reduced}
    with the same parameters, and for every matrix *mat* obtained from
    :func:`sp2gz_fundamental`, the determinant of the corresponding cocycle is
    at least `1-\eps`.

.. function:: void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)

    Sets *tau* to a random matrix in `\mathbb{H}_g`, possibly far from being
    reduced.

.. function:: void acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)

    Sets *tau* to a random reduced matrix in `\mathbb{H}_g` by calling
    :func:`acb_siegel_reduce` on a random matrix. The reduction may fail at low
    precisions for a given choice of *g* and *mag_bits*, in which case the
    output will be similar to :func:`acb_siegel_randtest`.

.. function:: void acb_siegel_randtest_nice(acb_mat_t tau, flint_rand_t state, slong prec)

    Sets *tau* to a random matrix that is well within the reduced domain in
    `\mathbb{H}_g`.

Theta characteristics
-------------------------------------------------------------------------------

.. function:: void acb_theta_char_get_slong(slong* n, ulong a, slong g)

    Sets each entry of *n* to the corresponding bit of *a*.

.. function:: ulong acb_theta_char_get_a(const slong* n, slong g)

    Returns the unique characteristic *a* such that `n\in 2\mathbb{Z}^g + a`.

.. function:: void acb_theta_char_get_arb(arb_ptr v, ulong a, slong g)

.. function:: void acb_theta_char_get_acb(acb_ptr v, ulong a, slong g)

    Sets *v* to `a/2` seen as an element of `\mathbb{R}^g` or `\mathbb{C}^g`
    respectively.

.. function:: slong acb_theta_char_dot(ulong a, ulong b, slong g)

    Returns `\sum_{i=0}^{g-1} a_i b_i` modulo 4 as an integer between 0 and 3,
    where `a_i, b_i` for `0\leq i < g` denote the bits of `a` and `b`
    respectively.

.. function:: slong acb_theta_char_dot_slong(ulong a, const slong* n, slong g)

    Returns `\sum_{i=0}^{g-1} a_i n_i` modulo 4 as an integer between 0 and 3.

.. function:: void acb_theta_char_dot_acb(acb_t x, ulong a, acb_srcptr z, slong g, slong prec)

    Sets *x* to `\sum_{i=0}^{g-1} a_i z_i`.

.. function:: int acb_theta_char_is_even(ulong ab, slong g)

    Returns true iff the characteristic `(a,b)` is even, i.e. `a^Tb` is divisible by 2.

.. function:: int acb_theta_char_is_goepel(ulong ch1, ulong ch2, ulong ch3, ulong ch4, slong g)

    Returns true iff the given characteristics define a Göpel quadruple,
    i.e. they are distinct even characteristics whose sum belongs to
    `2\mathbb{Z}^g`.

.. function:: int acb_theta_char_is_syzygous(ulong ch1, ulong ch2, ulong ch3, slong g)

    Returns true iff the given characteristics define a syzygous triple,
    i.e. they can be completed into a Göpel quadruple.


Ellipsoids: types and macros
-------------------------------------------------------------------------------

Following [DHBHS2004]_, we will consider partial sums of theta series over
points `n` in the lattice `\mathbb{Z}^g` contained in certain ellipsoids.

Fix `1\leq d\leq g`, an upper-triangular Cholesky matrix `C`, a radius `R\geq
0`, a vector `v\in \mathbb{R}^g`, and integers `n_{d},\ldots,
n_{g-1}`. Consider the ellipsoid `E` consisting of points `n =
(n_0,\ldots,n_{g-1})` satisfying `(v + Cn)^T(v + Cn)\leq R^2`. We encode `E` as
follows: we store the endpoints and midpoint of the interval of allowed values
for `n_{d-1}` as \type{slong}'s, and if `d\geq 1`, we also store a
`(d-1)`-dimensional ``child'' of `E` for each value of `n_{d-1}`. Children are
partitioned between left and right children depending on the position of
`n_{d-1}` relative to the midpoint (by convention, the midpoint is a right
child). When `d=g` and for a fixed Cholesky matrix `C`, this representation
uses `O(R^{g-1})` space for an ellipsoid of radius `R` containing approximately
`O(R^{g})` points.

.. type:: acb_theta_eld_struct

.. type:: acb_theta_eld_t

    An :type:`acb_theta_eld_t` is an array of length one of type
    :type:`acb_theta_eld_struct` encoding an ellipsoid as described above,
    permitting it to be passed by reference.

The following macros are available after `E` of type :func`acb_theta_eld_t` has
been initialized using :func:`acb_theta_eld_init` below:

.. macro:: acb_theta_eld_dim(E)

    Macro returning `d`.

.. macro:: acb_theta_eld_ambient_dim(E)

    Macro returning `g`.

The following macros are available after `E` has been initialized and then
computed using :func:`acb_theta_eld_fill` below:

.. macro:: acb_theta_eld_coord(E, k)

    Macro returning the common coordinate `n_k` of the points in *E*. This
    requires `d \leq k < g`.

.. function:: acb_theta_eld_min(E)

.. function:: acb_theta_eld_mid(E)

.. function:: acb_theta_eld_max(E)

    Macros returning the minimum, midpoint, and maximum of `n_{d-1}` in *E*
    respectively.

.. function:: acb_theta_eld_nr(E)

.. function:: acb_theta_eld_nl(E)

    Macros returning the number of right and left children of *E*
    respectively.

.. function:: acb_theta_eld_rchild(E, k)

.. function:: acb_theta_eld_lchild(E, k)

    Macros returning a pointer to the `k^{\text{th}}` right (resp. left) child
    of *E* as an :type:`acb_theta_eld_t`.

.. function:: acb_theta_eld_nb_pts(E)

    Macro returning the number of points contained in *E*.

.. function:: acb_theta_eld_nb_border(E)

    Macro returning the number of points in the border of *E*, defined as
    follows. If `d=1`, then it consists of the two points with `n_0` equal to
    :func:`acb_theta_eld_min(E)` - 1 and :func:`acb_theta_eld_max(E)` + 1
    respectively. If `d\geq 2`, then it is the reunion of the borders of all
    children of *E*. This is only used for testing.

.. function:: acb_theta_eld_box(E, k)

    Macro returning the smallest nonnegative integer `M_k` such that all the points
    in *E* satisfy `|n_k|\leq M_k`. This requires `0\leq k < d`.


Ellipsoids: memory management and computations
-------------------------------------------------------------------------------

.. function:: void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)

    Initializes *E* as a *d*-dimensional ellipsoid in ambient dimension *g*.

.. function:: void acb_theta_eld_clear(acb_theta_eld_t E)

    Clears *E* as well as any recursive data contained in it.

.. function:: void acb_theta_eld_interval(slong* min, slong* mid, slong* max, const arb_t ctr, const arf_t rad)

    Computes the minimum, midpoint, and maximum of a subinterval of
    `\mathbb{Z}` that is guaranteed to contain all points within a distance
    *rad* of the real number *ctr*. Both *ctr* and *rad* must be finite and the
    result must fit in :type:`slong`'s, otherwise an error is thrown.

.. function:: void acb_theta_eld_cho(arb_mat_t C, const acb_mat_t tau, slong prec)

Computes an upper-triangular Cholesky matrix *C} for the symmetric matrix
`\pi \mathrm{Im}(\tau)`. If one cannot determine that `\mathrm{Im}(\tau)` is
positive definite at the current working precision, *C} is set to an
indeterminate matrix.

\T check that `C^TC = \pi \mathrm{Im}(\tau)` on random input.

.. function:: void acb_theta_eld_fill(acb_theta_eld_t E, const arb_mat_t C, const arf_t R2,
  arb_srcptr v)

Sets *E} to represent an ellipsoid as defined above, where *R2}
indicates `R^2`. The matrix *C} must be an upper-triangular matrix with
positive diagonal entries, *R2} must be finite, and the coordinate of
ellipsoid points must fit in .. type:: slong}'s, otherwise an error is thrown.

\T see \func{acb_theta_eld_points}.

\subsubsection{Points in ellipsoids}

The following functions are available after \func{acb_theta_eld_fill} has been called.

.. function:: void acb_theta_eld_points(slong* pts, const acb_theta_eld_t E)

Sets *pts} to the list of all the points in *E}, as a
concatenation of vectors of length *g}.

\T generate a random ellipsoid *E}. Check that all the points of
*E} are inside the box. Then, generate random points: points inside the
ellipsoid according to \func{acb_theta_eld_contains} below must appear in the
list of points, and the norm of any point outside *E} must be at least
the radius of *E}.

.. function:: void acb_theta_eld_border(slong* pts, const acb_theta_eld_t E)

Sets *pts} to the list of all the points in the border of *E}.

\T check that the border points are not contained in the ellipsoid.

.. function:: int acb_theta_eld_contains(const acb_theta_eld_t E, slong* pt)

Returns true iff *pt} is contained in *E}. The vector *pt}
must be of length *g}.

\T see \func{acb_theta_eld_points} and \func{acb_theta_eld_border} above.

.. function:: void acb_theta_eld_print(const acb_theta_eld_t E)

Prints *E} to .. type:: stdout}. This describes *E} faithfully but may
be unwieldy in high dimensions.

\subsection{Precomputations in naive algorithms}

When running naive algorithms on an ellipsoid~`E` for a certain matrix
`\tau\in \mathbb{H}_g` and points `z^{(0),\ldots, z^{(n-1)\in \mathbb{C}^g`, we
precompute the following quantities:
\begin{itemize}
\item `\exp(i\pi (2 - \delta_{j,k})\tau_{j,k})` for `0\leq j\leq k < g`,
\item `\exp(i\pi j^2 \tau_{k,k})` for `0\leq k < g` and `j` between 0 and
  \func{acb_theta_eld_box(E,k),
\item `\exp(2 i\pi z^{(k)_j)` for `0\leq j < g` and `1\leq k\leq n`.
\end{itemize}
These complex numbers are stored in a structure of type
\func{acb_theta_precomp_t}. Considering several vectors `z` at the same time is
meant to accelerate the computation of `\theta_{a,b}(z,\tau)` for many values
of `z` and a fixed~`\tau`.

\subsubsection{Types and macros}

.. function:: acb_theta_precomp_struct}

.. function:: acb_theta_precomp_t}

An .. type:: acb_theta_precomp_t} is an array of length one of type
.. function:: acb_theta_precomp_struct} containing the above data, permitting it to be
passed by reference.

The following macros are available after calling \func{acb_theta_precomp_init}
and \func{acb_theta_precomp_set} below.

.. function:: acb_theta_precomp_dim(D)

Macro returning the ambient dimension `g`.

.. function:: acb_theta_precomp_nb(D)

Macro returning the number of vectors `z` stored in *D}.

.. function:: acb_theta_precomp_exp_mat(D)

Macro returning a pointer to an .. type:: acb_mat_t} whose entry `(j,k)` contains
`\exp(i \pi (2 - \delta_{j,k}) \tau_{j,k})` for every `0\leq j \leq k\leq g`.

.. function:: acb_theta_precomp_sqr_pow(D, k, j)

Macro returning a pointer to the complex number `\exp(i \pi j^2 \tau_{k,k})` as
an .. type:: acb_t}.

.. function:: acb_theta_precomp_exp_z(D, k, j)

Macro returning a pointer to the complex number `\exp(2\pi i z_k^{(j))` as an
.. type:: acb_t}.

\subsubsection{Memory management and basic manipulation}

.. function:: void acb_theta_precomp_init(acb_theta_precomp_t D, slong nb, slong g)

Initializes *D} for precomputations on *nb} vectors `z\in \mathbb{C}^g`.

.. function:: void acb_theta_precomp_clear(acb_theta_precomp_t D)

Clears *D}.

.. function:: void acb_theta_precomp_set(acb_theta_precomp_t D, acb_srcptr zs,
  const acb_mat_t tau, const acb_theta_eld_t E, slong prec)

Computes the above data for the provided matrix *tau}, vectors `zs` (a
concatenation of *nb} vectors of length `g`) and ellipsoid *E}. The
dimensions must match, in particular *E} must be an ellipsoid of
dimension `g`.

\T check that all entries are set to one on the phony input `z=0, \tau=0`.

\subsection{Naive algorithms}

Naive algorithms consist in summing terms of the theta series over a certain
ellipsoid and adding an error bound coming from the tail of the series. We will
compute the relevant ellipsoid using low-precision computations. When several
vectors `z` are present, we first reduce them to a common compact domain and
use only one ellipsoid, following~\cite{deconinck}. When `g = 1`, we call
functions from \myref{acb_modular.h}{acb_modular} instead.

\subsubsection{Ellipsoids and bounds}

By \cite{main}, for any `v\in \mathbb{R}^g` and any upper-triangular Cholesky
matrix `C`, and any `R` such that `R^2 \geq\max\{4,\mathit{ord}\}`, the sum
\[
  S = \sum_{n\in C\mathbb{Z}^g + v,\ \lVert n\rVert^2 \geq R^2} \lVert n\rVert^{\mathit{ord}} e^{-\lVert n\rVert^2}
\]
satisfies
\[
  S \leq 2^{2g+2} R^{g-1+p} e^{-R^2} \prod_{j=0}^{g-1} (1 + \gamma_j^{-1})
\]
where `\gamma_0,\ldots, \gamma_{g-1}` are the diagonal coefficients of~`C`.

.. function:: void acb_theta_naive_radius(arf_t R2, arf_t eps, const arb_mat_t C, slong ord, slong prec)

Sets *R2} and *eps} such that the above upper bound for *R2}
and the given *ord} is at most *eps}. We choose *eps} so that
the relative error on the output of the naive algorithm should be roughly
`2^{-\mathit{prec}}` if no cancellations occur in the sum, i.e.
`\mathit{eps} \simeq 2^{-\mathit{prec}} \prod_{j=0}^{g-1} (1 + \gamma_j^{-1})`.

\T evaluate the above upper bound on the tail for the computed *R2} on a
random Cholesky matrix *C} and check that it is not greater than *eps}.

.. function:: void acb_theta_naive_reduce(arb_ptr v, acb_ptr new_zs, acb_ptr cs,
  arb_ptr us, acb_srcptr zs, slong nb, const acb_mat_t tau, const arb_mat_t C,
  slong prec)

Performs the simultaneous reductions of the *nb} vectors stored in `zs`
with respect to the matrix `\tau`. This means the following. Let
`0\leq k\leq \mathit{nb}-1`, let `z` denote the `k^{\mathrm{th}}` vector stored
in `zs`, and let `X,Y` (resp. `x,y`) be the real and imaginary parts of `\tau`
(resp. `z`). Write `Y^{-1}y = r + a` where `a` is an even integral vector
and~`r` is bounded. Then
\[
  \begin{aligned}
  \theta_{0,b}(z,\tau) &= e^{\pi y^T Y^{-1} y} \sum_{n\in \mathbb{Z}^g}
                         e^{\pi i ((n - a)^T X (n - a) + 2(n - a)^T (x + \tfrac b2)) e^{-\pi (n + r)^T Y (n + r)\\
    &= e^{\pi y^T Y^{-1} y} e^{\pi i (a^T X a - 2a^T x + i r^T Y r) \theta_{0,b}((x - Xa) + iYr, \tau).
  \end{aligned}
\]
The reduction of `z` is defined as `(x - Xa) + i Y r`, which has a bounded
imaginary part, and this vector is stored as the `k^{\mathrm{th}}` vector of
*new_zs}. The quantity `u = \exp(\pi y^T Y^{-1} y)` is a multiplicative
factor for the error bound, and is stored as the `k^{\mathrm{th}}` entry of
*us}. the quantity `c = u \exp(\pi i (a^T X a - 2a^T x + i r^T Y r))` is
a multiplicative factor for the theta values, and is stored as the
`k^{\mathrm{th}}` entry of *cs}. The offset for the corresponding
ellipsoid is `v^{(k) = C r` which is also bounded independently of~`k`, and
the vector *v} is set to the \func{acb_union} of these vectors `v^{(k)`
for `0\leq k\leq \mathit{nb}-1`.

\T check that the results are sound or some special values, as follows. If
*zs} are real vectors, then *new_zs} must be equal to *zs},
the entries of *cs}, *us} must all be `1`, and *v} must be
zero. If `\mathrm{Im}(z) = -Yn + \varepsilon` where `n` is an even integral
vector and `\varepsilon` is small, then the result of
\func{acb_theta_naive_term} on `n` for `z` must overlap `c` times the term for
*new_z} attached to the lattice point `0`, and the computed offset `v`
must be small.

.. function:: void acb_theta_naive_ellipsoid(acb_theta_eld_t E, acb_ptr new_zs, acb_ptr cs,
  arb_ptr us, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

Sets the ellipsoid *E} and the vectors *new_zs}, *cs} and
*us} such that the following is satisfied: for each
`0\leq k\leq \mathit{nb}-1`, if `z` and `z'` denote the `k^{\mathrm{th}}`
vectors in `zs` and *new_zs} respectively, and `c, u` denote the
`k^{\mathrm{th}}` entries in *cs} and *us}, then summing
exponential terms involving `z'` over *E} and multiplying by `c` will
yield an approximation of theta values at `z` up to an error at most
`u`. Unless cancellations occur in the sum, the relative precision of the
resulting theta values should be roughly *prec}.

\T check that the sum of terms on the border of the ellipsoid *E} is at
most *u}.

.. function:: slong acb_theta_naive_fullprec(const acb_theta_eld_t E, slong prec)

Returns a good choice of full precision for the summation phase when working at
precision *prec}, which is at least `\mathit{prec} + \log_2(n)` where `n`
is the number of points in *E}.

\T no test, but used throughout in naive algorithms.

.. function:: acb_theta_naive_term(acb_t res, acb_srcptr z, const acb_mat_t tau,
  slong* tup, slong* n, slong prec)

Sets *res} to
`n_0^{k_0} \cdots n_{g-1}^{k_{g-1}}\exp(i\pi(n^T\tau n + 2 n^Tz))`, where the
`k_j` and `n_j` denotes the `j^{\mathrm{th}}` entry in *tup} and
*n} respectively. *tup} may be \func{NULL}, which is understood to
mean the zero tuple. This is used for testing and in
\func{acb_theta_naive_00} for very skewed ellipsoids.

\T if `g=1`, this should simply be `n^k\exp(i\pi (n^2\tau + 2nz))`.

\subsubsection{Workers in naive algorithms}

The main worker inside each version of the naive algorithm will process one
line inside the computed ellipsoid. Before calling this worker, for fixed
`\tau` and `z` and fixed coordinates `n_1,\ldots n_{g-1}` defining a line
inside the ellipsoid, if `n_{\mathrm{min}}` are `n_{\mathrm{max}}` are the
endpoints of the interval of allowed values for `n_0`, we (efficiently)
precompute:
\begin{itemize}
\item The vector `v_1` with entries `\exp(i \pi j^2 \tau_{0,0})` for
  `n_{\mathrm{min}}\leq j\leq n_{\mathrm{max}}`,
\item The vector `v_2` with entries `x^j` for `n_{\mathrm{min}}\leq j\leq n_{\mathrm{max}}` where
  \[
    x = \exp(2 \pi i z_0) \prod_{k = 1}^{g-1} \exp(2 i \pi n_k \tau_{0,k}),
  \]
\item The cofactor `c\in \mathbb{C}` given by
  \[
    c = \prod_{k = 1}^{g-1} \exp(2 i\pi n_k z_k) \cdot \prod_{1\leq j\leq k < g} \exp(\pi i (2 - \delta_{j,k}) n_j n_k \tau_{j,k}).
  \]
\end{itemize}
This allow us to use \func{acb_dot} in the workers while maintaining reasonable
memory costs, and to use an average of strictly less than two complex
multiplications per lattice point as `R\to \infty`. Moreover, these
multiplications are performed at only a fraction of the full precision for
lattice points far from the ellipsoid center.


.. function:: acb_theta_naive_worker_t}

A function pointer type. A function *worker} of this type has the
following signature:

.. function:: void worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2, const slong* precs, slong len,
  const acb_t c, const slong* coords, slong ord, slong g, slong prec, slong fullprec)

where
\begin{itemize}
\item *th} denotes the output vector of theta values to which terms will
  be added,
\item *v1}, *v2} and *c} are precomputed as above,
\item *precs} is a vector of working precisions for each term
  `n_{\mathrm{min}}\leq j\leq n_{\mathrm{max}}`,
\item *len} `= n_{\mathrm{max}} - n_{\mathrm{min}} + 1` is the common
  length of *v1}, *v2} and *precs},
\item *coords} is `(n_{\mathrm{min}}, n_1, \ldots, n_{g-1})`,
\item *ord} is the maximal derivation order (0 outside \func{acb_theta_jet_naive}
  functions),
\item *prec} is the working precision for this line inside the ellipsoid,
  and finally
\item *fullprec} is the working precision for summing into *th}.
\end{itemize}

.. function:: void acb_theta_naive_worker(acb_ptr th, slong len, const acb_t c, const arb_t u,
  const acb_theta_eld_t E, const acb_theta_precomp_t D, slong k, slong ord,
  slong prec, acb_theta_naive_worker_t worker)

Runs the naive algorithm on the ellipsoid *E} using precomputed data for
the `k^{\mathrm{th}}` vector stored in *D}. Here `c` and `u` are as
output by \func{acb_theta_naive_ellipsoid}, *ord} is passed as an
argument to \func{worker}, *prec} is the precision for summing into
the vector *th}, and *len} is the length of the output vector
*th}.

\T no test, but used throughout in naive algorithms.

\subsubsection{Main user functions}

.. function:: void acb_theta_naive_00(acb_ptr th, acb_srcptr zs, slong nb,
  const acb_mat_t tau, slong prec)

.. function:: void acb_theta_naive_0b(acb_ptr th, acb_srcptr zs, slong nb,
  const acb_mat_t tau, slong prec)

Evaluates either `\theta_{0,0}(z^{(k), \tau)`, or alternatively
`\theta_{0,b}(z^{(k), \tau)` for each `b\in \{0,1\}^g`, for each
`1\leq k \leq \mathit{nb}`. The associated worker performs one \func{acb_dot}
operation. The result *th} will be a concatenation of *nb} vectors of length
`1` or `2^g` respectively.

\T check that the result of \func{naive_00} overlaps the first entry of the
result of \func{naive_0b} on random input.

.. function:: void acb_theta_naive_fixed_ab(acb_ptr th, ulong ab, acb_srcptr zs, slong nb,
  const acb_mat_t tau, slong prec)

.. function:: void acb_theta_naive_fixed_a(acb_ptr th, ulong a, acb_srcptr zs, slong nb,
  const acb_mat_t tau, slong prec)

.. function:: void acb_theta_naive_all(acb_ptr th, acb_srcptr zs, slong nb,
  const acb_mat_t tau, slong prec)

Evaluates `\theta_{a,b}(z^{(k), \tau)` for, respectively: the given value of
`(a,b)`; all `(a,b)` for `b\in \{0,1\}^g` and the given value of~`a`; or all
`(a,b)\in\{0,1\}^{2g}`, for each `1\leq k\leq \mathit{nb}`. The result
*th} will be a concatenation of *nb} vectors of length `1`, `2^g`
or `2^{2g}` respectively. We reduce to calling \func{acb_theta_naive_00} or
\func{acb_theta_naive_0b} by writing
\[
\theta_{a,b}(z,\tau) = \exp(\pi i \tfrac{a^T}{2} \tau \tfrac a2) \exp(\pi i
a^T(z + \tfrac b 2)) \theta_{0,b}(z + \tau \tfrac{a}{2}, \tau).
\]

\T check that the results of \func{naive_fixed_ab} and \func{naive_fixed_a}
agree with the corresponding entries of the result of \func{naive_all}. Also
check that \func{naive_all} agrees with \func{acb_modular_theta} as follows: if
`\tau` is diagonal with coefficients `\tau_0,\ldots,\tau_{g-1}`, then we have
`\theta_{a,b}(z,\tau) = \prod_{j=0}^{g-1} \theta_{a_j,b_j}(z_j,\tau_j)`, so
these quantities must overlap.

\subsection{Naive algorithms for derivatives}

We only consider the successive partial derivatives of `\theta_{a,b}(z,\tau)`
with respect to the~`g` coordinates of~`z`, since derivatives with respect
to~`\tau` are accounted for by the heat equation
\[
  \frac{\partial\theta_{a,b}}{\partial \tau_{j,k}} = \frac{1}{2\pi i(1
    +\delta_{j,k}) \frac{\partial^2\theta_{a,b}}{\partial z_j \partial z_k}.
\]
We encode tuples of derivation orders, henceforth called ``derivation tuples'',
as vectors of type .. type:: slong} and length~`g`. In agreement with
\myref{acb_modular}{acb_modular}, we also normalize derivatives in the same way
as in the Taylor expansion, so that the tuple `(k_0,\ldots,k_{g-1})`
corresponds to the differential operator
\[
  \frac{1}{k_0!}\cdots\frac{1}{k_{g-1}!} \cdot
  \frac{\partial^{|k|}}{\partial z_0^{k_0}\cdots \partial z_{g-1}^{k_{g-1}}},
\]
where `|k|:=\sum k_i`.  We always consider all derivation tuples up to a total
order *ord}, and order them first by their total order, then
reverse-lexicographically. For example, in the case `g=2`, the sequence of
orders is `(0,0), (1,0), (0,1), (2,0), (1,1)`, etc.

The naive algorithms for derivatives will evaluate a partial sum of the
differentiated series:
\[
  \frac{\partial^{|k|}\theta_{a,b}}{\partial z_0^{k_0}\cdots \partial z_{g-1}^{k_{g-1}}}(z,\tau)
  = (2\pi i)^{|k|} \sum_{n\in \mathbb{Z}^g + \tfrac a2} n_0^{k_0} \cdots n_{g-1}^{k_{g-1}}
  \exp(\pi i n^T \tau n + 2\pi i n^T (z + \tfrac b2)).
\]

.. function:: slong acb_theta_jet_nb(slong ord, slong g)

Returns the number of derivation tuples with total order at most *ord}.

\T no test, but used e.g. in \func{acb_theta_jet_index}.

.. function:: slong acb_theta_jet_total_order(const slong* tup, slong g)

Returns the total derivation order for the given tuple *tup} of length *g}.

\T no test, but used e.g. in \func{acb_theta_jet_index}.

.. function:: void acb_theta_jet_tuples(slong* tups, slong ord, slong g)

Sets *tups} to the concatenation of all derivation tuples up to total
order *ord}.

\T generate the list of tuples, pick an index `i` at random, at check that the
result of \func{acb_theta_jet_index} below on the `i^{\mathrm{th}}` tuple is
indeed `i`.

.. function:: slong acb_theta_jet_index(const slong* tup, slong g)

Returns *n} such that *tup} is the `n^{\mathrm{th}}` derivation tuple of length *g}.

\T see \func{acb_theta_jet_tuples}.

.. function:: void acb_theta_jet_naive_radius(arf_t R2, arf_t eps, arb_srcptr v,
  const arb_mat_t C, slong ord, slong prec)

Assuming that *C} is the upper-triangular Cholesky matrix for
`\pi \mathrm{Im}(\tau)` and `v = C Y^{-1} y` where~`y, Y` are the imaginary
parts of `z` and `\tau` respectively, returns *R2} and *eps} so
that, when summing the above series on terms `n\in \mathbb{Z}^g` such that
`(v + C n)^T(v + C n)\leq \mathit{R2}`, the absolute value of the tail of the
series (ignoring leading multiplicative terms, see below) will be bounded above
by *eps}, for any derivation tuple `k` with `|k|\leq \mathit{ord}`.

We can rewrite the above sum as
\[
  \frac{\partial^{|k|}\theta_{a,b}}{\partial z_0^{k_0}\cdots \partial
    z_{g-1}^{k_{g-1}}}(z,\tau) = (2\pi i)^{|k|} e^{\pi y^T Y^{-1} y} \sum_{n\in
    \mathbb{Z}^g + \tfrac a2} n_0^{k_0} \cdots n_{g-1}^{k_{g-1}} e^{\pi
    i(\cdots) e^{-\pi (n + Y^{-1}y)^T Y (n + Y^{-1}y).
\]
Ignore the multiplicative factors in front of the sum. Writing
`m = C n + v`, we have
\[
  n_0^{k_0}\cdots n_{g-1}^{k_{g-1}}\leq
  (\lVert C^{-1}\rVert \lVert n\rVert + \lVert
  Y^{-1}y\rVert)^{|k|}.
\]
Here all norms are (induced) infinity norms, which for vectors are bounded
above by the `L^2` norm. Therefore, the absolute value of the tail of the
series is bounded above by
\[
  \biggl(\sum_{j=0}^{|k|} \binom{|k|}{j} \lVert C^{-1} \rVert^{j}
  R^j \lVert Y^{-1}y\rVert^{|k|-j}\biggr) \cdot 2^{2g+2} R^{g-1} e^{-R^2}
  \prod_{j=0}^{g-1} (1 + \gamma_j^{-1}).
\]
The inner sum is simply
`(\lVert C^{-1} \rVert R + \lVert Y^{-1}y \rVert)^{|k|}`.

Thus, we proceed as follows. We first compute *R2} and *eps} using
\func{acb_theta_naive_radius} with *ord = 0}. If
`R\leq \lVert Y^{-1}y\rVert/\lVert C^{-1}\rVert`, we simply multiply *eps} by
`\max\{1, 2 \lVert Y^{-1}y \rVert\}^{\mathit{ord}}`. Otherwise, we compute
*R2} and *eps} using \func{acb_theta_naive_radius} with the given
value of *ord}. We can then set *R2} to the maximum of *R2}
and `\lVert Y^{-1}y \rVert /\lVert C^{-1} \rVert`, and multiply *eps} by
`\max\{1, 2\lVert C^{-1}\rVert\}^{\mathit{ord}}`.

\T generate `C` and `v = Cy^{-1}y` randomly, compute *R2} and
*eps}, then check that
`\max\{1, \lVert C^{-1}\rVert R + \lVert Y^{-1}y\rVert\}^{\textit{ord}}\cdot
2^{2g+2}R^{g-1}e^{-R^2}\prod_{j=0}^{g-1} (1+\gamma_j^{-1})\leq \mathit{eps}`.

.. function:: void acb_theta_jet_ellipsoid(acb_theta_eld_t E, arb_t u, acb_srcptr z,
  const acb_mat_t tau, slong ord, slong prec)

Sets *E} and *u} so that summing over *E}
yields derivatives of theta functions up to an error of at most *u},
ignoring factorials and powers of `\pi`. After computing *R2} and
*eps} as in \func{acb_theta_jet_naive_radius}, we set the radius of
*E} to be *R2} and set
`u = e^{\pi y^T Y^{-1} y}\cdot \mathit{eps}`.

\T generate random *z} and *tau} and check that the sum of the
absolute values of terms on the border of the ellipsoid *E} is at most
*u}.

.. function:: void acb_theta_jet_naive_00(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
  slong ord, slong prec)

Sets *dth} to the vector of derivatives of `\theta_{0,0}` at the given
point `(z,\tau)` up to total order *ord}.

\T check that the values overlap with the result of \func{acb_theta_jet_naive_all} below.

.. function:: void acb_theta_jet_naive_fixed_ab(acb_ptr dth, ulong ab, acb_srcptr z, const acb_mat_t tau,
  slong ord, slong prec)

Sets *dth} to the vector of derivatives of `\theta_{a,b}` at the given
point `(z,\tau)` up to total order *ord}. We reduce to
\func{acb_theta_jet_naive_00} using the same formula as in
\func{acb_theta_naive_ind}, making suitable linear combinations of the
derivatives.

\T check that the values overlap with the result of \func{acb_theta_jet_naive_all} below.

.. function:: void acb_theta_jet_naive_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
  slong ord, slong prec)

Sets *dth} to the vector of derivatives of all the functions
`\theta_{a,b}` for `a,b\in \{0,1\}^g` up to total order *ord} at the
given point. The result will be a concatenation of `2^{2g}` vectors of length
\func{acb_theta_jet_nb(ord, g).

We use an ellipsoid to encode points in `\tfrac 12 \mathbb{Z}^g`, and divide
`\tau` by 4 and `z` by 2 to sum the correct terms. The bounds output by
\func{acb_theta_jet_naive_radius} are still valid, since this just has the
effect of multiplying `\lVert C^{-1} \rVert` and each
`\gamma_j^{-1}` by `2`.

\T check that for diagonal matrices, the results agree with
\func{acb_modular_theta_jet}.

.. function:: void acb_theta_jet_error_bounds(arb_ptr err, acb_srcptr z, const acb_mat_t tau,
  acb_srcptr dth, slong ord, slong prec)

Assuming that *dth} contains the derivatives of a function `\theta_{a,b}`
up to total order `\mathit{ord} + 2`, sets *err} to a vector with the
following property. Let `(z_0,\tau_0)` be the midpoint of `(z,\tau)`, and let
`(z_1,\tau_1)` be any point inside the ball specified by the given *z}
and *tau}. Then the vectors of derivatives of `\theta_{a,b}` at
`(z_0,\tau_0)` and `(z_1,\tau_1)` up to total order *ord} differ by at
most *err} elementwise.

\T generate two pairs `(z_1,\tau_1)` and `(z_2,\tau_2)` close to each other but
not overlapping, set `(z,\tau)` to be their reunion, and compute *err} on
`(z,\tau)`. The difference between the result of \func{acb_theta_jet_naive_all}
on `(z_1,\tau_1)` and `(z_2,\tau_2)` must be at most two times *err}.

\subsection{Quasi-linear algorithms on the reduced domain}

\subsubsection{A simple case: theta constants for `g=1`}

In this section, we present the quasi-linear algorithm in the simple case
`g=1` and `z=0`, with the hope that it will make the general case easier
to follow.

The algorithm is based on a \emph{duplication formula}: such formulas typically
related theta values at `\tau` and `2\tau`, and look like taking a step in
an arithmetic-geometric (AGM) sequence. For instance, one has for any `g`:
\begin{displaymath}
  \theta_{0,b}(0,2\tau)^2 = 2^{-g} \sum_{b'\in (\mathbb{Z}/2\mathbb{Z})^g}
  \theta_{0,b'}(0,\tau) \theta_{0, b+b'}(0,\tau).
\end{displaymath}
When `g=1`, this becomes
\begin{displaymath}
  \theta_0(0,2\tau)^2 = \tfrac12(\theta_0(0,\tau)^2 + \theta_{1}(0,\tau)^2),
  \quad \theta_1(0,2\tau)^2 = \theta_0(0,\tau)\theta_1(0,\tau).
\end{displaymath}
This is the formula that is used (in a convoluted way, using limits of
AGM sequences and a Newton scheme) to obtain quasi-linear algorithms for theta
values in low genera in \cite{dupont,labrande,kieffer}.

Here we describe a new algorithm which consists of using the duplication
formula without any Newton scheme. We just use duplications to transform `\tau`
into a point of `\mathbb{H}_g` where theta values are easier to evaluate using
the naive algorithm. Indeed, if `\lambda` denotes the smallest eigenvalue of
`\mathrm{Im}(\tau)`, then the number of lattice points to consider in the naive
algorithm to obtain theta values at absolute precision *prec} is
`O((*prec}/\lambda)^{g/2})`. Thus replacing `\tau` by `2\tau` divides
this number by `2^{g/2}`. This is already a huge gain, and the quasi-linear
algorithm will perform roughly `\log_2(*prec})` such duplication steps.

Clearly, to apply this strategy, the previous formula expressing theta values
at `2\tau` in terms of values at `\tau` goes in the wrong direction. Instead,
we use a formula relating the values `\theta_{a,0}` for
`a\in (\mathbb{Z}/2\mathbb{Z})^g`:
\begin{displaymath}
  \theta_{a,0}(0,\tau)^2 = \sum_{a'\in (\mathbb{Z}/2\mathbb{Z})^g}
  \theta_{a',0}(0,2\tau)\theta_{a+a',0}(0,2\tau).
\end{displaymath}
When `g=1`, this becomes
\begin{displaymath}
  \theta_0(0,\tau)^2 = \theta_0(0,2\tau)^2 + \theta_2(0,2\tau)^2,\quad
  \theta_2(0,\tau)^2 = 2\theta_0(0,2\tau)\theta_2(0,2\tau).
\end{displaymath}
We will explain later how to obtain all the theta values
`\theta_{a,b}(0,\tau)`, not just `\theta_{a,0}`. For now, we focus on how to
obtain a quasi-linear algorithm for computing `\theta_{a,0}` from this formula
in the case `g=1`. We can assume that `\tau` belongs to the usual fundamental
domain for the action of `\mathrm{SL}_2(\mathbb{Z})` thanks to the
transformation formula.

At each step, we will have to extract square roots to compute
`\theta_{a,0}(0,\tau)` from `\theta_{a,0}(0,\tau)^2`. We have to worry about 1)
the choice of sign, and 2) precision losses, since `\theta_{a,0}(0,\tau)` tends
to zero rapidly as the imaginary part of `\tau` gets large when `a\neq 0`. For
simplicity, we now restrict to `g=1`. We examine the series in terms of
`q = \exp(\pi i\tau)`:
\begin{displaymath}
  \theta_0(0,\tau) = 1 + 2q + 2q^4 + \cdots, \quad \theta_2(0,\tau) = 2q^{1/4} + 2q^{9/4} + \cdots
\end{displaymath}
Thus `\theta_{0}(0,2^n\tau)` tends to `1` as `n\to \infty`: this means that
taking square roots at each step loses just `O(1)` bits of absolute precision,
and that it is easy to determine the correct choice of sign, since we just have
to run the naive algorithm at precision `O(1)` at each step. (Indeed, we could
just say that the correct square root is the one closest to `1`, but we won't
have such a shortcut for general `g`).

What about `\theta_2`? If `y` denotes the imaginary part of `\tau`, then one
can expect (and it's easy to show in the `g=1` case) that
`|\theta_{2}(0,\tau)|` is roughly `2\exp(-\pi y/4)`. Taking a square root will
lead to a large precision loss in terms of absolute precision, so it is better
to think in terms of relative precision, or rather some kind of ``shifted
absolute precision'': if we compute `\theta_2(0,\tau)^2` up to an absolute
error of `\exp(-\pi (2y)/4) 2^{-\mathit{prec}}`, then we get `\theta_2(0,\tau)`
up to an absolute error of `\exp(-\pi y/4) 2^{-\mathit{prec}}`, with `O(1)`
bits of precision loss.

How do these shifted absolute precisions interact in the duplication formula?
For
\begin{displaymath}
  \theta_0(0,\tau)^2 = \theta_0(0,2\tau)^2 + \theta_2(0,2\tau)^2,
\end{displaymath}
we're fine: we know `\theta_2(0,2\tau)` to a larger precision than
necessary. For
\begin{displaymath}
  \theta_2(0,\tau)^2 = 2\theta_0(0,2\tau)\theta_2(0,2\tau),
\end{displaymath}
we're also fine: we get `2\theta_0(0,2\tau)\theta_2(0,2\tau)` up to an error of
`\exp(-\pi (2y)/4)2^{-\mathit{prec}}`, which is exactly what we want, with
`O(1)` bits of precision lost.  We note in passing that we can also determine
the correct choice of square root of each step for `\theta_2(0,\tau)` using the
naive algorithm on `O(1)` lattice points.

In the quasi-linear algorithm, we perform `n\simeq \log_2(\mathit{prec})` such
duplication steps. We initialize at `2^n\tau` using the naive algorithm. We
need `\theta_0(0,2^n\tau)` to an absolute precision *prec} (plus maybe a
logarithmic number of guard bits to account for precision losses at each step):
for this we need `O(1)` lattice points in the naive algorithm. We also need
`\theta_2(0,2^n\tau)` to an absolute precision
`\mathit{prec} + \frac{\pi}{\log(2)2^{n-2}`, which is still
`O(\mathit{prec})`. So we also need only `O(1)` lattice points to run the naive
algorithm and compute `\theta_2(0,2^n\tau)` to the right precision.

To conclude the `g=1` case, we explain how to recover all the theta values, not
just `\theta_{a,0}`. In fact, we have
\begin{displaymath}
  \theta_{a,b}(0,\tau)^2 = \sum_{a'\in (\mathbb{Z}/2\mathbb{Z})^g} (-1)^{a'^Tb}
  \theta_{a',0}(0,2\tau)\theta_{a+a',0}(0,2\tau).
\end{displaymath}
Thus, if we want the squared theta values `\theta_{a,b}`, it is enough to
compute `\theta_{a,0}` at `2\tau` for all `a`. If we want the actual values, we
add one last square-root step.

Finally, here are some indications of how the above strategy will
generalize to any `g`.
\begin{enumerate}
\item The concept of shifted absolute precision is important. We expect that
  `|\theta_{a,0}(0,\tau)|` is roughly `e^{-d^2}`, where `d` denotes the
  distance between `0` and `\mathbb{Z}^g + \tfrac a2` for the distance attached
  to the quadratic form `\mathrm{Im}(\tau)`. (There is a similar formula when
  `z\neq 0`.) These distances are computed using the \func{acb_theta_dist...}
  functions.
\item To avoid making `2^{2g}` multiplications in the duplication formula, we
  use the Hadamard matrix: then a duplication steps costs only `2^g`
  multiplications and square roots. In the `g=1` case, we would compute
  `x = (\theta_0 + \theta_2)^2` and `y = (\theta_0 - \theta_2)^2`, then write
  for instance `2\theta_0\theta_2 = \frac12(x - y)`. However, this would bring
  huge precision losses in terms of shifted absolute precisions! So we must
  compute the Hadamard products at a significantly higher precision and adjust
  the error bounds in the end. See \func{acb_theta_agm_mul_tight}.
\item When `g\geq 3`, some theta values `\theta_{a,0}(0,2^k\tau)` encountered
  in the algorithm may very well be much smaller than expected (in terms of
  lattice distances) or vanish altogether (This also happens for `g=2` and
  nonzero `z`, maybe even `g=1`). This is problematic for two reasons: we will
  lose precision when taking square roots, and perhaps more importantly, it
  won't be cheap anymore to compute the correct choice of sign with the naive
  algorithm. Luckily, we can circumvent this by introducing a random auxiliary
  real vector `t` and considering `\theta_{a,0}(2^kt, 2^k\tau)`: these will be
  large enough with overwhelming probability, and one can adapt the duplication
  formula to output `\theta_{a,0}(0,\tau)`. See \func{acb_theta_ql_roots},
  \func{acb_theta_ql_step_1} and \func{step_3}.
\item When `g\geq 2`, it may be the case that `\mathrm{Im}(\tau)` has
  eigenvalues of different orders of magnitude. In this case, the ellipsoids in
  the naive algorithms for `2^k\tau` as `k` grows will become very thin in some
  directions while still being thick in other directions. We can then do a few
  duplication steps and then fall back to computing theta values in smaller
  dimensions: this is implented in \func{acb_theta_ql_a0_split}.
\item The transformation formula has an analogue for any `g` using the action
  of `\mathrm{Sp}_{2g}(\mathbb{Z})`: see \func{acb_theta_transform_...} This is
  important because in order to determine the correct choice of square root at
  `\tau` using the naive algorithm, we want `\tau` to be reduced.
\end{enumerate}



\subsubsection{Distances}

.. function:: void acb_theta_dist_pt(arb_t d, arb_srcptr v, const arb_mat_t C, slong* n, slong prec)

Sets *d} to `\lVert v - Cn\rVert^2` for the
Euclidean norm.

\T check that the results for `v = Cn_1`, `n = n_2` and `v = Cn_2`, `n = n_1` overlap.

.. function:: void acb_theta_dist_lat(arb_t d, arb_srcptr v, const arb_mat_t C, slong prec)

Sets *d} to `\mathrm{Dist}(v, C \mathbb{Z}^g)^2` for the
Euclidean norm. We first compute an upper bound on the result by considering
the `2^g` vectors obtained by rounding the entries of `C^{-1}v` to
integers, up or down, then compute an ellipsoid to find the minimum distance.

\T for a random choice of `C` and `v`, compute the distance, then compute the
corresponding ellipsoid. Check that it has at least one point and that the
distance is correct.

.. function:: void acb_theta_dist_a0(arb_ptr d, acb_srcptr z, const acb_mat_t tau, slong prec)

Sets *d} to the vector containing
`\mathrm{Dist}(C \cdot(Y^{-1}y + \tfrac a2), C\cdot
\mathbb{Z}^g)^2` for `a\in \{0,1\}^g`, where `y, Y` are the imaginary parts of
`z, \tau` respectively and `C` is the upper-triangular Cholesky
matrix for `\pi \mathrm{Im}(\tau)`. The `a^{\mathrm{th}}` entry of *d}
is also `\mathrm{Dist}_\tau(-Y^{-1}y, \mathbb{Z}^g + \tfrac a2)^2`, where
`\mathrm{Dist}_\tau` denotes the distance attached to the quadratic form
`\mathrm{Im}(\tau)`.

\T when the imaginary part of `z` is `Y \tfrac{a}{2}` for some theta
characteristic `a`, check that the `a^{\mathrm{th}}` entry of the result of
\func{acb_theta_dist_a0} contains zero.

.. function:: slong acb_theta_dist_addprec(const arb_t d)

Returns an integer that is close to *d} divided by `\log(2)`. Requires
that *d} is finite and of reasonable size, otherwise an error is
thrown.

\T no test, but used throughout in quasi-linear algorithms.

\subsubsection{Duplication formulas}

.. function:: void acb_theta_agm_hadamard(acb_ptr res, acb_srcptr a, slong g, slong prec)

Sets *res} to the product of the Hadamard matrix
`\left(\begin{smallmatrix} 1 & 1 \\ 1 & -1\end{smallmatrix}\right)^{\otimes g}`
and the vector `a`. Both `r` and `a` must be vectors of length `2^g`. In other
words, for each `k\in \{0,1\}^g`, this sets the `k^{\mathrm{th}}` entry of
*res} to `\sum_{j\in \{0,1\}^g} (-1)^{k^T j} a_j`.

\T check that applying the Hadamard matrix twice is equivalent to mutiplying by `2^g`.

.. function:: void acb_theta_agm_sqrt(acb_ptr res, acb_srcptr a, acb_srcptr rts, slong nb, slong prec)

Sets the `k^{\mathrm{th}}` entry of *res} for `0\leq k < \mathit{nb}` to a square
root of the corresponding entry of `a`. The choice of sign is determined by
*rts}: each entry of `r` will overlap the corresponding entry of
*rts} but not its opposite. The result is indeterminate if both square
roots overlap, and an error is thrown if there is no overlap at all.

\T generate a random vector *t}, set *rts} to a low-precision
rounding of *t} and set *a} to the square of *t}
elementwise. The result of \func{acb_theta_agm_sqrt} must then overlap
*t} and the precision loss must be small (this is just checked on the
first entry).

.. function:: void acb_theta_agm_mul(acb_ptr res, acb_srcptr a1, acb_srcptr a2, slong g, slong prec)

For each `0\leq k < 2^g`, sets the `k^{\mathrm{th}}` entry of *res} to
`2^{-g}\sum_{b\in \{0,1\}^g} a_{1,b}\, a_{2, b + k}`, where addition is meant
in `(\mathbb{Z}/2\mathbb{Z}^g)` (a bitwise xor). Following \cite{labrande}, we
apply the Hadamard matrix twice with multiplications in-between. This causes
precision losses when the absolute values of the entries of *a1} and/or
*a2} are of different orders of magnitude. This function is faster when
*a1} and *a2} are equal as pointers, as we can use squarings
instead of multiplications.

\T check that the duplication formula holds: the result of
\func{acb_theta_agm_mul} on vectors containing `\theta_{0,b}(0,\tau)` and
`\theta_{0,b}(z,\tau)` for `b\in \{0,1\}^g` must contain
`\theta_{0,b}^2(2z,2\tau)`.

.. function:: void acb_theta_agm_rel_mag_err(arf_t m, arf_t eps, acb_srcptr a,
  arb_srcptr d, slong nb, slong prec)

Computes *m} and *eps} such that the following holds: for each
`0\leq k < \mathit{nb}`, if `d_k` (resp. `a_k`) denotes the `k^{\mathrm{th}}` entry of
*d} (resp. *a}), then the absolute value of `a_k` is at most
`m \cdot e^{-d_k}` and the radius of the complex ball `a_k` is at most
`\mathit{eps}\cdot e^{-d_k}`.

\T after choosing random *m}, *eps} and *d}, generate a
random vector *a} whose entries satisfy the corresponding inequalities,
making sure that equality holds for at least one entry of *a}. The result
`(m',\mathit{eps}')` of \func{agm_rel_mag_err} must then satisfy
`m'\geq m` and `\mathit{eps'}\geq \mathit{eps}`.

.. function:: void acb_theta_agm_mul_tight(acb_ptr res, acb_srcptr a0, acb_srcptr a,
  arb_srcptr d0, arb_srcptr d, slong g, slong prec)

Assuming that *d0} and *d} are obtained as the result of
\func{acb_theta_dist_a0} on `(0,\tau)` and `(z,\tau)` respectively, performs
the same computation as \func{acb_theta_agm_mul} on the vectors *a0} and
*a}, but manages the error bounds as follows. Let `m_0, \varepsilon_0`
(resp.~`m,\varepsilon`) be the result of \func{acb_theta_agm_rel_mag_err} on
`a_0,d_0` (resp. `a,d`). We call \func{acb_theta_agm_mul} on the midpoints of
*a0} and *a} at working precision
`\mathit{prec} + {}`\func{acb_theta_dist_addprec(dmax) where *dmax} is
the largest entry of `d`, then add an error bound to the `k^\mathrm{th}` entry
of *res} of the form
`e^{-d_k} (m_0 \varepsilon + m \varepsilon_0 + \varepsilon\varepsilon_0)`. The
resulting precision losses are very mild when `m_0` and `m` are relatively
small. The computation is valid for the following reason: for each
`b\in \{0,1\}^g`, we have (keeping notation from \func{acb_theta_dist_a0})
\[
  \mathrm{Dist}_\tau(-Y^{-1}y, \mathbb{Z}^g + \tfrac b2)^2 +
  \mathrm{Dist}_\tau(-Y^{-1} y, \mathbb{Z}^g + \tfrac{b + k}{2})^2 \leq
  \mathrm{Dist}_\tau(-Y^{-1}y, \mathbb{Z}^g + \tfrac{k}{2})^2
\]
by the parallelogram identity.

\T generate random `\tau` and `z` at precision *prec} and compute the
associated vectors *d0} and *d}. Set each entry of *a0}
(resp. *a}) to be of the form `z e^{-t}` where `z` is uniformly random
with `|z|\leq 1` and `t` is the corresponding entry of *d0}
(resp. *d}). Apply \func{agm_mul_tight}, then apply
\func{agm_rel_mag_err} on the result with respect to *d}. Check that the
resulting *m} satisfies `m \leq 1` and that *eps} is at most
`2^{-*prec} + \delta}` for some reasonable value of `\delta` (e.g.  25).

\subsubsection{AGM steps for `\theta_{a,0}`}

The first step in quasi-linear algorithms is to compute the quantities
`\theta_{a,0}(z,\tau)` for `a\in \{0,1\}^g` by repeated applications of the
duplication formula:
\[
  \theta_{a,0}(z,\tau) \theta_{a,0}(z',\tau) = \sum_{a'\in(\mathbb{Z}/2\mathbb{Z})^g}
  \theta_{a',0}(z+z',2\tau) \theta_{a+a',0}(z-z',2\tau).
\]
In particular,
\[
  \begin{aligned}
    \theta_{a,0}(z,\tau)^2 &= \sum_{a'\in (\mathbb{Z}/2\mathbb{Z})^g}
  \theta_{a',0}(2z,2\tau) \theta_{a+a',0}(0,2\tau),\\
  \theta_{a,0}(0,\tau)\theta_{a,0}(z,\tau) &= \sum_{a'\in(\mathbb{Z}/2\mathbb{Z})^g}
  \theta_{a',0}(z,2\tau) \theta_{a+a',0}(z,2\tau), \\
  \theta_{a,0}(0,\tau)^2 &= \sum_{a'\in (\mathbb{Z}/2\mathbb{Z})^g}
  \theta_{a',0}(0,2\tau) \theta_{a+a',0}(0,2\tau).
  \end{aligned}
\]
Say we wish to compute `\theta_{a,0}(0,\tau)` for all~`a\in
\{0,1\}^g`. Applying the last formula `n` times, we reduce to evaluating
`\theta_{a,0}(0,2^n\tau)`, and we expect that its absolute value is roughly
`\exp(-2^n\mathrm{Dist}_\tau(0, \mathbb{Z}^g + \tfrac a2))`. Provided that
`n \simeq \log(\mathit{prec})`, we have to sum only `O_g(1)` terms in the naive
algorithm to evaluate `\theta_{a,0}(0,2^n\tau)` at ``shifted absolute
precision'' *prec}, i.e. absolute precision *prec} +
\func{acb_theta_dist_addprec}`(2^n \mathrm{Dist}_\tau(0, \mathbb{Z}^g + \tfrac
a2))`. In order to recover `\theta_{a,0}(0,\tau)`, we then perform `n` AGM
steps. The precision loss when applying \func{acb_theta_agm_mul_tight} is
`O_g(1)` bits in terms of shifted absolute precision. One also has to take
square roots at each step. For this, we assume that each
`|\theta_{a,0}(0, 2^k\tau)|` is indeed of the expected order of
magnitude. Then, using the naive algorithm with `O_g(1)` terms will be
sufficient to determine the correct choices of square roots at each step, with
a precision loss of `O(1)` bits as well. At the end of this algorithm, we
indeed obtain `\theta_{a,0}(0,\tau)` at shifted absolute precision
`\mathit{prec - O_g(n)` bits for each~`a`.

We make the following adjustments to make the algorithm work in general:
\begin{itemize}
\item If we see (after applying the naive algorithm) that some value
  `\theta_{a,0}(0,2^k\tau)` is too small, we introduce an auxiliary real
  vector~`t`. At each step, starting from `\theta_{a,0}(0,2^{k+1}\tau)`,
  `\theta_{a,0}(2^{k+1}t, 2^{k+1}\tau)` and
  `\theta_{a,0}(2^{k+2}t, 2^{k+1}\tau)`, we compute
  `\theta_{a,0}(2^{k}t, 2^k\tau)` and `\theta_{a,0}(2^{k+1}t, 2^k\tau)` using
  square roots (second formula), then `\theta_{a,0}(0, 2^k\tau)` using a
  division (third formula). For a huge majority of such `t`, none of the theta
  values `\theta_{a,0}(2^kt, 2^k\tau)` and~`\theta_{a,0}(2^{k+1}t, 2^k\tau)`
  will be too small \cite{main}. In practice, we choose `t` at random and obtain a
  probabilistic algorithm with a negligible failure probability.
\item When computing `\theta_{a,0}(z,\tau)` for a nonzero~`z`, we compute
  `\theta_{a,0}(0, 2^k\tau)` and `\theta_{a,0}(2^k z, 2^k\tau)` using the
  second and fourth formulas at each step.
\item Finally, these two techniques can be combined by evaluating theta values
  at the six vectors `2^k v` for `v\in\{0, t, 2t, z, z + t, z + 2t\}`. Note
  that we only have to compute `\theta_{a,0}(2^kz, 2^k\tau)` at the very last
  step `k=0`.
\end{itemize}

We use an additional improvement when the eigenvalues of `\mathrm{Im}(\tau)`
have different sizes. Let `\gamma_i` for `0\leq i < g` be the diagonal
coefficients of a Cholesky matrix for `\pi\mathrm{Im}(\tau)`. Let
`1\leq s < g`, and assume that `\gamma_s` is significantly bigger than
`\gamma_{s-1}`, so that one can find `n` such that
`2^{n}\gamma_s^2 \simeq \mathit{prec}` while `2^n \gamma_{s-1}^2` is much
smaller. One can then split the theta series for `\theta_{a,0}(z, 2^n\tau)` and
reduce to computing `O_g(1)` theta values in dimension~`s`. See
\func{acb_theta_ql_a0_split} below for more details.

Finally, we note that the formulas above still hold after replacing each
occurrence of `\theta_{a,0}(z,\tau)` by
`e^{-\pi y^T Y^{-1} y}\theta_{a,0}(z,\tau)`. We use the latter quantities
instead for convenience, since their magnitude does not increase as `y` gets
farther from zero, and their expected absolute values are easily expressed in
terms of lattice distances.

The functions in this section will work best when `\tau` lies in the reduced
domain and the eigenvalues of `\mathrm{Im}(\tau)` are not too large, say in
`O(\mathit{prec})`.

.. function:: slong acb_theta_ql_nb_steps(const arb_mat_t C, slong s, slong prec)

Returns an integer `n` such that `2^n \gamma_s^2 \simeq \mathit{prec}` in the
above notation, meant to be the number of steps to use in the quasi-linear
algorithm for `\theta_{a,0}` (before applying the splitting strategy, in the
case `s > 0`). The precise value of `n` is chosen to optimize performance: see
\func{acb_theta/profile/p-ql_a0_steps}.

\T no test, but used in \func{acb_theta_ql_a0}.

.. function:: void acb_theta_ql_log_rescale(acb_t res, acb_srcptr z, const acb_mat_t tau, slong prec)

Sets *res} to `i y^T Y^{-1} y`. This is used to rescale theta values as explained above.

\T generate `z` and `x` such that `y = C^Tx` where `C` is obtained from
\func{acb_theta_eld_cho}, and check that the result is `i\pi\lVert x\rVert^2`
(for the `L^2` norm).

.. function:: int acb_theta_ql_roots(acb_ptr rts, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
  arb_srcptr d, const acb_mat_t tau, slong nb_steps, slong guard, slong prec)

Attempts to set *rts} to the collection of low-precision roots for the
given choice of `z` and `t`. It is assumed that *d0} (resp. *d})
contains the result of \func{acb_theta_dist_a0} on `(0,\tau)`
(resp. `(z,\tau)`), and that `t` is a real vector.

More precisely, for each `0\leq k < n`, each `v\in \{t, 2t, z + t, z + 2t\}`,
and each `a\in \{0,1\}^g`, we run \func{acb_theta_naive_ind} to evaluate
`\theta_{a,0}(2^kv, 2^k\tau)` at working precision *guard} +
\func{acb_theta_dist_addprec}`(2^k d_k)`, where `d_k` denotes the
`k^{\mathrm{th}}` entry of *d0} or *d}, according to the
imaginary part of `v`. If none of these complex balls contains zero, returns
`1` and sets *rts} to the resulting vector of length `4 \times n \times 2^g`;
otherwise, returns `0` and leaves *rts} undefined. The number of output values is
reduced to `2\times n\times 2^g` or `n\times 2^g` when `z = 0`, `t = 0`, or
both.

\T when `g = 2`, `z = t = 0`, and `\tau` is inside the Siegel fundamental
domain, it is known that the theta values are bounded away from zero, and thus
the return value must be 1.

.. function:: void acb_theta_ql_step_1(acb_ptr res, acb_srcptr th0, acb_srcptr th,
  acb_srcptr rts, arb_srcptr d0, arb_srcptr d, slong g, slong prec)

Given `\theta_{a,0}(0, 2\tau)` (stored in *th0}) and
`\theta_{a,0}(2z, 2\tau)` (stored in *th}), sets *res} to the values
`\theta_{a,0}(z,\tau)` for `a\in \{0,1\}^g`. We assume that *d0}
(resp. *d}) contains the result of \func{acb_theta_dist_a0} on
`(0,2\tau)` (resp. `(2z, 2\tau)`), and that *rts} contains
low-precision approximations of `\theta_{a,0}(z,\tau)`. We call
\func{acb_theta_agm_mul_tight} and \func{acb_theta_agm_sqrt} once each.

\T working at low precision, check that the duplication formula holds by
generating input using \func{acb_theta_naive_fixed_ab}, applying
\func{ql_step_1}, and checking the output against
\func{acb_theta_naive_fixed_ab} as well.

.. function:: void acb_theta_ql_step_3(acb_ptr res, acb_srcptr th0, acb_srcptr th,
  acb_srcptr rts, arb_srcptr d0, arb_srcptr d, slong g, slong prec)

Given `\theta_{a,0}(2v, 2\tau)` for `v\in\{0, t, 2t\}` (stored in *th0}
as a vector of length `3\times 2^g`) and for `v\in\{z, z + t, z + 2t\}` (stored
in *th}), sets *res} to the vector of length `3\times 2^g` containing
`\theta_{a,0}(v,\tau)` for `v\in\{ z, z + t, z + 2t\}` and `a\in
\{0,1\}^g`. The assumptions on *d0} and *d} are as above, and
*rts} must contain low-precision approximations of `\theta(v,\tau)` for
`v\in \{z+t, z+ 2t\}`. We make three calls to \func{acb_theta_agm_mul}, take
`2^{g+1}` square roots, and make `2^g` divisions.

\T check against the naive algorithm as in \func{acb_theta_ql_step_1}.

.. function:: void acb_theta_ql_step_2(acb_ptr res, acb_srcptr th0, acb_srcptr th,
  acb_srcptr rts, arb_srcptr d0, arb_srcptr d, slong g, slong prec)

Same as \func{acb_theta_ql_step_3}, but does not perform the divisions. The first
`2^g` entries of *res} are set to zero.

\T check against the naive algorithm as in \func{acb_theta_ql_step_1}.

.. function:: void acb_theta_ql_dupl(acb_ptr th2, acb_srcptr th0, acb_srcptr th,
  arb_srcptr d0, arb_srcptr d, slong g, slong prec)

Given input as in \func{acb_theta_ql_step_1} (except that *rts} is not
needed), sets `r` to the vector of squared theta values
`\theta_{a,b}(z,\tau)^2` for all `a,b\in \{0,1\}^g`. We use the following
version of the duplication formula:
\[
  \theta_{a,b}(z,\tau)^2 = \sum_{a'\in (\mathbb{Z}/2\mathbb{Z})^g}
  (-1)^{a'^Tb} \theta_{a',0}(2z,2\tau) \theta_{a+a',0}(0,2\tau),
\]
making `2^g` calls to \func{acb_theta_agm_mul_tight}.

\T check against the naive algorithm as in \func{acb_theta_ql_step_1}.

\subsubsection{Quasi-linear algorithms for `\theta_{a,0}`}

The functions in this section will work best when `\tau` lies in the reduced
domain and the eigenvalues of `\mathrm{Im}(\tau)` are not too large, say in
`O(\mathit{prec})`.

.. function:: acb_theta_ql_worker_t}

A function pointer type. A function *worker} of this type has the
following signature:

.. function:: int worker(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_scptr d0,
  arb_srcptr d, const acb_mat_t tau, slong guard, slong prec)

Such a worker will attempt to set *th} to the values `\theta_{a,0}(v,\tau)` for
`v\in \{0,t,2t,z,z+t,z+2t\}` and `a\in \{0,1\}^g` at shifted absolute
precision *prec}, return `1` on success and `0` on failure. The vectors
*d0} and *d} must contain the result of
\func{acb_theta_dist_a0} on `(0,\tau)` and `(z,\tau)`. If `z = 0`, `t = 0`, or
both, we only compute `3`, `2`, or `1` vectors of `2^g` values
respectively. Two functions of this type are available:
\func{acb_theta_ql_a0_naive} and the main function
\func{acb_theta_ql_a0}. Using function pointers allows us to write independent
test code for the main workhorses \func{acb_theta_ql_a0_steps} and
\func{acb_theta_ql_a0_split} below.

.. function:: int acb_theta_ql_a0_naive(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
  arb_srcptr d, const acb_mat_t tau, slong guard, slong prec)

Follows the specifications of a function of type .. type:: acb_theta_ql_worker_t}
using the naive algorithm only. The return value is always `1`.

\T no test, but used throughout in quasi-linear algorithms.

.. function:: int acb_theta_ql_a0_split(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d,
  const acb_mat_t tau, slong s, slong guard, slong prec, acb_theta_ql_worker_t worker)

Follows the specifications of a function of type .. type:: acb_theta_ql_worker_t},
except for the additional arguments *s} and *worker}. We split the
theta series according to the first `s` coordinates of `n\in \mathbb{Z}^g`,
writing `n = (n_0,n_1)` where `n_0\in \mathbb{Z}^s` and
`n_1\in \mathbb{Z}^{g - s}`. Then *worker} is called to evaluate each
term in the split theta series. We must have `1\leq s\leq g -1`.

More precisely, for each `0\leq a < 2^g`, we compute *R2} and *eps}
as in \func{acb_theta_naive_radius} at precision
*prec}`{} + {}`\func{acb_theta_dist_addprec}`(d_a)`. Note that
`n^T \mathrm{Im}(\tau) n\geq \lVert C_1 n_1\rVert^2`, where `C_1` denotes the
lower-right block of `C` of dimensions `(g-s)\times(g-s)`. Thus, in order to
compute `\theta_{a,0}(z, 2^n\tau)` at shifted absolute precision *prec},
it is enough to consider those `n_1\in \mathbb{Z}^{g - s}` that lie in a
certain ellipsoid of radius *R2} for the Cholesky matrix `C_1`. This
ellipsoid is meant to contain very few points, and we list all of them. Then,
for a given choice of `n_1`, the sum of the corresponding terms in the theta
series is
\[
  \begin{aligned}
    &c \sum_{n_0\in \mathbb{Z}^s} \exp(\pi i ((n_0 + \tfrac{a_0}{2})^T\tau_0
    (n_0 + \tfrac{a_0}{2})
    + 2 (n_0 + \tfrac{a_0}{2})^T x (n_1 + \tfrac{a_1}{2}) + 2(n_0 + \tfrac{a_0}{2}) z_0))\\
    &\qquad\qquad = c\, \theta_{a_0,0}(z_0 + x (n_1 + \tfrac{a_1}{2}), \tau_0)
  \end{aligned}
\]
with
\[
  c = \exp(\pi i ((n_1 + \tfrac{a_1}{2})\tau_1 (n_1 + \tfrac{a_1}{2}) + 2 (n_1
  + \tfrac{a_1}{2}) z_1)).
\]
where `\tau = (\begin{smallmatrix} \tau_0 & x\\x^T & \tau_1\end{smallmatrix})`
and `z = (z_0,z_1)`. For each `n_1`, we adjust the shifted absolute precision
for the corresponding term according to the distance between `n_1` and the
center of the above ellipsoid. The return value is 1 iff *worker}
succeeds for each `n_1`.

\T check that the result agrees with \func{acb_theta_ql_a0_naive} on random
input in case of success, using \func{acb_theta_ql_a0_naive} as *worker}.

.. function:: int acb_theta_ql_a0_steps(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
  arb_srcptr d, const acb_mat_t tau, slong nb_steps, slong s,
  slong guard, slong prec, acb_theta_ql_worker_t worker)

Follows the specifications of a function of type \func{acb_theta_ql_worker_t},
except for the additional arguments *nb_steps}, *s} and
*worker}, by performing `k := {}`*nb_steps} AGM steps. We first
call \func{acb_theta_ql_roots} with *guard} bits of shifted absolute
precision, then call \func{acb_theta_ql_a0_naive} or
\func{acb_theta_ql_a0_split} on `(2^k t, 2^k z, 2^k\tau)` depending on whether
*s} is zero or not, and finally we perform the AGM steps. If any
subprocedure fails, we end the computation and return 0, and otherwise return
1.

\T same as \func{acb_theta_ql_a0_split}.

.. function:: int acb_theta_ql_a0(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
  arb_srcptr d, const acb_mat_t tau, slong guard, slong prec)

Follows the specifications of a function of type
\func{acb_theta_ql_worker_t}. We first decide how many AGM steps we should use
and whether we should use the splitting strategy. Then we run
\func{acb_theta_ql_a0_steps} on the midpoints of `t,z` and `\tau` at a slightly
higher precision to account for precision losses in the duplication formulas,
using a recursive call to \func{acb_theta_ql_a0} as *worker}. If the
return value is 1, we finally compute provable error bounds on the result using
\func{acb_theta_jet_naive_ind} and \func{acb_theta_jet_error_bounds}.

\T check that the result agrees with \func{acb_theta_ql_a0_naive} on random
input if successful.

\subsubsection{Quasi-linear algorithms for `\theta_{a,b}`}

The function \func{acb_theta_ql_a0} may fail for an unlucky choice of auxiliary
vector `t` or when *guard} is too small. Thus, we implement a
probabilistic algorithm where we gradually increase *guard} and choose
first `t = 0`, then a random choice of `t` at each step. The following
functions will still work best when `\tau` is reduced, however `\mathrm{Im}(\tau)`
may have large eigenvalues.

.. function:: ACB_THETA_QL_TRY}

The number of times that a new `t` should be picked before abandoning and
setting the result to indeterminate. This is currently set to 100 for
a negligible failure probability.

.. function:: slong acb_theta_ql_reduce(acb_ptr new_z, acb_t c, arb_t u, slong* n1,
  acb_srcptr z, const acb_mat_t tau, slong prec)

Sets *new_z}, *c}, *u}, *n1} and returns
`-1\leq s\leq g` such that the following holds. When `s\geq 0`,
`z':=*new_z}` is a vector of length `s` and `n_1` is a vector of length
`g-s`, and for each characteristic `(a,b)`, we have (borrowing notation from
\func{acb_theta_ql_a0_split}): eitner
\[
  |\theta_{a,b}(z,\tau) - c i^{\,n_1^Tb_1} \theta_{a_0,b_0}(z', \tau_0)| \leq u
\]
when the last `g-s` coordinates of `a` equal
`a_1 =`\func{acb_theta_char_get_a}`(n_1)`, or
\[
  |\theta_{a,b}(z,\tau)|\leq u
\]
otherwise. When `s=-1`, *n1}, *c} and *new_z} are left
undefined and we have `\theta_{a,b}(z,\tau)\leq u` for all characteristics
`(a,b)`. This filters out very large eigenvalues of `\mathrm{Im}(\tau)` that
have a negligible impact on theta values but would give rise to unreasonable
choices of precisions in the duplication formula.

This works as follows. We first compute *R2} and *eps} as in
\func{acb_theta_naive_radius}, then set *c}, *u} and *new_z}
as in \func{acb_theta_naive_reduce} in dimension `g`. We set `s` such that for
each `s\leq j < g`, we have `\gamma_j^2 > 4R^2`, where `\gamma_j` is the
`j^{\mathrm{th}}` diagonal coefficient of the Cholesky matrix `C` for
`\pi\mathrm{Im}(\tau)`. We may assume that `s< g`, otherwise there is nothing
to be done. Then the ellipsoid `E` of radius `R^2` for `C` that we are
interested in, when intersected with `\frac12\mathbb{Z}^g`, is either empty or
consists of points whose last `g-s` coordinates are fixed. In the latter case,
we return `s = -1`. Now assume that `E` is not empty, let `n_1` be the vector
of these fixed last `g-s` coordinates, and let `a_1\in \{0,1\}^{g-s}` be the
corresponding characteristic. We can then write the sum defining `\theta_{a,b}`
over `E` as
\[
  e^{i\pi (\tfrac{n_1^T}{2} \tau_1 \tfrac{n_1}{2} + n_1^T(z_1 + \tfrac{b_1}{2}))
  \sum_{n_0\in E_0 \cap (\mathbb{Z}^s + \tfrac{a_0}{2}) e^{i\pi(n_0^T \tau_0 n_0 + 2n_0^T(z_0 + x \tfrac{n_1}{2} + \tfrac{b_0}{2}))
\]
if the last `g-s` coordinates of `a` are equal to `a_1` on the nose; the sum is
zero otherwise. Thus we can set `z'` to `z_0 + x\tfrac{n_1}{2}` and multiply
`c` by `\exp(i\pi (\tfrac{n_1^T}{2}\tau_1\tfrac{n_1}{2} + n_1^Tz_1))`.

\T generate random *tau} and *z} that are likely to lead to
*s} being positive or `-1` and *n1} being nonzero, and check that
the above conditions hold when computing theta values with the naive algorithm.

.. function:: void acb_theta_ql_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec)

Sets *th} to the collection of `\theta_{a,b}(z,\tau)` for all
`a,b\in \{0,1\}^g`. After calling \func{acb_theta_ql_reduce}, we generally use
the duplication formula on the result of \func{acb_theta_ql_a0} at `2\tau` and
a final square-root step. At low precisions, we call \func{acb_theta_naive_all}
instead.

\T check that the result agrees with \func{acb_theta_naive_all} on random input.

.. function:: void acb_theta_ql_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)

Sets *th2} to the collection of `\theta_{a,b}(z,\tau)^2` for all
`a,b\in \{0,1\}^g`. After calling \func{acb_theta_ql_reduce}, we use
the duplication formula on the result of \func{acb_theta_ql_a0} at `2\tau`.

\T check that the result agrees with \func{acb_theta_naive_all} on random input.

\subsection{The transformation formula}

The functions in this section implement the theta transformation formula
\cite[p.\,176]{igusa}.

.. function:: ulong acb_theta_transform_char(slong* e, const fmpz_mat_t mat, ulong ab)

Returns the theta characteristic `(a',b')` and sets `0\leq e < 8` such that the
transformation formula reads: for every `\tau\in \mathbb{H}_g`,
\[
  \theta_{a,b}(0,\mathit{mat}\cdot \tau) = c\,\zeta_8^e\, \theta_{a',b'}(0,\tau)
\]
where `c` depends only on *mat} and `\tau` and `\zeta_8=\exp(i\pi/4)`. In
Igusa's notation, *e} is `\phi_m(\mathit{mat})`.

\T check that the `a` component of any characteristic remains the same when
*mat} is a trigonal symplectic matrix as in \func{sp2gz_trig}.

.. function:: slong acb_theta_transform_kappa(const fmpz_mat_t mat)

Returns an integer `0\leq e < 8` such that in the transformation formula, we
have `\kappa(\mathit{mat}) = \zeta_8^e`. The sign of `\kappa(\mathit{mat})` is
fixed by making an arbitrary choice of `\det(c \tau + d)` when `\tau` is `i`
times the identity matrix.

\T check that for a block-diagonal symplectic matrix *mat} constructed
from `U\in \mathrm{GL}_g(\mathbb{Z})`, `\kappa^2` is `\det(U)` in agreement
with Igusa's results.

.. function:: void acb_theta_transform_sqrtdet(acb_t res, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)

Sets *res} to `\sqrt{\det(\gamma\tau + \delta)` where `\gamma,\delta`
denote the lower `g\times g` blocks of *mat}. The choice of square root
is made so that the transformation formula holds, and is determined by
computing theta values at low precision.

\T check that the result squares to the determinant of \func{acb_siegel_cocycle}.

.. function:: void acb_theta_transform_proj(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th,
  int sqr, slong prec)

Assuming that *sqr} is 0 (false) and that *th} contains
`\theta_{a,b}(z,\tau)` for some `z\in \mathbb{C}^g` and `\tau\in \mathbb{H}_g`,
sets *res} to contain the values
`\theta_{a,b}(\mathit{mat}\cdot (z,\tau))` (where *mat} acts as in
\func{acb_theta_transform_z}) up to a common scalar factor in
`\mathbb{C}^\times`. This only permutes the theta values and multiplies them by
a suitable eighth root of unity. If *sqr} is nonzero (true), does the
same computation for squared theta values `\theta_{a,b}(z,\tau)^2` instead.

\T check that applying \func{acb_theta_transform_proj} using a random
*mat} then its inverse gives back the initial projective point.

.. function:: void acb_theta_transform(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th,
  acb_srcptr z, const acb_mat_t tau, slong kappa, int sqr, slong prec)

Assuming that *sqr} is 0, that *kappa} is precomputed as in
\func{acb_theta_transform_kappa}, and that *th} contains
`\theta_{a,b}(z,\tau)`, sets *res} to vector of values
`\theta_{a,b}(\mathit{mat}\cdot(z,\tau))` for `a,b\in\{0,1\}^g`. If *sqr}
is nonzero, does the same computation for squared theta values instead.

\T check that the result agrees with \func{acb_modular_theta} when `g=1` on
random input. We restrict to `g=1` to avoid calling the naive algorithm on a
matrix that is far from the fundamental domain.

.. function:: void acb_theta_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, int sqr, slong prec)

Sets *th} to the vector of theta values `\theta_{a,b}(z,\tau)` or
`\theta_{a,b}(z,\tau)^2` for `a,b\in \{0,1\}^g`, depending on whether
*sqr} is 0 (false) or not. We reduce `\tau` using
\func{acb_theta_siegel_reduce}, call \func{acb_theta_ql_all} or
\func{acb_theta_ql_all_sqr}, and then apply the transformation formula. If the
reduction is not successful, we call the naive algorithm at a lower precision
instead.

\T check that the result agrees with \func{acb_theta_naive_all} on random
input. The matrix *tau} is chosen to be a priori non-reduced but
reasonably close to the fundamental domain.

\subsection{Quasi-linear algorithms for derivatives}

We implement an algorithm for derivatives of theta functions based on finite
differences. It is quasi-linear in terms of the precision and the number of
derivatives to be computed.

Consider the Fourier expansion:
\[
  \begin{aligned}
  \theta_{a,b}(z + h, \tau) &= \sum_{k\in \mathbb{Z}^g,\ k\geq 0}
  \frac{1}{k_0!}\cdots \frac{1}{k_{g-1}!}
  \frac{\partial^{|k|}\theta_{a,b}}{\partial z_0^{k_0}\cdots \partial
    z_{g-1}^{k_{g-1}}}(z,\tau)\cdot h_0^{k_0}\cdots
  h_{g-1}^{k_{g-1}}\\
  &=: \sum_{k\in \mathbb{Z}^g,\ k\geq 0} a_k\, h_0^{k_0}\cdots h_{g-1}^{k_{g-1}}.
  \end{aligned}
\]
The basic observation is that if one chooses
`h = h_n = (\varepsilon \zeta^{n_0},\ldots, \varepsilon \zeta^{n_{g-1}})` where
`\varepsilon > 0` and `\zeta` is a primitive `j^{\mathrm{th}}` root of unity
and lets `n` run through all vectors in `\{0,\ldots, j - 1\}^g`, then taking a
discrete Fourier transform of the resulting values will compute the individual
Taylor coefficient for each derivation tuple that is bounded by `j`
elementwise. A constant proportion, for fixed `g`, of this set consists of all
tuples of total order at most `j`. More precisely, fix `p\in
\mathbb{Z}^g`. Then
\[
  \sum_{n\in \{0,\ldots,m-1\}^g} \zeta^{-p^T n} \theta_{a,b}(z + h_n, \tau) =
  m^g \sum_{\substack{k\in \mathbb{Z}^g,\ k\geq 0,\\\ k = p\ (\text{mod } m)}
  a_k\,\varepsilon^{|k|} =: m^g (a_p\,\varepsilon^{|p|} + T).
\]
Observe that the magnitude gap between the leading term in
the latter sum and the next ones is `\varepsilon^m`, not `\varepsilon`. This is
is crucial for the algorithm to remain quasi-linear in the precision and the
number of derivatives to be computed.


In order to give an upper bound on `T`, we use the Cauchy integration
formula. Assume that `|\theta_{a,b}(z,\tau)|\leq c` uniformly on a ball of
radius `\rho` centered in `z` for `\lVert\cdot\rVert_\infty`. Then
`|a_k|\leq c/\rho^{|k|}`, so, assuming that `(2g)^{1/m}\varepsilon \leq \rho`,
\[
  |T|\leq c \left(\frac{\varepsilon}{\rho}\right)^{|p|} \sum_{j\geq 1} \binom{g
    - 1 + j}{j} \left(\frac{\varepsilon}{\rho}\right)^{mj} \leq 2 c
  \left(\frac{\varepsilon}{\rho}\right)^{|p|}
  \frac{(g\varepsilon/\rho)^m}{1 - (g\varepsilon/\rho)^m} \leq 2c g\,\frac{\varepsilon^{|p|+m}}{\rho^m}.
\]
Since we divide by `\varepsilon^{|p|}` to get `a_p`, we will add an error of
`2c g (\varepsilon/\rho)^m`.

.. function:: void acb_theta_jet_bounds(arb_t c, arb_t rho, acb_srcptr z, const
  acb_mat_t tau, slong ord)

Sets *c} and *rho} such that on every ball centered at (a point
contained in) *z} of radius *rho}, the functions `|\theta_{a,b}|`
for all characteristics `(a,b)` are uniformly bounded by `c`. The choice of
*rho} is tuned to get interesting upper bounds on derivatives of
`\theta_{a,b}` up to order *ord}.

We proceed as follows. First, we compute `c_0`, `c_1`, `c_2` such that for any
choice of `\rho`, one can take `c = c_0\exp((c_1 + c_2\rho)^2)`
above. Following \func{acb_theta_naive_reduce}, we get
\[
  |\theta_{a,b}(z',\tau)| \leq c_0\exp(\pi y'^T Y^{-1} y')
\]
where `c_0 = 2^g \prod_{j=0}^{g-1} (1 + 2\gamma_j^{-1})` \cite{main}. In turn,
if `\lVert z' - z\rVert_\infty\leq \rho`, then
\[
  \pi y'^T Y^{-1} y' \leq \bigl(\sqrt{\pi y^T Y^{-1} y} + \rho \sup_{\lVert x
    \rVert_\infty\leq 1} \sqrt{\pi x^T Y^{-1} x}\bigr)^2
\]
by the triangle inequality for the quadratic form `\pi Y^{-1}`. An upper bound
`c_2` on the sup is easily computed from the Cholesky matrix for `\pi Y^{-1}`,
while we can take `c_1 = \sqrt{\pi y^T Y^{-1} y}`.

Once `c_0, c_1, c_2` are computed, we look for a value of `\rho` that minimizes
`\exp((c_1 + c_2\rho)^2)/\rho^{m}` where `m = \mathit{ord}+1`, so we set `\rho`
to the positive root of `2c_2\rho (c_1 + c_2\rho) = m`.

\T on reasonable input, check that *c} and *rho} are finite, and
check that their definition is satisfied by sampling theta values on the
corresponding ball at low precisions.

.. function:: void acb_theta_jet_fd_radius(arf_t eps, arf_t err, const arb_t c, const arb_t rho,
  slong ord, slong g, slong prec)

Sets *eps} and *err} to be a suitable radius and error bound for
computing derivatives up to total order *ord} at precision *prec},
given *c} and *rho} as above. We want
`(2 g)^{1/m} \varepsilon \leq \rho` and
`2 c g (\varepsilon/\rho)^{m} \leq 2^{-\mathit{prec}}` where
`m = \mathit{ord} + 1`, so we set `\varepsilon` to a lower bound for
`\rho \cdot (\min\{2^{-\mathit{prec}}/c, 1\}/2g)^{1/m}`. We also set
*err} to `2^{-\mathit{prec}}`.

\T check that these inequalities are satisfied on random choices of *c} and *rho}.

.. function:: void acb_theta_jet_fd(acb_ptr dth, const arf_t eps, const arf_t err, acb_srcptr val,
  slong ord, slong g, slong prec)

Assuming that *val} contains the values `\theta_{a,b}(z + h_n,\tau)`
where `h_n = (\varepsilon \zeta^{n_0},\ldots, \varepsilon \zeta^{n_{g-1}})` for
a root of unity `\zeta` of order `\mathit{ord} + 1`, and assuming that
*eps} and *err} has been computed as in
\func{acb_theta_jet_fd_radius}, sets *dth} to the vector of partial
derivatives of `\theta_{a,b}` at `(z,\tau)` up to total order *ord}. The
vector *val} should be indexed in lexicographic order as in
\func{acb_dft}, i.e. writing `j = \overline{a_{g-1}\cdots a_0}` in basis `m`,
the `j^{\mathrm{th}}` entry of *val} corresponds to
`n = (a_0,\ldots, a_{g-1})`. The output derivatives are normalized as in the
Taylor expansion.

\T check that this computes the correct Fourier coefficients for the
exponential function `\exp(z_0+\cdots+z_{g-1})`, setting *c} and
*rho} by hand instead of calling \func{acb_theta_jet_bounds}.

.. function:: void acb_theta_jet_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec)

Sets *dth} to the derivatives of all functions `\theta_{a,b}` for
`a,b\in \{0,1\}^g` at `(z,\tau)`, as a concatenation of `2^{2g}` vectors of
length \func{acb_theta_jet_nb(ord, g). This algorithm runs in quasi-linear
time in `\mathit{prec}\cdot \mathit{ord}^g` for any fixed `g`.

We first compute *c}, *rho}, *err} and *eps} as above,
then compute theta values `\theta_{a,b}(z + h_n,\tau)` at a higher precision to
account for division by `\varepsilon^{\mathit{ord}}\cdot
(\mathit{ord}+1)^g`. For this, we first extract the midpoint of `z` and
`\tau`. Finally, we adjust the error bounds using
\func{acb_theta_jet_error_bounds} and the naive algorithm for derivatives of
order `\mathit{ord} + 2`.

\T check that the output agrees with \func{acb_theta_jet_naive_all} on random input.

\subsection{Dimension~`2` specifics}

In the `g=2` case, one can use theta functions to evaluate many fundamental
Siegel modular forms. This section contains functions to do so, in analogy with
\func{acb_modular_delta}, \func{acb_modular_eisenstein}, etc. when `g=1`.

We use the following notation. For `k,j\geq 0`, a Siegel modular form of weight
`\det^k\otimes \mathrm{Sym}^j` is by definition an analytic function
`f: \mathbb{H}_g\to \mathbb{C}_j[X]` (the vector space of polynomials of degree
at most~`j`) such that for any `\tau\in \mathbb{H}_g` and
`m\in \mathrm{Sp}_4(\mathbb{Z})`, we have
\[
  f((\alpha\tau + \beta)(\gamma\tau + \delta)^{-1}) = \det(\gamma\tau +
  \delta)^k\cdot \mathrm{Sym}^j(\gamma\tau + \delta)(f(\tau)).
\]
Here `\alpha,\beta,\gamma,\delta` are the `g\times g` blocks of `m`, and
`\mathrm{Sym}^j(r)` for
`r = \smash{(\begin{smallmatrix} a & b\\ c & d\end{smallmatrix})\in
\mathrm{GL}_2(\mathbb{C})` is the map
`P(X) \mapsto (b X + d)^j P(\tfrac{a X + c}{b X + d})`. For a nonzero `f` to
exist, `j` must be even.

Siegel modular forms generate a bi-graded ring which is not finitely
generated. However, if we relax the definition of a Siegel modular form and
allow them to have a pole along the diagonal
`\mathbb{H}_1^2 = \{(\begin{smallmatrix} \tau_1 & 0 \\ 0 &
  \tau_2\end{smallmatrix})\}\subset \mathbb{H}_2` of a certain order (depending
on the weight), we indeed find a finitely generated ring with 26 generators
corresponding to classical \emph{covariants}, studied e.g. by Clebsch
\cite{clebsch}. For historical reasons, covariants are classified in terms of
their degree `k` and index `j`, corresponding to Siegel modular functions of
weight `\det^{k - j/2}\otimes \mathrm{Sym}^j`.

.. function:: ACB_THETA_G2_COV_NB}

Macro giving the number of generators of the ring of covariants, equal to 26.

.. function:: void acb_theta_g2_jet_naive_1(acb_ptr dth, const acb_mat_t tau, slong prec)

Sets *dth} in the same way as .. function:: acb_theta_jet_naive_all(dth, z, tau,
  1, prec) for `z = 0`, but works more efficiently, since the value
(resp. gradients) of `\theta_{a,b}(z,\tau)` at `z = 0` vanish if `(a,b)` is odd
(resp. even). The attached worker uses one of two available strategies (doing
multiplications and then summing, or calling \func{acb_dot} twice) depending on
*prec}.

\T check that the output agrees with \func{acb_theta_jet_naive_all}.

.. function:: void acb_theta_g2_detk_symj(acb_poly_t res, const acb_mat_t m, const acb_poly_t f,
  slong k, slong j, slong prec)

Sets *res} to `\det(m)^k \mathrm{Sym}^j(m)(f)`. The polynomial `f` should
be of degree at most `j` (any coefficients of larger degree are ignored).

\T check that the chain rule holds when `m` is obtained as a product of two matrices.

.. function:: void acb_theta_g2_transvectant(acb_poly_t res, const acb_poly_t g, const acb_poly_t h,
  slong m, slong n, slong k, slong prec)

Sets *res} to the `k^{\mathrm{th}}` transvectant of the polynomials `g`
and `h` of degrees `m` and `n`: considering `g` and `h` as homogeneous
polynomials of degree `m` (resp. `n`) in `x_1,x_2`, this sets *res} to
\[
  (g,h)_k := \frac{(m-k)!(n-k)!}{m!n!}  \sum_{j=0}^{k} (-1)^{k-j} \binom{k}{j}
  \frac{\partial^k g}{\partial x_1^{k-j}\partial x_2^j} \frac{\partial^k
    h}{\partial x_1^{j}\partial x_2^{k-j}}.
\]
Any coefficients of `g` or `h` of larger degree than `m` (resp. `n`) are
ignored.

\T check that for `f = \sum_{j=0}^6 a_jx^{6-j}`, we have
`(f,f)_6 = -3a_3^2 + 8a_2a_4 - 20a_1a_5 + 120a_0a_6`.

.. function:: void acb_theta_g2_transvectant_lead(acb_t res, const acb_poly_t g, const acb_poly_t h,
  slong m, slong n, slong k, slong prec)

Sets *res} to the leading coefficient of `(g,h)_k` in `x_1`, with the same
conventions as above.

\T check that we indeed get the leading term of the transvectant computed using
\func{acb_theta_g2_transvectant}.

.. function:: void acb_theta_g2_psi4(acb_t res, acb_srcptr th2, slong prec)

.. function:: void acb_theta_g2_psi6(acb_t res, acb_srcptr th2, slong prec)

.. function:: void acb_theta_g2_chi10(acb_t res, acb_srcptr th2, slong prec)

.. function:: void acb_theta_g2_chi12(acb_t res, acb_srcptr th2, slong prec)

Sets *res} to the value of the Eisenstein series `\psi_4`, `\psi_6` or
the cusp forms `\chi_{10}, \chi_{12}` corresponding to the given vector of
squared theta values. We use the formulas from \cite[§7.1]{streng}, with the
following normalizations: `\psi_4 = h_4/4`, `\psi_6 = h_6/4`,
`\chi_{10} = -2^{-12} h_{10}`, `\chi_{12} = 2^{-15}h_{12}`. We warn that these
differ from the classical notation of Igusa. Writing
`\tau = (\begin{smallmatrix} \tau_1 & \tau_2 \\ \tau_2 &
  \tau_3\end{smallmatrix})` and `q_j = \exp2(\pi i \tau_j)`, the Fourier
expansions of these modular forms begin as follows:
\[\begin{aligned}
    \psi_4(\tau) &= 1 + 240(q_1 + q_3) + \cdots\\
    \psi_6(\tau) &= 1 - 504(q_1 + q_3) + \cdots\\
    \chi_{10}(\tau) &= (q_2 - 2 + q_2^{-1}) q_1 q_3 + \cdots\\
    \chi_{12}(\tau) &= (q_2 + 10 + q_2^{-1}) q_1 q_3 + \cdots.
  \end{aligned}
\]

\T check that the values transform as they should under \func{acb_theta_transform_proj}.

.. function:: void acb_theta_g2_chi5(acb_t res, acb_srcptr th, slong prec)

Sets *res} to the value of
`\chi_5 = - 2^{-6} \prod_{(a,b)\text{ even}} \theta_{a,b}` corresponding to the
given vector of theta values. The form `\chi_5` is a Siegel cusp form with
character: see \cite{clery} for more details.

\T check that `\chi_5^2=\chi_{10}`.

.. function:: void acb_theta_g2_chi35(acb_t res, acb_srcptr th, slong prec)

Sets *res} to the value of the cusp form `\chi_{35}` corresponding to the vector
of theta values. The form `\chi_{35}` is the unique scalar-valued Siegel
modular form of weight `\det^{35}\otimes \mathrm{Sym}^0` up to scalars, and is
normalized as follows:
\[
  \chi_{35}(\tau) = q_1^2 q_3^2 (q_1 - q_3 )(q_2 - q_2^{-1}) + \cdots
\]

\T check that the values transform as they should under \func{acb_theta_transform_proj}.

.. function:: void acb_theta_g2_chi3_6(acb_poly_t res, acb_srcptr dth, slong prec)

Sets *res} to the value of the vector-valued cusp form with character
`\chi_{6,3}` of weight `\det^3\otimes \mathrm{Sym}^6` corresponding to the
given values of *dth}, computed as in e.g.
\func{acb_theta_g2_jet_naive_1}. We have by \cite{clery}:
\[
  \chi_{3,6}(\tau) = \frac{1}{64\pi^6} \prod_{(a,b) \text{ odd}}
  \left(\frac{\partial \theta_{a,b}}{\partial z_1}(0,\tau) x_1 +
    \frac{\partial\theta_{a,b}}{\partial z_2}(0,\tau) x_2\right).
\]

\T check that `\chi_5\chi_{3,6}` defines a modular form of weight `\det^{8}\mathrm{Sym}^6`.

.. function:: void acb_theta_g2_sextic(acb_poly_t res, const acb_mat_t tau, slong prec)

Sets *res} to the value of `\chi_{-2,6}:=\chi_{3,6}/\chi_5` at `\tau`. We
reduce `\tau` to the Siegel fundamental domain and call either
\func{acb_theta_g2_jet_naive_1} or \func{acb_theta_jet_all} to compute theta
gradients, depending on *prec}. Under the correspondence between Siegel
modular functions and covariants of binary sextics, `\chi_{-2,6}` corresponds
to the binary sextic itself, hence the name.

\T check that the discriminant of \func{acb_theta_g2_sextic} is `2^{12}\chi_{10}`.

.. function:: void acb_theta_g2_covariants(acb_poly_struct* res, const acb_poly_t f, slong prec)

Sets *res} to the vector of 26 generators of the ring of covariants
evaluated at the given sextic *f} (any terms of degree `>6` are ignored),
in the following order:
\begin{multicols}{2}
\begin{enumerate}
  \setcounter{enumi}{-1}
\item `C_{1,6}=f`
\item `C_{2,0}= 60(f,f)_6`
\item `C_{2,4}= 75(f,f)_4`
\item `C_{2,8}= 90(f,f)_2`
\item `C_{3,2}= 30(f,C_{2,4})_4`
\item `C_{3,6}= 30(f,C_{2,4})_2`
\item `C_{3,8}= 6(f,C_{2,4})_1`
\item `C_{3,12}= 6 (f,C_{2,8})_1`
\item `C_{4,0}= 2 (C_{2,4},C_{2,4})_4`
\item `C_{4,4}= 30 (f,C_{3,2})_2`
\item `C_{4,6}= 6(f,C_{3,2})_1`
\item `C_{4,10}= 2(C_{2,8},C_{2,4})_1`
\item `C_{5,2}=(C_{2,4},C_{3,2})_2`
\item `C_{5,4}=\frac 25 (C_{2,4},C_{3,2})_1`
\item `C_{5,8}= 2(C_{2,8},C_{3,2})_1`
\item `C_{6,0}= 2(C_{3,2},C_{3,2})_2`
\item `C_{6,6}^{(1)= \frac 25(C_{3,6},C_{3,2})_1`
\item `C_{6,6}^{(2)= \frac 83(C_{3,8},C_{3,2})_2`
\item `C_{7,2}= 30(f,C_{3,2}^2)_4`
\item `C_{7,4}= 12(f,C_{3,2}^2)_3`
\item `C_{8,2}= \frac 25(C_{2,4},C_{3,2}^2)_3`
\item `C_{9,4}= 4(C_{3,8},C_{3,2}^2)_4`
\item `C_{10,0}= 20(f,C_{3,2}^3)_6`
\item `C_{10,2}= \frac 65(f,C_{3,2}^3)_5`
\item `C_{12,2}= \frac 85(C_{3,8},C_{3,2}^3)_6`
\item `C_{15,0}= \frac{1}{30000} (C_{3,8},C_{3,2}^4)_8`.
\end{enumerate}
\end{multicols}
The scalar factors are chosen so that when evaluated at a formal sextic
`f = \sum a_i x_1^{6-i}x_2^i`, the covariants are integral and primitive as
multivariate polynomials in `a_0,\ldots,a_6,x_1,x_2`.

\T check that the output agrees with \func{acb_theta_g2_psi4} using the
relation `\psi_4 = -(C_{2,0} - 3C_{4,0})/20`. Also check that covariants
transform as they should under the action of `\mathrm{Sp}_4(\mathbb{Z})`, and
that covariants take integral values on integral polynomials.

.. function:: void acb_theta_g2_covariants_lead(acb_ptr res, const acb_poly_t f, slong prec)

Sets *res} to the vector of leading coefficients in `x_1` of the 26
covariants evaluated at *f}. This is more efficient than taking leading
coefficients of \func{acb_theta_g2_covariants}, since we can use
\func{acb_theta_g2_transvectant_lead} instead of
\func{acb_theta_g2_transvectant}.

\T check that the result agrees with taking leading coefficients of
\func{acb_theta_g2_covariants}.
