.. _acb-theta:

**acb_theta.h** -- Riemann theta functions
===============================================================================

This module provides methods for the numerical evaluation of theta functions in
any dimension `g\geq 1`. The algorithms will be detailed in the forthcoming paper
[EK2025]_. In the case `g=1`, we rely on, but also improve on functionality
from :ref:`acb_modular.h <acb-modular>`.

In the context of this module, *tau* or `\tau` always denotes an element of the
Siegel upper half-space `\mathbb{H}_g`, which consists of all symmetric
`g\times g` complex matrices with positive definite imaginary part. The letter
`z` denotes an element of `\mathbb{C}^g`. For each `a,b\in \{0,1\}^g`, the
Riemann theta function (of level 2) of characteristic `(a,b)` is the following
analytic function in `\tau\in \mathbb{H}_g` and `z\in \mathbb{C}^g`:

    .. math::

        \theta_{a,b}(z,\tau) = \sum_{n\in \mathbb{Z}^{g} + \tfrac a2} \exp(\pi i n^T\tau n + 2\pi i n^T (z + \tfrac b2)),

considering `a`, `b` and `z` as column vectors.

We encode a theta characteristic `a\in \{0,1\}^g` as the :type:`ulong` between
`0` and `2^g-1` that has the corresponding expansion in base 2: thus `a =
(1,0,0)` for `g = 3` will be numbered `4`. We also use this encoding to order
vectors of theta values throughout. Similarly, a pair of characteristics
`(a,b)` is encoded as an :type:`ulong` between `0` and `2^{2g}-1`, where `a`
corresponds to the `g` more significant bits. With these conventions, the
output of :func:`acb_modular_theta` is
`(-\theta_3,\theta_2,\theta_0,\theta_1)` in dimension 1.

The main method to evaluate theta functions is

.. function:: void acb_theta_all(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, int sqr, slong prec)

Here *zs* should be a vector of length *nb* `\times g`, and encodes a tuple of
*nb* elements `z_1,\ldots,z_n\in \mathbb{C}^g`. The output, placed in *th*, is
the concatenation of the vector containing `\theta_{a,b}(z_1,\tau)` or
`\theta_{a,b}(z_1,\tau)^2` for all `a,b\in \{0,1\}^g`, depending on whether
*sqr* is 0 (false) or nonzero (true), followed by the same for `z_2`, etc., up
to `z_n`. (This function is faster when *sqr* is nonzero.) These intput/output
conventions hold for several methods in this module, even if we don't restate
them explicitly.

:func:`acb_theta_all` first reduces the input using the action of the Siegel
modular group `\mathrm{Sp}_{2g}(\mathbb{Z})` (the symplectic group) on
`\mathbb{C}^g\times \mathbb{H}_g`, evaluates theta functions on the reduced
arguments, and finally applies the transformation formula for theta functions
under `\mathrm{Sp}_{2g}(\mathbb{Z})`. The second step (evaluating theta
functions) uses an advanced algorithm based on duplication forumas that has a
uniform, quasi-linear complexity in terms of the required precision.

This module also provides functionality to evaluate derivatives of theta
functions, and to evaluate Siegel modular forms in terms of theta functions
when `g=2`.

We conclude this presentation with a word on error bounds and on the argument
*prec*. Barring unexpected cancellations, the "expected" absolute value of
`\theta_{a,b}(z,\tau)` is

    .. math::

        \left|\theta_{a,b}(z,\tau)\right| \approx \exp(\pi y^T Y^{-1} y) \exp(- d^2)

where

- `Y` and `y` denote the imaginary parts of `\tau` and `z` respectively (we
  keep this notation throughout);
- `d` denotes the distance between the point `v = -Y^{-1}y \in \R^g` and the
  shifted lattice `\Z^g + \tfrac{a}{2} \subset \R^g` for the Euclidean norm
  given by the Gram matrix `\pi Y`.

This leads us to define the normalized functions

    .. math::

        \widetilde{\theta}_{a,b}(z,\tau) = \exp(-\pi y^T Y^{-1} y) \theta_{a,b}(z,\tau)

which are no longer holomorphic, but are uniformly bounded on `\mathbb{C}^g`
for a fixed `\tau`. (Adding an element of `\mathbb{Z}^g + \tau\mathbb{Z}^g` to
`z` multiplies the value by a complex number of absolute value one.)

In most cases, when evaluating theta functions `\theta_{a,b}`
(resp. `\widetilde{\theta}_{a,b}`) to precision *prec*, the expected error
bound on the output should be of the order of `\exp(\pi y^Y Y^{-1} y) \cdot
2^{-\mathit{prec}}` (resp. `2^{-\mathit{prec}}`) to prevent unreasonable
computations when `y` is very far from zero. Some internal functions also take
the factor `\exp(-d^2)` into account, and are documented as such.

In any case, the numerical functions in this module always compute certified
error bounds: for instance, if `\tau` is represented by an :type:`acb_mat_t`
which is not certainly positive definite at the chosen working precision, the
output will have an infinite radius.

Main user functions
-------------------------------------------------------------------------------

In addition to :func:`acb_theta_all`, the following functions are the main
outcome of this module. Note that when the inputs `(z,\tau)` are already known
to be reduced, one can bypass the transformation formula using
:func:`acb_theta_all_notransform` and similar functions documented later on.

.. function:: void acb_theta_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

    Sets *th* to the vector containing `\theta_{0,0}(z,\tau)` for all input
    values of `z`.

    We stress that the intermediate evaluation step in this case is not
    quasi-linear in terms of the required precision: rather, it is polynomial
    with exponent `1+g/2`. Depending on `g`, the required precision, and the
    shape of the matrix `tau`, calling `acb_theta_all` then extracting the
    desired entries might be faster.

.. function:: void acb_theta_jet_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, slong prec)

    Sets *th* to the vector or partial derivatives of `\theta_{0,0}` respect to
    `z` up to total order *ord* at the given points `(z,\tau)`.

    The same warning as in :func:`acb_theta_00` applies.

.. function:: void acb_theta_jet_all(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, slong prec)

    Sets *th* to the vector of partial derivatives with respect to `z` up to
    total order *ord* of all functions `\theta_{a,b}` for `a,b\in \{0,1\}^g` at
    the given points `(z,\tau)`. This is a concatenation of *nb* vectors, each
    of them being a concatenation of `2^{2g}` vectors of derivatives of an
    individual `\theta_{a,b}`. (See below for conventions on the numbering and
    normalization of derivatives.)

    The algorithm we use is based on finite differences on the output of
    :func:`acb_theta_all`. At low precisions, using a direct summation
    algorithm will probably be faster: see :func:`acb_theta_sum_jet_all` below.

Example of usage
-------------------------------------------------------------------------------

The following code snippet constructs the period matrix `\tau = iI_2` for `g =
2`, computes the associated theta values at `z = 0` at 10000 bits of precision,
and prints them.

.. code-block:: c

    #include "acb_theta.h"

    int main()
    {
        acb_mat_t tau;
        acb_ptr th, z;
        slong prec = 10000;

        acb_mat_init(tau, 2, 2);
        z = _acb_vec_init(2);
        th = _acb_vec_init(16);

        acb_mat_onei(tau);
        acb_theta_all(th, z, 1, tau, 0, prec);
        _acb_vec_printd(th, 16, 5);

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2);
        _acb_vec_clear(th, 16);
        flint_cleanup();
        return 0;
    }

::

       (1.1803 + 0j)  +/-  (2.23e-3010, 1.23e-3010j), (0.99254 + 0j)  +/-  (1.73e-3010, 1.23e-3010j), (0.99254 + 0j)  +/-  (1.73e-3010, 1.23e-3010j), (0.83463 + 0j)  +/-  (1.73e-3010, 1.23e-3010j), (0.99254 + 0j)  +/-  (1.73e-3010, 1.23e-3010j), (0 + 0j)  +/-  (1.23e-3010, 1.23e-3010j), (0.83463 + 0j)  +/-  (1.73e-3010, 1.23e-3010j), (0 + 0j)  +/-  (1.23e-3010, 1.23e-3010j), (0.99254 + 0j)  +/-  (1.73e-3010, 1.23e-3010j), (0.83463 + 0j)  +/-  (1.73e-3010, 1.23e-3010j), (0 + 0j)  +/-  (1.23e-3010, 1.23e-3010j), (0 + 0j)  +/-  (1.23e-3010, 1.23e-3010j), (0.83463 + 0j)  +/-  (1.73e-3010, 1.23e-3010j), (0 + 0j)  +/-  (1.23e-3010, 1.23e-3010j), (0 + 0j)  +/-  (1.23e-3010, 1.23e-3010j), (0 + 0j)  +/-  (1.23e-3010, 1.23e-3010j)

The Siegel modular group
-------------------------------------------------------------------------------

We use the type :type:`fmpz_mat_t` to handle matrices in
`\operatorname{Sp}_{2g}(\mathbb{Z})`. In addition to the functions in this
section, methods from :ref:`fmpz_mat.h <fmpz-mat>` such as
:func:`fmpz_mat_equal` can thus be used on symplectic matrices directly.

In the following functions (with the exception of :func:`sp2gz_is_correct`) we
always assume that the input matrix *mat* is square of even size `2g`, and
write it as

    .. math::

        m = \begin{pmatrix} \alpha&\beta\\ \gamma&\delta \end{pmatrix}

where `\alpha,\beta,\gamma,\delta` are `g\times g` blocks.

.. function:: slong sp2gz_dim(const fmpz_mat_t mat)

    Returns `g`, which is half the number of rows (or columns) of *mat*.

.. function:: void sp2gz_set_blocks(fmpz_mat_t mat, const fmpz_mat_t alpha, const fmpz_mat_t beta, const fmpz_mat_t gamma, const fmpz_mat_t delta)

    Sets *mat* to `\bigl(\begin{smallmatrix} \alpha&\beta\\ \gamma&\delta
    \end{smallmatrix}\bigr)`. The dimensions must match.

.. function:: void sp2gz_j(fmpz_mat_t mat)

    Sets *mat* to the symplectic matrix `J = \Bigl(\begin{smallmatrix}
    0&I_g\\-I_g&0 \end{smallmatrix}\Bigr)`.

.. function:: void sp2gz_block_diag(fmpz_mat_t mat, const fmpz_mat_t U)

    Sets *mat* to the symplectic matrix `\Bigl(\begin{smallmatrix}
    U&0\\0&U^{-T} \end{smallmatrix}\Bigr)`. We require that `U\in
    \operatorname{GL}_g(\mathbb{Z})`.

.. function:: void sp2gz_trig(fmpz_mat_t mat, const fmpz_mat_t S)

    Sets *mat* to `\Bigl(\begin{smallmatrix} I_g&S\\0&I_g
    \end{smallmatrix}\Bigr)`, where *S* is a symmetric `g\times g` matrix.

.. function:: void sp2gz_embed(fmpz_mat_t res, const fmpz_mat_t mat)

    Assuming that *mat* is a symplectic matrix of size `2r\times 2r` and *res*
    is square of size `2g\times 2g` for some `g\geq r`, sets *res* to the symplectic matrix

        .. math::

            \begin{pmatrix} \alpha && \beta & \\ & I_{g-r} && 0_{g-r} \\ \gamma &&\delta &\\ & 0_{g-r} && I_{g-r} \end{pmatrix}

    where `\alpha,\beta,\gamma,\delta` are the `r\times r` blocks of *mat*.

.. function:: void sp2gz_restrict(fmpz_mat_t res, const fmpz_mat_t mat)

    Assuming that *mat* is a symplectic matrix of size `2g\times 2g` and *res*
    is square of size `2r\times 2r` for some `r\leq g`, sets *res* to the
    matrix whose `r\times r` blocks are the upper left corners of the
    corresponding `g\times g` block of *mat*. The result may not be a
    symplectic matrix in general.

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
    symplectic matrix, and store this matrix in *res*. Otherwise, returns
    false (0) and leaves *res* undefined.

.. function:: void sp2gz_inv(fmpz_mat_t inv, const fmpz_mat_t mat)

    Sets *inv* to the inverse of the symplectic matrix *mat*. In contrast with
    :func:`fmpz_mat_inv`, this involves no computation.

.. function:: fmpz_mat_struct * sp2gz_decompose(slong * nb, const fmpz_mat_t mat)

    Returns a vector *res* of symplectic matrices and store its length in *nb*
    such that the following holds: *mat* is the product of the elements of
    *res* from left to right, and each element of *res* is block-diagonal,
    trigonal, the `J` matrix, an embedded `J` matrix from a lower dimension, or
    an embedded matrix from dimension 1. The output vector *res* will need to
    be freed by the user as follows:

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

.. function:: void acb_siegel_cho_yinv(arb_mat_t cho, arb_mat_t yinv, const acb_mat_t tau, slong prec)

    Sets *yinv* to the inverse of the imaginary part `Y` of *tau*, and sets
    *cho* to an upper-triangular Cholesky matrix for `\pi Y`, i.e. to the
    upper-triangular matrix `C` with positive diagonal entries such that `\pi Y
    = C^T C`. If one cannot determine that `Y` is positive definite at the
    current working precision, *yinv* and *cho* are set to indeterminate
    matrices.

.. function:: void acb_siegel_reduce(fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Sets *mat* to a symplectic matrix such that `\mathit{mat}\cdot\tau` is as
    reduced as possible, repeatedly reducing the imaginary and real parts of
    `\tau` and applying fundamental symplectic matrices. If the coefficients of
    `\tau` do not have a reasonable size or if `\det \mathrm{Im}(\tau)` is
    vanishingly small, we simply set *mat* to the identity.

.. function:: int acb_siegel_is_reduced(const acb_mat_t tau, slong tol_exp, slong prec)

    Returns true (nonzero) iff it is certainly true that `\tau` belongs to the
    reduced domain defined by the tolerance parameter `\varepsilon =
    2^{\mathit{tol\_exp}}`. This means the following:
    `|\mathrm{Re}(\tau_{j,k})| < \frac12 + \varepsilon` for all `0\leq j,k <
    g`; the imaginary part of `\tau` passes :func:`arb_mat_spd_is_lll_reduced`
    with the same parameters; and for every matrix obtained from
    :func:`sp2gz_fundamental`, the determinant of the corresponding cocycle is
    at least `1-\varepsilon`.

.. function:: slong acb_siegel_kappa(acb_t sqrtdet, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)

    Returns `0\leq r < 8` such that `\kappa(\mathit{mat}) = \zeta_8^r` and sets
    *sqrtdet* to the corresponding square root of `\det(\gamma\tau + \delta)`
    in the theta transformation formula of [Igu1972]_, p. 176 and
    [Mum1983]_, p. 189: for any symplectic matrix `m`, any `(z,\tau)\in
    \mathbb{C}^g\times \mathbb{H}_g`, and any characteristic `(a,b)`, we have

        .. math::

            \theta_{a,b}(m\cdot(z,\tau)) = \kappa(m) \zeta_8^{e(m, a, b)} \det(\gamma\tau + \delta)^{1/2} e^{\pi i z^T (\gamma\tau + \delta)^{-1} \gamma z} \theta_{a',b'}(z,\tau)

    where

    - `\gamma,\delta` are the lower `g\times g` blocks of `m`,
    - `a',b'` is another characteristic depending on `m,a,b`,
    - `\zeta_8=\exp(i\pi/4)`,
    - `e(m,a,b)` is an integer given by an explicit formula in terms of `m,a,b`
      (this is `\phi_m` in Igusa's notation), and
    - `\kappa(m)` is an `8^{\mathrm{th}}` root of unity, only well-defined up
      to sign unless we choose a particular branch of `\det(\gamma\tau +
      \delta)^{1/2}` on `\mathbb{H}_g`.

    We proceed as follows. After applying :func:`sp2gz_decompose`, we only have
    to consider four special cases for *mat*. If *mat* is trigonal or
    block-diagonal, one can compute its action on `\theta_{0,0}` directly. If
    *mat* is an embedded matrix from `\mathrm{SL}_2(\mathbb{Z})`, we rely on
    :func:`acb_modular_theta_transform`. Finally, if *mat* is an embedded `J`
    matrix from dimension `0\leq r\leq g`, then `\kappa(m) \zeta_8^{e(m,0,0)}
    i^{r/2} \det(\tau_0)^{1/2} = 1`, where `\tau_0` denotes the upper left
    `r\times r` submatrix of `\tau` and the branch of the square root is chosen
    such that the result is `i^{g/2}\det(Y)` when `\tau = iY` is purely
    imaginary. To compute `\det(\tau_0)^{1/2}`, we pick a purely imaginary
    matrix *A* and consider the polynomial `P(t) = \det(A + \tfrac{t+1}{2}
    (\tau_0 - A))`. Up to choosing another `A`, we may assume that it has
    degree `g` and that its roots (as complex balls) do not intersect the
    segment `[-1,1]\subset \mathbb{C}`. We then find the correct branch of
    `P(t)^{1/2}` between `t=-1` and `t=1` following [MN2019]_.

.. function:: slong acb_siegel_kappa2(const fmpz_mat_t mat)

    Returns `0\leq r < 3` such that `\kappa(m)^2 = i^r`, which makes sense
    without reference to a branch of `\det(\gamma\tau + \delta)^{1/2}`.

    We adopt a similar strategy to :func:`acb_theta_transform_kappa`, but do
    not have to choose the correct square root of `\det(\tau_0)`.

.. function:: void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)

    Sets *tau* to a random matrix in `\mathbb{H}_g`, possibly far from being
    reduced.

.. function:: void acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)

    Sets *tau* to a random reduced matrix in `\mathbb{H}_g` whose imaginary
    part possibly has large entries.

.. function:: void acb_siegel_randtest_compact(acb_mat_t tau, flint_rand_t state, int exact, slong prec)

    Sets *tau* to a random reduced matrix in `\mathbb{H}_g` whose imaginary
    part has bounded entries. If *exact* is nonzero (true), then the entries of
    *tau* are set to exact (dyadic) complex numbers.

.. function:: void acb_siegel_randtest_vec(acb_ptr z, flint_rand_t state, slong g, slong prec)

    Sets *z* to a random vector of length *g*, possibly with large entries.

.. function:: void acb_siegel_randtest_vec_reduced(acb_ptr zs, flint_rand_t state, slong nb, const acb_mat_t tau, int exact, slong prec)

    Sets *zs* to the concatenation of *nb* random vectors *z* sampled from
    `[-1,1]^g + \tau[-1,1]^g`, i.e.~close to being reduced with respect to
    `tau`. If *exact* is nonzero (true), then the entries of *zs* are set to
    exact (dyadic) complex numbers.

Theta characteristics
-------------------------------------------------------------------------------

.. function:: void acb_theta_char_get_arb(arb_ptr v, ulong a, slong g)

.. function:: void acb_theta_char_get_acb(acb_ptr v, ulong a, slong g)

    Sets *v* to `a/2` seen as an element of `\mathbb{R}^g` or `\mathbb{C}^g`
    respectively.

.. function:: slong acb_theta_char_dot(ulong a, ulong b, slong g)

    Returns `\sum_{i=0}^{g-1} a_i b_i` modulo 4 as an integer between 0 and 3,
    where `a_i, b_i` for `0\leq i < g` denote the bits of `a` and `b`
    respectively.

.. function:: slong acb_theta_char_dot_slong(ulong a, const slong * n, slong g)

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

.. function:: void acb_theta_char_table(ulong * ch, slong * e, const fmpz_mat_t mat, slong ab)

    If *ab* encodes a valid characteristic, sets *ch* to the theta
    characteristic `(a',b')` and sets *e* to `e(\mathit{mat},a,b)` as in the
    transformation formula. If *ab* is negative, then sets *ch* and *e* to
    vectors of length `2^{2g}` containing this output for all characteristics
    from 0 to `2^{2g}-1`.

.. function:: void acb_theta_char_shuffle(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th, int sqr, slong prec)

    Partially applies the theta transformation formula to the given vector *th*
    for the symplectic matrix *mat* and stores the output in *res*. This omits
    the `\kappa`, determinant, and exponential factors from the formula. If
    *sqr* is nonzero (true), then replaces `\zeta_8` in the formula by `i` to
    mimic the transformation formula on squared theta values.

Ellipsoids
-------------------------------------------------------------------------------

The most direct way of evaluating Riemann theta functions consists in
evaluating a partial sum of the exponential series defining them, then adding
an error bound coming from the tail of the series. We refer to this strategy as
the *summation algorithms*.

The upper bound on the tail will be obtained from the triangle inequality. To
analyze the absolute value of each term in the sum defining
`\theta_{0,b}(z,\tau)`, we write:

    .. math::

        \bigl| \exp(i\pi n^T\tau n + 2n^T (z + \tfrac{b}{2}) \bigr| = \exp(\pi y^T Y^{-1} y) \exp (-\lVert n + Y^{-1}y \rVert_\tau^2)

where `Y` and `y` denote the imaginary parts of `tau` and `z` respectively, and
`\lVert\cdot\rVert_\tau` denotes the Euclidean norm attached to the Gram matrix
`\pi Y`. Thus, the exponential terms whose absolute values are less than a
given threshold correspond to lattice points `n\in \mathbb{Z}^g` lying outside
a certain ball centered in `v = -Y^{-1}y` for `\lVert\cdot\rVert_\tau`; in
other words, we should be computing partial sums over points `n\in
\mathbb{Z}^g` lying in certain ellipsoids, as in [DHBHS2004]_. (We use the
relation

    .. math::

        \theta_{a,b}(z,\tau) = \exp(\pi i a^T (z + \tfrac b2) + \pi i a^T\tau a/4) \theta_{0,b}(z + \tau\tfrac{a}{2},\tau)

to avoid summing over `\mathbb{Z}^g + \tfrac{a}{2}` with a nonzero `a`.)

This section gathers methods to manipulate such ellipsoids directly. Fix an
upper-triangular matrix `C` with positive diagonal entries (henceforth called a
"Cholesky matrix"), a radius `R\geq 0`, a vector `v\in \mathbb{R}^g`, and
`1\leq d\leq g`. Consider the ellipsoid `E` consisting of points `n =
(n_0,\ldots,n_{g-1})` satisfying `(v + Cn)^T(v + Cn)\leq R^2` and such that
their last coordinates `n_{d},\ldots, n_{g-1}` are fixed. We encode `E` as
follows: we store the endpoints and midpoint of the interval of allowed values
for `n_{d-1}` as :type:`slong`'s, and if `d\geq 1`, we store a
`(d-1)`-dimensional "child" of `E` for each value of `n_{d-1}` as another
ellipsoid in a recursive way. Children are partitioned between left and right
children depending on the position of `n_{d-1}` relative to the midpoint (by
convention, the midpoint is a right child). When `d=g` and for a fixed Cholesky
matrix `C`, this representation uses `O(R^{g-1})` space for an ellipsoid of
radius `R` containing approximately `O(R^{g})` points.

.. type:: acb_theta_eld_struct

.. type:: acb_theta_eld_t

    An :type:`acb_theta_eld_t` is an array of length one of type
    :type:`acb_theta_eld_struct` encoding an ellipsoid as described above,
    alllowing it to be passed by reference.

.. function:: void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)

    Initializes *E* as a *d*-dimensional ellipsoid in ambient dimension *g*. We
    require `1\leq d\leq g`.

.. function:: void acb_theta_eld_clear(acb_theta_eld_t E)

    Clears *E* as well as any recursive data contained in it.

.. function:: int acb_theta_eld_set(acb_theta_eld_t E, const arb_mat_t C, const arf_t R2, arb_srcptr v)

    Assuming that *C* is upper-triangular with positive diagonal entries,
    attempts to set *E* to represent an ellipsoid as defined above, where *R2*
    indicates `R^2`, and returns 1 upon success. If the ellipsoid points do not
    fit in :type:`slong`'s or if the ellipsoid is unreasonably large, returns 0
    instead and leaves *E* undefined.

The following functions are available after *E* has been initialized and then
computed using :func:`acb_theta_eld_init` and :func:`acb_theta_eld_set`.

.. function:: acb_theta_eld_nb_pts(E)

    Returns the number of points contained in `E`, which is stored in the data
    structure.

.. function:: void acb_theta_eld_points(slong * pts, const acb_theta_eld_t E)

    Sets *pts* to the list of all the points in `E`, as a concatenation of
    vectors of length *g*.

.. function:: acb_theta_eld_nb_border(E)

    Macro returning the number of points in the border of `E`, defined as
    follows. If `d=1`, then it consists of the two points with `n_0` equal to
    `m - 1` and `M + 1`, where `m` and `M` are the result of
    :macro:`acb_theta_eld_max` and :macro:`acb_theta_eld_min` respectively. If
    `d\geq 2`, then it is the reunion of the borders of all children of
    `E`. This number is stored in the data structure.

.. function:: void acb_theta_eld_border(slong * pts, const acb_theta_eld_t E)

    Sets *pts* to the list of all the points in the border of `E`. This is only
    used for testing.

.. function:: int acb_theta_eld_contains(const acb_theta_eld_t E, slong * pt)

    Returns true (nonzero) iff *pt* is contained in `E`. The vector *pt* must
    be of length *g*.

.. function:: void acb_theta_eld_print(const acb_theta_eld_t E)

    Prints a faithful description of `E`. This may be unwieldy in high
    dimensions.

.. function:: void acb_theta_agm_distances(arb_ptr ds, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

    Sets *ds* to the concatenation of the following *nb* vectors of length
    `2^g`: for each input vector `z`, we compute `\mathrm{Dist}_\tau(-Y^{-1}y,
    \mathbb{Z}^g + \tfrac a2)^2` for `a\in \{0,1\}^g`, where `y, Y` are the
    imaginary parts of `z, \tau` respectively and `\mathrm{Dist}_\tau` denotes
    the distance attached to `\lVert \cdot \rVert_\tau`.

Toolbox for derivatives
-------------------------------------------------------------------------------

In this module, we only consider the successive partial derivatives of
`\theta_{a,b}(z,\tau)` with respect to the `g` coordinates of `z`, because
derivatives with respect to `\tau` are accounted for by the heat equation

    .. math::

        \frac{\partial\theta_{a,b}}{\partial \tau_{j,k}} = \frac{1}{2\pi i(1 +\delta_{j,k})}
        \frac{\partial^2\theta_{a,b}}{\partial z_j \partial z_k}.

We encode tuples of derivation orders, henceforth called "derivation tuples",
as vectors of type :type:`slong` and length `g`. In agreement with
:ref:`acb_modular.h <acb-modular>`, we also normalize derivatives in the same way
as in the Taylor expansion, so that the tuple `(k_0,\ldots,k_{g-1})`
corresponds to the differential operator

    .. math::

        \frac{1}{k_0!}\cdots\frac{1}{k_{g-1}!} \cdot \frac{\partial^{|k|}}{\partial z_0^{k_0}\cdots \partial z_{g-1}^{k_{g-1}}},

where `|k|:=\sum k_i`. We always consider all derivation tuples up to a total
order *ord*, and order them first by their total order, then
reverse-lexicographically. For example, in the case `g=2`, the sequence of
orders is `(0,0)`, `(1,0)`, `(0,1)`, `(2,0)`, `(1,1)`, etc.

This sections gathers methods to work with vectors of partial derivatives of
holomorphic functions in general.

.. function:: slong acb_theta_jet_nb(slong ord, slong g)

    Returns the number of derivation tuples with total order at most *ord*. The
    result will be zero if *ord* is negative.

.. function:: slong acb_theta_jet_total_order(const slong * tup, slong g)

    Returns the total derivation order for the given tuple *tup* of length *g*.

.. function:: void acb_theta_jet_tuples(slong * tups, slong ord, slong g)

    Sets *tups* to the concatenation of all derivation tuples up to total order
    *ord*.

.. function:: slong acb_theta_jet_index(const slong * tup, slong g)

    Returns *n* such that *tup* is the `n^{\mathrm{th}}` derivation tuple of
    length *g*.

.. function:: void acb_theta_jet_mul(acb_ptr res, acb_srcptr v1, acb_srcptr v2, slong ord, slong g, slong prec)

    Sets *res* to the vector of derivatives of the product `fg`, assuming that
    *v1* and *v2* contains the derivatives of `f` and `g` respectively.

.. function:: void acb_theta_jet_compose(acb_ptr res, acb_srcptr v, const acb_mat_t N, slong ord, slong prec)

    Sets *res* to the vector of derivatives of the composition `f(Nz)`,
    assuming that *v* contains the derivatives of *f* at the point `Nz`.

.. function:: void acb_theta_jet_exp_pi_i(acb_ptr res, arb_srcptr a, slong ord, slong g, slong prec)

    Sets *res* to the vector of derivatives of the function `\exp(\pi i (a_0
    z_1 + \cdots + a_{g-1} z_{g-1}))` at `z = 0`, where `a_0,\ldots a_{g-1}` are
    the entries of *a*.

.. function:: void acb_theta_jet_exp_qf(acb_ptr res, acb_srcptr z, const acb_mat_t N, slong ord, slong prec)

    Sets *res* to the vector of derivatives of the function `\exp(\pi i z^T N
    z)` at the chosen point `z`.

Context structures in summation algorithms
-------------------------------------------------------------------------------

After the relevant ellipsoid has been computed, summation algorithms only
involve exponential terms in `\tau` and `z`. Sometimes, especially in the
setting of the quasi-linear algorithms below, these exponentials can be
computed once and for all at high precision, then used for several calls to the
summation functions. This section introduces context structures to make these
manipulations easier.

.. type:: acb_theta_ctx_tau_struct

.. type:: acb_theta_ctx_tau_t

    An :type:`acb_theta_ctx_tau_t` is an array of length one of type
    :type:`acb_theta_ctx_tau_struct` containing all the necessary data to run
    the summation algorithm on a given matrix `\tau\in\mathcal{H}_g`. In
    particular, it contains a matrix `exp_tau_div_4` whose `(j,k)` entry (when
    `j\leq k`) is `\exp(\pi i (2 - \delta_{j,k}) \tau_{j,k}/4)`, where `\delta`
    denotes the Kronecker symbol, as well as a Cholesky matrix for `\pi Y` if
    `g>1`.

.. type:: acb_theta_ctx_z_struct

.. type:: acb_theta_ctx_z_t

    An :type:`acb_theta_ctx_z_t` is an array of length one of type
    :type:`acb_theta_ctx_tau_struct` containing all the necessary data to run
    the summation algorithm on a given vector `z` (provided that an element of
    type :type:`acb_theta_ctx_tau_t` is also given.) In particular, it contains
    the values `\exp(2\pi i z_j)` for all `1\leq j\leq g`.

.. function:: void acb_theta_ctx_tau_init(acb_theta_ctx_tau_t ctx, int allow_shift, slong g)

    Initializes *ctx* for use in dimension *g*. If *allow_shift* is nonzero
    (true), then additional fields in *ctx* are initialized to allow for the
    evaluation of theta functions `\theta_{a,0}` for nonzero `a`.

.. function:: void acb_theta_ctx_tau_clear(acb_theta_ctx_tau_t ctx)

    Clears *ctx*.

.. function:: void acb_theta_ctx_z_init(acb_theta_ctx_z_t ctx, slong g)

    Initializes *ctx* for use in dimension *g*.

.. function:: void acb_theta_ctx_z_clear(acb_theta_ctx_z_t ctx)

    Clears *ctx*.

.. function:: acb_theta_ctx_z_struct * acb_theta_ctx_z_vec_init(slong nb, slong g)

    Returns a pointer to a vector of *nb* initialized elements of type
    :type:`acb_theta_ctx_z_struct`.

.. function:: void acb_theta_ctx_z_vec_clear(acb_theta_ctx_z_struct * vec, slong nb)

    Clears the elements of type :type:`acb_theta_ctx_z_struct` pointed to by
    *vec* as well as the pointer itself.

.. function:: void acb_theta_ctx_exp_inv(acb_t exp_inv, const acb_t exp, const acb_t x, int is_real, slong prec)

    Given a complex value *x* and given *exp* containing `\exp(\pi i x)`, sets
    *exp_inv* to `\exp(-\pi i x)`. This is computed by complex conjugation from
    *exp* if *is_real* is nonzero (true). Otherwise, it is computed by
    inverting *exp*, except if the result is indeterminate, in which case we
    recompute *exp_inv* from *x* directly.

.. function:: void acb_theta_ctx_sqr_inv(acb_t sqr_inv, const acb_t inv, const acb_t sqr, int is_real, slong prec)

    Given *inv* and *sqr* containing complex values `\exp(-\pi i x)` and
    `\exp(2\pi i x)` respectively, sets *sqr_inv* to `\exp(-2\pi i x)`. This
    uses complex conjugation from *sqr* if *is_real* is nonzero (true), and
    otherwise a complex squaring from *inv*.

.. function:: void acb_theta_ctx_tau_set(acb_theta_ctx_tau_t ctx, const acb_mat_t tau, slong prec)

    Computes and stores in *ctx* the required data for the input matrix
    `\tau`. The dimensions must match.

.. function:: void acb_theta_ctx_tau_dupl(acb_theta_ctx_tau_t ctx, slong prec)

    Modifies *ctx* in place to correspond to the matrix `2\tau` instead of
    `\tau`. This is much cheaper than calling :func:`acb_theta_ctx_tau_set`
    again.

.. function:: void acb_theta_ctx_z_set(acb_theta_ctx_z_t ctx, acb_srcptr z, const acb_theta_ctx_tau_t ctx_tau, slong prec)

    Computes and stores in *ctx* the required data for the complex vector
    *z*. Here *ctx_tau* should contain context data for the matrix `\tau`. The
    dimensions must match.

.. function:: void acb_theta_ctx_z_dupl(acb_theta_ctx_z_t ctx, slong prec)

    Modifies *ctx* in place to correspond to the pair `(2z,2\tau)` instead of
    `(z,\tau)`. This is much cheaper than calling :func:`acb_theta_ctx_z_set`
    again.

.. function:: void acb_theta_ctx_z_add_real(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx, const acb_theta_ctx_z_t ctx_real, slong prec)

    Assuming that *ctx* and *ctx_real* correspond to pairs `(z,\tau)` and `(t,
    \tau)` respectively where `t` is a real vector, sets *res* to a valid
    context for the pair `(z + t,\tau)`.

.. function:: void acb_theta_ctx_z_shift_a0(acb_theta_ctx_z_t res, acb_t c, const acb_theta_ctx_z_t ctx, const acb_theta_ctx_tau_t ctx_tau, ulong a, slong prec)

    Assuming that *ctx* and *ctx_tau* correspond to a pair `(z,\tau)`, and that
    *allow_shift* was set to true when computing *ctx_tau*, sets *res* to a
    valid context for the pair `(z + \tau \tfrac{a}{2},\tau)` and sets `c` to
    the complex value such that for all `0\leq b\leq 2^g-1`,

    .. math::

        \theta_{a,b}(z,\tau) = c \exp(\pi i a^T b/2) \theta_{0,b}(z + \tau\tfrac{a}{2},\tau).

    We have `c = \exp(\pi i a^T z + \pi i a^T\tau a/4)`.

.. function:: void acb_theta_ctx_z_common_v(arb_ptr v, const acb_theta_ctx_z_struct * vec, slong nb, slong prec)

    Given a vector *vec* of valid contexts for pairs
    `(z_1,\tau),\ldots,(z_n,\tau)`, sets *v* to a valid ellipsoid center for
    use in :func:`acb_theta_eld_set` when running the summation algorithm for
    all these pairs, constructed using :func:`arb_union`.

.. function:: int acb_theta_ctx_z_overlaps(const acb_theta_ctx_z_t ctx1, const acb_theta_ctx_z_t ctx2)

    Returns true iff the data contained in *ctx1* and *ctx2* overlap in the
    sense of :func:`acb_overlaps`. This is only used for testing.

Summation algorithms
-------------------------------------------------------------------------------

In this module, summation algorithms are mainly used for low to moderate
precisions due to their higher complexity (except in special cases such as
:func:`acb_theta_00`). Nevertheless, running summation algorithms at low
precisions is a key step in the quasi-linear algorithms. They have been
optimized in many ways and should already compare very favorably to other
software packages that evaluate theta functions.

We always assume in this section that the vectors `z` have been reduced. This
allows us to use only one ellipsoid when several vectors `z` are given.

In order to compute the correct ellipsoids for a target working precision, we
use the following upper bound on the tail of the series: by [EK2025]_, for any
`v\in \mathbb{R}^g` and any upper-triangular Cholesky matrix `C`, any
`\mathit{ord}\geq 0`, and any `R` such that `R^2 \geq\max(4,\mathit{ord})`, we
have

    .. math::

        \sum_{n\in C\mathbb{Z}^g + Cv,\ \lVert n\rVert^2 \geq R^2} \lVert n\rVert^{\mathit{ord}} e^{-\lVert n\rVert^2}
        \leq 2^{2g+2} R^{g-1+p} e^{-R^2} \prod_{j=0}^{g-1} (1 + \gamma_j^{-1})

where `\gamma_0,\ldots, \gamma_{g-1}` are the diagonal coefficients of
`C`.

After the relevant ellipsoid *E* has been computed, the main worker inside each
version of the summation algorithm will process one line (i.e. 1-dimensional
ellipsoid) in *E*. Before calling this worker, for fixed `\tau` and `z` and
fixed coordinates `n_1,\ldots n_{g-1}` defining a line inside the ellipsoid, if
`n_{\mathrm{min}}` are `n_{\mathrm{max}}` are the endpoints of the interval of
allowed values for `n_0`, we (efficiently) compute:

- the vector `v_1` with entries `\exp(\pi i j^2 \tau_{0,0})` for
  `n_{\mathrm{min}}\leq j\leq n_{\mathrm{max}}`,
- the vector `v_2` with entries `x^j` for `n_{\mathrm{min}}\leq j\leq
  n_{\mathrm{max}}`, where

    .. math::

        x = \exp(2 \pi i z_0) \prod_{k = 1}^{g-1} \exp(2 \pi i n_k \tau_{0,k}),

- the cofactor `c\in \mathbb{C}` given by

    .. math::

        c = \prod_{k = 1}^{g-1} \exp(2 \pi i n_k z_k) \cdot
        \prod_{1\leq j\leq k < g} \exp(\pi i (2 - \delta_{j,k}) n_j n_k \tau_{j,k}).

This allow us to use :func:`acb_dot` in the workers while maintaining
reasonable memory costs, and to use an average of strictly less than two
complex multiplications per lattice point as `R\to \infty`. Moreover, these
multiplications are performed at only a fraction of the full precision for
lattice points far from the ellipsoid center. Different versions of the
summation algorithm will rely on slightly different workers, so introducing a
function pointer type is helpful to avoid code duplication.

When `g=1`, the code does not rely on these worker functions, and calls
:func:`acb_modular_theta_sum` from :ref:`acb_modular.h <acb-modular>` instead.

.. type:: acb_theta_sum_worker_t

    A function pointer type. A function *worker* of this type has the
    following signature:

    .. function:: void worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2, const slong * precs, slong len, const acb_t c, const slong * coords, slong ord, slong g, slong prec, slong fullprec)

    where:

    - *th* denotes the output vector of theta values to which terms will be added,
    - *v1*, *v2* and *c* are precomputed as above,
    - *precs* contains working precisions for each term `n_{\mathrm{min}}\leq
      j\leq n_{\mathrm{max}}`,
    - *len* `= n_{\mathrm{max}} - n_{\mathrm{min}} + 1` is the common length of
      *v1*, *v2* and *precs*,
    - *coords* is `(n_{\mathrm{min}}, n_1, \ldots, n_{g-1})`,
    - *ord* is the maximal derivation order,
    - *prec* is the working precision for this line inside the ellipsoid, and
      finally
    - *fullprec* is the working precision for summing into *th*.

.. function:: void acb_theta_sum_radius(arf_t R2, arf_t eps, const arb_mat_t cho, slong ord, slong prec)

    Sets *R2* and *eps* such that the above upper bound for *R2* and the given
    *ord* is at most *eps*, where `C` is *cho*. We choose *eps* so that the
    relative error on the output of the summation algorithm should be roughly
    `2^{-\mathit{prec}}` if no cancellations occur in the sum, i.e.
    `\mathit{eps} \simeq 2^{-\mathit{prec}} \prod_{j=0}^{g-1} (1 +
    \gamma_j^{-1})`.

.. function:: void acb_theta_sum_jet_radius(arf_t R2, arf_t eps, const arb_mat_t cho, arb_srcptr v, slong ord, slong prec)

    Assuming that `C = \mathit{cho}` is the upper-triangular Cholesky matrix
    for `\pi Y` and `v = C Y^{-1} y`, returns *R2* and *eps* so that, when
    summing the above series on terms `n\in \mathbb{Z}^g` such that `(v + C
    n)^T(v + C n)\leq R^2`, the absolute value of the tail of the series
    (before multiplying by the leading factor `(2\pi i)^{|k|} e^{\pi y^T Y^{-1}
    y}`, see below) will be bounded above by *eps*, for any derivation tuple
    `k` with `|k|\leq \mathit{ord}`.

    We can rewrite the the differentiated series

        .. math::

            \frac{\partial^{|k|}\theta_{a,b}}{\partial z_0^{k_0}\cdots \partial z_{g-1}^{k_{g-1}}}(z,\tau) = (2\pi i)^{|k|} \sum_{n\in \mathbb{Z}^g + \tfrac a2} n_0^{k_0} \cdots n_{g-1}^{k_{g-1}}
            e^{\pi i n^T \tau n + 2\pi i n^T (z + \tfrac b2)}.

    as

        .. math::

            (2\pi i)^{|k|} e^{\pi y^T Y^{-1} y} \sum_{n\in \mathbb{Z}^g + \tfrac a2} n_0^{k_0} \cdots n_{g-1}^{k_{g-1}} e^{\pi i(\cdots)} e^{-\pi (n + Y^{-1}y)^T Y (n + Y^{-1}y)}.

    We ignore the leading multiplicative factor. Writing `m = C n + v`, we have

        .. math::

            n_0^{k_0}\cdots n_{g-1}^{k_{g-1}}\leq
            (\lVert C^{-1}\rVert_\infty \lVert n\rVert_2 + \lVert Y^{-1}y\rVert_\infty)^{|k|}.

    Using the upper bound from :func:`acb_theta_sum_radius`, we see that the
    absolute value of the tail of the series is bounded above by

        .. math::

            (\lVert C^{-1} \rVert_\infty R + \lVert Y^{-1}y \rVert_\infty)^{|k|}
             2^{2g+2} R^{g-1} e^{-R^2} \prod_{j=0}^{g-1} (1 + \gamma_j^{-1}).

    Thus, we proceed as follows. We first compute *R2* and *eps* using
    :func:`acb_theta_sum_radius` with *ord* = 0. If `R\leq \lVert
    Y^{-1}y\rVert_\infty/\lVert C^{-1}\rVert_\infty`, we simply multiply *eps*
    by `\max\{1, 2 \lVert Y^{-1}y \rVert\}^{\mathit{ord}}`. Otherwise, we
    compute *R2* and *eps* using :func:`acb_theta_sum_radius` with the given
    value of *ord*. We can then set *R2* to the maximum of *R2* and `\lVert
    Y^{-1}y \rVert_\infty /\lVert C^{-1} \rVert_\infty`, and multiply *eps* by
    `\max\{1, 2\lVert C^{-1}\rVert\}^{\mathit{ord}}`.

.. function:: void acb_theta_sum_term(acb_t res, acb_srcptr z, const acb_mat_t tau, slong * tup, slong * n, slong prec)

    Sets *res* to `n_0^{k_0} \cdots n_{g-1}^{k_{g-1}}\exp(\pi i(n^T\tau n + 2
    n^Tz))`, where the `k_j` and `n_j` denotes the `j^{\mathrm{th}}` entry in
    *tup* and *n* respectively. The vector *tup* may be *NULL*, which is
    understood to mean the zero tuple. This is only used for testing.

.. function:: void acb_theta_sum_work(acb_ptr th, slong len, acb_srcptr exp_zs, acb_srcptr exp_zs_inv, slong nb, const acb_mat_t exp_tau, const acb_mat_t exp_tau_inv, const acb_theta_eld_t E, slong ord, slong prec, acb_theta_sum_worker_t worker)

    Runs the summation algorithm on the ellipsoid *E*, assuming `g\geq 2`. The
    input related to `\tau` is the matrices *exp_tau* and *exp_tau_inv* whose
    `(j,k)` entries for each `1\leq j\leq k\leq g` should contain `\exp(\pi i
    (2 - \delta_{j,k}) \tau_{j,k})` and its inverse, respectively. The input
    related to the *nb* vectors `z` is the vectors `exp_zs` and `exp_zs_inv`
    which should contain the values `\exp(\pi i z_j)` (for each `z` and each
    `1\leq j\leq g`) and their inverses, respectively. The parameters *len*,
    *ord* and the output vector *th* are passed to *worker* when processing
    each individual line in *E*. The data associated with *zs* and `\tau` is
    typically stored in contexts of type :type:`acb_theta_ctx_tau_t` and
    :type:`acb_theta_ctx_z_t` respectively. No error bound coming from the tail
    of the theta series is added.

.. function:: void acb_theta_sum_00(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb, const acb_theta_ctx_tau_t ctx_tau, slong prec)

.. function:: void acb_theta_sum_0b(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb, const acb_theta_ctx_tau_t ctx_tau, slong prec)

    Evaluates either `\theta_{0,0}(z, \tau)`, or alternatively `\theta_{0,b}(z,
    \tau)` for each `b\in \{0,1\}^g`, for each pair `(z,\tau)` corresponding to
    a context stored in *vec* together with *ctx_tau*. The result *th* will be
    a concatenation of *nb* vectors of length `1` or `2^g` respectively. The
    associated worker performs one :func:`acb_dot` operation.

.. function:: void acb_theta_sum_a0_tilde(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb, const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, slong prec)

    Evaluates `\widetilde{\theta_{a,0}}(z,\tau)` for each `a\in \{0,1\]^g` at
    each of the *nb* pairs `(z,\tau)` corresponding to a context stored in
    *vec* together with *ctx_tau*. The result *th* will be a concatenation of
    *nb* vectors of length `2^g`.

    In this function, the absolute error radius we add on
    `\widetilde{\theta_{a,0}}(z,\tau)` from the tail of the exponential series
    depend on `a`. The amount of precision added is controlled by *distances*,
    which could be computed as in :func:`acb_theta_agm_distances` (although
    other values sometimes make sense, such as 0.) Since this vector is the
    same for all vectors *z*, this function makes the most sense when the
    different values of *z* differ by real vectors.

.. function:: void acb_theta_sum_all_tilde(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb, const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, slong prec)

    Evaluates `\widetilde{\theta}_{a,b}(z,\tau)` for all characteristics `a,b`
    and all *nb* pairs `(z,\tau)` specified by the contexts. The precision used
    depends on *a* and the vector *distances* exactly as in
    :func:`acb_theta_sum_a0_tilde`.

.. function:: void acb_theta_sum_jet_00(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb, const acb_theta_ctx_tau_t ctx_tau, slong ord, slong prec)

    Sets *th* to the vector of derivatives of `\theta_{0,0}` at the given *nb*
    points `(z,\tau)` up to total order *ord*.

.. function:: void acb_theta_sum_jet_all(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb, const acb_theta_ctx_tau_t ctx_tau, slong ord, slong prec)

    Sets *th* to the vector of derivatives of `\theta_{a,b}` for all
    characteristics `a,b`, at all the given *nb* points `(z,\tau)`, up to total
    order *ord*. The result will be a concatenation of *nb* vectors (one for
    each *z*), each piece being the concatenation of `2^{2g}` vectors of
    derivatives.

AGM steps
-------------------------------------------------------------------------------

The quasi-linear algorithm to evaluate theta functions uses the following
*duplication formula*: for all `z,z'\in \mathbb{C}^g` and `\tau\in
\mathbb{H}_g`,

    .. math::

        \theta_{a,0}(z,\tau) \theta_{a,0}(z',\tau) = \sum_{a'\in(\mathbb{Z}/2\mathbb{Z})^g}
        \theta_{a',0}(z+z',2\tau) \theta_{a+a',0}(z-z',2\tau).

In particular,

    .. math::

        \begin{aligned}
        \theta_{a,0}(z,\tau)^2 &= \sum_{a'\in (\mathbb{Z}/2\mathbb{Z})^g}
        \theta_{a',0}(2z,2\tau) \theta_{a+a',0}(0,2\tau),\\
        \theta_{a,0}(0,\tau)\theta_{a,0}(z,\tau) &= \sum_{a'\in(\mathbb{Z}/2\mathbb{Z})^g}
        \theta_{a',0}(z,2\tau) \theta_{a+a',0}(z,2\tau), \\
        \theta_{a,0}(0,\tau)^2 &= \sum_{a'\in (\mathbb{Z}/2\mathbb{Z})^g}
        \theta_{a',0}(0,2\tau) \theta_{a+a',0}(0,2\tau).
        \end{aligned}

Applying one of these duplication formulas amounts to taking a step in a
(generalized) AGM sequence. Note that the duplication formula still holds if we
replace `\theta_{a,0}` by the normalized version
`\widetilde{\theta}_{a,0}`.

This section gathers methods to apply such duplication formulas efficiently
while minimizing precision losses. In the case `z = z'`, the duplication
formula is typically followed by extracting square roots, using low-precision
computations as a guide to make the correct choice.

.. function:: void acb_theta_agm_sqrt(acb_ptr res, acb_srcptr a, acb_srcptr roots, slong nb, slong prec)

    Sets the `k^{\mathrm{th}}` entry `r_k` of *res* for `0\leq k < \mathit{nb}`
    to a square root of the corresponding entry `a_k` of `a`. The choice of
    sign is determined by *rts*: each `r_k` will overlap the corresponding
    entry of *rts* but not its opposite. Exceptional cases are handled as
    follows: if both square roots of `a_k` overlap *rts*, then `r_k` is set to
    their :func:`acb_union`; if none ovelaps *rts*, then `r_k` is set to an
    indeterminate value.

.. function:: void acb_theta_agm_mul(acb_ptr res, acb_srcptr a1, acb_srcptr a2, slong g, slong prec)

    For each `0\leq k < 2^g`, sets the `k^{\mathrm{th}}` entry of *res* to
    `2^{-g}\sum_{b\in \{0,1\}^g} a_{1,b}\, a_{2, b + k}`, where addition is
    meant in `(\mathbb{Z}/2\mathbb{Z}^g)` (a bitwise xor).

    Following [LT2016]_, we apply the Hadamard matrix twice with
    multiplications in-between. This causes precision losses when the absolute
    values of the entries in *a1* and/or *a2* are of different orders of
    magnitude. This function is faster when *a1* and *a2* are equal as
    pointers, as we can use squarings instead of multiplications.

.. function:: slong acb_theta_agm_addprec(const arb_t d)

    Returns an integer that is close to `d/\log(2)` if *d* is
    finite and of reasonable size, and otherwise returns 0.

.. function:: void acb_theta_agm_mul_tight(acb_ptr res, acb_srcptr a0, acb_srcptr a, arb_srcptr d0, arb_srcptr d, slong g, slong prec)

    Assuming that *d0* and *d* are obtained as the result of
    :func:`acb_theta_agm_distances` on `(0,\tau)` and `(z,\tau)` respectively,
    performs the same computation as :func:`acb_theta_agm_mul` on the vectors
    *a0* and *a* with a different management of error bounds. The resulting
    error bounds on *res* will be tighter when the absolute value of `a_k` is
    roughly `e^{-d_k}` for each `0\leq k < 2^g`, and similarly for *a0* and
    *d0*, for instance when applying the duplication formula on modified theta
    values.

    When `g>1`, we manage error bounds as follows. We compute `m, \varepsilon`
    such that the following holds: for each `0\leq k < \mathit{nb}`, if `d_k`
    (resp. `a_k`) denotes the `k^{\mathrm{th}}` entry of *d* (resp. *a*), then
    the absolute value of `a_k` is at most `m \cdot e^{-d_k}` and the radius of
    the complex ball `a_k` is at most `\mathit{eps}\cdot e^{-d_k}`. We proceed
    similarly on *a0* and *d0* to obtain `m_0, \varepsilon_0`. Then we call
    :func:`acb_theta_agm_mul` on the midpoints of *a0* and *a* at a higher
    working precision, and finally add `e^{-d_k} (m_0 \varepsilon + m
    \varepsilon_0 + \varepsilon\varepsilon_0)` to the error bound on the
    `k^\mathrm{th}` entry of *res*. This is valid for the following reason:
    keeping notation from :func:`acb_theta_agm_distances`, for each `b\in
    \{0,1\}^g`, the sum

        .. math::

            \mathrm{Dist}_\tau(-Y^{-1}y, \mathbb{Z}^g + \tfrac b2)^2
            + \mathrm{Dist}_\tau(-Y^{-1} y, \mathbb{Z}^g + \tfrac{b + k}{2})^2

    is at most `\mathrm{Dist}_\tau(-Y^{-1}y, \mathbb{Z}^g + \tfrac{k}{2})^2` by
    the parallelogram identity.

Quasilinear algorithms on exact, reduced input
-------------------------------------------------------------------------------

Suppose that we wish to compute `\theta_{a,0}(0,\tau)` for all `a\in \{0,1\}^g`
and a reduced matrix `\tau\in \mathbb{H}_g`. Applying the last of the above
duplication formulas `n` times, we reduce to evaluating
`\theta_{a,0}(0,2^n\tau)`. We expect that the absolute value of this complex
number is roughly `\exp(-d^2)` for `d = 2^n\mathrm{Dist}_\tau(0, \mathbb{Z}^g +
\tfrac a2))`, where `\mathrm{Dist}_\tau` denotes the distance in `\mathbb{R}^g`
attached to the quadratic form `\mathrm{Im}(\tau)`. Provided that `n \simeq
\log_2(\mathit{prec})`, we have to sum only `O_g(1)` terms in the summation
algorithm to evaluate `\theta_{a,0}(0,2^n\tau)` at "shifted absolute precision"
*prec*, i.e. absolute precision `\mathit{prec} + d^2/\log(2)`.

In order to recover `\theta_{a,0}(0,\tau)`, we then perform `n` AGM
steps. *Assuming* that each `|\theta_{a,0}(0, 2^k\tau)|` is indeed of the
expected order of magnitude, we can ensure that the precision loss is `O_g(1)`
bits at each step in terms of shifted absolute precision, and we can make the
correct choices of square roots at each step by computing low-precision
approximations with the summation algorithm. However, depending on the choice
of `\tau`, this assumption may not always hold.

We make the following adjustments to make the algorithm work for all `\tau`, as
well as for theta values at `z\neq 0`:

- If we discover (after applying the summation algorithm) that some value
  `\theta_{a,0}(0,2^k\tau)` is too small, we introduce an auxiliary real vector
  `t`. At each step, starting from `\theta_{a,0}(0,2^{k+1}\tau)`,
  `\theta_{a,0}(2^{k+1}t, 2^{k+1}\tau)` and `\theta_{a,0}(2^{k+2}t,
  2^{k+1}\tau)`, we compute `\theta_{a,0}(2^{k}t, 2^k\tau)` and
  `\theta_{a,0}(2^{k+1}t, 2^k\tau)` using square roots (second formula above),
  then `\theta_{a,0}(0, 2^k\tau)` using divisions (third formula). For a huge
  majority of such `t`, none of the values `\theta_{a,0}(2^kt, 2^k\tau)` and
  `\theta_{a,0}(2^{k+1}t, 2^k\tau)` will be too small. In practice, we choose
  `t` at random and obtain a probabilistic algorithm with a negligible failure
  probability.

- When computing `\theta_{a,0}(z,\tau)` for a nonzero `z`, we compute
  `\theta_{a,0}(0, 2^k\tau)` and `\theta_{a,0}(2^k z, 2^k\tau)` using the
  second and fourth formulas at each step. We actually replace each occurrence
  of `\theta_{a,0}(z,\tau)` by `e^{-\pi y^T Y^{-1} y}\theta_{a,0}(z,\tau)`, as
  the absolute values of the latter quantities do not increase as `y` gets
  farther from zero, and they still satisfy the duplication formulas.

- These two techniques can be combined by evaluating theta values at the six
  vectors `2^k v` for `v = 0, t, 2t, z, z + t, z + 2t`. Note that we only have
  to compute `\theta_{a,0}(2^kz, 2^k\tau)` at the last step `k=0`.

- Finally, if the eigenvalues of `\mathrm{Im}(\tau)` have different orders of
  magnitude, then the ellipsoids from the summation algorithm will become very
  thin in one direction while still being thick in other directions. In such a
  case, we can rewrite the theta values as a sum of `O(1)` theta values in
  lower dimensions. This increases the efficiency of the algorithm while
  ensuring that the absolute precisions we consider are always in
  `O(\mathit{prec})`.

Finally, we note that the duplication formulas also have analogues for all
theta values, not just `\theta_{a,0}`: for instance, we have

    .. math::

        \theta_{a,b}(0,\tau)^2 = \sum_{a'\in (\mathbb{Z}/2\mathbb{Z})^g} (-1)^{a'^Tb}
        \theta_{a',0}(0,2\tau)\theta_{a+a',0}(0,2\tau).

We use those generalized formulas for the very last duplication step when
needed.

We always assume in this section that the inputs `(z,\tau)` are provided as
exact dyadic numbers and that the pairs `(z,\tau)` have been reduced.

.. function:: int acb_theta_ql_nb_steps(slong * pattern, const acb_mat_t tau, slong prec)

    Determines how many duplication steps we should perform to evaluate theta
    functions at `\tau` at precision *prec*, and at which steps we should fall
    back to lower dimensions, if any.

    The output is stored in *pattern*, a vector of length `g`. Roughly
    speaking, the `j`-th entry of *pattern* is a nonnegative integer `m` such
    that `2^m c_j^2` is of the order of *prec*, where `c_j` denotes the `j`-th
    diagonal coefficient of the Cholesky matrix attached to `pi Y`. In other
    words, the ellipsoid we need to consider in the summation algorithm at
    `2^m\tau` has width `O(1)` in the direction of the `j`-th canonical basis
    vector in `\mathbb{R}^g`. Because `\tau` is assumed to be reduced, we
    expect *pattern* to be a roughly decreasing vector.

    If some entries of the Cholesky matrix are interminate or too extreme for a
    reasonable `m` to be computed, then the output is 0, and otherwise 1.

    Modifying this function is the main way to tune the behavior and
    performance of the quasi-linear algorithms to evaluate theta functions.

.. function:: int acb_theta_ql_lower_dim(acb_ptr * new_zs, acb_ptr * cofactors, slong ** pts, slong * nb, arf_t err, slong * fullprec, acb_srcptr z, const acb_mat_t tau, arb_srcptr distances, slong s, ulong a, slong prec)

    Implements the dimension-lowering strategy for evaluating theta functions.

    The input is as follows:
    - `(z,\tau)` should be an exact element of `\mathbb{C}^g\times
      \mathcal{H}_g` (ideally reduced)
    - *distances* should be the output of :func:`acb_theta_agm_distances` on
      this pair
    - *s* should be an integer between `1` and `g-1`; we will reduce the
      evaluation of theta functions from dimension `g` to dimension `s`
    - *a* should be an integer between `0` and `2^{g-s}-1` included; we will
      only decompose `\theta_{a',0}(z,\tau)` when the last `g - s` bits of `a'`
      correspond to those of *a*.

    We then proceed as follows:
    - *fullprec* is set to the binary precision at which those theta values
      `\theta_{a',0}(z,\tau)` should be computed. We take distances into
      account, so *fullprec* is *prec* plus additional guard bits derived from
      the maximum of the entries in *distances* corresponding to the
      characteristics *a'*.
    - *R2* and *err* are set as in :func:`acb_theta_sum_radius` for this choice
      of *fullprec*. (Note that *R2* is not part of the output.) Thus,
      `\theta_{a',0}(z,\tau)` can be obtain by summing over an ellipsoid of
      squared radius *R2* and adding an error *err* coming from the tail.
    - Let `n\in \mathbb{Z}^g + \tfrac{a'}{2}` be a point in that
      ellipsoid. Write `a' = (a_0,a)` and `n = (n_0,n_1)` where `n_0\in
      \mathbb{Z}^s + \tfrac{a_0}{2}` and `n_1\in \mathbb{Z}^{g - s} +
      \tfrac{a}{2}`. By the Pythagorean theorem, the possible values for `n_1`
      all lie in an ellipsoid of radius *R2* in dimension `g-s`, whose Cholesky
      matrix is the lower-right part of a Cholesky matrix for `\pi Y`. This
      ellipsoid is meant to contain very few points. We list all possible
      values for `n_1 - \tfrac{a}{2}` (which lies in `\mathbb{Z}^g`) in *pts*,
      and set *nb* to the number of those points. Note that *pts* will have to
      be freed by the user afterwards.
    - For each `n_1 - \tfrac{a}{2}` listed in *pts*, then the sum of the
      corresponding terms in the theta series is

        .. math::

            e^{\pi i \bigl(n_1^T \tau_1 n_1 + 2 n_1^T z_1\bigr)}
            \theta_{a_0,0}(z_0 + x n_1, \tau_0).

      where `\tau = \Bigl(\begin{smallmatrix} \tau_0 & x\\x^T &
      \tau_1\end{smallmatrix}\Bigr)` and `z = (z_0,z_1)`. Thus, we allocate
      *new_zs* to contain *nb* vectors of length `g` and set the corresponding
      entry to `z_0 + x n_1` (which is still exact). We also allocate
      *cofactors* to be a vector of length *nb* and set the corresponding entry
      to the above exponential. Both *new_zs* and *cofactors* will have to be
      freed by the user.

.. function:: void acb_theta_ql_recombine(acb_ptr th, acb_srcptr th0, acb_srcptr cofactors, const slong * pts, slong nb, const arf_t err, slong fullprec, slong s, ulong a, int all, slong g, slong prec)

    Performs the converse to :func:`acb_theta_ql_lower_dim`, namely recovers
    theta values `\theta_{a',0}(z,\tau)` from the output of
    :func:`acb_theta_ql_lower_dim` and theta values in dimension `s`.

    The input is as follows:
    - *cofactors*, *pts*, *nb*, *err*, *fullprec*, *s*, *a*, *g*, *prec* should
      be as output by :func:`acb_theta_ql_lower_dim`.
    - If *all* is true (nonzero), then `th0` should be a concatenation of *nb*
      vectors of length `2^{2s}` containing `\theta_{a_0,b_0}(z_0,\tau_0)` for
      all characteristics `(a_0,b_0)` in dimension `s`, where `z_0` runs
      through *new_zs* as output by :func:`acb_theta_ql_lower_dim`, and
      `\tau_0` is defined as above. If *all* is false (zero), then *th0* should
      be a concatenation of *nb* vectors of length `2^{s}` containing
      `\theta_{a_0,0}(z_0,\tau_0)` only.

    The output, stored in *th*, is either the vector containing
    `\theta_{a,b}(z,\tau)` for all `g`-dimensional characteristics `(a,b)` (if
    *all* is true) or only `\theta_{a,0}(z,\tau)` for all `a` (if *all* is
    false), where `(z,\tau)` were the initial input to
    :func:`acb_theta_ql_lower_dim`.

.. function:: int acb_theta_ql_setup(acb_ptr rts, acb_ptr rts_all, acb_ptr t, slong * guard, slong * easy_steps, acb_srcptr zs, slong nb, const acb_mat_t tau, arb_srcptr distances, slong nb_steps, int all, slong prec)

    Sets up the structure of AGM steps to evaluate theta functions at the given
    *nb* pairs `(z,\tau)` where `z` runs through *zs*, which are assumed to be
    exact and reduced, using *nb_steps* duplication steps. The parameters *nb*
    and *nb_steps* must be at least one, and *zs* must begin with the zero
    vector in `\mathbb{C}^g`.

    The rest of the input is as follows:
    - *distances* should be the concatenation of *nb* vectors of length `2^g`
      computed by :func:`acb_theta_agm_distances` for each pair `(z,\tau)`.
    - *nb_steps* should be the number of times we wish to apply the duplication
      formulas before falling back to either the summation algorithms or the
      dimension-lowering strategy.
    - if *all* is nonzero (true), then we will compute `\theta_{a,b}(z,\tau)`
      for all characteristics `(a,b)`, and otherwise only
      `\theta_{a,0}(z,\tau)`.

    The vectors *rts*, *rts_all*, *t*, and *easy_steps* should be
    preinitialized with lengths `2^g \times 3\times
    \mathit{nb}\times\mathit{nb_steps}`, `2^{2g}\times\mathit{nb}` (only used
    if *all* is true), `g` and *nb* respectively, while *guard* is only a
    pointer to an integer.

    We proceed as follows. Initially, *guard* is set to a small value such as 16.
    1. For each `z`, we use the summation algorithms to obtain approximations
      of `\theta_{a,b}(z,\tau)` (if *all* is true) or `\theta_{a,0}(z,\tau)`
      (if *all* is false), and store them in *rts_all* or *rts*
      respectively. We adjust the choice of precision in terms of *distances*
      and add *guard* bits, so that the computed approximations do not contain
      zero with a good probability. If none of the computed approximations
      contains zero, it means that we can successfully apply the last (and
      simplest) duplication formula for the last step of the quasi-linear
      algorithm. In that case, we go on and compute approximations of
      `\theta_{a,0}(2^k z,2^k\tau)`, for `k = 1,2,` etc., up to *nb_steps*-1 or
      until one of the approximations we compute contains zero, taking
      distances into account at each step. We store the computed values in
      *rts*, and set the corresponding entry of *easy_steps* to be the number
      of steps for which the simplest duplication formula can be applied.
    2. At that point, if the entries of *easy_steps* are all equal to
      *nb_steps*, we are done. Otherwise, we pick an auxiliary vector `t` at
      random. The 1st entry of *easy_steps*, corresponding to `z=0`, is set to
      the minimal value in *easy_steps* (this is necessary to be able to apply
      the duplication formulas.)
    3. For each `z`, if `m` denotes the corresponding entry of *easy_steps*, we
      use the summation algorithms to compute approximations of
      `\theta_{a,0}(2^k (z + t), 2^k\tau)` and `\theta_{a,0}(2^k(z + 2t),
      2^k\tau)` for each `k` between `m` and *nb_steps*-1 included at low
      precision. If one of these values contains zero, we restart step 3 with
      another `t` (we allow a small number of retries, such as 4). We store
      those approximations in *rts*. If `k=0` and *all* is true, then the
      values we need are `\theta_{a,b}(z+2t,\tau)` for all `(a,b)` instead;
      those are stored in *rts_all*.
    4. If no suitable `t` was found in step 4, then we double *guard* and go
      back to step 1. We allow this until *guard* reaches *prec*. After that,
      if `t` still cannot be found, then we declare failure and output 0. This
      should only happen with negligible probability for well-formed input. The
      output value if 1 if a suitable `t` was found.

.. function:: void acb_theta_ql_steps(acb_ptr th, acb_ptr th_init, acb_srcptr rts, acb_srcptr rts_all, slong nb, slong nb_steps, arb_srcptr distances, const slong * easy_steps, int all, slong g, slong prec)

    Performs AGM steps in the context of the quasi-linear algorithm for theta
    functions. The parameters *nb* and *nb_steps* must be at least one, and
    *zs* must begin with the zero vector.

    The input is as follows:
    - *rts*, *rts_all*, *nb*, *nb_steps*, *distances*, *easy_steps*, *all*, *g*
      and *prec* are as input to/output by :func:`acb_theta_ql_setup`
    - *th_init* should be a vector of length *nb* times 3 times `2^{2g}`
      containing the initial theta values needed in the duplication
      algorithm. (These may differ depending on the content of *easy_steps*.)

    The output is either the collection of theta values
    `\widetilde{\theta}_{a,b}(z,\tau)` for all `a,b` or
    `\widetilde{\theta}_{a,0}` for all `a` (depending on whether *all* is true
    or not) for each vector `z` in *zs*, and is stored in *th*. The precise
    duplication steps taken depends on the contents of *easy_steps*.

.. function:: int acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, const slong * pattern, int all, int shifted_prec, slong prec)

    Runs the full quasi-linear algorithm to evaluate theta functions at the
    given *nb* pairs `(z,\tau)` where `z` runs through *zs*, which are assumed
    to be exact and reduced.

    The output is either the collection of theta values
    `\widetilde{\theta}_{a,b}(z,\tau)` for all `a,b` or
    `\widetilde{\theta}_{a,0}` for all `a` (depending on whether *all* is true
    or not) for each vector `z` in *zs*, and is stored in *th*. If
    *shifted_prec* is nonzero (true), then the precision to which these values
    are computed will take distances into account similarly to
    :func:`acb_theta_sum_a0_tilde`.

    The input *pattern* conditions how many duplication steps will be performed
    and whether to apply the dimension-lowering strategy. If no duplication
    steps are needed at all, we call :func:`acb_theta_sum_a0_tilde` or
    :func:`acb_theta_sum_all_tilde` directly (using phony 0 values as
    *distances* if *shifted_prec* is false). Otherwise, we call
    :func:`acb_theta_ql_setup`, which we expect to succeed with overwhelming
    probability. The initial theta values required for the final call to
    :func:`acb_theta_ql_steps` are computed either by the summation algorithms
    or, if the dimension-lowering strategy is used, using
    :func:`acb_theta_ql_lower_dim`, making a recursive call to
    :func:`acb_theta_ql_exact` in a lower dimension but (possibly) a longer
    list of vectors *zs*, and finally recombining the values with
    :func:`acb_theta_ql_recombine`.

Main functions on reduced input
-------------------------------------------------------------------------------

This section wraps up the quasi-linear and summation algorithms on inputs that
are reduced, but not necessarily exact.

.. function:: void acb_theta_sum_bound(arb_t c, arb_t rho, acb_srcptr z, const acb_mat_t tau, slong ord)

    Sets *c* and *rho* such that on every ball centered at (a point contained
    in) *z* of radius *rho*, the functions `|\theta_{a,b}|` for all
    characteristics `(a,b)` are uniformly bounded by `c`. The choice of *rho*
    is tuned to get interesting upper bounds on derivatives of `\theta_{a,b}`
    up to order *ord*.

    We proceed as follows. First, we compute `c_0`, `c_1`, `c_2` such that for
    any choice of `\rho`, one can take `c = c_0\exp((c_1 + c_2\rho)^2)`
    above. We can take

        .. math::

            c_0 = 2^g \prod_{j=0}^{g-1} (1 + 2\gamma_j^{-1}),

        .. math::

            c_1 = \sqrt{\pi y^T Y^{-1} y},

        .. math::

            c_2 = \sup_{\lVert x \rVert_\infty\leq 1} \sqrt{\pi x^T Y^{-1} x}.

    One can easily compute an upper bound on `c_2` from the Cholesky
    decomposition of `\pi Y^{-1}`. We then look for a value of `\rho` that
    minimizes `\exp((c_1 + c_2\rho)^2)/\rho^{2m-1}` where `m = \mathit{ord}+1`,
    i.e. we set `\rho` to the positive root of `2c_2\rho (c_1 + c_2\rho) =
    2m-1`.

.. function:: void acb_theta_jet_error(arb_ptr err, acb_srcptr z, const acb_mat_t tau, acb_srcptr dth, slong ord, slong prec)

    Assuming that *dth* contains the derivatives of a theta function
    `\theta_{a,b}` up to total order `\mathit{ord} + 2`, sets *err* to a vector
    with the following property. Let `(z_0,\tau_0)` be the midpoint of
    `(z,\tau)`, and let `(z_1,\tau_1)` be any point inside the ball specified
    by the given *z* and *tau*. Then the vectors of derivatives of
    `\theta_{a,b}` at `(z_0,\tau_0)` and `(z_1,\tau_1)` up to total order *ord*
    differ by at most *err* elementwise.

.. function:: void acb_theta_00_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

    Same as :func:`acb_theta_00`, but does not attempt to reduce the input
    under the Siegel modular group.

.. function:: void acb_theta_one_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, ulong ab, slong prec);

    Same as :func:`acb_theta_00_notransform`, but evalues `\theta_{a,b}` for a
    fixed characteristic instead of `\theta_{0,0}`.

.. function:: void acb_theta_all_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, int sqr, slong prec)

    Same as :func:`acb_theta_all`, but does not attempt to reduce the input
    under the Siegel modular group.

.. function:: void acb_theta_jet_00_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, slong prec)

    Same as :func:`acb_theta_jet_00`, but does not attempt to reduce the input
    under the Siegel modular group.

.. function:: void acb_theta_jet_one_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, ulong ab, slong prec)

    Same as :func:`acb_theta_jet_00_notransform`, but evaluates the derivatives
    of `\theta_{a,b}` for a fixed characteristic instead of derivatives of
    `\theta_{0,0}`.

.. function:: void acb_theta_jet_all_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, slong prec)

    Same as :func:`acb_theta_jet_all`, but does not attempt to reduce the input
    under the Siegel modular group.

Reduction and main functions
-------------------------------------------------------------------------------

The only building blocks we still need to assemble :func:`acb_theta_all`,
:func:`acb_theta_all`, :func:`acb_theta_jet_00` and :func:`acb_theta_jet_all`
are the following reduction functions.

.. function:: int acb_theta_reduce_tau(acb_ptr new_zs, acb_mat_t new_tau, fmpz_mat_t mat, acb_mat_t N, acb_mat_t ct, acb_ptr exps, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

    Reduces the input matrix `\tau` under the action of the modular group
    `\mathrm{Sp}_{2g}(\mathbb{Z})` and modifies the *nb* input vectors stored
    in *zs* according to the theta transformation formula.

    The output is as follows:
    - *new_tau* is the reduced matrix,
    - *mat* is the symplectic matrix such that *new_tau* is the result of *mat*
      acting on *tau*,
    - *ct* is the transpose of `(\gamma \tau + \delta)^{-1}` where
      `\gamma,\delta` are the lower `g\times g` blocks of *mat*,
    - *N* is the matrix `i \pi \gamma \cdot \mathit{ct}` that appears in the
      transformation formula,
    - and finally *new_zs* is the list of *nb* vectors in `\mathbb{C}^g`
      corresponding to the elements in *zs* multiplied by *ct* on the left.

    If reduction was unsuccessful (usually indicating that the input is
    malformed or that the working precision is insufficient to detect that `Y`
    is positive definite), the return value is 0 and the above output is left
    undefined. Otherwise, the return value is 1.

.. function:: int acb_theta_reduce_z(acb_ptr new_zs, arb_ptr rs, acb_ptr cs, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

    Reduces the *nb* vectors stored in *zs* in the context of evaluating theta
    functions with *tau* as a second argument.

    The output vectors *new_zs*, *rs* and *cs* should have lengths *nb* times
    `g`, *nb* times `g`, and *nb* respectively. For a given vector *z*
    appearing in *zs*, we round `Y^{-1}y` to the nearest even, integral vector
    `r`, and store it in *rs*. Then, we consider the vector `z - \tau r`,
    substract the nearest even integral vector from its real part, and store
    the result `z'` into *new_zs*. For all characteristics `a,b`, we have

        .. math::

            \theta_{a,b}(z',\tau) &= e^{- i \pi r^T (z + z')} \theta_{a,b}(z,\tau).

    Finally, we store this exponential as the corresponding entry of *cs*.

    If rounding the imaginary part to integers does not succeed due to extreme
    values, then the return value is 0 and the output vectors are left
    undefined. Otherwise, the return value is 1.

Let us now detail how the algorithm for :func:`acb_theta_jet_all_notransform`
based on finite differences works. Consider the Taylor expansion:

    .. math::

        \theta_{a,b}(z + h, \tau)
        = \sum_{k\in \mathbb{Z}^g,\ k\geq 0} a_k\, h_0^{k_0}\cdots h_{g-1}^{k_{g-1}}.

If one chooses `h = h_n = (\varepsilon \zeta^{n_0},\ldots, \varepsilon
\zeta^{n_{g-1}})` where `\varepsilon > 0` and `\zeta` is a primitive
`m^{\mathrm{th}}` root of unity and lets `n` run through all vectors in
`\{0,\ldots, m - 1\}^g`, then taking a discrete Fourier transform of the
resulting values will compute the individual Taylor coefficient for each
derivation tuple that is bounded by `m-1` elementwise. (A constant proportion,
for fixed `g`, of this set consists of all tuples of total order at most
`m-1`.) More precisely, fix `p\in \mathbb{Z}^g`. Then

    .. math::

        \sum_{n\in \{0,\ldots,m-1\}^g} \zeta^{-p^T n} \theta_{a,b}(z + h_n, \tau)
        = m^g \sum_{\substack{k\in \mathbb{Z}^g,\ k\geq 0,\\ k = p\ (\text{mod } m)}}
        a_k\,\varepsilon^{|k|}.

We obtain an upper bound on the tail of this series from the Cauchy integration
formula: if `|\theta_{a,b}(z,\tau)|\leq c` uniformly on a ball of radius `\rho`
centered in `z` for `\lVert\cdot\rVert_\infty`, then the sum is `m^g
(a_p\,\varepsilon^{|p|} + T)` with

    .. math::

        |T|\leq 2c g\,\frac{\varepsilon^{|p|+m}}{\rho^m}.

Since we divide by `\varepsilon^{|p|}` to get `a_p`, we will add an error of
`2c g \varepsilon^m/\rho^{m+|p|}` to the result of the discrete Fourier
transform.

The algorithm first decides whether the summation algorithm
:func:`acb_theta_sum_jet_all` should be used directly. (Right now, we use
direct summation when less than 3 duplication steps would be used; this could
be optimized more carefully in the future.) If not, it then computes `c` and
`\rho` using :func:`acb_theta_sum_bound`, chooses a suitable `\varepsilon`,
strips `(z,\tau)` of their error bounds, increases the working precision to
account for division by `\varepsilon^{\mathit{ord}}\cdot (\mathit{ord}+1)^g`,
calls :func:`acb_theta_all_notransform` on all the auxiliary points, performs
the discrete Fourier transforms, and finally restores provably correct error
bounds on the results using :func:`acb_theta_jet_error` and derivatives to
order *ord* + 2 computed at low precision. This algorithm runs in quasi-linear
time in `\mathit{prec}\cdot \mathit{ord}^{\,g}` for any fixed `g`.

Dimension 2 specifics
-------------------------------------------------------------------------------

In the `g=2` case, one can use theta functions to evaluate many fundamental
Siegel modular forms. This section contains methods to do so, in analogy with
:func:`acb_modular_delta` or :func:`acb_modular_eisenstein` when `g=1`.

We use the following notation. Fix `k,j\geq 0`. A Siegel modular form of weight
`\det^k\otimes \mathrm{Sym}^j` is by definition an analytic function
`f: \mathbb{H}_g\to \mathbb{C}_j[X]` (the vector space of polynomials of degree
at most `j`) such that for any `\tau\in \mathbb{H}_g` and
`m\in \mathrm{Sp}_4(\mathbb{Z})`, we have

    .. math::

        f((\alpha\tau + \beta)(\gamma\tau + \delta)^{-1}) = \det(\gamma\tau +
        \delta)^k\cdot \mathrm{Sym}^j(\gamma\tau + \delta)(f(\tau)).

Here `\alpha,\beta,\gamma,\delta` are the `g\times g` blocks of `m`, and the
notation `\mathrm{Sym}^j(r)` where `r = \bigl(\begin{smallmatrix} a & b\\ c &
d\end{smallmatrix}\bigr)\in \mathrm{GL}_2(\mathbb{C})` stands for the map

    .. math::

        P(X) \mapsto (b X + d)^j P\bigl(\tfrac{a X + c}{b X + d}\bigr).

For a nonzero `f` to exist, `j` must be even.

Siegel modular forms generate a bi-graded ring which is not finitely
generated. However, if we relax the definition of a Siegel modular form and
allow them to have a pole along the diagonal `\mathbb{H}_1^2 =
\bigl\{\bigl(\begin{smallmatrix} \tau_1 & 0 \\ 0 &
\tau_2\end{smallmatrix}\bigr)\bigr\}\subset \mathbb{H}_2` of a certain order
(depending on the weight), we indeed find a finitely generated ring
corresponding to classical "covariants" of a binary sextic. Historically,
covariants are classified in terms of their degree `k` and index `j`,
corresponding to Siegel modular functions of weight `\det^{k - j/2}\otimes
\mathrm{Sym}^j`. See [CFG2017]_ for more details on the correspondence between
modular forms and covariants.

.. macro:: ACB_THETA_G2_COV_NB

    Macro giving the number of generators of the ring of covariants, equal to 26.

.. function:: void acb_theta_g2_detk_symj(acb_poly_t res, const acb_mat_t m, const acb_poly_t f, slong k, slong j, slong prec)

    Sets *res* to `\det(m)^k \mathrm{Sym}^j(m)(f)`. The polynomial `f` should
    be of degree at most `j` (any coefficients of larger degree are ignored).

.. function:: void acb_theta_g2_transvectant(acb_poly_t res, const acb_poly_t g, const acb_poly_t h, slong m, slong n, slong k, slong prec)

    Sets *res* to the `k^{\mathrm{th}}` transvectant of the polynomials `g` and
    `h` of degrees `m` and `n`: considering `g` and `h` as homogeneous
    polynomials of degree `m` (resp. `n`) in `x_1,x_2`, this sets *res* to

        .. math::

            (g,h)_k := \frac{(m-k)!(n-k)!}{m!n!}  \sum_{j=0}^{k} (-1)^{k-j} \binom{k}{j}
            \frac{\partial^k g}{\partial x_1^{k-j}\partial x_2^j}
            \frac{\partial^k h}{\partial x_1^{j}\partial x_2^{k-j}}.

    Any coefficients of `g` or `h` of larger degree than `m` (resp. `n`) are
    ignored.

.. function:: void acb_theta_g2_transvectant_lead(acb_t res, const acb_poly_t g, const acb_poly_t h, slong m, slong n, slong k, slong prec)

    Sets *res* to the leading coefficient of `(g,h)_k` in `x_1`, with the same
    conventions as in :func:`acb_theta_g2_transvectant`.

.. function:: slong acb_theta_g2_character(const fmpz_mat_t mat)

    Returns the value in `\mathbb{Z}/2\mathbb{Z}` (0 or 1) of the unique
    nontrivial character of `\mathrm{Sp}_4(\mathbb{Z})` at *mat*, following
    [CFG2019]_, §12.

.. function:: void acb_theta_g2_psi4(acb_t res, acb_srcptr th2, slong prec)

.. function:: void acb_theta_g2_psi6(acb_t res, acb_srcptr th2, slong prec)

.. function:: void acb_theta_g2_chi10(acb_t res, acb_srcptr th2, slong prec)

.. function:: void acb_theta_g2_chi12(acb_t res, acb_srcptr th2, slong prec)

    Sets *res* to the value of the Eisenstein series `\psi_4`, `\psi_6` or the
    cusp forms `\chi_{10}, \chi_{12}` corresponding to the given vector *th2* of
    squared theta values (of length 16).

    We use the formulas from §7.1 in [Str2014]_, with the following normalizations:

        .. math::

            \psi_4 = h_4/4, \quad \psi_6 = h_6/4,\quad \chi_{10} = -2^{-12} h_{10},
            \quad \chi_{12} = 2^{-15}h_{12}.

    We warn that `\chi_{10}` and `\chi_{12}` differ from the classical notation
    of Igusa [Igu1979]_ by scalar factors. Writing `\tau =
    \bigl(\begin{smallmatrix} \tau_1 & \tau_2 \\ \tau_2 &
    \tau_3\end{smallmatrix}\bigr)` and `q_j = \exp(2\pi i \tau_j)`, the Fourier
    expansions of these modular forms begin as follows:

        .. math::

            \begin{aligned} \psi_4(\tau) &= 1 + 240(q_1 + q_3) + \cdots\\
            \psi_6(\tau) &= 1 - 504(q_1 + q_3) + \cdots\\
            \chi_{10}(\tau) &= (q_2 - 2 + q_2^{-1}) q_1 q_3 + \cdots\\
            \chi_{12}(\tau) &= (q_2 + 10 + q_2^{-1}) q_1 q_3 + \cdots.
            \end{aligned}

.. function:: void acb_theta_g2_chi5(acb_t res, acb_srcptr th, slong prec)

    Sets *res* to the value of `\chi_5 = - 2^{-6} \prod_{(a,b)\text{ even}}
    \theta_{a,b}` corresponding to the given theta values *th*. The form
    `\chi_5` is a Siegel cusp form with character: see [CFG2019]_ for more
    details.

.. function:: void acb_theta_g2_chi35(acb_t res, acb_srcptr th, slong prec)

    Sets *res* to the value of the cusp form `\chi_{35}` corresponding to the vector
    of theta values *th*. The form `\chi_{35}` is the unique scalar-valued Siegel
    modular form of weight `\det^{35}\otimes \mathrm{Sym}^0` up to scalars, and is
    normalized as follows:

        .. math::

            \chi_{35}(\tau) = q_1^2 q_3^2 (q_1 - q_3 )(q_2 - q_2^{-1}) + \cdots

    An explicit formula for `\chi_{35}` in terms of theta values is given in
    [Bol1887]_. See also [Mum1984]_, Prop. 6.2 p. 98 for how to translate
    Bolza's notation in terms of theta characteristics.

.. function:: void acb_theta_g2_chi3_6(acb_poly_t res, acb_srcptr dth, slong prec)

    Sets *res* to the value of the vector-valued cusp form with character
    `\chi_{6,3}` of weight `\det^3\otimes \mathrm{Sym}^6` corresponding to the
    given values of *dth*, computed as in e.g. :func:`acb_theta_jet_all` with
    `\mathit{ord}=1`. We have by [CFG2017]_:

        .. math::

            \chi_{3,6}(\tau) = \frac{1}{64\pi^6} \prod_{(a,b) \text{ odd}}
            \left(\frac{\partial \theta_{a,b}}{\partial z_1}(0,\tau) x_1 +
            \frac{\partial\theta_{a,b}}{\partial z_2}(0,\tau) x_2\right).

.. function:: void acb_theta_g2_sextic_chi5(acb_poly_t res, acb_t chi5, const acb_mat_t tau, slong prec)

    Sets *res* and *chi5* to the values of `\chi_{-2,6}` and `\chi_5` at
    `\tau`. We reduce `\tau` to the Siegel fundamental domain, call
    :func:`acb_theta_jet_all_notransform`, then apply the transformation
    formula for Siegel modular forms (which is simpler than the transformation
    formula on derivatives of theta functions.) Under the correspondence
    between Siegel modular functions and covariants of binary sextics,
    `\chi_{-2,6}` corresponds to the binary sextic itself, hence the name.

.. function:: void acb_theta_g2_sextic(acb_poly_t res, const acb_mat_t tau, slong prec)

    Same as :func:`acb_theta_g2_sextic_chi5`, but does not output *chi5*.

.. function:: void acb_theta_g2_covariants(acb_poly_struct * res, const acb_poly_t f, slong prec)

    Sets *res* to the vector of 26 generators of the ring of covariants
    evaluated at the sextic *f* (any terms of degree `>6` are ignored), in the
    following order:

    0. `C_{1,6}=f`
    1. `C_{2,0}= 60(f,f)_6`
    2. `C_{2,4}= 75(f,f)_4`
    3. `C_{2,8}= 90(f,f)_2`
    4. `C_{3,2}= 30(f,C_{2,4})_4`
    5. `C_{3,6}= 30(f,C_{2,4})_2`
    6. `C_{3,8}= 6(f,C_{2,4})_1`
    7. `C_{3,12}= 6 (f,C_{2,8})_1`
    8. `C_{4,0}= 2 (C_{2,4},C_{2,4})_4`
    9. `C_{4,4}= 30 (f,C_{3,2})_2`
    10. `C_{4,6}= 6(f,C_{3,2})_1`
    11. `C_{4,10}= 2(C_{2,8},C_{2,4})_1`
    12. `C_{5,2}=(C_{2,4},C_{3,2})_2`
    13. `C_{5,4}=\frac 25 (C_{2,4},C_{3,2})_1`
    14. `C_{5,8}= 2(C_{2,8},C_{3,2})_1`
    15. `C_{6,0}= 2(C_{3,2},C_{3,2})_2`
    16. `C_{6,6}^{(1)}= \frac 25(C_{3,6},C_{3,2})_1`
    17. `C_{6,6}^{(2)}= \frac 83(C_{3,8},C_{3,2})_2`
    18. `C_{7,2}= 30(f,C_{3,2}^2)_4`
    19. `C_{7,4}= 12(f,C_{3,2}^2)_3`
    20. `C_{8,2}= \frac 25(C_{2,4},C_{3,2}^2)_3`
    21. `C_{9,4}= 4(C_{3,8},C_{3,2}^2)_4`
    22. `C_{10,0}= 20(f,C_{3,2}^3)_6`
    23. `C_{10,2}= \frac 65(f,C_{3,2}^3)_5`
    24. `C_{12,2}= \frac 85(C_{3,8},C_{3,2}^3)_6`
    25. `C_{15,0}= \frac{1}{30000} (C_{3,8},C_{3,2}^4)_8`.

    The scalar factors are chosen so that when evaluated at a formal sextic `f
    = \sum a_i x_1^{6-i}x_2^i`, the covariants are integral and primitive as
    multivariate polynomials in `a_0,\ldots,a_6,x_1,x_2`.

.. function:: void acb_theta_g2_covariants_lead(acb_ptr res, const acb_poly_t f, slong prec)

    Sets *res* to the vector of leading coefficients in `x_1` of the 26
    covariants evaluated at *f*. This is more efficient than taking leading
    coefficients of :func:`acb_theta_g2_covariants`, since we can use
    :func:`acb_theta_g2_transvectant_lead` instead of
    :func:`acb_theta_g2_transvectant`.

Tests
-------------------------------------------------------------------------------

.. code-block:: bash

    ./build/acb_theta/test/main sp2gz_set_blocks

Generates a random `2g\times 2g` matrix, calls :func:`sp2gz_set_blocks` on its
four `g\times g` windows, and checks that the result equals the original
matrix.

.. code-block:: bash

    ./build/acb_theta/test/main sp2gz_is_correct

Checks that the return value of :func:`sp2gz_is_correct` is 1 on matrices
generated by :func:`sp2gz_j`, :func:`sp2gz_block_diag`, :func:`sp2gz_trig` and
:func:`sp2gz_fundamental`, and 0 on the identity matrix if it is not square of
even size.

.. code-block:: bash

    ./build/acb_theta/test/main sp2gz_inv

Checks that the result of :func:`sp2gz_inv` agrees with :func:`fmpz_mat_inv` on
random input.

.. code-block:: bash

    ./build/acb_theta/test/main sp2gz_decompose

Checks that the result of :func:`sp2gz_decompose` on random input only consists
of symplectic matrices of the allowed types, and that their product equals the
original matrix.

.. code-block:: bash

    ./build/acb_theta/test/main acb_siegel_cocycle

Checks that the chain rule holds: if `m'' = m'm` is a product of two symplectic
matrices and `\tau\in \mathbb{H}_g`, then `\gamma''\tau + \delta'' =
(\gamma'\tau' + \delta')(\gamma\tau+\delta)` where `\tau' = m\tau`. These
quantities are computed using :func:`acb_siegel_cocycle` and
:func:`acb_siegel_transform`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_siegel_transform

Checks that the chain rule holds, i.e. :func:`acb_siegel_transform` defines an
action of the group `\mathrm{Sp}_{2g}(\mathbb{Z})` on `\mathbb{H}_g`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_siegel_transform_z

Checks that :func:`acb_siegel_transform` and :func:`acb_siegel_transform_z`
agree on random input, and that :func:`acb_siegel_transform_z` on the inverse
of any matrix yields the inverse transformation.

.. code-block:: bash

    ./build/acb_theta/test/main acb_siegel_reduce

Generates an input matrix `\tau` at a working precision that is not too low
compared to the size of its coefficients, and calls :func:`acb_siegel_reduce`.
Checks that the resulting matrix `m` is symplectic and that `m\tau` is reduced
with a tolerance of `2^{-10}` using :func:`acb_siegel_is_reduced`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_siegel_is_reduced

Checks that :func:`acb_siegel_is_reduced` returns 1 on the matrix `i I_g`, but
0 on other matrices specially constructed to not be reduced.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_char_get_a

Generates a random characteristic *a*, sets *n* to the result of
:func:`acb_theta_char_get_slong` on *a*, and checks that the result of
:func:`acb_theta_char_get_a` on *n* gives back *a*.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_char_dot

Checks that dot products computed by :func:`acb_theta_char_dot`,
:func:`acb_theta_char_dot_slong` and :func:`acb_theta_char_dot_acb` agree on
random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_char_is_even

Checks that the 10 even theta characteristics for `g=2` are 0, 1, 2, 3, 4, 6,
8, 9, 12, 15.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_char_is_goepel

Checks that there are exactly 15 Göpel quadruples for `g=2`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_char_is_syzygous

Checks that there are exactly 60 syzygous triples for `g=2`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_eld_points

Generates a random ellipsoid *E* using :func:`acb_theta_eld_set`, computes its
points using :func:`acb_theta_eld_points`, and checks that each of these points
lies within the box specified by :macro:`acb_theta_eld_box`. Then, generates
random points *pt*: if *pt* is in *E* according to
:func:`acb_theta_eld_contains`, then *pt* must appear in the list of points,
otherwise the norm of *pt* according to the chosen Cholesky matrix must be at
least the radius of *E*.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_eld_border

Generates a random ellipsoid *E*, computes its border using
:func:`acb_theta_eld_border`, and checks that none of these border points lie
in *E* nor any of its children.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_naive_radius

Generates a reduced matrix `\tau` in `\mathbb{H}_g` and vector `z\in
\mathbb{C}^g`, calls :func:`acb_theta_naive_radius`, constructs the associated
ellipsoid *E*, and checks that the sums of absolute values of terms of the
theta series on the border of *E* is at most the specified bound.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_naive_reduce

Checks that the results of :func:`acb_theta_naive_reduce` are sound on some
special values of the input, namely when *zs* has only real entries and when
`\mathrm{Im}(z) = -\mathrm{Im}(\tau) n + \varepsilon` where *n* is an even
integral vector and `\varepsilon` is small.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_naive_term

Checks that the result of :func:`acb_theta_naive_term` is `n^k
\exp(i\pi(n^2\tau + 2nz))` in the `g=1` case.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_naive_00

Checks that the ouput of :func:`acb_theta_naive_00` overlaps the first entry of
the output of :func:`acb_theta_naive_0b`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_naive_all

Checks that the results of :func:`acb_theta_naive_all` agree with
:func:`acb_modular_theta` as follows: if the input matrix `\tau` is diagonal
with coefficients `\tau_0,\ldots, \tau_{g-1}`, then for all characteristics
`(a,b)` and vectors `z`, we have

    .. math::

       \theta_{a,b}(z,\tau) = \prod_{j=0}^{g-1} \theta_{a_j,b_j}(z_j,\tau_j).

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_naive_fixed_a

Checks that the output of :func:`acb_theta_naive_fixed_a` overlaps the relevant
entries of :func:`acb_theta_naive_all` on random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_naive_fixed_ab

Checks that the output of :func:`acb_theta_naive_fixed_ab` overlaps the relevant
entries of :func:`acb_theta_naive_all` on random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_tuples

For random *g* and *ord*, generates the list of derivation tuples using
:func:`acb_theta_jet_tuples`, picks an index `i` at random, and checks that the
result of :func:`acb_theta_jet_index` on the `i^{\mathrm{th}}` tuple is indeed
`i`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_mul

Checks that the results of :func:`acb_theta_jet_mul` agrees with the result of
:func:`fmpz_mpoly_mul` on any input with integral entries.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_compose

Checks that the chain rule holds: if `N_3 = N_2 N_1`, then applying
:func:`acb_theta_jet_compose` with `N_2`, then `N_1` corresponds to applying
:func:`acb_theta_jet_compose` with `N_3` directly.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_naive_radius

Generates a reduced matrix `\tau` in `\mathbb{H}_g` and vector `z\in
\mathbb{C}^g`, chooses a random order of derivation, calls
:func:`acb_theta_jet_naive_radius`, constructs the associated ellipsoid *E*,
and checks that the sums of absolute values of terms of the differentiated
theta series on the border of *E* is at most the specified bound.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_naive_all

Checks that the results of :func:`acb_theta_jet_naive_all` agree with
:func:`acb_modular_theta_jet` as follows: if the input matrix `\tau` is
diagonal with coefficients `\tau_0,\ldots, \tau_{g-1}`, then for all
characteristics `(a,b)`, any vector `z`, and any derivation tuple
`(k_0,\ldots,k_{g-1})`, we have

    .. math::

       \frac{\partial^{|k|} \theta_{a,b}} {\partial z_0^{k_0}\cdots \partial
       z_{g-1}^{k-1}}(z,\tau) = \prod_{j=0}^{g-1}
       \frac{\partial^{k_j}\theta_{a_j,b_j}}{\partial z^{k_j}}(z_j,\tau_j).

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_naive_00

Checks that the output of :func:`acb_theta_jet_naive_00` agrees with the
relevant entries of :func:`acb_theta_jet_naive_all` on random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_naive_fixed_ab

Checks that the output of :func:`acb_theta_jet_naive_fixed_ab` agrees with the
relevant entries of :func:`acb_theta_jet_naive_all` on random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_error_bounds

Generates two pairs `(z_1,\tau_1)` and `(z_2,\tau_2)` close to each other but
not overlapping, sets `(z,\tau)` to be their reunion (as complex balls on each
coefficient), and calls :func:`acb_theta_jet_error_bounds` on `(z,\tau)` for
some choice of derivation order. The difference between the results of
:func:`acb_theta_jet_naive_all` on `(z_1,\tau_1)` and `(z_2,\tau_2)` must then
be at most two times the computed error.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_dist_pt

Checks that for a random Cholesky matrix `C` and integral vectors `n_1,n_2`,
the results of :func:`acb_theta_dist_pt` on `(v,n) = (Cn_1, n_2)` and `(Cn_2,
n_1)` agree.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_dist_lat

Picks a random Cholesky matrix `C` and vector `v`, calls
:func:`acb_theta_dist_lat`, and computes the ellipsoid *E* whose radius is the
computed distance. Checks that *E* contains at least one point and that the
minimum distance is correct by looping over all the points in *E*.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_dist_a0

Checks that when `z = \mathrm{Im}(\tau) \tfrac{a}{2}` for some theta
characteristic `a`, the result of :func:`acb_theta_dist_a0` on `(z,\tau)`
contains zero in its `a^{\mathrm{th}}` entry.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_agm_hadamard

Checks that calling :func:`acb_theta_agm_hadamard` twice on random input is
equivalent to multiplying by `2^g`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_agm_sqrt

Generates a random complex number *t*, sets *rts* to a low-precision rounding
of *t* (possibly containing zero), and sets *a* to the square of *t*. Checks
that the result of :func:`acb_theta_agm_sqrt` on this input is finite, contains
*t*, and that the precision loss is small when *rts* does not contain zero.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_agm_mul

Checks that the duplication formula holds: the result of
:func:`acb_theta_agm_mul` on vectors containing `\theta_{0,b}(0,\tau)` and
`\theta_{0,b}(z,\tau)` for all `b\in\{0,1\}^g` and any choice of `(z,\tau)`
contains the squared theta values `\theta_{0,b}^2(2z,2\tau)`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_agm_mul_tight

Generates random `\tau` and `z` at working precision *prec*, computes the
associated vectors of distances *d0* and *d* using :func:`acb_theta_dist_a0`,
and constructs vectors *a0* and *a* with entries of the form `x e^{-t}` where
`x` is uniformly random with `|x|\leq 1` (generated by :func:`acb_urandom`) and
*t* is the corresponding entry of *d0* (resp. *d*). Calls
:func:`acb_theta_agm_mul_tight` at a lower precision *mprec*. For each `0\leq
k< 2^g`, checks that the absolute value of `k^{\mathrm{th}}` entry of the
result *res* is at most `e^{-d_k}`, and that the error bound on that entry is
at most `2^{-\mathit{mprec} + \delta} e^{-d_k}` for a reasonable value of
`\delta` (e.g. 25).

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_ql_a0_split

Checks that the result of :func:`acb_theta_ql_a0_split` (using
:func:`acb_theta_ql_a0_naive` as *worker*) agrees with that of
:func:`acb_theta_ql_a0_naive` in case of success.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_ql_a0_steps

Checks that the result of :func:`acb_theta_ql_a0_steps` (using
:func:`acb_theta_ql_a0_naive` as *worker*) agrees with that of
:func:`acb_theta_ql_a0_naive` in case of success.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_ql_a0

Checks that :func:`acb_theta_ql_a0`, if successful, agrees with
:func:`acb_theta_ql_a0_naive` on random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_ql_reduce

Generates random values `\tau` and `z` in such a way that
:func:`acb_theta_ql_reduce` is likely to output `s < g` and a nonzero *n1*, and
checks that the claimed inequalities in that function's documentation hold when
computing theta values using :func:`acb_theta_naive_all`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_ql_all

Checks that :func:`acb_theta_ql_all` agrees with :func:`acb_theta_naive_all` on
random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_ql_bounds

Generates random `(z,\tau)` at a working precision that is not too low and
calls :func:`acb_theta_jet_ql_bounds` to compute the bounds *c* and
*rho*. Checks that they are finite and that their definition is satisfied by
sampling theta values on the corresponding neighborhood of `z` at low
precisions with :func:`acb_theta_naive_all`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_ql_radius

Checks that the result of :func:`acb_theta_jet_ql_radius` on random input
satisfies the required inequalities.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_ql_finite_diff

Checks that :func:`acb_theta_jet_ql_finite_diff` computes the correct Taylor
coefficients for the function `\exp(z_0+\cdots+z_{g-1})` at zero. Correct input
can be generated by :func:`acb_theta_jet_ql_radius`, as the bounds *c* and
*rho* can be computed directly for this function.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_ql_all

Checks that :func:`acb_theta_jet_ql_all` agrees with
:func:`acb_theta_jet_naive_all` on random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_transform_char

Checks that the `a` component of any theta characteristic remains the same
after applying :func:`acb_theta_transform_char` when the symplectic matrix is
trigonal as in :func:`sp2gz_trig`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_transform_sqrtdet

Checks that the result of :func:`acb_theta_transform_sqrtdet` on any input
`\tau\in \mathbb{H}_g` squares to `\det(\tau)`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_transform_kappa

Checks that :func:`acb_theta_transform_kappa` and
:func:`acb_theta_transform_kappa2` agree on random input (i.e. they are
congruent modulo 4).

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_transform_proj

Checks that applying :func:`acb_theta_transform_proj` with a random symplectic
matrix, then its inverse gives back the initial vector up to scaling.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_all

Checks that :func:`acb_theta_all` agrees with :func:`acb_theta_naive_all` on
random input. The matrix `\tau` is chosen to be a priori non-reduced but still
reasonably close to the reduced domain.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_jet_all

Checks that :func:`acb_theta_jet_all` agrees with
:func:`acb_theta_jet_naive_all` on random input. The matrix `\tau` is chosen to
be a priori non-reduced but still reasonably close to the reduced domain.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_jet_naive_1

Checks that :func:`acb_theta_g2_jet_naive_1` agrees with
:func:`acb_theta_jet_naive_all` with `g=2`, `z = 0` and `\mathit{ord} = 1` on a
random matrix `\tau`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_detk_symj

Checks that the chain rule holds for the representation `\det^k \mathrm{Sym}^j`
of `\mathrm{GL}_2(\mathbb{C})` as computed by :func:`acb_theta_g2_detk_symj`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_transvectant

Checks that on any sextic polynomial `f = \sum_{j=0}^6 a_j x^{6-j}`, the
transvectant `(f,f)_6` as computed by :func:`acb_theta_g2_transvectant` is
`-3a_2^3 + 8a_2 a_4 - 20a_1 a_5 + 120a_0 a_6`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_transvectant_lead

Checks that the result of :func:`acb_theta_g2_transvectant_lead` is indeed the
leading term of the result of :func:`acb_theta_g2_transvectant` on random
input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_character

Checks that the results of :func:`acb_theta_g2_character` and
:func:`acb_theta_transform_kappa2` for `g=2` are compatible, using the fact
that the product `\chi_5` of the ten even theta constants is a Siegel modular
form with character.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_psi4

Checks that the result of :func:`acb_theta_g2_psi4` is invariant when applying
:func:`acb_theta_transform_proj` on any input vector.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_psi6

Checks that the result of :func:`acb_theta_g2_psi6` is multiplied by `\pm 1`
when applying :func:`acb_theta_transform_proj` on any input vector. The correct
sign is given by :func:`acb_theta_transform_kappa2`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_chi10

Checks that the result of :func:`acb_theta_g2_chi10` is multiplied by `\pm 1`
when applying :func:`acb_theta_transform_proj` on any input vector. The correct
sign is given by :func:`acb_theta_transform_kappa2`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_chi12

Checks that the result of :func:`acb_theta_g2_chi12` is invariant when applying
:func:`acb_theta_transform_proj` on any input vector.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_chi5

Checks that the result of :func:`acb_theta_g2_chi5` squares to the result of
:func:`acb_theta_g2_chi10` on any input vector.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_chi35

Checks that the result of :func:`acb_theta_g2_chi35` is multiplied by `i^k`
when applying :func:`acb_theta_transform_proj` on an input vector of theta
values. The exponent `k` is given by :func:`acb_theta_transform_kappa2`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_chi3_6

Checks that the product `\chi_{8,6} = \chi_{5}\chi_{3,6}`, computed using
:func:`acb_theta_g2_chi5` and :func:`acb_theta_g2_chi3_6`, indeed defines a
modular form of weight `\det^8\mathrm{Sym}^6` by evaluating both sides of the
transformation law on random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_sextic

Checks that the discriminant of the result of :func:`acb_theta_g2_sextic` on a
random matrix `\tau` is `2^{12}\chi_{10}(\tau)`, as computed by
:func:`acb_theta_g2_chi10`.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_sextic_chi5

Checks that the results of :func:`acb_theta_g2_sextic_chi5` agree with those of
:func:`acb_theta_g2_sextic` and :func:`acb_theta_g2_chi5` on random input.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_covariants

Checks that the output of :func:`acb_theta_g2_covariants` agrees with that of
:func:`acb_theta_g2_psi4` using the relation `20\psi_4 = - C_{2,0} + 3
C_{4,0})`. Also checks that each covariant, when evaluated on the result of
:func:`acb_theta_g2_sextic`, defines a Siegel modular function of the correct
weight by evaluating the transformation law, and that covariants take integral
values when the input polynomial is integral.

.. code-block:: bash

    ./build/acb_theta/test/main acb_theta_g2_covariants_lead

Checks that the results of :func:`acb_theta_g2_covariants_lead` are indeed the
leading terms of the results of :func:`acb_theta_g2_covariants` on random
input.

Profiling
-------------------------------------------------------------------------------

.. code-block:: bash

    ./build/acb_theta/profile/p-siegel_reduce g pstep pmax dstep dmax

Prints quick performance measurements for :func:`acb_siegel_reduce`: for the
given `g`, for *d* `\leq` *dmax* by steps of *dstep*, and *prec* `\leq` *pmax*
by steps of *pstep*, constructs an input matrix `w` as `\tau/d` where `\tau` is
generated by :func:`acb_siegel_randtest_reduced` and runs
:func:`acb_siegel_reduce` on `w` at working precision *prec*.

This is meant to show that reduction is generally not a critical step when
evaluating theta functions.

.. code-block:: bash

    ./build/acb_theta/profile/p-ql_a0_split g prec cstep cmax

Prints quick performance measurements for :func:`acb_theta_ql_a0_split`: for
the given `g` and at the given working precision *prec*, generates an input
matrix `\tau` as in :func:`acb_siegel_randtest_reduced`, but whose lower right
`(g-s)\times (g-s)` submatrix is subsequently multiplied by `c`, where `s` runs
between `1` and `g-1` and `c\leq` *cmax* is increased by steps of *cstep*. The
running times of :func:`acb_theta_ql_a0_steps` with or without splitting at `s`
are then compared on each of these matrices, as well as the running time of
:func:`acb_theta_ql_a0`.

This is meant to provide information on how the choice of splitting in
:func:`acb_theta_ql_a0` should be made.

.. code-block:: bash

    ./build/acb_theta/profile/p-ql_a0_steps g pstep pmax

Prints quick performance measurements for :func:`acb_theta_ql_a0_steps`: for
the given `g` and for a working precision *prec* `\leq` *pmax* increasing by
steps of *pstep*, generates a random matrix `\tau` in the reduced domain and
compares the running time of :func:`acb_theta_ql_a0_steps` with different
parameters *nb_steps*.

This is meant to provide information on the correct value to return in
:func:`acb_theta_ql_a0_nb_steps`.

.. code-block:: bash

    ./build/acb_theta/profile/p-all g nb_steps hasz

Prints quick performance measurements for the functions :func:`acb_theta_all`,
:func:`acb_theta_ql_a0`, :func:`acb_theta_ql_all` and
:func:`acb_theta_naive_all` at different precisions on a specific input matrix
of the specified dimension *g*. We start at precision 32, then double it
*nb_steps* times. The parameter *hasz* should be either 0 (theta constants) or
1 (theta values at a nonzero point).

This is meant to show whether the main user function is slower than naive
algorithms at low precisions. (This is currently the case.)

.. code-block:: bash

    ./build/acb_theta/profile/p-jet_all

Prints quick performance measurements for the functions
:func:`acb_theta_jet_all` and :func:`acb_theta_jet_naive_all` at different
precisions and order 1 on a specific input matrix for `g=2`.

This is meant to show whether the main user function is slower than naive
algorithms at low precisions. (This is currently the case.)

