
Computing isomorphisms and embeddings of finite fields
--------------------------------------------------------------------------------


.. function:: void fq_nmod_embed_gens(fq_nmod_t gen_sub, fq_nmod_t gen_sup, nmod_poly_t minpoly, const fq_nmod_ctx_t sub_ctx, const fq_nmod_ctx_t sup_ctx)

    Given two contexts ``sub_ctx`` and ``sup_ctx``, such that
    ``degree(sub_ctx)`` divides ``degree(sup_ctx)``, compute:

    \begin{itemize}
    \item an element ``gen_sub`` in ``sub_ctx`` such that
      ``gen_sub`` generates the finite field defined by
      ``sub_ctx``,
    \item its minimal polynomial ``minpoly``,
    \item a root ``gen_sup`` of ``minpoly`` inside the field
      defined by ``sup_ctx``.
    \end{itemize}
    
    These data uniquely define an embedding of ``sub_ctx`` into
    ``sup_ctx``.

    
.. function:: void _fq_nmod_embed_gens_naive(fq_nmod_t gen_sub, fq_nmod_t gen_sup, nmod_poly_t minpoly, const fq_nmod_ctx_t sub_ctx, const fq_nmod_ctx_t sup_ctx)

    Given two contexts ``sub_ctx`` and ``sup_ctx``, such that
    ``degree(sub_ctx)`` divides ``degree(sup_ctx)``, compute an
    embedding of ``sub_ctx`` into ``sup_ctx`` defined as follows:

    \begin{itemize}
    \item ``gen_sub`` is the canonical generator of ``sup_ctx``
      (i.e., the class of `X`),
    \item ``minpoly`` is the defining polynomial of ``sub_ctx``,
    \item ``gen_sup`` is a root of ``minpoly`` inside the field
      defined by ``sup_ctx``.
    \end{itemize}

.. function:: void fq_nmod_embed_matrices(nmod_mat_t embed, nmod_mat_t project, const fq_nmod_t gen_sub, const fq_nmod_ctx_t sub_ctx, const fq_nmod_t gen_sup, const fq_nmod_ctx_t sup_ctx, const nmod_poly_t gen_minpoly)

    Given:

    \begin{itemize}
    \item two contexts ``sub_ctx`` and ``sup_ctx``, of
      respective degrees `m` and `n`, such that `m` divides `n`;
    \item a generator ``gen_sub`` of ``sub_ctx``, its minimal
      polynomial ``gen_minpoly``, and a root ``gen_sup`` of
      ``gen_minpoly`` in ``sup_ctx``, as returned by
      ``fq_nmod_embed_gens``;
    \end{itemize}
    
    Compute:

    \begin{itemize}
    \item the `n\times m` matrix ``embed`` mapping ``gen_sub``
      to ``gen_sup``, and all their powers accordingly;
    \item an `m\times n` matrix ``project`` such that
      ``project`` `\times` ``embed`` is the `m\times m` identity
      matrix.
    \end{itemize}

.. function:: void fq_nmod_embed_trace_matrix(nmod_mat_t res, const nmod_mat_t basis, const fq_nmod_ctx_t sub_ctx, const fq_nmod_ctx_t sup_ctx)

    Given:

    \begin{itemize}
    \item two contexts ``sub_ctx`` and ``sup_ctx``, of degrees
      `m` and `n`, such that `m` divides `n`;
    \item an `n\times m` matrix ``basis`` that maps ``sub_ctx``
      to an isomorphic subfield in ``sup_ctx``;
    \end{itemize}

    Compute the `m\times n` matrix of the trace from ``sup_ctx`` to
    ``sub_ctx``.

    This matrix is computed as
    \[
      \texttt{embed\_dual\_to\_mono\_matrix(\_, sub\_ctx)}
      \times \texttt{basis}^t \times
      \texttt{embed\_mono\_to\_dual\_matrix(\_, sup\_ctx)}.
    \]
    \textbf{Note:} if
    `m=n`, ``basis`` represents a Frobenius, and the result is its
    inverse matrix.

.. function:: void fq_nmod_embed_composition_matrix(nmod_mat_t matrix, const fq_nmod_t gen, const fq_nmod_ctx_t ctx)

    Compute the \emph{composition matrix} of ``gen``.

    For an element `a\in\mathbf{F}_{p^n}`, its composition matrix is the
    matrix whose columns are `a^0, a^1, \ldots, a^{n-1}`.

.. function:: void fq_nmod_embed_composition_matrix_sub(nmod_mat_t matrix, const fq_nmod_t gen, const fq_nmod_ctx_t ctx, slong trunc)

    Compute the \emph{composition matrix} of ``gen``, truncated to
    ``trunc`` columns.

.. function:: void fq_nmod_embed_mul_matrix(nmod_mat_t matrix, const fq_nmod_t gen, const fq_nmod_ctx_t ctx)

    Compute the \emph{multiplication matrix} of ``gen``.

    For an element `a` in `\mathbf{F}_{p^n}=\mathbf{F}_p[x]`, its
    multiplication matrix is the matrix whose columns are `a, ax,
    \dots, ax^{n-1}`.

.. function:: void fq_nmod_embed_mono_to_dual_matrix(nmod_mat_t res, const fq_nmod_ctx_t ctx)

    Compute the change of basis matrix from the monomial basis of
    ``ctx`` to its dual basis.

.. function:: void fq_nmod_embed_dual_to_mono_matrix(nmod_mat_t res, const fq_nmod_ctx_t ctx)

    Compute the change of basis matrix from the dual basis of
    ``ctx`` to its monomial basis.

.. function:: void fq_nmod_modulus_pow_series_inv(nmod_poly_t res, const fq_nmod_ctx_t ctx, slong trunc)

    Compute the power series inverse of the reverse of the modulus of
    ``ctx`` up to `O(x^\texttt{trunc})`.

.. function:: void fq_nmod_modulus_derivative_inv(fq_nmod_t m_prime, fq_nmod_t m_prime_inv, const fq_nmod_ctx_t ctx)

    Compute the derivative ``m_prime`` of the modulus of ``ctx``
    as an element of ``ctx``, and its inverse ``m_prime_inv``.
