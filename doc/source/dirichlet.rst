.. _dirichlet:

**dirichlet.h** -- Dirichlet characters
===================================================================================

*Warning: the interfaces in this module are experimental and may change
without notice.*

This module allows working with Dirichlet characters algebraically.
For evaluations of characters as complex numbers, see :ref:`acb-dirichlet`.

Dirichlet characters
-------------------------------------------------------------------------------

Working with Dirichlet characters mod *q* consists mainly
in going from residue classes mod *q* to exponents on a set
of generators of the group.

This implementation relies on the Conrey numbering scheme
introduced in the
`L-functions and Modular Forms DataBase <http://www.lmfdb.org/Character/Dirichlet>`_,
which is an explicit choice of generators
allowing to represent Dirichlet characters via the pairing

.. math::

   \begin{array}{ccccc}
   (\mathbb Z/q\mathbb Z)^\times \times (\mathbb Z/q\mathbb Z)^\times & \to & \bigoplus_i \mathbb Z/\phi_i\mathbb Z \times \mathbb Z/\phi_i\mathbb Z & \to &\mathbb C \\
   (m,n) & \mapsto& (a_i,b_i) &\mapsto& \chi_q(m,n) = \exp(2i\pi\sum \frac{a_ib_i}{\phi_i} )
   \end{array}

We call *number* a residue class `m` modulo *q*, and *log* the
corresponding vector `(a_i)` of exponents of Conrey generators.

Going from a *log* to the corresponding *number* is a cheap
operation we call exponential, while the converse requires computing discrete
logarithms.

Multiplicative group modulo *q*
-------------------------------------------------------------------------------

.. type:: dirichlet_group_struct

.. type:: dirichlet_group_t

    Represents the group of Dirichlet characters mod *q*.

    An *dirichlet_group_t* is defined as an array of *dirichlet_group_struct*
    of length 1, permitting it to be passed by reference.

.. function:: int dirichlet_group_init(dirichlet_group_t G, ulong q)

    Initializes *G* to the group of Dirichlet characters mod *q*.

    This method computes a canonical decomposition of *G* in terms of cyclic
    groups, which are the mod `p^e` subgroups for `p^e\|q`, plus
    the specific generator described by Conrey for each subgroup.

    In particular *G* contains:

    - the number *num* of components

    - the generators

    - the exponent *expo* of the group

    It does *not* automatically precompute lookup tables
    of discrete logarithms or numerical roots of unity, and can therefore
    safely be called even with large *q*.

    For implementation reasons, the largest prime factor of *q* must not
    exceed `10^{16}`. This restriction could
    be removed in the future. The function returns 1 on success and 0
    if a factor is too large.

.. function:: void dirichlet_subgroup_init(dirichlet_group_t H, const dirichlet_group_t G, ulong h)

   Given an already computed group *G* mod `q`, initialize its subgroup *H*
   defined mod `h\mid q`. Precomputed discrete log tables are inherited.

.. function:: void dirichlet_group_clear(dirichlet_group_t G)

    Clears *G*. Remark this function does *not* clear the discrete logarithm
    tables stored in *G* (which may be shared with another group).

.. function:: ulong dirichlet_group_size(const dirichlet_group_t G)

   Returns the number of elements in *G*, i.e. `\varphi(q)`.

.. function:: ulong dirichlet_group_num_primitive(const dirichlet_group_t G)

   Returns the number of primitive elements in *G*.

.. function:: void dirichlet_group_dlog_precompute(dirichlet_group_t G, ulong num)

    Precompute decomposition and tables for discrete log computations in *G*,
    so as to minimize the complexity of *num* calls to discrete logarithms.

    If *num* gets very large, the entire group may be indexed.

.. function:: void dirichlet_group_dlog_clear(dirichlet_group_t G)

   Clear discrete logarithm tables in *G*. When discrete logarithm tables are
   shared with subgroups, those subgroups must be cleared before clearing the
   tables.

Character type
-------------------------------------------------------------------------------

.. type:: dirichlet_char_struct

.. type:: dirichlet_char_t

    Represents a Dirichlet character.
    This structure contains both a *number* (residue class) and
    the corresponding *log* (exponents on the group generators).

    An *dirichlet_char_t* is defined as an array of *dirichlet_char_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void dirichlet_char_init(dirichlet_char_t chi, const dirichlet_group_t G)

    Initializes *chi* to an element of the group *G* and sets its value
    to the principal character.

.. function:: void dirichlet_char_clear(dirichlet_char_t chi)

    Clears *chi*.

.. function:: void dirichlet_char_print(const dirichlet_group_t G, const dirichlet_char_t chi)

    Prints the array of exponents representing this character.

.. function:: void dirichlet_char_log(dirichlet_char_t x, const dirichlet_group_t G, ulong m)

    Sets *x* to the character of number *m*, computing its log using discrete
    logarithm in *G*.

.. function:: ulong dirichlet_char_exp(const dirichlet_group_t G, const dirichlet_char_t x)

    Returns the number *m* corresponding to exponents in *x*.

.. function:: ulong _dirichlet_char_exp(dirichlet_char_t x, const dirichlet_group_t G)

    Computes and returns the number *m* corresponding to exponents in *x*.
    This function is for internal use.

.. function:: void dirichlet_char_one(dirichlet_char_t x, const dirichlet_group_t G)

    Sets *x* to the principal character in *G*, having *log* `[0,\dots 0]`.

.. function:: void dirichlet_char_first_primitive(dirichlet_char_t x, const dirichlet_group_t G)

    Sets *x* to the first primitive character of *G*, having *log* `[1,\dots 1]`,
    or `[0, 1, \dots 1]` if `8\mid q`.

.. function:: void dirichlet_char_set(dirichlet_char_t x, const dirichlet_group_t G, const dirichlet_char_t y)

    Sets *x* to the element *y*.

.. function:: int dirichlet_char_next(dirichlet_char_t x, const dirichlet_group_t G)

    Sets *x* to the next character in *G* according to lexicographic ordering
    of *log*.

    The return value
    is the index of the last updated exponent of *x*, or *-1* if the last
    element has been reached.

    This function allows to iterate on all elements of *G* looping on their *log*.
    Note that it produces elements in seemingly random *number* order.

    The following template can be used for such a loop::

        dirichlet_char_one(chi, G);
        do {
            /* use character chi */
        } while (dirichlet_char_next(chi, G) >= 0);

.. function:: int dirichlet_char_next_primitive(dirichlet_char_t x, const dirichlet_group_t G)

    Same as :func:`dirichlet_char_next`, but jumps to the next primitive character of *G*.

.. function:: ulong dirichlet_index_char(const dirichlet_group_t G, const dirichlet_char_t x)

    Returns the lexicographic index of the *log* of *x* as an integer in `0\dots \varphi(q)`.

.. function:: void dirichlet_char_index(dirichlet_char_t x, const dirichlet_group_t G, ulong j)

    Sets *x* to the character whose *log* has lexicographic index *j*.

.. function:: int dirichlet_char_eq(const dirichlet_char_t x, const dirichlet_char_t y)

.. function:: int dirichlet_char_eq_deep(const dirichlet_group_t G, const dirichlet_char_t x, const dirichlet_char_t y)

   Return 1 if *x* equals *y*.

   The second version checks every byte of the representation and is intended for testing only.

Character properties
-------------------------------------------------------------------------------

As a consequence of the Conrey numbering, all these numbers are available at the
level of *number* and *char* object. Both case require no discrete log computation.

.. function:: int dirichlet_char_is_principal(const dirichlet_group_t G, const dirichlet_char_t chi)

   Returns 1 if *chi* is the principal character mod *q*.

.. function:: ulong dirichlet_conductor_ui(const dirichlet_group_t G, ulong a)

.. function:: ulong dirichlet_conductor_char(const dirichlet_group_t G, const dirichlet_char_t x)

   Returns the *conductor* of `\chi_q(a,\cdot)`, that is the smallest `r` dividing `q`
   such `\chi_q(a,\cdot)` can be obtained as a character mod `r`.

.. function:: int dirichlet_parity_ui(const dirichlet_group_t G, ulong a)

.. function:: int dirichlet_parity_char(const dirichlet_group_t G, const dirichlet_char_t x)

   Returns the *parity* `\lambda` in `\{0,1\}` of `\chi_q(a,\cdot)`, such that
   `\chi_q(a,-1)=(-1)^\lambda`.

.. function:: ulong dirichlet_order_ui(const dirichlet_group_t G, ulong a)

.. function:: ulong dirichlet_order_char(const dirichlet_group_t G, const dirichlet_char_t x)

   Returns the order of `\chi_q(a,\cdot)` which is the order of `a\bmod q`.

.. function:: int dirichlet_char_is_real(const dirichlet_group_t G, const dirichlet_char_t chi)

   Returns 1 if *chi* is a real character (iff it has order `\leq 2`).

.. function:: int dirichlet_char_is_primitive(const dirichlet_group_t G, const dirichlet_char_t chi)

   Returns 1 if *chi* is primitive (iff its conductor is exactly *q*).

Character evaluation
-------------------------------------------------------------------------------

Dirichlet characters take value in a finite cyclic group of roots of unity plus zero.

Evaluation functions return a *ulong*, this number corresponds to the
power of a primitive root of unity, the special value *DIRICHLET_CHI_NULL*
encoding the zero value.

.. function:: ulong dirichlet_pairing(const dirichlet_group_t G, ulong m, ulong n)

.. function:: ulong dirichlet_pairing_char(const dirichlet_group_t G, const dirichlet_char_t chi, const dirichlet_char_t psi)

   Compute the value of the Dirichlet pairing on numbers *m* and *n*, as
   exponent modulo *G->expo*.

   The *char* variant takes as input two characters, so that no discrete
   logarithm is computed.

   The returned value is the numerator of the actual value exponent mod the group exponent *G->expo*.

.. function:: ulong dirichlet_chi(const dirichlet_group_t G, const dirichlet_char_t chi, ulong n)

   Compute the value `\chi(n)` as the exponent modulo *G->expo*.

.. function:: void dirichlet_chi_vec(ulong * v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)

   Compute the list of exponent values *v[k]* for `0\leq k < nv`, as exponents
   modulo *G->expo*.

.. function:: void dirichlet_chi_vec_order(ulong * v, const dirichlet_group_t G, const dirichlet_char_t chi, ulong order, slong nv)

   Compute the list of exponent values *v[k]* for `0\leq k < nv`, as exponents
   modulo *order*, which is assumed to be a multiple of the order of *chi*.

Character operations
-------------------------------------------------------------------------------

.. function:: void dirichlet_char_mul(dirichlet_char_t chi12, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2)

   Multiply two characters of the same group *G*.

.. function:: void dirichlet_char_pow(dirichlet_char_t c, const dirichlet_group_t G, const dirichlet_char_t a, ulong n)

   Take the power of a character.

.. function:: void dirichlet_char_lift(dirichlet_char_t chi_G, const dirichlet_group_t G, const dirichlet_char_t chi_H, const dirichlet_group_t H)

    If *H* is a subgroup of *G*, computes the character in *G* corresponding to
    *chi_H* in *H*.

.. function:: void dirichlet_char_lower(dirichlet_char_t chi_H, const dirichlet_group_t H, const dirichlet_char_t chi_G, const dirichlet_group_t G)

    If *chi_G* is a character of *G* which factors through *H*, sets *chi_H* to
    the corresponding restriction in *H*.

    This requires `c(\chi_G)\mid q_H\mid q_G`, where `c(\chi_G)` is the
    conductor of `\chi_G` and `q_G, q_H` are the moduli of G and H.

