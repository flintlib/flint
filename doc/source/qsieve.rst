.. _qsieve:

**qsieve.h** -- Quadratic sieve
================================================================================

.. function:: ulong qsieve_knuth_schroeppel(qs_t qs_inf)

    Return the Knuth-Schroeppel multiplier for the `n`, integer to be factored
    based upon the Knuth-Schroeppel function.

.. function:: ulong qsieve_primes_init(qs_t qs_inf)

    Compute the factor base prime along with there inverse for `kn`, where `k` 
    is Knuth-Schroeppel multiplier and `n` is the integer to be factored. It
    also computes the square root of `kn` modulo factor base primes.

.. function:: ulong qsieve_primes_increment(qs_t qs_inf, ulong delta)

    It increase the number of factor base primes by amount 'delta' and 
    calculate inverse of those primes along with the square root of `kn` modulo
    those primes.

.. function:: void qsieve_init_A0(qs_t qs_inf)

    First it chooses the possible range of factor of `A _0`, based on the number
    of bits in optimal value of `A _0`. It tries to select range such that we have
    plenty of primes to choose from as well as number of factor in `A _0` are 
    sufficient. For input of size less than 130 bit, this selection method doesn't
    work therefore we randomly generate 2 or 3-subset of all the factor base prime
    as the factor of `A _0`.
    Otherwise, if we have to select `s` factor for `A _0`, we generate `s - 1`-
    subset from odd indices of the possible range of factor and then search last 
    factor using binary search from the even indices of possible range of factor
    such that value of `A _0` is close to it's optimal value.

.. function:: int qsieve_next_A(qs_t qs_inf)

    Find the next value of `A`. If `A` has `s > 3` factors, the first `s - 1`
    are chosen at random among the smallest primes of the allowed range whose
    index is not divisible by 4 (the number of candidates being adjusted as we
    go so that their products are about right), and the last factor is found
    by binary search among those with index divisible by 4, so that `A` is
    within a factor of two of its optimal value. Values of `A` used before are
    rejected. With at most 3 factors, the next subset in lexicographic order
    is used. Return 0 if no new value is found.

.. function:: void qsieve_compute_pre_data(qs_t qs_inf)
    
    Precompute all the data associated with factor's of `A _0`, since `A _0` is going
    to be fixed for several `A`.

.. function:: void qsieve_init_poly_first(qs_t qs_inf)
    
    Initializes the value of `A = q _0 * A _0`, where `q _0` is non-factor base prime.
    precompute the data necessary for generating different `B` value using grey code
    formula. Combine the data calculated for the factor of `A _0` along with the 
    parameter `q _0` to obtain data as for factor of `A`. It also calculates the sieve
    offset for all the factor base prime, for first polynomial.

.. function:: void qsieve_init_poly_next(qs_t qs_inf, slong i)

    Generate next polynomial or next `B` value for particular `A` and also updates the
    sieve offsets for all the factor base prime, for this `B` value.

.. function:: void qsieve_compute_C(fmpz_t C, qs_t qs_inf, qs_poly_t poly)

    Given `A` and `B`, calculate `C = (B ^2 - A) / N`.

.. function:: void qsieve_do_sieving(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)

    First initialize the sieve array to zero, then for each `p \in` ``factor base``, add
    `\log_2(p)` to the locations `\operatorname{soln1} _p + i * p` and `\operatorname{soln2} _p + i * p` for 
    `i = 0, 1, 2,\dots`, where `\operatorname{soln1} _p` and `\operatorname{soln2} _p` are the sieve offsets calculated
    for `p`.

.. function:: void qsieve_do_sieving2(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)

    Perform the same task as above but instead of sieving over whole array at once divide
    the array in blocks and then sieve over each block for all the primes in factor base.
    The primes below ``QS_SMALL_BOUND`` are sieved over sub-blocks of ``QS_L1_SIZE``
    bytes, to make better use of the cache.

.. function:: void qsieve_do_sieving3(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)

    Perform the same task as :func:`qsieve_do_sieving2`, but for the primes
    larger than ``BLOCK_SIZE``, use a bucket sieve: each such prime is visited
    once, and its hits are recorded in buckets, one per sub-block, which are
    then applied to the sub-blocks in turn. The buckets are also used by
    :func:`qsieve_evaluate_candidate` to find the large factor base primes
    dividing a candidate without trial division. This is the routine used
    unless the sieve is small.

.. function:: slong qsieve_evaluate_candidate(qs_t qs_inf, ulong i, unsigned char * sieve, qs_poly_t poly)

    For location `i` in sieve array value at which, is greater than sieve threshold, check
    the value of `Q(x)` at position `i` for smoothness. If value is found to be smooth then
    store it for later processing, else check the residue for the partial if it is found to
    be partial then store it for late processing.

.. function:: slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)

    Scan the sieve array for location at, which accumulated value is greater than sieve
    threshold.
    
.. function:: slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve)

    Call for initialization of polynomial, sieving, and scanning of sieve
    for all the possible polynomials for particular hypercube i.e. `A`.

.. function:: void qsieve_relations_init(qs_t qs_inf)

    Initialise the store which holds the relations found by sieving. Called by
    :func:`qsieve_init_with_tune`.

.. function:: void qsieve_relations_clear(qs_t qs_inf)

    Free the relation store. Called by :func:`qsieve_clear`.

.. function:: void qsieve_relations_reset(qs_t qs_inf)

    Discard all the relations found so far, retaining the allocations made for
    them. Used when the linear algebra fails to split `n` and sieving has to
    start again with a larger factor base.

.. function:: void qsieve_add_relation(qs_t qs_inf, ulong prime, ulong prime2, const fmpz_t Y, const qs_poly_t poly)

    Append a relation to the store, where *prime* and *prime2* are its large
    primes, *prime2* being 1 if it has at most one and both being 1 if the
    relation is full, and *poly* supplies the small prime exponents and the
    remaining factors. The relations are held in memory for the duration of the
    run. As this is called from the sieving threads, the caller must hold
    ``qs_inf->mutex``.

.. function:: void qsieve_add_partial_edge(qs_t qs_inf, ulong prime, ulong prime2)

    Add the edge *prime* -- *prime2* (with *prime2* equal to 1 for a single
    large prime) to the graph of large primes of the partial relations,
    incrementing ``qs_inf->num_cycles`` if it closes a new independent cycle,
    which is tracked with a union-find structure. As this is called from the
    sieving threads, the caller must hold ``qs_inf->mutex``.

.. function:: hash_t * qsieve_get_table_entry(qs_t qs_inf, ulong prime)

    Return the pointer to the location of 'prime' is hash table if it exist, else
    create and entry for it in hash table and return pointer to that.

.. function:: void qsieve_add_to_hashtable(qs_t qs_inf, ulong prime)
    
    Add 'prime' to the hast table.

.. function:: relation_t qsieve_get_relation(qs_t qs_inf, slong i)

    Return a copy of the *i*-th relation in the store. The store retains its
    own copy; the caller owns, and must free, the one returned.

.. function:: relation_t qsieve_merge_relation(qs_t qs_inf, relation_t  a, relation_t  b)

    Given two partial relation having same large prime, merge them to obtain a full
    relation.

.. function:: int qsieve_compare_relation(const void * a, const void * b)

    Compare two relation based on, first large prime, then number of factor and then
    offsets of factor in factor base.

.. function:: int qsieve_remove_duplicates(relation_t * rel_list, slong num_relations)

    Remove duplicate from given list of relations by sorting relations in the list.

.. function:: void qsieve_insert_relation2(qs_t qs_inf, relation_t * rel_list, slong num_relations)

    Given a list of relations, insert each relation from the list into the matrix for
    further processing. 

.. function:: int qsieve_process_relation(qs_t qs_inf)

    After we have accumulated required number of relations, process the
    relations found so far: remove duplicates, and combine the partial
    relations along the cycles of the graph of their large primes (found
    from a spanning forest) to obtain full relations. If there are enough,
    insert them in the matrix and return 1. Otherwise lower the estimate
    ``qs_inf->num_cycles`` by the shortfall and return 0, so that sieving
    continues. Return `-1` if a factor of `n` was found among the large primes.

.. function:: void qsieve_factor(fmpz_factor_t factors, const fmpz_t n)

    Factor `n` using the quadratic sieve method. It is required that `n` is not a
    prime and not a perfect power. There is no guarantee that the factors found will
    be prime, or distinct.

.. function:: void qsieve_factor_with_tune(fmpz_factor_t factors, const fmpz_t n, ulong ks_primes, slong fb_primes, slong small_primes, slong sieve_size, ulong sieve_bits)

    As :func:`qsieve_factor`, but with the provided tuning parameters
    instead of default values.

