---
title: Lock-free shared caches
summary: Read-mostly tables shared by all threads, read with one acquire load and grown by publishing a new version under a lock; retired versions are kept, and geometric growth bounds the memory they use.
complexity: A read costs one atomic load. A rebuild costs a lock plus the build, and happens $O(\log N)$ times for data growing geometrically to size $N$.
order: 3
---

Many computations depend on tables that are expensive to compute but read
very often: powers of a constant, reduction tables, precomputed primes. A
shared cache stores one *version* of such a table for the whole process.
The design aims at making the common case, a read of a version that is
large enough, as cheap as possible, at the price of some extra memory.

## The protocol {#protocol}

A cache holds a pointer to the current version, a lock and a list of
retired versions.

::: algorithm flint_shared_cache_ensure {#alg-ensure}
**Input:** a cache, a predicate `sufficient` and a builder `build`.

1. Load the current version $v$ with acquire semantics. If $v$ is not
   `NULL` and `sufficient(v)` holds, return $v$.
2. Take the lock. Load $v$ again; another thread may have published a
   larger version since step 1.
3. If $v$ is still not sufficient, compute $v' = $ `build(v)`, append $v$
   to the retired list, store $v'$ as the current version with release
   semantics, and let $v = v'$.
4. Release the lock and return $v$.
:::

This is double-checked locking, made correct by the memory ordering of
step 1 and step 3: the release store in step 3 is ordered after all writes
that built $v'$, and an acquire load that observes the new pointer also
observes those writes [@BA2008]. A reader can therefore never see a
partially built version. Builds are serialised by the lock, so two threads
needing a larger table at the same time build it once.

## Why versions are never freed early {#retire}

A reader holds no lock while it uses a version, so the cache cannot know
when the last reader of a superseded version has finished. Instead, old
versions are only retired: they stay allocated until
`flint_cleanup_master`, which runs when no other thread is computing. The
memory this costs is bounded if versions grow geometrically.

::: proposition {#prop-memory}
If each new version is at least $r > 1$ times as large as the version it
replaces, then the retired versions together use less than $1/(r-1)$ times
the memory of the current version.
:::

::: proof
If the current version has size $s$, the retired ones have sizes at most
$s/r, s/r^2, \ldots$, which sum to less than $s/(r - 1)$.
:::

With doubling, $r = 2$, the cache never uses more than twice the memory of
its current version. The builder receives the current version so that it
can extend it, reusing the work already done, rather than start from
scratch.

## Use in FLINT {#use}

The `mp_real` module keeps its tables in shared caches, for example those
of the bitwise exponential and trigonometric functions. The pattern
applies to any lazily extended read-only data, and avoids both global
locks on the read path and per-thread copies of large tables.

::: functions flint_shared_cache_ensure flint_shared_cache_get flint_shared_cache_publish :::
