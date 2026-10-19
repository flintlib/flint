"""
Parallel map over independent work items in forked worker processes.

The workers inherit the parent's memory (the loaded site, the rendered
pages), so only the items and results are sent between processes. Where
fork is not available (Windows) or only one job is requested, the work is
done in the current process.
"""

import multiprocessing
import os

JOBS = 0        # 0: the number of CPUs; set by the -j option

_state = {}     # data for the workers, inherited through fork


def jobs():
    return JOBS or os.cpu_count() or 1


def can_fork():
    return "fork" in multiprocessing.get_all_start_methods()


def pmap(func, items, progress=None, chunksize=None, state=None):
    """[func(x) for x in items], computed by worker processes; func must be a
    module-level function. `state` (any object) is available to the workers
    as parallel.state() and is not pickled. Updates `progress` as results
    arrive; the results are returned in the order of `items`."""
    items = list(items)
    n = jobs()
    _state["value"] = state
    try:
        if n <= 1 or len(items) < 2 or not can_fork():
            out = []
            for x in items:
                out.append(func(x))
                if progress:
                    progress.step()
            return out
        if chunksize is None:
            chunksize = max(1, min(64, len(items) // (4 * n)))
        ctx = multiprocessing.get_context("fork")
        out = []
        with ctx.Pool(n) as pool:
            for r in pool.imap(func, items, chunksize):
                out.append(r)
                if progress:
                    progress.step()
        return out
    finally:
        _state.pop("value", None)


def state():
    return _state.get("value")
