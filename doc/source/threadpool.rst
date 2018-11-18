.. _threadpool:

**threadpool.h** -- thread pool
===============================================================================


Threadpool
--------------------------------------------------------------------------------


.. function:: void threadpool_init(threadpool_t T, slong l)

    Initialise `T` and create `l` sleeping threads that are available to work.
    If `l \le 0` no threads are created and future calls to
    ``threadpool_request`` will return `0`. 

.. function:: slong threadpool_size(threadpool_t T)

    Return the number of threads in `T`.

.. function:: slong threadpool_request(threadpool_t T, threadpool_threadhandle * out, slong requested)

    Put at most ``requested`` threads in the unavailable state and return
    their handles. The handles are written to ``out`` and the number of
    handles written is returned. These threads must be released by a call to
    ``threadpool_giveback``.

.. function:: void threadpool_wake(threadpool_t T, threadpool_threadhandle i, void (*f)(void*), void * a)

    Wake up a sleeping thread `i` and have it work on ``f(a)``.

.. function:: void threadpool_wait(threadpool_t T, threadpool_threadhandle i)

    Wait for thread `i` to finish working and go back to sleep.

.. function:: void threadpool_giveback(threadpool_t T, threadpool_threadhandle i)

    Put thread `i` back in the available state. This thread should be sleeping
    when this function is called.

.. function:: void threadpool_clear(threadpool_t T)

    Release any resources used by `T`. All threads should be given back before
    this function is called.
