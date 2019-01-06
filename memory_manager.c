/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2016 Claus Fieker
    Copyright (C) 2016 William Hart.
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef __unix__
#include <unistd.h> /* sysconf */
#endif

#if defined(_WIN32) || defined(WIN32)
#include <windows.h> /* GetSystemInfo */
#endif

#include <string.h> /* memset */
#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "thread_pool.h"

#if HAVE_GC
#include "gc.h"
#endif

#define PAGES_PER_BLOCK 16

FLINT_TLS_PREFIX slong flint_page_size;
FLINT_TLS_PREFIX slong flint_page_mask;
FLINT_TLS_PREFIX int flint_pools_initialised = 0;
FLINT_TLS_PREFIX void ** flint_pool_ptr[16]; /* pools 2^i limbs */
FLINT_TLS_PREFIX slong flint_pool_num_free[16]; /* no. free in pool */
FLINT_TLS_PREFIX slong flint_pool_free_alloc[16]; /* alloc size of free lists*/

static void * _flint_malloc(size_t);
static void * _flint_calloc(size_t, size_t);
static void * _flint_realloc(void *, size_t);
static void _flint_free(void *);

static void  *(*__flint_allocate_func) (size_t) = _flint_malloc;
static void  *(*__flint_callocate_func) (size_t, size_t) = _flint_calloc;
static void  *(*__flint_reallocate_func) (void *, size_t) = _flint_realloc;
static void  (*__flint_free_func) (void *) = _flint_free;

#if FLINT_REENTRANT && !HAVE_TLS
#include <pthread.h>

static pthread_once_t register_initialised = PTHREAD_ONCE_INIT;
pthread_mutex_t register_lock;

static pthread_once_t alloc_func_init = PTHREAD_ONCE_INIT;
pthread_mutex_t alloc_func_lock;

void __flint_set_memory_functions_init()
{
   pthread_mutex_init(&alloc_func_lock, NULL);
}
#endif

static void flint_memory_error(size_t size)
{
    flint_printf("Exception (FLINT memory_manager). Unable to allocate memory (%ld).\n", size);
    flint_abort();
}

void __flint_set_memory_functions(void *(*alloc_func) (size_t),
                             void *(*calloc_func) (size_t, size_t),
                             void *(*realloc_func) (void *, size_t),
                             void (*free_func) (void *))
{  
#if FLINT_REENTRANT && !HAVE_TLS
    pthread_once(&alloc_func_init, __flint_set_memory_functions_init);
    pthread_mutex_lock(&alloc_func_lock);
#endif

  __flint_allocate_func = alloc_func;
  __flint_callocate_func = calloc_func;
  __flint_reallocate_func = realloc_func;
  __flint_free_func = free_func;

#if FLINT_REENTRANT && !HAVE_TLS
    pthread_mutex_unlock(&alloc_func_lock);
#endif
}

void * _flint_malloc(size_t size)
{
   void * ptr;

#if HAVE_GC
   ptr = GC_malloc(size);
#else
   ptr = malloc(size);
#endif

   return ptr;
}

void * flint_malloc(size_t size)
{
   void * ptr = (*__flint_allocate_func)(size);

   if (ptr == NULL)
        flint_memory_error(size);

   return ptr;
}

void * _flint_realloc(void * ptr, size_t size)
{
    void * ptr2;

#if HAVE_GC
    ptr2 = GC_realloc(ptr, size);
#else
    ptr2 = realloc(ptr, size);
#endif

    return ptr2;
}


void * flint_realloc(void * ptr, size_t size)
{
    void * ptr2;
  
    if (ptr)
      ptr2 = (*__flint_reallocate_func)(ptr, size);
    else
      ptr2 = (*__flint_allocate_func)(size);

    if (ptr2 == NULL)
        flint_memory_error(size);

    return ptr2;
}

void * _flint_calloc(size_t num, size_t size)
{
   void * ptr;

#if HAVE_GC
    ptr = GC_malloc(num*size);
#else
    ptr = calloc(num, size);
#endif

    return ptr;
}

void * flint_calloc(size_t num, size_t size)
{
   void * ptr;

    ptr = (*__flint_callocate_func)(num, size);

    if (ptr == NULL)
        flint_memory_error(size);

    return ptr;
}

void _flint_free(void * ptr)
{
#if !HAVE_GC
    free(ptr);
#endif
}

void flint_free(void * ptr)
{
   (*__flint_free_func)(ptr);
}


FLINT_TLS_PREFIX size_t flint_num_cleanup_functions = 0;

FLINT_TLS_PREFIX flint_cleanup_function_t * flint_cleanup_functions = NULL;

#pragma omp threadprivate(flint_num_cleanup_functions, flint_cleanup_functions)

#if FLINT_REENTRANT && !HAVE_TLS
void register_init()
{
   pthread_mutex_init(&register_lock, NULL);
}
#endif

void flint_register_cleanup_function(flint_cleanup_function_t cleanup_function)
{
#if FLINT_REENTRANT && !HAVE_TLS
    pthread_once(&register_initialised, register_init);
    pthread_mutex_lock(&register_lock);
#endif

    flint_cleanup_functions = flint_realloc(flint_cleanup_functions,
        (flint_num_cleanup_functions + 1) * sizeof(flint_cleanup_function_t));

    flint_cleanup_functions[flint_num_cleanup_functions] = cleanup_function;

    flint_num_cleanup_functions++;

#if FLINT_REENTRANT && !HAVE_TLS
    pthread_mutex_unlock(&register_lock);
#endif
}

void _fmpz_cleanup();

void flint_cleanup()
{
    size_t i;

#if FLINT_REENTRANT && !HAVE_TLS
    pthread_mutex_lock(&register_lock);
#endif

    for (i = 0; i < flint_num_cleanup_functions; i++)
        flint_cleanup_functions[i]();

    flint_free(flint_cleanup_functions);
    flint_cleanup_functions = NULL;
    flint_num_cleanup_functions = 0;

    mpfr_free_cache();
    _fmpz_cleanup();

    if (flint_pools_initialised)
       flint_pool_cleanup();
 
#if FLINT_REENTRANT && !HAVE_TLS
    pthread_mutex_unlock(&register_lock);
#endif

}

void flint_cleanup_master()
{
    if (global_thread_pool_initialized)
    {
        thread_pool_clear(global_thread_pool);
        global_thread_pool_initialized = 0;
    }
    flint_cleanup();
}

/******************************************************************************

   Threaded memory manager

*******************************************************************************/

/*
   This memory manager improves threaded performance by using a pooled
   allocator. To use it, one must set both the Flint and GMP memory allocation
   functions to use the _flint_pooled_malloc/realloc/free below.

   The allocator requires thread local storage and a pthreads implementation.
*/

slong flint_get_page_size()
{
#if defined(__unix__)
   return sysconf(_SC_PAGESIZE);
#elif defined(_WIN32) || defined(WIN32)
   SYSTEM_INFO si;
   GetSystemInfo(&si);
   return si.dwPageSize;
#else
   return 4096;
#endif
}

void * flint_align_ptr(void * ptr, slong size)
{
    slong mask = ~(size - 1);

    return (void *)((mask & (slong) ptr) + size);
}

void * flint_pooled_malloc(size_t n)
{
   void * aligned_ptr, * ptr;
   slong i, skip, num, nind;

   /* round up to multiple of 64 bit limb size and set n to number of limbs */
   n = (n + 7) >> 3;
   if (n < 1) n = 1;

   if (!flint_pools_initialised) /* check if pools are initialised */
   {
      /* initialise data for pages and blocks */
      flint_page_size = flint_get_page_size();
      flint_page_mask = ~(flint_page_size - 1);

      /* initialise data for pools */
      for (i = 0; i < 16; i++)
      {
         flint_pool_ptr[i] = NULL;
         flint_pool_num_free[i] = 0;
	 flint_pool_free_alloc[i] = 0;
      }

      flint_pools_initialised = 1;
   }

   if (n > (flint_page_size >> 4)) /* large allocation, single block per alloc */
   {
      /* allocate new block */
      ptr = malloc(n*8 + sizeof(flint_pool_header_s) + flint_page_size);

      /* align to page boundary */
      aligned_ptr = flint_align_ptr(ptr, flint_page_size);

      /* set alloc count to 1, must count down to zero before freeing block */
      ((flint_pool_header_s *) ptr)->count = 1;

      /* set pointer in page to start of entire block */
      ((flint_pool_header_s *) aligned_ptr)->address = ptr;

      /* set n in header to number of words per alloc */
      ((flint_pool_header_s *) aligned_ptr)->n = n;

      return (void *) ((slong) aligned_ptr + sizeof(flint_pool_header_s));
   }
   
   /* pooled allocation of n words */
      
   /* round up to power of two limbs, i.e. set n = 2^nind */
   nind = FLINT_BIT_COUNT(n - 1);
   n = (WORD(1) << nind);
      
   if (flint_pool_num_free[nind] == 0) /* no allocs left in pool */
   {
      void * ptri;
      slong j, k;
      slong block_size = PAGES_PER_BLOCK*flint_page_size;

      /* allocate new block */
      ptr = malloc(block_size + flint_page_size);

      /* align to page boundary */
      aligned_ptr = flint_align_ptr(ptr, flint_page_size);

      /* no. n word blocks dedicated to header, per page */
      skip = (sizeof(flint_pool_header_s) - 1)/(8*n) + 1;

      /* total number of n word allocs worth per page */
      num = flint_page_size/(8*n);

      /* set no. allocs, must count down to 0 before freeing block */
      ((flint_pool_header_s *) ptr)->count = PAGES_PER_BLOCK*(num - skip);
#if HAVE_PTHREAD
      ((flint_pool_header_s *) ptr)->thread = pthread_self();
#endif

      /* no. allocs in block */
      flint_pool_num_free[nind] = PAGES_PER_BLOCK*(num - skip);
           
      /* allocate free list for pool */
      if (flint_pool_free_alloc[nind] == 0)
      {
         flint_pool_ptr[nind] = malloc(flint_pool_num_free[nind]*sizeof(void *));
         flint_pool_free_alloc[nind] = flint_pool_num_free[nind];
      }

      /* set up headers for pages in block */
      for (i = 0; i < PAGES_PER_BLOCK; i++)
      {
          void * page_ptr = (void *)((slong) aligned_ptr + i*flint_page_size);

          /* set pointer in each page to start of entire block */
          ((flint_pool_header_s *) page_ptr)->address = ptr;

          /* set n in each page to number of words per alloc */
          ((flint_pool_header_s *) page_ptr)->n = n;
      }

      /* set pointer to end of block */
      ptri = (void *) ((slong) aligned_ptr + block_size);
      
      j = 0;
      
      for (i = 0; i < PAGES_PER_BLOCK; i++)
      {
         for (k = 0; k < num - skip; j++, k++)
	 {
	    ptri = (void *) ((slong) ptri - 8*n);
	    
	    flint_pool_ptr[nind][j] = ptri; 
         }

         /* skip over header to end of next page of block */
         ptri = (void *) ((slong) ptri & flint_page_mask);
      }
   }

   /* get next entry from pool */
   flint_pool_num_free[nind]--;
      
   return flint_pool_ptr[nind][flint_pool_num_free[nind]];
}

void * flint_pooled_calloc(size_t num, size_t size)
{
   size_t n = num*size;
   void * ptr = flint_pooled_malloc(n);

   memset(ptr, 0, n);

   return ptr;
}

void flint_pooled_free(void * ptr)
{
   if (ptr != NULL)
   {
      /* get block header pointer */
      flint_pool_header_s * header_ptr = (flint_pool_header_s *)((slong) ptr & flint_page_mask);
      flint_pool_header_s * header_ptr1 = (flint_pool_header_s *) header_ptr->address;
      int n = header_ptr->n;
      
      if (n > (flint_page_size >> 4)) /* large allocation */
         free(header_ptr1);
      else /* pooled allocation */
      {
         int nind = FLINT_BIT_COUNT(n - 1);

         /* clean up if left from another thread */
#if HAVE_PTHREAD
         if (!pthread_equal(header_ptr1->thread, pthread_self()))
         {
            if (--header_ptr1->count == 0)
               free(header_ptr1);
         } else
#endif
         {
            /* check free list has enough space alloc'd */
            if (flint_pool_num_free[nind] == flint_pool_free_alloc[nind])
            {
               flint_pool_free_alloc[nind] *= 2;

               flint_pool_ptr[nind] = realloc(flint_pool_ptr[nind],
   	                                 2*flint_pool_free_alloc[nind]*sizeof(void *));
            }

            flint_pool_ptr[nind][flint_pool_num_free[nind]++] = ptr;
         }
      }
   }
}

void flint_pooled_free_with_size(void * ptr, size_t size)
{
   flint_pooled_free(ptr);
}

void * flint_pooled_realloc(void * ptr, size_t n)
{
   void * new_ptr;
   slong i;

   /* get header pointer */
   flint_pool_header_s * header_ptr;

   /* round up to multiple of 64 bit limb size and set n to number of limbs */
   n = (n + 7) >> 3;
   if (n < 1) n = 1;

   if (ptr != NULL)
   {
      header_ptr = (flint_pool_header_s *)((slong) ptr & flint_page_mask);

      /* if current alloc is big enough, and new alloc not small, return old one */
      if (n <= header_ptr->n && 2*n >= header_ptr->n)
         return ptr;

      /* else create new allocation, copy, free old alloc and return new one */
      new_ptr = flint_pooled_malloc(n*8);

      for (i = 0; i < FLINT_MIN(header_ptr->n, n); i++)
         ((slong *) new_ptr)[i] = ((slong *) ptr)[i];

      flint_pooled_free(ptr);
   } else
      new_ptr = flint_pooled_malloc(n*8);

   return new_ptr;
}

void * flint_pooled_realloc_with_old_size(void * ptr, size_t old_n, size_t n)
{
   return flint_pooled_realloc(ptr, n);
}

void flint_pool_cleanup()
{
   slong i, j;

   /* for each pool */
   for (i = 0; i < 16; i++)
   {
      for (j = flint_pool_num_free[i] - 1; j >= 0; j--)
      {
         void * ptr = flint_pool_ptr[i][j];

	 flint_pool_header_s * header_ptr = (flint_pool_header_s *)((slong) ptr & flint_page_mask);

	 header_ptr = header_ptr->address;

         if (--header_ptr->count == 0)
            free(header_ptr);
      }

      free(flint_pool_ptr[i]);

      flint_pool_num_free[i] = 0;
      flint_pool_free_alloc[i] = 0;
      flint_pool_ptr[i] = NULL;
   }
}

void flint_pooled_init()
{
   __gmp_set_memory_functions(flint_pooled_malloc,
        flint_pooled_realloc_with_old_size, flint_pooled_free_with_size);
   __flint_set_memory_functions(flint_pooled_malloc, flint_pooled_calloc,
                            flint_pooled_realloc, flint_pooled_free);
}

