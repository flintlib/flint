/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "hashmap.h"

void hashmap1_init(hashmap1_t h)
{
   h->data = (hashmap1_elem_s *) flint_calloc(HASHMAP_START_SIZE, 
                                                       sizeof(hashmap1_elem_s));
   h->alloc = HASHMAP_START_SIZE;
   h->mask  = HASHMAP_START_MASK;
   h->num_used = 0;
}

void hashmap1_init2(hashmap1_t h, slong size)
{
   slong bits = 10;

   if (2*size >= 0)
      size *= 2;

   while ((WORD(1) << bits) < size)
      bits++;
 
   h->alloc = (WORD(1) << bits);
   h->mask  = h->alloc - 1;
   h->num_used = 0;

   h->data = (hashmap1_elem_s *) flint_calloc(h->alloc, 
                                                       sizeof(hashmap1_elem_s));
}

void hashmap1_clear(hashmap1_t h)
{
   flint_free(h->data);
}

/* find location in data to store value with 1 word key */
slong hashmap1_hash(ulong key, hashmap1_t h)
{
   slong loc, i;

   if (h->num_used == h->alloc/2)
      return -WORD(1); /* hashmap is full */

   loc = (slong) hashmap1_hash_key(key, h);

   for (i = 0; i < h->alloc; i++)
   {
      if (h->data[loc].in_use == 0 || h->data[loc].key == key)
         return loc;

      loc++;
      if (loc == h->alloc)
         loc = 0;
   }

   return -WORD(1); /* map needs rehashing */
}

/* rehash a full hashmap to twice current size */
void hashmap1_rehash(hashmap1_t h)
{
   slong i;
   hashmap1_elem_s * tmp;

   tmp = h->data;
   h->data = (hashmap1_elem_s *) flint_calloc(2*h->alloc, sizeof(hashmap1_elem_s));

   h->alloc = 2*h->alloc;
   h->mask = h->alloc - 1;
   h->num_used = 0;

   for (i = 0; i < h->alloc/2; i++)
   {
      if (tmp[i].in_use == 1)
         hashmap1_insert(tmp[i].key, tmp[i].value, h);
   }
   
   flint_free(tmp);
}

/* insert key, value pair into hashmap */
void hashmap1_insert(ulong key, void * value, hashmap1_t h)
{
   slong loc;

   loc = hashmap1_hash(key, h);
   if (loc == -WORD(1))
   {
      hashmap1_rehash(h);
      loc = hashmap1_hash(key, h);

      if (loc == -WORD(1))
      {
         /* should never be reached */
         flint_printf("Rehashing unsuccessful\n");
         flint_abort();
      }
   }

   h->data[loc].value = value;
   h->data[loc].key = key;
   h->data[loc].in_use = 1;
   h->num_used += 1;
}

/*
   set *ptr to location of value corresponding to key in hashmap
   return 1 if found, otherwise return 0 (in which case *ptr = NULL)
*/
int hashmap1_find(void ** ptr, ulong key, hashmap1_t h)
{
   slong i, loc;

   loc = hashmap1_hash_key(key, h);

   for (i = 0; i < h->alloc; i++)
   {
      if (h->data[loc].in_use == 0)
      {
          (*ptr) = NULL;
          return 0;
      }

      if (h->data[loc].key == key)
      {
         (*ptr) = h->data[loc].value;
         return 1;
      }

      loc++;
      if (loc == h->alloc)
         loc = 0;
   }

   (*ptr) = NULL;

   return 0;
}
