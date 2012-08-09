#include <stdlib.h>
#include <stdarg.h>
#include <gmp.h>

#include "ecc_point.h"

/* copied from lib/nettle */

int mp_init_multi(mpz_t *a, ...)
{
   mpz_t    *cur = a;
   int       np  = 0;
   va_list   args;

   va_start(args, a);
   while (cur != NULL) {
       mpz_init(*cur);
       ++np;
       cur = va_arg(args, mpz_t*);
   }
   va_end(args);
   return 0;
}

void mp_clear_multi(mpz_t *a, ...)
{
   mpz_t    *cur = a;
   va_list   args;

   va_start(args, a);
   while (cur != NULL) {
       mpz_clear(*cur);
       cur = va_arg(args, mpz_t*);
   }
   va_end(args);
}

/*
   Allocate a new ECC point
   @return A newly allocated point or NULL on error
*/
ecc_point *
ecc_new_point (void)
{
  ecc_point *p;
  p = calloc (1, sizeof (*p));
  if (p == NULL)
    {
      return NULL;
    }
  if (mp_init_multi (&p->x, &p->y, &p->z, NULL) != 0)
    {
      free (p);
      return NULL;
    }
  return p;
}

/* Free an ECC point from memory
  @param p   The point to free
*/
void
ecc_del_point (ecc_point * p)
{
  /* prevents free'ing null arguments */
  if (p != NULL)
    {
      mp_clear_multi (&p->x, &p->y, &p->z, NULL);       /* note: p->z may be NULL but that's ok with this function anyways */
      free (p);
    }
}

/* end of copied */
