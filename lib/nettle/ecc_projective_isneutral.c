/*
 * Copyright (C) 2011-2012 Free Software Foundation, Inc.
 *
 * This file is part of GNUTLS.
 *
 * The GNUTLS library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 */

#include "ecc.h"

/*
   Check if the given point is the neutral point.
   @param P        The point to check
   @return  0 if given point == neutral point
   @return  1 if given point != neutral point
   @return -1 otherwise
*/
int
ecc_projective_isneutral (ecc_point * P)
{
  if (P == NULL)
    return -1;

  /* neutral point is a point with projective
   * coordinates (k,k,0) where k is any real number
   * excluding point (0,0,0)
   */
  if ( (mpz_cmp_ui(P->z, 0)) || (mpz_cmp(P->x, P->y)) )
    return 1;

  if (!(mpz_cmp_ui(P->x, 0)))
    /* we have excluded (0,0,0) point */
    return -1;

  return 0;
}
