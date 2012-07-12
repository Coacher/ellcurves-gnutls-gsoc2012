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
   Negate an ECC point
   @param P        The point to negate
   @param R        [out] The destination of the negate
   @return 0 on success
*/
int
ecc_projective_negate_point (ecc_point * P, ecc_point * R)
{

  if (P == NULL || R == NULL)
    return -1;

  /* Set R == -P
   * Since we are dealing only with curves of type
   * y^2 = x^3 + ax + b negation of point is
   * the same as negation of its y coordinate
   */
  mpz_set(R->x, P->x);
  mpz_neg(R->y, P->y);
  mpz_set(R->z, P->z);

  return 0;
}
