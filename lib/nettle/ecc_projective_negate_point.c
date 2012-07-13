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
   @param modulus  The modulus of the field the ECC curve is in
   @return 0 on success
*/
int
ecc_projective_negate_point (ecc_point * P, ecc_point * R, mpz_t modulus)
{

  if (P == NULL || R == NULL)
    return -1;

  /* we set R.y to modulus - P.y to avoid negative coordinates */
  mpz_set(R->x, P->x);
  mpz_sub(R->y, modulus, P->y);
  mpz_set(R->z, P->z);

  return 0;
}
