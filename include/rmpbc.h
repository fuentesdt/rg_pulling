/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _rmpbc_h
#define _rmpbc_h
#include "visibility.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gmx_rmpbc *gmx_rmpbc_t;

GMX_LIBGMX_EXPORT
gmx_rmpbc_t gmx_rmpbc_init(t_idef *idef, int ePBC, int natoms,
                           matrix box);

GMX_LIBGMX_EXPORT
void gmx_rmpbc_done(gmx_rmpbc_t gpbc);

GMX_LIBGMX_EXPORT
void gmx_rmpbc(gmx_rmpbc_t gpbc, int natoms, matrix box, rvec x[]);
/* Correct coordinates x for atoms within every molecule for the periodic
 * boundary conditions such that every molecule is whole.
 * natoms is the size x and can be smaller than the number
 * of atoms in idef, but should only contain complete molecules.
 * When ePBC=-1, the type of pbc is guessed from the box matrix.
 */

GMX_LIBGMX_EXPORT
void gmx_rmpbc_copy(gmx_rmpbc_t gpbc, int natoms, matrix box, rvec x[],
                    rvec x_s[]);
/* As gmx_rmpbc, but outputs in x_s and does not modify x. */

GMX_LIBGMX_EXPORT
void gmx_rmpbc_trxfr(gmx_rmpbc_t gpbc, t_trxframe *fr);
/* As gmx_rmpbc but operates on a t_trxframe data structure. */

/*void rm_pbc(t_idef *idef,int ePBC,int natoms,
   matrix box,rvec x[],rvec x_s[]);*/
/* Convenience function that still holds a static variable. */

GMX_LIBGMX_EXPORT
void rm_gropbc(t_atoms *atoms, rvec x[], matrix box);
/* Simple routine for use in analysis tools that just have a pdb or
 * similar file.
 */

#ifdef __cplusplus
}
#endif

#endif  /* _rmpbc_h */
