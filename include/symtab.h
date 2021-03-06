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

#ifndef _symtab_h
#define _symtab_h

#include <stdio.h>
#include "visibility.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This module handles symbol table manipulation. All text strings
 * needed by an application are allocated only once. All references
 * to these text strings use handles returned from the put_symtab()
 * routine. These handles can easily be converted to address independent
 * values by invoking lookup_symtab(). So when writing structures to
 * a file which contains text strings, this value can be written in stead
 * of the text string or its address. This value can easily be converted
 * back to a text string handle by get_symtab_handle().
 */

GMX_LIBGMX_EXPORT
void open_symtab(t_symtab *symtab);
/* Initialises the symbol table symtab.
 */

GMX_LIBGMX_EXPORT
void close_symtab(t_symtab *symtab);
/* Undoes the effect of open_symtab(), after invoking this function,
 * no value can be added to the symbol table, only values can be
 * retrieved using get_symtab().
 */

void free_symtab(t_symtab *symtab);
/* Frees the space allocated by the symbol table itself */

GMX_LIBGMX_EXPORT
void done_symtab(t_symtab *symtab);
/* Frees the space allocated by the symbol table, including all
 * entries in it */

GMX_LIBGMX_EXPORT
char **put_symtab(t_symtab *symtab, const char *name);
/* Enters a string into the symbol table symtab, if it was not
 * available, a reference to a copy is returned else a reference
 * to the earlier entered value is returned. Strings are trimmed
 * of spaces.
 */

int lookup_symtab(t_symtab *symtab, char **name);
/* Returns a unique handle for **name, without a memory reference.
 * It is a failure when name cannot be found in the symbol table,
 * it should be entered before with put_symtab().
 */

char **get_symtab_handle(t_symtab *symtab, int name);
/* Returns a text string handle for name. Name should be a value
 * returned from lookup_symtab(). So get_symtab_handle() and
 * lookup_symtab() are inverse functions.
 */

long wr_symtab(FILE *fp, t_symtab *symtab);
/* Writes the symbol table symtab to the file, specified by fp.
 * The function returns the number of bytes written.
 */

long rd_symtab(FILE *fp, t_symtab *symtab);
/* Reads the symbol table symtab from the file, specified by fp.
 * This will include allocating the needed space. The function
 * returns the number of bytes read. The symtab is in the closed
 * state afterwards, so no strings can be added to it.
 */

GMX_LIBGMX_EXPORT
void pr_symtab(FILE *fp, int indent, const char *title, t_symtab *symtab);
/* This routine prints out a (human) readable representation of
 * the symbol table symtab to the file fp. Ident specifies the
 * number of spaces the text should be indented. Title is used
 * to print a header text.
 */

#ifdef __cplusplus
}
#endif

#endif  /* _symtab_h */
