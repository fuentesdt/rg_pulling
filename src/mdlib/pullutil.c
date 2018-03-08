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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include "sysstuff.h"
#include "princ.h"
#include "futil.h"
#include "statutil.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "names.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "network.h"
#include "pbc.h"
#include "pull.h"
#include "gmx_ga2la.h"

static void pull_set_pbcatom(t_commrec *cr, t_pullgrp *pg,
                             t_mdatoms *md, rvec *x,
                             rvec x_pbc)
{
    int a, m;

    if (cr && PAR(cr))
    {
        if (DOMAINDECOMP(cr))
        {
            if (!ga2la_get_home(cr->dd->ga2la, pg->pbcatom, &a))
            {
                a = -1;
            }
        }
        else
        {
            a = pg->pbcatom;
        }

        if (a >= md->start && a < md->start+md->homenr)
        {
            copy_rvec(x[a], x_pbc);
        }
        else
        {
            clear_rvec(x_pbc);
        }
    }
    else
    {
        copy_rvec(x[pg->pbcatom], x_pbc);
    }
}

static void pull_set_pbcatoms(t_commrec *cr, t_pull *pull,
                              t_mdatoms *md, rvec *x,
                              rvec *x_pbc)
{
    int g, n, m;

    n = 0;
    for (g = 0; g < 1+pull->ngrp; g++)
    {
        if ((g == 0 && PULL_CYL(pull)) || pull->grp[g].pbcatom == -1)
        {
            clear_rvec(x_pbc[g]);
        }
        else
        {
            pull_set_pbcatom(cr, &pull->grp[g], md, x, x_pbc[g]);
            for (m = 0; m < DIM; m++)
            {
                if (pull->dim[m] == 0)
                {
                    x_pbc[g][m] = 0.0;
                }
            }
            n++;
        }
    }

    if (cr && PAR(cr) && n > 0)
    {
        /* Sum over the nodes to get x_pbc from the home node of pbcatom */
        gmx_sum((1+pull->ngrp)*DIM, x_pbc[0], cr);
    }
}

/* switch function, x between r and w */
static real get_weight(real x, real r1, real r0)
{
    real weight;

    if (x >= r0)
    {
        weight = 0;
    }
    else if (x <= r1)
    {
        weight = 1;
    }
    else
    {
        weight = (r0 - x)/(r0 - r1);
    }

    return weight;
}

static void make_cyl_refgrps(t_commrec *cr, t_pull *pull, t_mdatoms *md,
                             t_pbc *pbc, double t, rvec *x, rvec *xp)
{
    int         g, i, ii, m, start, end;
    rvec        g_x, dx, dir;
    double      r0_2, sum_a, sum_ap, dr2, mass, weight, wmass, wwmass, inp;
    t_pullgrp  *pref, *pgrp, *pdyna;
    gmx_ga2la_t ga2la = NULL;

    if (pull->dbuf_cyl == NULL)
    {
        snew(pull->dbuf_cyl, pull->ngrp*4);
    }

    if (cr && DOMAINDECOMP(cr))
    {
        ga2la = cr->dd->ga2la;
    }

    start = md->start;
    end   = md->homenr + start;

    r0_2 = dsqr(pull->cyl_r0);

    /* loop over all groups to make a reference group for each*/
    pref = &pull->grp[0];
    for (g = 1; g < 1+pull->ngrp; g++)
    {
        pgrp  = &pull->grp[g];
        pdyna = &pull->dyna[g];
        copy_rvec(pgrp->vec, dir);
        sum_a          = 0;
        sum_ap         = 0;
        wmass          = 0;
        wwmass         = 0;
        pdyna->nat_loc = 0;

        for (m = 0; m < DIM; m++)
        {
            g_x[m] = pgrp->x[m] - pgrp->vec[m]*(pgrp->init[0] + pgrp->rate*t);
        }

        /* loop over all atoms in the main ref group */
        for (i = 0; i < pref->nat; i++)
        {
            ii = pull->grp[0].ind[i];
            if (ga2la)
            {
                if (!ga2la_get_home(ga2la, pref->ind[i], &ii))
                {
                    ii = -1;
                }
            }
            if (ii >= start && ii < end)
            {
                pbc_dx_aiuc(pbc, x[ii], g_x, dx);
                inp = iprod(dir, dx);
                dr2 = 0;
                for (m = 0; m < DIM; m++)
                {
                    dr2 += dsqr(dx[m] - inp*dir[m]);
                }

                if (dr2 < r0_2)
                {
                    /* add to index, to sum of COM, to weight array */
                    if (pdyna->nat_loc >= pdyna->nalloc_loc)
                    {
                        pdyna->nalloc_loc = over_alloc_large(pdyna->nat_loc+1);
                        srenew(pdyna->ind_loc, pdyna->nalloc_loc);
                        srenew(pdyna->weight_loc, pdyna->nalloc_loc);
                    }
                    pdyna->ind_loc[pdyna->nat_loc] = ii;
                    mass   = md->massT[ii];
                    weight = get_weight(sqrt(dr2), pull->cyl_r1, pull->cyl_r0);
                    pdyna->weight_loc[pdyna->nat_loc] = weight;
                    sum_a += mass*weight*inp;
                    if (xp)
                    {
                        pbc_dx_aiuc(pbc, xp[ii], g_x, dx);
                        inp     = iprod(dir, dx);
                        sum_ap += mass*weight*inp;
                    }
                    wmass  += mass*weight;
                    wwmass += mass*sqr(weight);
                    pdyna->nat_loc++;
                }
            }
        }
        pull->dbuf_cyl[(g-1)*4+0] = wmass;
        pull->dbuf_cyl[(g-1)*4+1] = wwmass;
        pull->dbuf_cyl[(g-1)*4+2] = sum_a;
        pull->dbuf_cyl[(g-1)*4+3] = sum_ap;
    }

    if (cr && PAR(cr))
    {
        /* Sum the contributions over the nodes */
        gmx_sumd(pull->ngrp*4, pull->dbuf_cyl, cr);
    }

    for (g = 1; g < 1+pull->ngrp; g++)
    {
        pgrp  = &pull->grp[g];
        pdyna = &pull->dyna[g];

        wmass         = pull->dbuf_cyl[(g-1)*4+0];
        wwmass        = pull->dbuf_cyl[(g-1)*4+1];
        pdyna->wscale = wmass/wwmass;
        pdyna->invtm  = 1.0/(pdyna->wscale*wmass);

        for (m = 0; m < DIM; m++)
        {
            g_x[m]         = pgrp->x[m] - pgrp->vec[m]*(pgrp->init[0] + pgrp->rate*t);
            pdyna->x[m]    = g_x[m] + pgrp->vec[m]*pull->dbuf_cyl[(g-1)*4+2]/wmass;
            if (xp)
            {
                pdyna->xp[m] = g_x[m] + pgrp->vec[m]*pull->dbuf_cyl[(g-1)*4+3]/wmass;
            }
        }

        if (debug)
        {
            fprintf(debug, "Pull cylinder group %d:%8.3f%8.3f%8.3f m:%8.3f\n",
                    g, pdyna->x[0], pdyna->x[1],
                    pdyna->x[2], 1.0/pdyna->invtm);
        }
    }
}

static double atan2_0_2pi(double y, double x)
{
    double a;

    a = atan2(y, x);
    if (a < 0)
    {
        a += 2.0*M_PI;
    }
    return a;
}

/* calculates center of mass of selection index from all coordinates x */
void pull_calc_coms(t_commrec *cr,
                    t_pull *pull, t_mdatoms *md, t_pbc *pbc, double t,
                    rvec x[], rvec *xp)
{
    int        g, i, ii, m;
    real       mass, w, wm, twopi_box = 0;
    double     wmass, wwmass, invwmass;
    dvec       com, comp;
    double     cm, sm, cmp, smp, ccm, csm, ssm, csw, snw;
    rvec      *xx[2], x_pbc = {0, 0, 0}, dx;
    t_pullgrp *pgrp;

    if (pull->rbuf == NULL)
    {
        snew(pull->rbuf, 1+pull->ngrp);
    }
    if (pull->dbuf == NULL)
    {
        snew(pull->dbuf, 3*(1+pull->ngrp));
    }

    if (pull->bRefAt)
    {
        pull_set_pbcatoms(cr, pull, md, x, pull->rbuf);
    }

    if (pull->cosdim >= 0)
    {
        for (m = pull->cosdim+1; m < pull->npbcdim; m++)
        {
            if (pbc->box[m][pull->cosdim] != 0)
            {
                gmx_fatal(FARGS, "Can not do cosine weighting for trilinic dimensions");
            }
        }
        twopi_box = 2.0*M_PI/pbc->box[pull->cosdim][pull->cosdim];
    }

    for (g = 0; g < 1+pull->ngrp; g++)
    {
        pgrp = &pull->grp[g];
        clear_dvec(com);
        clear_dvec(comp);
        wmass  = 0;
        wwmass = 0;
        cm     = 0;
        sm     = 0;
        cmp    = 0;
        smp    = 0;
        ccm    = 0;
        csm    = 0;
        ssm    = 0;
        if (!(g == 0 && PULL_CYL(pull)))
        {
            if (pgrp->epgrppbc == epgrppbcREFAT)
            {
                /* Set the pbc atom */
                copy_rvec(pull->rbuf[g], x_pbc);
            }
            w = 1;
            for (i = 0; i < pgrp->nat_loc; i++)
            {
                ii   = pgrp->ind_loc[i];
                mass = md->massT[ii];
                if (pgrp->epgrppbc != epgrppbcCOS)
                {
                    if (pgrp->weight_loc)
                    {
                        w = pgrp->weight_loc[i];
                    }
                    wm      = w*mass;
                    wmass  += wm;
                    wwmass += wm*w;
                    if (pgrp->epgrppbc == epgrppbcNONE)
                    {
                        /* Plain COM: sum the coordinates */
                        for (m = 0; m < DIM; m++)
                        {
                            com[m]    += wm*x[ii][m];
                        }
                        if (xp)
                        {
                            gmx_fatal(FARGS, "For Rg pulling, should have only pull = umbrella\n");
                            for (m = 0; m < DIM; m++)
                            {
                                comp[m] += wm*xp[ii][m];
                            }
                        }
                    }
                    else
                    {
                        /* Sum the difference with the reference atom */
                        pbc_dx(pbc, x[ii], x_pbc, dx);
                        for (m = 0; m < DIM; m++)
                        {
                            com[m]    += wm*dx[m];
                        }
                        if (xp)
                        {
                            gmx_fatal(FARGS, "For Rg pulling, should have only pull = umbrella\n");
                            /* For xp add the difference between xp and x to dx,
                             * such that we use the same periodic image,
                             * also when xp has a large displacement.
                             */
                            for (m = 0; m < DIM; m++)
                            {
                                comp[m] += wm*(dx[m] + xp[ii][m] - x[ii][m]);
                            }
                        }
                    }
                }
                else
                {
                    /* Determine cos and sin sums */
                    csw  = cos(x[ii][pull->cosdim]*twopi_box);
                    snw  = sin(x[ii][pull->cosdim]*twopi_box);
                    cm  += csw*mass;
                    sm  += snw*mass;
                    ccm += csw*csw*mass;
                    csm += csw*snw*mass;
                    ssm += snw*snw*mass;

                    if (xp)
                    {
                        csw  = cos(xp[ii][pull->cosdim]*twopi_box);
                        snw  = sin(xp[ii][pull->cosdim]*twopi_box);
                        cmp += csw*mass;
                        smp += snw*mass;
                    }
                }
            }
        }

        /* Copy local sums to a buffer for global summing */
        switch (pgrp->epgrppbc)
        {
            case epgrppbcNONE:
            case epgrppbcREFAT:
                copy_dvec(com, pull->dbuf[g*3]);
                copy_dvec(comp, pull->dbuf[g*3+1]);
                pull->dbuf[g*3+2][0] = wmass;
                pull->dbuf[g*3+2][1] = wwmass;
                pull->dbuf[g*3+2][2] = 0;
                break;
            case epgrppbcCOS:
                pull->dbuf[g*3  ][0] = cm;
                pull->dbuf[g*3  ][1] = sm;
                pull->dbuf[g*3  ][2] = 0;
                pull->dbuf[g*3+1][0] = ccm;
                pull->dbuf[g*3+1][1] = csm;
                pull->dbuf[g*3+1][2] = ssm;
                pull->dbuf[g*3+2][0] = cmp;
                pull->dbuf[g*3+2][1] = smp;
                pull->dbuf[g*3+2][2] = 0;
                break;
        }
    }

    if (cr && PAR(cr))
    {
        /* Sum the contributions over the nodes */
        gmx_sumd((1+pull->ngrp)*3*DIM, pull->dbuf[0], cr);
    }

    for (g = 0; g < 1+pull->ngrp; g++)
    {
        pgrp = &pull->grp[g];
        if (pgrp->nat > 0 && !(g == 0 && PULL_CYL(pull)))
        {
            if (pgrp->epgrppbc != epgrppbcCOS)
            {
                /* Determine the inverse mass */
                wmass    = pull->dbuf[g*3+2][0];
                wwmass   = pull->dbuf[g*3+2][1];
                invwmass = 1/wmass;
                /* invtm==0 signals a frozen group, so then we should keep it zero */
                if (pgrp->invtm > 0)
                {
                    pgrp->wscale = wmass/wwmass;
                    pgrp->invtm  = 1.0/(pgrp->wscale*wmass);
                }
                /* Divide by the total mass */
                for (m = 0; m < DIM; m++)
                {
                    pgrp->x[m]    = pull->dbuf[g*3  ][m]*invwmass;
                    if (xp)
                    {
                        pgrp->xp[m] = pull->dbuf[g*3+1][m]*invwmass;
                    }
                    if (pgrp->epgrppbc == epgrppbcREFAT)
                    {
                        pgrp->x[m]    += pull->rbuf[g][m];
                        if (xp)
                        {
                            pgrp->xp[m] += pull->rbuf[g][m];
                        }
                    }
                }
            }
            else
            {
                /* Determine the optimal location of the cosine weight */
                csw                   = pull->dbuf[g*3][0];
                snw                   = pull->dbuf[g*3][1];
                pgrp->x[pull->cosdim] = atan2_0_2pi(snw, csw)/twopi_box;
                /* Set the weights for the local atoms */
                wmass  = sqrt(csw*csw + snw*snw);
                wwmass = (pull->dbuf[g*3+1][0]*csw*csw +
                          pull->dbuf[g*3+1][1]*csw*snw +
                          pull->dbuf[g*3+1][2]*snw*snw)/(wmass*wmass);
                pgrp->wscale = wmass/wwmass;
                pgrp->invtm  = 1.0/(pgrp->wscale*wmass);
                /* Set the weights for the local atoms */
                csw *= pgrp->invtm;
                snw *= pgrp->invtm;
                for (i = 0; i < pgrp->nat_loc; i++)
                {
                    ii                  = pgrp->ind_loc[i];
                    pgrp->weight_loc[i] = csw*cos(twopi_box*x[ii][pull->cosdim]) +
                        snw*sin(twopi_box*x[ii][pull->cosdim]);
                }
                if (xp)
                {
                    csw                    = pull->dbuf[g*3+2][0];
                    snw                    = pull->dbuf[g*3+2][1];
                    pgrp->xp[pull->cosdim] = atan2_0_2pi(snw, csw)/twopi_box;
                }
            }
            if (debug)
            {
                fprintf(debug, "Pull group %d wmass %f wwmass %f invtm %f\n",
                        g, wmass, wwmass, pgrp->invtm);
            }
        }
    }

    if (PULL_CYL(pull))
    {
        /* Calculate the COMs for the cyclinder reference groups */
        make_cyl_refgrps(cr, pull, md, pbc, t, x, xp);
    }
}

/* Compute radii of gyration or the linearized radius of gyration

   With epullgRADIUS_GYR:  Rg = [ M^[-1] * sum_i m_i (r_i - Rs)^2 ]^(1/2)

   With epullgRADIUS_GYR1: Rg = M^[-1] * sum_i m_i |r_i - Rs|
 */
void pull_calc_rg(t_commrec *cr,
                    t_pull *pull, t_mdatoms *md, t_pbc *pbc, double t,
                    rvec x[], rvec *xp)
{
    int        g, i, ii, m;
    real       mass, w, wm, twopi_box = 0;
    double     wmass, wwmass, invwmass;
    dvec       com, comp, rg_sum, rg_sum_tmp, dvec_tmp; // Jochen: use dvec for rg_sum, so I can use the same global sum routine as for COM above
    double     cm, sm, cmp, smp, ccm, csm, ssm, csw, snw, dist2;
    rvec      *xx[2], x_pbc = {0, 0, 0}, dx, com_rvec;
    t_pullgrp *pgrp;

    // Jochen: for Rg pulling, xp == NULL since I do not have 
    
    if (pull->rbuf == NULL)
    {
        snew(pull->rbuf, 1+pull->ngrp);
    }
    if (pull->dbuf == NULL)
    {
        snew(pull->dbuf, 3*(1+pull->ngrp));
    }

    if (pull->bRefAt)
    {
        pull_set_pbcatoms(cr, pull, md, x, pull->rbuf);
    }

    if (pull->cosdim >= 0)
    {
        for (m = pull->cosdim+1; m < pull->npbcdim; m++)
        {
            if (pbc->box[m][pull->cosdim] != 0)
            {
                gmx_fatal(FARGS, "Can not do cosine weighting for trilinic dimensions");
            }
        }
        twopi_box = 2.0*M_PI/pbc->box[pull->cosdim][pull->cosdim];
    }

    for (g = 0; g < 1+pull->ngrp; g++)
    {
        pgrp = &pull->grp[g];
        clear_dvec(com);
        clear_dvec(comp);
        wmass  = 0;
        wwmass = 0;
        cm     = 0;
        sm     = 0;
        cmp    = 0;
        smp    = 0;
        ccm    = 0;
        csm    = 0;
        ssm    = 0;
        clear_dvec(rg_sum);
        clear_dvec(rg_sum_tmp);
        if (!(g == 0 && PULL_CYL(pull)))
        {
            if (pgrp->epgrppbc == epgrppbcREFAT)
            {
                /* Set the pbc atom */
                copy_rvec(pull->rbuf[g], x_pbc);
            }
            w = 1;
            for (i = 0; i < pgrp->nat_loc; i++)
            {
                ii   = pgrp->ind_loc[i];
                mass = md->massT[ii];
                if (pgrp->epgrppbc != epgrppbcCOS)
                {
                    if (pgrp->weight_loc)
                    {
                        w = pgrp->weight_loc[i];
                    }
                    wm      = w*mass;
                    wmass  += wm;
                    wwmass += wm*w;
                    if (pgrp->epgrppbc == epgrppbcNONE)
                    {
                        gmx_fatal(FARGS, "Rg pulling, should not be here (epgrppbcNONE)\n"); 
                        /* Plain COM: sum the coordinates */
                    }
                    else
                    {
                        /* JOCHEN: get difference wrt. center of mass of this pull group */
                        for (m = 0; m < DIM; m++)
                        {
                            /* Get rvec version of center of mass, to be used in pbc_dx */
                            com_rvec[m] = pgrp->x[m];
                        }
                        pbc_dx(pbc, x[ii], com_rvec, dx);
                        dist2 = 0;
                        clear_dvec(dvec_tmp);
                        for (m = 0; m < DIM; m++)
                        {
                            /* sum up m_i(r_i - Rs)^2 for requested dimensions */
                            if (pull->dim[m])
                            {
                                /* Get square distance of atom to centerf of mass */
                                dist2 += dsqr(dx[m]);

                                /* Keep vector r_j - R, if dimension m is used */
                                dvec_tmp[m] = dx[m];
                            }
                        }

                        if (pull->eGeom == epullgRADIUS_GYR)
                        {
                            /* For normal distance-square-weighted radius of gyration */
                            rg_sum[0] += wm*dist2;
                        }
                        else if (pull->eGeom == epullgRADIUS_GYR1)
                        {
                            /* For distance-linear-weighted radius of gyration */ 
                            rg_sum[0]     += wm*sqrt(dist2);

                            /* Keep m_j * (vec{r}_j - vec{R})/|vec{r}_j - vec{R}| */
                            dsvmul(wm/dnorm(dvec_tmp), dvec_tmp, dvec_tmp);
                            dvec_inc(rg_sum_tmp, dvec_tmp);
                            /*fprintf(stderr,"i = %d, dx       = %8g %8g %8g -- rg_sum_tmp = %8g %8g %8g\n", i, 
                                    dx[0], dx[1], dx[2], rg_sum_tmp[0], rg_sum_tmp[1], rg_sum_tmp[2]);
                            fprintf(stderr,"i = %d, dvec_tmp = %8g %8g %8g -- rg_sum_tmp = %8g %8g %8g\n", i,
                            dvec_tmp[0], dvec_tmp[1], dvec_tmp[2], rg_sum_tmp[0], rg_sum_tmp[1], rg_sum_tmp[2]); */

                        }
                        else
                        {
                            gmx_fatal(FARGS, "Expected pull geometry radius of gyration in pull_calc_rg()\n");
                        }
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "Rg pulling, should not be here\n");
                }
            }
        }

        /* Copy local sums to a buffer for global summing */
        switch (pgrp->epgrppbc)
        {
            case epgrppbcNONE:
            case epgrppbcREFAT:
                copy_dvec(rg_sum,     pull->dbuf[g*3]);   // Changed by Jochen
                copy_dvec(rg_sum_tmp, pull->dbuf[g*3+1]); // Changed by Jochen
                pull->dbuf[g*3+2][0] = wmass;
                pull->dbuf[g*3+2][1] = wwmass;
                pull->dbuf[g*3+2][2] = 0;
                break;
            case epgrppbcCOS:
                pull->dbuf[g*3  ][0] = cm;
                pull->dbuf[g*3  ][1] = sm;
                pull->dbuf[g*3  ][2] = 0;
                pull->dbuf[g*3+1][0] = ccm;
                pull->dbuf[g*3+1][1] = csm;
                pull->dbuf[g*3+1][2] = ssm;
                pull->dbuf[g*3+2][0] = cmp;
                pull->dbuf[g*3+2][1] = smp;
                pull->dbuf[g*3+2][2] = 0;
                break;
        }
    }

    if (cr && PAR(cr))
    {
        /* Sum the contributions over the nodes */
        gmx_sumd((1+pull->ngrp)*3*DIM, pull->dbuf[0], cr);
    }

    for (g = 0; g < 1+pull->ngrp; g++)
    {
        pgrp = &pull->grp[g];
        if (pgrp->nat > 0 && !(g == 0 && PULL_CYL(pull)))
        {
            if (pgrp->epgrppbc != epgrppbcCOS)
            {
                /* Determine the inverse mass */
                wmass    = pull->dbuf[g*3+2][0];
                wwmass   = pull->dbuf[g*3+2][1];
                invwmass = 1/wmass;
                /* invtm==0 signals a frozen group, so then we should keep it zero */
                if (pgrp->invtm > 0)
                {
                    pgrp->wscale = wmass/wwmass;
                    pgrp->invtm  = 1.0/(pgrp->wscale*wmass);
                }

                // Added by Jochen
                /* Divide com_rvec by the total mass */
                if (pull->eGeom == epullgRADIUS_GYR)
                {
                    /* For normal distance-square-weighted radius of gyration, do sqrt */
                    pgrp->radius_gyr = sqrt(pull->dbuf[g*3  ][0]*invwmass);
                }
                else if (pull->eGeom == epullgRADIUS_GYR1)
                {
                    /* For distance-linear-weighted radius of gyration */ 
                    pgrp->radius_gyr = pull->dbuf[g*3  ][0]*invwmass;

                    /* Keep: sum_i m_i/M (vec{r}_j - vec{R})/|vec{r}_j - vec{R}| */
                    dsvmul(invwmass, pull->dbuf[g*3+1], pgrp->aver_dr_norm);
                }
                // printf("Rg: %8g %8g aver_dr_norm = %8g %8g %8g - invwmass = %g\n", t, pgrp->radius_gyr, pgrp->aver_dr_norm[XX], pgrp->aver_dr_norm[YY], pgrp->aver_dr_norm[ZZ], invwmass); fflush(stdout);
            }
            else
            {
                gmx_fatal(FARGS, "Radii of gyration pulling, found strange pgrp->epgrppbc\n");
            }
            if (debug)
            {
                fprintf(debug, "Pull group %d wmass %f wwmass %f invtm %f\n",
                        g, wmass, wwmass, pgrp->invtm);
            }
        }
    }

    if (PULL_CYL(pull))
    {
        /* Calculate the COMs for the cyclinder reference groups */
        make_cyl_refgrps(cr, pull, md, pbc, t, x, xp);
    }
}
