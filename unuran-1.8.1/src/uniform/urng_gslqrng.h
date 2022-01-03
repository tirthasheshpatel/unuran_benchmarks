/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_gslqmc.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *     Function prototypes for using URNG of type GSL-:                       *
 *     Function prototypes for using uniform of type GSL:                    *
 *     uniform random number from GSL (GNU Scientific Library),              *
 *     see http://www.gnu.org/software/gsl/.                                 *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifndef URNG_GSLQRNG_H_SEEN
#define URNG_GSLQRNG_H_SEEN
/*---------------------------------------------------------------------------*/
#include <gsl/gsl_qrng.h>
/*---------------------------------------------------------------------------*/

/* 
   =NODE  URNG-GSLQRNG  Interface to GSL generators for quasi-random points

   =UP URNG [30]

   =DESCRIPTION
      Interface to the generators for quasi-random points (also called
      low discrepancy point sets) from the GNU Scientific Library (GSL).
      Documentation and source code of this library is available from
      @uref{http://www.gnu.org/software/gsl/}.

      The interface to the GSL must be compiled into UNU.RAN using the
      configure flag @code{--with-urng-gsl}.
      Notice that the GSL has to be installed before running
      @code{./configure}.

   =HOWTOUSE
      When using this interface @file{unuran_urng_gsl.h} must be included
      in the corresponding C file, i.e., one must add the line
      @smallexample
      #include <unuran_urng_gsl.h>
      @end smallexample

      Moreover, one must not forget to link the executable against
      @file{libgsl}.

      The following routines are supported for URNG objects of this
      type:

      @itemize @minus
      @item unur_urng_sample()
      @item unur_urng_sample_array()
      @item unur_urng_reset() 
      @item unur_urng_sync() 
      @item unur_urng_free()
      @end itemize

      unur_urng_sync() is used to jump to the first coordinate of
      the next point generated by the generator.

   =END

*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

/*---------------------------------------------------------------------------*/

UNUR_URNG *unur_urng_gslqrng_new( const gsl_qrng_type *qrngtype, unsigned int dim );
/*
   Make object for quasi-random point generators for dimension
   @var{dim} from the @file{GSL} (GNU Scientific Library). 
   @var{qrngtype} is the type of the chosen generator as described in
   the GSL manual (see section Quasi-Random Sequences). 
   This library is available from @uref{http://www.gnu.org/software/gsl/}.
*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* URNG_GSLQRNG_H_SEEN */
/*---------------------------------------------------------------------------*/