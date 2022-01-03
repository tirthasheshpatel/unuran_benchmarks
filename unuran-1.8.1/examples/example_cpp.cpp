/* ------------------------------------------------------------- */
/* File: example1.cpp                                            */
/*                                                               */
/* This example simply shows that UNU.RAN can directly be used   */
/* with C++.                                                     */
/* (Of course the example could be much less "C-stylish".)       */
/* ------------------------------------------------------------- */

/* Include UNURAN header file.                                   */
#include <unuran.h>

#include <iostream>

/* ------------------------------------------------------------- */

int main(void)
{
  int    i;     /* loop variable                                 */
  double x;     /* will hold the random number                   */

  /* Declare the three UNURAN objects.                           */
  UNUR_DISTR *distr;    /* distribution object                   */
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */

  /* Use a predefined standard distribution:                     */
  /*   Gaussian with mean zero and standard deviation 1.         */
  /*   Since this is the standard form of the distribution,      */
  /*   we need not give these parameters.                        */
  distr = unur_distr_normal(NULL, 0);

  /* Choose a method: AROU.                                      */
  /* For other (suitable) methods replace "arou" with the        */
  /* respective name (in lower case letters).                    */
  par = unur_arou_new(distr);

  /* Now you can change some of the default settings for the     */
  /* parameters of the chosen method. We don't do it here.       */

  /* Create the generator object.                                */
  gen = unur_init(par);

  /* Notice that this call has also destroyed the parameter      */
  /* object `par' as a side effect.                              */

  /* It is important to check if the creation of the generator   */
  /* object was successful. Otherwise `gen' is the NULL pointer  */ 
  /* and would cause a segmentation fault if used for sampling.  */
  if (gen == NULL) {
    std::cerr << "ERROR: cannot create generator object" << std::endl;
    exit (EXIT_FAILURE);
  }

  /* It is possible to reuse the distribution object to create   */
  /* another generator object. If you do not need it any more,   */
  /* it should be destroyed to free memory.                      */
  unur_distr_free(distr);

  /* Now you can use the generator object `gen' to sample from   */
  /* the standard Gaussian distribution.                         */
  /* Eg.:                                                        */
  for (i=0; i<10; i++) {
    x = unur_sample_cont(gen);
    std::cout << "x(" << i << ") = " << x << std::endl;
  }

  /* When you do not need the generator object any more, you     */
  /* can destroy it.                                             */
  unur_free(gen);

  exit (EXIT_SUCCESS);

} /* end of main() */

/* ------------------------------------------------------------- */

