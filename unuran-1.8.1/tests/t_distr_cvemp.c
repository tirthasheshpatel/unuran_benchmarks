/*
	file automatically generated by make_test_files.pl
	Tue Apr 19 14:01:02 2011
*/

/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************/
    
/**
 ** Tests for DISTR_CVEMP
 **/
    
/*---------------------------------------------------------------------------*/
#include "testunuran.h"

#ifdef UNUR_URNG_DEFAULT_RNGSTREAM
#include <RngStream.h>
#endif
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

static FILE *TESTLOG;               /* test log file                         */
static FILE *UNURANLOG;             /* unuran log file                       */

static int test_ok = TRUE;          /* all tests ok (boolean)                */
static int fullcheck = FALSE;       /* whether all checks are performed      */ 

static TIMER watch;                 /* stop watch                            */

/*---------------------------------------------------------------------------*/

void run_verify_generator( FILE *LOG, int line, UNUR_PAR *par );

int unur_distr_cvemp_set_verify( UNUR_PAR *par, int verify);


/*---------------------------------------------------------------------------*/

void test_new (void);
void test_set (void);
void test_get (void);
void test_chg (void);
void test_init (void);
void test_reinit (void);
void test_sample (void);
void test_validate (void);
void test_special(void);

/*---------------------------------------------------------------------------*/



/* prototypes */

#define COMPARE_SAMPLE_SIZE  (500)
#define VIOLATE_SAMPLE_SIZE   (20)




/*---------------------------------------------------------------------------*/

#ifndef CHI2_FAILURES_TOLERATED
#  define CHI2_FAILURES_TOLERATED DEFAULT_CHI2_FAILURES_TOLERATED
#endif

/*---------------------------------------------------------------------------*/
/* [verbatim] */



/*---------------------------------------------------------------------------*/
/* [new] */

void test_new (void)
{
        int n_tests_failed;          /* number of failed tests */

	/* start test */
	printf("[new "); fflush(stdout);
	fprintf(TESTLOG,"\n[new]\n");

	/* reset counter */
	n_tests_failed = 0;

	/* set stop watch */
	stopwatch_lap(&watch);
  
{ /* invalid data */
UNUR_DISTR *distr = NULL;
	distr = NULL; 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,26,(unur_distr_cvemp_new(1)))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,26,UNUR_ERR_DISTR_SET)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_new() */

/*---------------------------------------------------------------------------*/
/* [set] */

void test_set (void)
{
        int n_tests_failed;          /* number of failed tests */

	/* start test */
	printf("[set "); fflush(stdout);
	fprintf(TESTLOG,"\n[set]\n");

	/* reset counter */
	n_tests_failed = 0;

	/* set stop watch */
	stopwatch_lap(&watch);
  
{ /* invalid NULL ptr */
UNUR_DISTR *distr = NULL;
   distr = NULL; 


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,37,(unur_distr_cvemp_set_data( distr, NULL, 0 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,37,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,40,(unur_distr_cvemp_read_data( distr, "junk" )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,40,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid distribution type */
UNUR_DISTR *distr = NULL;
	distr = unur_distr_discr_new(); 


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,47,(unur_distr_cvemp_set_data( distr, NULL, 0 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,47,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,50,(unur_distr_cvemp_read_data( distr, "junk" )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,50,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid NULL ptr */
UNUR_DISTR *distr = NULL;
   distr = unur_distr_cvemp_new(2); 


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,57,(unur_distr_cvemp_set_data( distr, NULL, 0 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,57,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid parameters */
UNUR_DISTR *distr = NULL;
   double data[] = {1.,2.,3.,4.};
   distr = unur_distr_cvemp_new(2); 


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,65,(unur_distr_cvemp_set_data( distr, data, 0 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,65,UNUR_ERR_DISTR_SET)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,68,(unur_distr_cvemp_read_data( distr, "there-should-be-no-such-file.junk" )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,68,UNUR_ERR_GENERIC)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_set() */

/*---------------------------------------------------------------------------*/
/* [get] */

void test_get (void)
{
        int n_tests_failed;          /* number of failed tests */

	/* start test */
	printf("[get "); fflush(stdout);
	fprintf(TESTLOG,"\n[get]\n");

	/* reset counter */
	n_tests_failed = 0;

	/* set stop watch */
	stopwatch_lap(&watch);
  
{ /* invalid NULL ptr */
UNUR_DISTR *distr = NULL;
   const double *sample;
   distr = NULL; 


unur_reset_errno();
n_tests_failed += (check_expected_zero(TESTLOG,80,(unur_distr_cvemp_get_data( distr, &sample )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,80,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid distribution type */
UNUR_DISTR *distr = NULL;
   const double *sample;
   distr = unur_distr_discr_new(); 


unur_reset_errno();
n_tests_failed += (check_expected_zero(TESTLOG,88,(unur_distr_cvemp_get_data( distr, &sample )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,88,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_get() */


/*---------------------------------------------------------------------------*/
/* run generator in verifying mode */

void run_verify_generator( FILE *LOG, int line, UNUR_PAR *par )
{
	UNUR_GEN *gen;
	int i;

	/* switch to verifying mode */
	unur_distr_cvemp_set_verify(par,1);

	/* initialize generator */
	gen = unur_init( par ); abort_if_NULL(LOG, line, gen);

	/* run generator */
	for (i=0; i<VIOLATE_SAMPLE_SIZE; i++)
		unur_sample_cont(gen);

	/* destroy generator */
	unur_free(gen); 

} /* end of run_verify_generator() */

int unur_distr_cvemp_set_verify(UNUR_PAR *par ATTRIBUTE__UNUSED, int verify ATTRIBUTE__UNUSED) {return 0;}

/*---------------------------------------------------------------------------*/

int main(void)
{ 
        unsigned long seed;
	char *str_seed, *str_tail;

	/* start stop watch */
	stopwatch_init();
	stopwatch_start(&watch);

        /* open log file for unuran and set output stream for unuran messages */
        UNURANLOG = fopen( "t_distr_cvemp_unuran.log","w" );
        abort_if_NULL( stderr,-1, UNURANLOG );
        unur_set_stream( UNURANLOG );

        /* open log file for testing */
	TESTLOG = fopen( "t_distr_cvemp_test.log","w" );
	abort_if_NULL( stderr,-1, TESTLOG );

        /* seed for uniform generators */

	/* seed set by environment */
	str_seed = getenv("SEED");

	if (str_seed != NULL) {
	    seed = strtol(str_seed, &str_tail, 10);
	    if (seed == 0u) 
		seed = 7294;
	}
	else {
#ifdef SEED
	    seed = SEED;
#else
	    seed = 7294;
#endif
	}

        /* seed build-in uniform generators */
        unur_urng_MRG31k3p_seed(NULL,seed);
        unur_urng_fish_seed(NULL,seed);
	unur_urng_mstd_seed(NULL,seed);

	/* seed uniform random number generator */
#ifdef UNUR_URNG_UNURAN
#  ifdef UNUR_URNG_DEFAULT_RNGSTREAM
	{
	        unsigned long sa[6];
	        int i;
	        for (i=0; i<6; i++) sa[i] = seed;
                RngStream_SetPackageSeed(sa);
        }
#  else
	if (unur_urng_seed(NULL,seed) != UNUR_SUCCESS) {
	        fprintf(stderr,"WARNING: Seed could not be set at random\n");
                seed = ~0u;
	}
#  endif  /* UNUR_URNG_DEFAULT_RNGSTREAM */
#endif  /* UNUR_URNG_UNURAN */
 
	/* set default debugging flag */
	unur_set_default_debug(UNUR_DEBUG_ALL);

        /* detect required check mode */
        fullcheck = (getenv("UNURANFULLCHECK")==NULL) ? FALSE : TRUE;

	/* write header into log file */
        print_test_log_header( TESTLOG, seed, fullcheck );

	/* set timer for sending SIGALRM signal */
	set_alarm(TESTLOG);

	/* start test */
	printf("distr_cvemp: ");

	/* run tests */
test_new();
test_set();
test_get();


	/* test finished */
	printf("\n");  fflush(stdout);

	/* close log files */
	fprintf(TESTLOG,"\n====================================================\n\n");
	if (test_ok)
		fprintf(TESTLOG,"All tests PASSED.\n");
	else
		fprintf(TESTLOG,"Test(s) FAILED.\n");

	/* timing */
	stopwatch_print(TESTLOG,"\n<*>total time = %.0f ms\n\n", stopwatch_stop(&watch));

	fclose(UNURANLOG);
	fclose(TESTLOG);

	/* free memory */
	compare_free_memory();
	unur_urng_free(unur_get_default_urng());
	unur_urng_free(unur_get_default_urng_aux());

	/* exit */
	exit( (test_ok) ? EXIT_SUCCESS : EXIT_FAILURE );

} /* end of main */

