/*
	file automatically generated by make_test_files.pl
	Tue Apr 19 14:01:03 2011
*/

/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************/
    
/**
 ** Tests for UNIF
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

int unur_unif_set_verify( UNUR_PAR *par, int verify);


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

#define COMPARE_SAMPLE_SIZE   (10000)
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
  

	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_new() */

/*---------------------------------------------------------------------------*/
/* [init] */

void test_init (void)
{
        int n_tests_failed;          /* number of failed tests */

	/* start test */
	printf("[init "); fflush(stdout);
	fprintf(TESTLOG,"\n[init]\n");

	/* reset counter */
	n_tests_failed = 0;

	/* set stop watch */
	stopwatch_lap(&watch);
  
{ /* invalid NULL ptr */
UNUR_PAR   *par = NULL;
   par = NULL; 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,42,(unur_init( par )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,42,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_init() */

/*---------------------------------------------------------------------------*/
/* [reinit] */

void test_reinit (void)
{
        int n_tests_failed;          /* number of failed tests */

	/* start test */
	printf("[reinit "); fflush(stdout);
	fprintf(TESTLOG,"\n[reinit]\n");

	/* reset counter */
	n_tests_failed = 0;

	/* set stop watch */
	stopwatch_lap(&watch);
  
{ /* exist */
UNUR_PAR   *par = NULL;
UNUR_GEN   *gen = NULL;
   par = unur_unif_new(NULL);
   gen = unur_init( par ); 
abort_if_NULL(TESTLOG, 50,    gen );


unur_reset_errno();
n_tests_failed += (check_expected_reinit(TESTLOG,54,(unur_reinit( gen )))==UNUR_SUCCESS)?0:1;
unur_free(gen);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_reinit() */

/*---------------------------------------------------------------------------*/
/* [sample] */

void test_sample (void)
{
        int n_tests_failed;          /* number of failed tests */

	/* start test */
	printf("[sample "); fflush(stdout);
	fprintf(TESTLOG,"\n[sample]\n");

	/* reset counter */
	n_tests_failed = 0;

	/* set stop watch */
	stopwatch_lap(&watch);
  
{ /* compare */
UNUR_PAR   *par = NULL;
   par = NULL; 


unur_reset_errno();
;
compare_sequence_urng_start(TESTLOG,64,COMPARE_SAMPLE_SIZE);

unur_reset_errno();
par = unur_unif_new(NULL);
n_tests_failed += (compare_sequence_par(TESTLOG,67,par,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;
}

{ /* compare clone */
UNUR_PAR   *par = NULL;
UNUR_GEN   *gen = NULL;
   UNUR_GEN *clone;
   par = NULL;
   gen = NULL; 


unur_reset_errno();
/* original generator object */
par = unur_unif_new(NULL);
gen = unur_init(par);
n_tests_failed += (compare_sequence_gen_start(TESTLOG,79,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
/* clone */
clone = unur_gen_clone(gen);
unur_free(gen);
gen = clone;
n_tests_failed += (compare_sequence_gen(TESTLOG,85,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;
unur_free(gen);
}

{ /* compare reinit */
UNUR_PAR   *par = NULL;
UNUR_GEN   *gen = NULL;
   par = NULL;
   gen = NULL; 


unur_reset_errno();
/* original generator object */
par = unur_unif_new(NULL);
gen = unur_init(par);
n_tests_failed += (compare_sequence_gen_start(TESTLOG,96,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
/* reinit */
unur_reinit(gen);
n_tests_failed += (compare_sequence_gen(TESTLOG,100,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;
unur_free(gen);
}

{ /* compare stringparser */
UNUR_PAR   *par = NULL;
UNUR_GEN   *gen = NULL;
   par = NULL;
   gen = NULL; 


unur_reset_errno();
par = unur_unif_new(NULL);
gen = unur_init(par);
n_tests_failed += (compare_sequence_gen_start(TESTLOG,110,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
unur_free(gen); gen = NULL;
gen = unur_str2gen( "uniform & method = unif" );
n_tests_failed += (compare_sequence_gen(TESTLOG,114,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;
unur_free(gen);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_sample() */


/*---------------------------------------------------------------------------*/
/* run generator in verifying mode */

void run_verify_generator( FILE *LOG, int line, UNUR_PAR *par )
{
	UNUR_GEN *gen;
	int i;

	/* switch to verifying mode */
	unur_unif_set_verify(par,1);

	/* initialize generator */
	gen = unur_init( par ); abort_if_NULL(LOG, line, gen);

	/* run generator */
	for (i=0; i<VIOLATE_SAMPLE_SIZE; i++)
		unur_sample_cont(gen);

	/* destroy generator */
	unur_free(gen); 

} /* end of run_verify_generator() */

int unur_unif_set_verify(UNUR_PAR *par ATTRIBUTE__UNUSED, int verify ATTRIBUTE__UNUSED) {return 0;}

/*---------------------------------------------------------------------------*/

int main(void)
{ 
        unsigned long seed;
	char *str_seed, *str_tail;

	/* start stop watch */
	stopwatch_init();
	stopwatch_start(&watch);

        /* open log file for unuran and set output stream for unuran messages */
        UNURANLOG = fopen( "t_unif_unuran.log","w" );
        abort_if_NULL( stderr,-1, UNURANLOG );
        unur_set_stream( UNURANLOG );

        /* open log file for testing */
	TESTLOG = fopen( "t_unif_test.log","w" );
	abort_if_NULL( stderr,-1, TESTLOG );

        /* seed for uniform generators */

	/* seed set by environment */
	str_seed = getenv("SEED");

	if (str_seed != NULL) {
	    seed = strtol(str_seed, &str_tail, 10);
	    if (seed == 0u) 
		seed = 773681;
	}
	else {
#ifdef SEED
	    seed = SEED;
#else
	    seed = 773681;
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
	printf("unif: ");

	/* run tests */
test_new();
test_init();
test_reinit();
test_sample();


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

