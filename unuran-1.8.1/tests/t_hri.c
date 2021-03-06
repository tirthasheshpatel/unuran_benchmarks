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
 ** Tests for HRI
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

double HR_invalid(double x, const UNUR_DISTR *distr);

double HR_decreasing(double x, const UNUR_DISTR *distr);
double CDF_decreasing(double x, const UNUR_DISTR *distr);

double HR_constant(double x, const UNUR_DISTR *distr);
double CDF_constant(double x, const UNUR_DISTR *distr);

double HR_increasing_wb31(double x, const UNUR_DISTR *distr);
double CDF_increasing_wb31(double x, const UNUR_DISTR *distr);

double HR_increasing_gm31(double x, const UNUR_DISTR *distr);
double CDF_increasing_gm31(double x, const UNUR_DISTR *distr);

double HR_increasing_gm51(double x, const UNUR_DISTR *distr);
double CDF_increasing_gm51(double x, const UNUR_DISTR *distr);

int unur_hri_set_pedantic( UNUR_PAR *par, int pedantic );

#define COMPARE_SAMPLE_SIZE   (10000)
#define VIOLATE_SAMPLE_SIZE   (20)




/*---------------------------------------------------------------------------*/

#ifndef CHI2_FAILURES_TOLERATED
#  define CHI2_FAILURES_TOLERATED DEFAULT_CHI2_FAILURES_TOLERATED
#endif

/*---------------------------------------------------------------------------*/
/* [verbatim] */



double HR_invalid(double x ATTRIBUTE__UNUSED, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* example with invalid hazard rate */
{ return 0.; }


double HR_decreasing(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* example with decreasing hazard rate */
{ return (1./(1.+x)); }

double CDF_decreasing(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* CDF for example with decreasing hazard rate */
{ return (x/(1. + x)); }


double HR_constant(double x ATTRIBUTE__UNUSED, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* example with constant hazard rate: standard expontential */
{ return 1.; }

double CDF_constant(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* CDF for example with constant hazard rate: standard expontential */
{ return (1.-exp(-x)); }


double HR_increasing_wb31(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* example with increasing hazard rate: Weibull (alpha=3,beta=1) */
{ return (3*x*x); }

double CDF_increasing_wb31(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* CDF for example with increasing hazard rate: Weibull (alpha=3,beta=1) */
{ return (1.-exp(-x*x*x)); }


double HR_increasing_gm31(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* example with increasing hazard rate: gamma (alpha=3,beta=1) */
{ return (x*x)/(x*x+2*x+2.); }

double CDF_increasing_gm31(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* CDF for example with increasing hazard rate: gamma (alpha=3,beta=1) */
{ return 1.-(x*x*exp(-x))/2.-(x*exp(-x))-exp(-x); }


double HR_increasing_gm51(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* example with increasing hazard rate: gamma (alpha=5,beta=1) */
{ return pow(x,4.)/(pow(x,4.)+4*pow(x,3.)+12.*x*x+24.*x+24.); }

double CDF_increasing_gm51(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
     /* CDF for example with increasing hazard rate: gamma (alpha=5,beta=1) */
{ return (-exp(-x)/24.)*((pow(x,4.)+4*pow(x,3.)+12.*x*x+24.*x+24.))+1.; }


/* dummy function */
int unur_hri_set_pedantic( UNUR_PAR *par ATTRIBUTE__UNUSED, int pedantic ATTRIBUTE__UNUSED)
{ return 1; }

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
  
{ /* invalid NULL ptr */
UNUR_DISTR *distr = NULL;
   distr = NULL; 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,45,(unur_hri_new( distr )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,45,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid distribution type */
UNUR_DISTR *distr = NULL;
   distr = unur_distr_discr_new(); 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,51,(unur_hri_new( distr )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,51,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* data missing in distribution object */
UNUR_DISTR *distr = NULL;
   distr = unur_distr_cont_new(); 


unur_reset_errno();
/* hazard rate */
n_tests_failed += (check_expected_NULL(TESTLOG,58,(unur_hri_new( distr )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,58,UNUR_ERR_DISTR_REQUIRED)==UNUR_SUCCESS)?0:1;
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
UNUR_PAR   *par = NULL;
   par = NULL; 


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,68,(unur_hri_set_verify(par,1)))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,68,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,71,(unur_hri_set_p0(par,1.)))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,71,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
}

{ /* invalid parameter object */
UNUR_DISTR *distr = NULL;
UNUR_PAR   *par = NULL;
   distr = unur_distr_normal(NULL,0);
   par = unur_arou_new(distr); 


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,79,(unur_hri_set_verify(par,1)))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,79,UNUR_ERR_PAR_INVALID)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,82,(unur_hri_set_p0(par,1.)))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,82,UNUR_ERR_PAR_INVALID)==UNUR_SUCCESS)?0:1;
unur_par_free(par);
unur_distr_free(distr);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_set() */

/*---------------------------------------------------------------------------*/
/* [chg] */

void test_chg (void)
{
        int n_tests_failed;          /* number of failed tests */

	/* start test */
	printf("[chg "); fflush(stdout);
	fprintf(TESTLOG,"\n[chg]\n");

	/* reset counter */
	n_tests_failed = 0;

	/* set stop watch */
	stopwatch_lap(&watch);
  
{ /* invalid generator object */
UNUR_DISTR *distr = NULL;
UNUR_PAR   *par = NULL;
UNUR_GEN   *gen = NULL;
   distr = unur_distr_normal(NULL,0);
   par = unur_arou_new(distr);
   unur_set_debug(par,0);
   gen = unur_init( par ); 
abort_if_NULL(TESTLOG, 98,    gen );


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,102,(unur_hri_chg_verify(gen,1)))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,102,UNUR_ERR_GEN_INVALID)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
unur_free(gen);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_chg() */

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
  
{ /* left border no valid upper bound */
UNUR_GEN   *gen = NULL;
UNUR_DISTR *distr = NULL;
UNUR_PAR   *par = NULL;
   gen = NULL;
   distr = unur_distr_cont_new();
   unur_distr_cont_set_hr(distr,HR_invalid);
   par = unur_hri_new(distr); 
abort_if_NULL(TESTLOG, 112,    par );


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,116,(gen = unur_init( par )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,116,UNUR_ERR_GEN_CONDITION)==UNUR_SUCCESS)?0:1;
unur_free(gen);
unur_distr_free(distr);
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
  
{ /* does not exist */
UNUR_DISTR *distr = NULL;
UNUR_PAR   *par = NULL;
UNUR_GEN   *gen = NULL;
   distr = unur_distr_cont_new();
   unur_distr_cont_set_hr(distr,HR_increasing_wb31);
   par = unur_hri_new(distr);
   gen = unur_init( par ); 
abort_if_NULL(TESTLOG, 126,    gen );


unur_reset_errno();
n_tests_failed += (check_expected_reinit(TESTLOG,130,(unur_reinit( gen )))==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
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
UNUR_DISTR *distr = NULL;
UNUR_PAR   *par = NULL;
    distr = unur_distr_cont_new();
    unur_distr_cont_set_hr(distr,HR_increasing_wb31);
    par = NULL; 


unur_reset_errno();
/* default algorithm */
par = unur_hri_new(distr);
n_tests_failed += (compare_sequence_par_start(TESTLOG,143,par,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
/* default algorithm - verifying mode */
par = unur_hri_new(distr);
unur_hri_set_verify(par,1);
n_tests_failed += (compare_sequence_par(TESTLOG,148,par,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
/* default algorithm - ignore invalid boundary */
unur_distr_cont_set_domain(distr,-10.,10.);
par = unur_hri_new(distr);
n_tests_failed += (compare_sequence_par_start(TESTLOG,153,par,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* compare reinit */
UNUR_DISTR *distr = NULL;
UNUR_PAR   *par = NULL;
UNUR_GEN   *gen = NULL;
   distr = unur_distr_cont_new();
   unur_distr_cont_set_hr(distr,HR_increasing_wb31);
   par = NULL;
   gen = NULL; 


unur_reset_errno();
/* original generator object */
par = unur_hri_new(distr);
gen = unur_init(par);
n_tests_failed += (compare_sequence_gen_start(TESTLOG,166,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
/* reinit */
unur_reinit(gen);
n_tests_failed += (compare_sequence_gen(TESTLOG,170,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
unur_free(gen);
}

{ /* compare stringparser */
UNUR_DISTR *distr = NULL;
UNUR_PAR   *par = NULL;
UNUR_GEN   *gen = NULL;
   distr = NULL;
   par = NULL;
   gen = NULL; 


unur_reset_errno();
distr = unur_distr_cont_new();
unur_distr_cont_set_hrstr(distr,"3*x^2");
par = unur_hri_new(distr);
gen = unur_init(par);
n_tests_failed += (compare_sequence_gen_start(TESTLOG,183,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
unur_free(gen); gen = NULL;
unur_distr_free(distr); distr = NULL;
gen = unur_str2gen( "cont; hr = \"3*x^2\"; domain = (0,inf)& \
	             method = hri" );
n_tests_failed += (compare_sequence_gen(TESTLOG,189,gen,COMPARE_SAMPLE_SIZE)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
unur_free(gen);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_sample() */

/*---------------------------------------------------------------------------*/
/* [validate] */

/*---------------------------------------------------------------------------*/

/* [validate] */

void test_validate (void)
{
	UNUR_DISTR *distr[5];
	UNUR_PAR *par;
	UNUR_GEN *gen;
	int n_tests_failed;
	int rcode;
	double *darray;
	double fpm[10];

	rcode = 0;

	/* start test */
	printf("[validate "); fflush(stdout);
	fprintf(TESTLOG,"\n[validate]\n");

	/* reset counter */
	n_tests_failed = 0;

	/* set stop watch */
	stopwatch_lap(&watch);


/* distributions: 5 */
{
distr[0] = unur_distr_cont_new();
unur_distr_set_name(distr[0],"decreasing HR");
unur_distr_cont_set_hr(distr[0],HR_decreasing);
unur_distr_cont_set_cdf(distr[0],CDF_decreasing);
}

{
distr[1] = unur_distr_cont_new();
unur_distr_set_name(distr[1],"constant HR");
unur_distr_cont_set_hr(distr[1],HR_constant);
unur_distr_cont_set_cdf(distr[1],CDF_constant);
}

{
distr[2] = unur_distr_cont_new();
unur_distr_set_name(distr[2],"increasing HR wb31");
unur_distr_cont_set_hr(distr[2],HR_increasing_wb31);
unur_distr_cont_set_cdf(distr[2],CDF_increasing_wb31);
}

{
distr[3] = unur_distr_cont_new();
unur_distr_set_name(distr[3],"increasing HR gm31");
unur_distr_cont_set_hr(distr[3],HR_increasing_gm31);
unur_distr_cont_set_cdf(distr[3],CDF_increasing_gm31);
}

{
distr[4] = unur_distr_cont_new();
unur_distr_set_name(distr[4],"increasing HR gm51");
unur_distr_cont_set_hr(distr[4],HR_increasing_gm51);
unur_distr_cont_set_cdf(distr[4],CDF_increasing_gm51);
}

	/* timing */
	stopwatch_print(TESTLOG,"\n<*>setup time = %.3f ms\n", stopwatch_lap(&watch));

	printf("\n(chi^2) "); fflush(stdout);

/* chi^2 tests: 10 */

	unur_set_default_debug(~UNUR_DEBUG_SAMPLE);
	fprintf( TESTLOG,"\nChi^2 Test:\n");

/* distribution [0] */

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[0]);
par = unur_hri_new(distr_localcopy);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[0],'-');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[0]);
par = unur_hri_new(distr_localcopy);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[0],'-');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [1] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_hri_new(distr_localcopy);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[1],'+');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_hri_new(distr_localcopy);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[1],'+');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [2] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_hri_new(distr_localcopy);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[2],'+');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_hri_new(distr_localcopy);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[2],'+');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [3] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_hri_new(distr_localcopy);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[3],'+');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_hri_new(distr_localcopy);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[3],'+');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [4] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_hri_new(distr_localcopy);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[4],'+');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_hri_new(distr_localcopy);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	rcode = run_validate_chi2(TESTLOG,0,gen,distr[4],'+');
	n_tests_failed += (rcode==UNUR_SUCCESS)?0:1;
	n_tests_failed += (rcode==UNUR_FAILURE)?1000:0;
	unur_free(gen);
	unur_distr_free(distr_localcopy);
	} while (0);
	}

	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.0f ms\n", stopwatch_lap(&watch));

	printf("\n(verify hat) "); fflush(stdout);

/* verify hat tests: 10 */

	unur_set_default_debug(~UNUR_DEBUG_SAMPLE);
	fprintf( TESTLOG,"\nVerify Hat Test (squeeze <= PDF <= hat):\n");

/* distribution [0] */

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[0]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[0],'-')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(fullcheck) {
	printf("."); fflush(stdout);
	}

/* distribution [1] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[1],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[1],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [2] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[2],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[2],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [3] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[3],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[3],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [4] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[4],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_hri_new(distr_localcopy);
	unur_hri_set_pedantic(par,0);
unur_hri_set_p0(par,2.);
	gen = unur_init(par);
	if (gen) unur_hri_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[4],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.0f ms\n", stopwatch_lap(&watch));


/* free distributions */
	unur_distr_free(distr[0]);
	unur_distr_free(distr[1]);
	unur_distr_free(distr[2]);
	unur_distr_free(distr[3]);
	unur_distr_free(distr[4]);

	/* test finished */
	test_ok &= (n_tests_failed>CHI2_FAILURES_TOLERATED) ? 0 : 1;
	/* we accept CHI2_FAILURES_TOLERATED failures */
	(n_tests_failed>CHI2_FAILURES_TOLERATED) ? printf(" ==> failed] ") : printf(" ==> ok] ");

	/* prevent compiler from making useless annoying warnings */
	distr[0] = NULL;
	par = NULL;
	gen = NULL;
	darray = NULL;
	fpm[0] = 0.;

} /* end of test_validate */


/*---------------------------------------------------------------------------*/
/* run generator in verifying mode */

void run_verify_generator( FILE *LOG, int line, UNUR_PAR *par )
{
	UNUR_GEN *gen;
	int i;

	/* switch to verifying mode */
	unur_hri_set_verify(par,1);

	/* initialize generator */
	gen = unur_init( par ); abort_if_NULL(LOG, line, gen);

	/* run generator */
	for (i=0; i<VIOLATE_SAMPLE_SIZE; i++)
		unur_sample_cont(gen);

	/* destroy generator */
	unur_free(gen); 

} /* end of run_verify_generator() */


/*---------------------------------------------------------------------------*/

int main(void)
{ 
        unsigned long seed;
	char *str_seed, *str_tail;

	/* start stop watch */
	stopwatch_init();
	stopwatch_start(&watch);

        /* open log file for unuran and set output stream for unuran messages */
        UNURANLOG = fopen( "t_hri_unuran.log","w" );
        abort_if_NULL( stderr,-1, UNURANLOG );
        unur_set_stream( UNURANLOG );

        /* open log file for testing */
	TESTLOG = fopen( "t_hri_test.log","w" );
	abort_if_NULL( stderr,-1, TESTLOG );

        /* seed for uniform generators */

	/* seed set by environment */
	str_seed = getenv("SEED");

	if (str_seed != NULL) {
	    seed = strtol(str_seed, &str_tail, 10);
	    if (seed == 0u) 
		seed = 106801;
	}
	else {
#ifdef SEED
	    seed = SEED;
#else
	    seed = 106801;
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
	printf("hri: ");

	/* run tests */
test_new();
test_set();
test_chg();
test_init();
test_reinit();
test_sample();
test_validate();


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

