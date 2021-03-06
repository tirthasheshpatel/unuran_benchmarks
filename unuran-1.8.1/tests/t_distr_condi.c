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
 ** Tests for DISTR_CONDI
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

int unur_distr_condi_set_verify( UNUR_PAR *par, int verify);


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

int unur_distr_condi_set_pedantic(UNUR_PAR *par, int pedantic);
int unur_distr_condi_chg_verify(UNUR_GEN *gen, int verify);





/*---------------------------------------------------------------------------*/

#ifndef CHI2_FAILURES_TOLERATED
#  define CHI2_FAILURES_TOLERATED DEFAULT_CHI2_FAILURES_TOLERATED
#endif

/*---------------------------------------------------------------------------*/
/* [verbatim] */



int unur_distr_condi_set_pedantic(UNUR_PAR *par ATTRIBUTE__UNUSED, int pedantic ATTRIBUTE__UNUSED)  { return UNUR_FAILURE; }

int unur_distr_condi_chg_verify(UNUR_GEN *gen, int verify)
{
  if (unur_arou_chg_verify(gen,verify)==UNUR_SUCCESS) return UNUR_SUCCESS;
  if (unur_srou_chg_verify(gen,verify)==UNUR_SUCCESS) return UNUR_SUCCESS;
  if (unur_tabl_chg_verify(gen,verify)==UNUR_SUCCESS) return UNUR_SUCCESS;
  if (unur_tdr_chg_verify(gen,verify)==UNUR_SUCCESS) return UNUR_SUCCESS;
  if (unur_ars_chg_verify(gen,verify)==UNUR_SUCCESS) return UNUR_SUCCESS;
  return UNUR_FAILURE;
}	

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
n_tests_failed += (check_expected_NULL(TESTLOG,31,(unur_distr_condi_new( distr, NULL, NULL, 0 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,31,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid NULL ptr */
UNUR_DISTR *distr = NULL;
   distr = unur_distr_multinormal(3,NULL,NULL); 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,38,(unur_distr_condi_new( distr, NULL, NULL, 0 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,38,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid distribution type */
UNUR_DISTR *distr = NULL;
   double p[] = { 1., 2., 3., 4.}; 
   distr = unur_distr_discr_new(); 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,46,(unur_distr_condi_new( distr, p, NULL, 4 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,46,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid data */
UNUR_DISTR *distr = NULL;
   double p[] = { 1., 2., 3., 4.};
   distr = unur_distr_multinormal(3,NULL,NULL); 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,54,(unur_distr_condi_new( distr, p, NULL, -1 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,54,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,57,(unur_distr_condi_new( distr, p, NULL, 3 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,57,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;
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
n_tests_failed += (check_expected_setfailed(TESTLOG,68,(unur_distr_condi_set_condition( distr, NULL, NULL, 0 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,68,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid distribution type */
UNUR_DISTR *distr = NULL;
	distr = unur_distr_discr_new(); 


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,75,(unur_distr_condi_set_condition( distr, NULL, NULL, 0 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,75,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid data */
UNUR_DISTR *distr = NULL;
   double p[] = { 1., 2., 3., 4.};
   UNUR_DISTR *condi;
   distr = unur_distr_multinormal(3,NULL,NULL); 
   condi = unur_distr_condi_new( distr, p, NULL, 0 ); 


unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,85,(unur_distr_condi_set_condition( condi, p, NULL, -1 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,85,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,88,(unur_distr_condi_set_condition( condi, p, NULL, 3 )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,88,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;
unur_distr_free(condi);
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
   const double *pos, *dir;
   int k;
   distr = NULL; 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,103,(unur_distr_condi_get_distribution( distr )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,103,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,106,(unur_distr_condi_get_condition( distr, &pos, &dir, &k )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,106,UNUR_ERR_NULL)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}

{ /* invalid distribution type */
UNUR_DISTR *distr = NULL;
   const double *pos, *dir;
   int k;
   distr = unur_distr_cont_new(); 


unur_reset_errno();
n_tests_failed += (check_expected_NULL(TESTLOG,115,(unur_distr_condi_get_distribution( distr )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,115,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
n_tests_failed += (check_expected_setfailed(TESTLOG,118,(unur_distr_condi_get_condition( distr, &pos, &dir, &k )))==UNUR_SUCCESS)?0:1;
n_tests_failed += (check_errorcode(TESTLOG,118,UNUR_ERR_DISTR_INVALID)==UNUR_SUCCESS)?0:1;
unur_distr_free(distr);
}


	/* timing */
	stopwatch_print(TESTLOG,"\n<*>time = %.3f ms\n\n", stopwatch_lap(&watch));

	/* test finished */
	test_ok &= (n_tests_failed) ? 0 : 1;
	(n_tests_failed) ? printf(" ==> failed] ") : printf(" ==> ok] ");

} /* end of test_get() */

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
  
{ /* check for memory leaks */
UNUR_DISTR *distr = NULL;
   double p[] = { 1., 2., 3., 4.};
   double dir[] = { 1., -1., 2., -2.};
   UNUR_DISTR *condi;
   distr = unur_distr_multinormal(3,NULL,NULL); 
   condi = unur_distr_condi_new( distr, p, NULL, 0 ); 


unur_reset_errno();
unur_distr_condi_set_condition( condi, p, NULL, 1 );
n_tests_failed += (check_errorcode(TESTLOG,145,UNUR_SUCCESS)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
unur_distr_condi_set_condition( condi, p, dir, 1 );
n_tests_failed += (check_errorcode(TESTLOG,148,UNUR_SUCCESS)==UNUR_SUCCESS)?0:1;
unur_distr_free(condi);
unur_distr_free(distr);
}

{ /* check for memory leaks */
UNUR_DISTR *distr = NULL;
   double p[] = { 1., 2., 3., 4.};
   double dir[] = { 1., -1., 2., -2.};
   UNUR_DISTR *condi;
   distr = unur_distr_multinormal(4,NULL,NULL); 
   condi = unur_distr_condi_new( distr, p, dir, 0 ); 


unur_reset_errno();
unur_distr_condi_set_condition( condi, p, dir, 1 );
n_tests_failed += (check_errorcode(TESTLOG,161,UNUR_SUCCESS)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
unur_distr_condi_set_condition( condi, p, NULL, 1 );
n_tests_failed += (check_errorcode(TESTLOG,164,UNUR_SUCCESS)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
unur_distr_condi_set_condition( condi, p, dir, 3 );
n_tests_failed += (check_errorcode(TESTLOG,167,UNUR_SUCCESS)==UNUR_SUCCESS)?0:1;
unur_distr_free(condi);
unur_distr_free(distr);
}

{ /* check for memory leaks */
UNUR_DISTR *distr = NULL;
   double p[] = { 1., 2., 3., 4.};
   double dir[] = { 1., -1., 2., -2.};
   UNUR_DISTR *condi;
   distr = unur_distr_multinormal(4,NULL,NULL); 
   condi = unur_distr_condi_new( distr, p, NULL, 0 ); 


unur_reset_errno();
unur_distr_condi_set_condition( condi, p, NULL, 1 );
n_tests_failed += (check_errorcode(TESTLOG,180,UNUR_SUCCESS)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
unur_distr_condi_set_condition( condi, p, dir, 1 );
n_tests_failed += (check_errorcode(TESTLOG,183,UNUR_SUCCESS)==UNUR_SUCCESS)?0:1;

unur_reset_errno();
unur_distr_condi_set_condition( condi, p, NULL, 3 );
n_tests_failed += (check_errorcode(TESTLOG,186,UNUR_SUCCESS)==UNUR_SUCCESS)?0:1;
unur_distr_free(condi);
unur_distr_free(distr);
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
	UNUR_DISTR *distr[7];
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


/* distributions: 7 */
{
#define dim (3)
UNUR_DISTR *normal = unur_distr_multinormal(dim,NULL,NULL);
int i;
double p[dim]; 
for(i=0;i<dim;i++) p[i]=3.*unur_urng_sample(NULL);  
distr[0] = unur_distr_condi_new( normal, p, NULL, 0 );
unur_distr_set_name(distr[0],"condi_standardmultinormal_3");
unur_distr_cont_get_mode(distr[0]);
unur_distr_free(normal);
#undef dim
}

{
#define dim (3)
UNUR_DISTR *normal = unur_distr_multinormal(dim,NULL,NULL);
int i;
double p[dim]; 
for(i=0;i<dim;i++) p[i]=3.*unur_urng_sample(NULL);  
distr[1] = unur_distr_condi_new( normal, p, NULL, 2 );
unur_distr_set_name(distr[1],"condi_standardmultinormal_3");
unur_distr_cont_get_mode(distr[1]);
unur_distr_free(normal);
#undef dim
}

{
#define dim (4)
UNUR_DISTR *normal = unur_distr_multinormal(dim,NULL,NULL);
int i; 
double p[dim], dir[dim]; 
for(i=0;i<dim;i++) p[i]=3.*unur_urng_sample(NULL);  
for(i=0;i<dim;i++) dir[i]=0.5+unur_urng_sample(NULL);  
distr[2] = unur_distr_condi_new( normal, p, dir, 0 );
unur_distr_set_name(distr[2],"condi_standardmultinormal_4");
unur_distr_cont_get_mode(distr[2]);
unur_distr_free(normal);
#undef dim
}

{
#define dim (3)
int i; 
double p[dim], dir[dim]; 
double mean[dim], covar[dim*dim];
UNUR_DISTR *normal;
UNUR_DISTR *covar_distr;
UNUR_GEN *covar_gen;
UNUR_GEN *mean_gen;
for(i=0;i<dim;i++) p[i]=3.*unur_urng_sample(NULL);  
for(i=0;i<dim;i++) dir[i]=0.5+unur_urng_sample(NULL);  
mean_gen = unur_str2gen("normal(0,3)");
for (i=0; i<dim; i++) mean[i] = unur_sample_cont(mean_gen);
unur_free(mean_gen); 
covar_distr = unur_distr_correlation(dim);
covar_gen = unur_init(unur_mcorr_new(covar_distr));
do { unur_sample_matr(covar_gen,covar); 
   normal = unur_distr_multinormal(dim,mean,covar); 
} while (normal==NULL);
unur_distr_free(covar_distr);
unur_free(covar_gen);
distr[3] = unur_distr_condi_new( normal, p, dir, 0 );
unur_distr_set_name(distr[3],"condi_multinormal_random");
unur_distr_cont_get_mode(distr[3]);
unur_distr_free(normal);
#undef dim
}

{
#define dim (3)
int i; 
double p[dim], dir[dim]; 
double ll[3] = {0.,0.,0.};
double ru[3] = {UNUR_INFINITY,UNUR_INFINITY,UNUR_INFINITY};
UNUR_DISTR *normal = unur_distr_multinormal(dim,NULL,NULL);
unur_distr_cvec_set_domain_rect(normal,ll,ru);
for(i=0;i<dim;i++) p[i]=0.01 + 3.*unur_urng_sample(NULL);  
for(i=0;i<dim;i++) dir[i]=0.5+unur_urng_sample(NULL);  
distr[4] = unur_distr_condi_new( normal, p, dir, 0 );
unur_distr_set_name(distr[4],"condi_standardmultinormal_domain");
unur_distr_cont_get_mode(distr[4]);
unur_distr_free(normal);
#undef dim
}

{
#define dim (3)
int i; 
double p[dim]; 
double ll[3] = {0.,0.,0.};
double ru[3] = {UNUR_INFINITY,UNUR_INFINITY,UNUR_INFINITY};
UNUR_DISTR *normal = unur_distr_multinormal(dim,NULL,NULL);
unur_distr_cvec_set_domain_rect(normal,ll,ru);
for(i=0;i<dim;i++) p[i]=0.01 + 3.*unur_urng_sample(NULL);  
distr[5] = unur_distr_condi_new( normal, p, NULL, 2 );
unur_distr_set_name(distr[5],"condi_standardmultinormal_domain");
unur_distr_cont_get_mode(distr[5]);
unur_distr_free(normal);
#undef dim
}

{
#define dim (3)
double p[dim]; 
double ll[3] = {-1.,0.,1.};
double ru[3] = {1.,1.,2.};
UNUR_DISTR *normal = unur_distr_multinormal(dim,NULL,NULL);
unur_distr_cvec_set_domain_rect(normal,ll,ru);
p[0] = -0.5 + unur_urng_sample(NULL);  
p[1] = unur_urng_sample(NULL);  
p[2] = 1.01 + unur_urng_sample(NULL);  
distr[6] = unur_distr_condi_new( normal, p, NULL, 2 );
unur_distr_set_name(distr[6],"condi_standardmultinormal_domain");
unur_distr_cont_get_mode(distr[6]);
unur_distr_free(normal);
#undef dim
}

	/* timing */
	stopwatch_print(TESTLOG,"\n<*>setup time = %.3f ms\n", stopwatch_lap(&watch));

	printf("\n(verify hat) "); fflush(stdout);

/* verify hat tests: 35 */

	unur_set_default_debug(~UNUR_DEBUG_SAMPLE);
	fprintf( TESTLOG,"\nVerify Hat Test (squeeze <= PDF <= hat):\n");

/* distribution [0] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[0]);
par = unur_arou_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[0],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[0]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[0],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[0]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
unur_tdr_set_c(par,0.);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[0],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[0]);
par = unur_ars_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[0],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[0]);
par = unur_tabl_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[0],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [1] */

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_arou_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[1],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[1],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
unur_tdr_set_c(par,0.);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[1],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_ars_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[1],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[1]);
par = unur_tabl_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[1],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [2] */

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_arou_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[2],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[2],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
unur_tdr_set_c(par,0.);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[2],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_ars_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[2],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[2]);
par = unur_tabl_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[2],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [3] */

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_arou_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[3],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[3],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
unur_tdr_set_c(par,0.);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[3],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_ars_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[3],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[3]);
par = unur_tabl_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[3],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [4] */

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_arou_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[4],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[4],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
unur_tdr_set_c(par,0.);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[4],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_ars_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[4],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[4]);
par = unur_tabl_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[4],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

/* distribution [5] */

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[5]);
par = unur_arou_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[5],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[5]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[5],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[5]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
unur_tdr_set_c(par,0.);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[5],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[5]);
par = unur_ars_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[5],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	printf("."); fflush(stdout);
	}

/* distribution [6] */

	if(fullcheck) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[6]);
par = unur_arou_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[6],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[6]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[6],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[6]);
par = unur_tdr_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
unur_tdr_set_c(par,0.);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[6],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[6]);
par = unur_ars_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[6],'~')==UNUR_SUCCESS)?0:1000;
	unur_free(gen);

	unur_distr_free(distr_localcopy);
	} while (0);
	}

	if(TRUE) {
	unur_reset_errno();
	do {
	UNUR_DISTR *distr_localcopy = unur_distr_clone(distr[6]);
par = unur_tabl_new(distr_localcopy);
	unur_distr_condi_set_pedantic(par,0);
	gen = unur_init(par);
	if (gen) unur_distr_condi_chg_verify(gen,1);
	n_tests_failed += (run_validate_verifyhat(TESTLOG,0,gen,distr[6],'~')==UNUR_SUCCESS)?0:1000;
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
	unur_distr_free(distr[5]);
	unur_distr_free(distr[6]);

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
	unur_distr_condi_set_verify(par,1);

	/* initialize generator */
	gen = unur_init( par ); abort_if_NULL(LOG, line, gen);

	/* run generator */
	for (i=0; i<VIOLATE_SAMPLE_SIZE; i++)
		unur_sample_cont(gen);

	/* destroy generator */
	unur_free(gen); 

} /* end of run_verify_generator() */

int unur_distr_condi_set_verify(UNUR_PAR *par ATTRIBUTE__UNUSED, int verify ATTRIBUTE__UNUSED) {return 0;}

/*---------------------------------------------------------------------------*/

int main(void)
{ 
        unsigned long seed;
	char *str_seed, *str_tail;

	/* start stop watch */
	stopwatch_init();
	stopwatch_start(&watch);

        /* open log file for unuran and set output stream for unuran messages */
        UNURANLOG = fopen( "t_distr_condi_unuran.log","w" );
        abort_if_NULL( stderr,-1, UNURANLOG );
        unur_set_stream( UNURANLOG );

        /* open log file for testing */
	TESTLOG = fopen( "t_distr_condi_test.log","w" );
	abort_if_NULL( stderr,-1, TESTLOG );

        /* seed for uniform generators */

	/* seed set by environment */
	str_seed = getenv("SEED");

	if (str_seed != NULL) {
	    seed = strtol(str_seed, &str_tail, 10);
	    if (seed == 0u) 
		seed = 679209;
	}
	else {
#ifdef SEED
	    seed = SEED;
#else
	    seed = 679209;
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
	printf("distr_condi: ");

	/* run tests */
test_new();
test_set();
test_get();
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

