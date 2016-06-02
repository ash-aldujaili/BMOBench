#ifndef GLOBALDECLARE_H
#define	GLOBALDECLARE_H
#include <stdlib.h>
#include <stdbool.h> // to use the bool datatype, required C99
#include <string.h>
//#ifndef _AVL_H
//#include "avl.h"
//#endif

#ifdef	__cplusplus
extern "C" {
#endif


#ifndef NaN
#define NaN    999999999
#endif

#ifndef INF
#define INF    1.0e14
#endif

#ifndef LENFULLFNAME
#define LENFULLFNAME    300  //length of the full file name with location
#endif

#ifndef LENFNAME
#define LENFNAME        1024  //predefined maximum length of a file name
#endif


/*
=============================================
 Set the number of runs and budget multiplier
 ============================================
 */
#define BUDGET_MULTIPLIER 100 /*this makes the function evaluation budget = BUDGET_MULTIPLIER * dimension of the problem */
#define NUM_RUNS 10


// Global declarations
extern size_t nPop;
extern size_t nVars;
extern size_t nObjs;
extern size_t maxFuncEvals;
extern bool lastRecord;
//extern size_t nRuns;
extern size_t iRun;
extern size_t problem_id;
extern size_t *func_id;
extern size_t *inst_id;
extern char *algoName;  // name of the running algorithm
extern struct SolutionsArchive *myarchive;
//extern avl_tree_t *avltree_archive;
extern char outputDir[LENFULLFNAME]; // string storing the output directory

extern double lobound;
extern double upbound;
extern double pcross;
extern double eta_cross;
extern double eta_mut;



/**
 * Create a duplicate copy of ${string} and return a pointer to it
 * The caller is responsible for free()ing the memory allocated.
 */
char *my_strdup(const char *string);
void instance_constructor(void);
void instance_destructor(void);
void initInstance(void);
void initInstanceFuzzy(void);
void finiInstance(void);
void finiInstanceFuzzy(void);
void freeStarStarDbl(double ** A, int r);
void freeStarStarInt(int ** A, int r);
void freeStarStarChar(char ** A, int r);
void free3StarsDbl(double *** A, int r, int c);
void free3StarsInt(int *** A, int r, int c);
void init_decodingGT(void);
void fini_decodingGT(void);
double decodingGT(int ** sched, int * chrom);
void init_decodingGTcrisp(void);
void fini_decodingGTcrisp(void);
double decodingGTcrisp(int ** sched, int * chrom);
void init_decodingGTfuzzy(void);
void fini_decodingGTfuzzy(void);
void decodingGTfuzzy(double * makespan, int ** sched, int * chrom);
void induce_fuzzy_makespan(double * makespan, int ** sched);
void JOX(int **pop, int nPop, int nVar, double pcross);
void GOX(int **pop, int nPop, int nVar, double pcross);
void GPMX(int **pop, int nPop, int nVar, double pcross);
void PPX(int **pop, int nPop, int nVar, double pcross, int fashion);
void exchangemutation(int **pop, int nPop, int nVar, double pmut);
void swapmutation(int **pop, int nPop, int nVar, double pmut);
void swapmutation2(int **pop, int nPop, int nVar, double pmut);
void swapmutation3(int **pop, int nPop, int nVar, double pmut);
void initializationPR(int **pop, int nPop, int nVar);
void tournamentselect(int **offspring, int **pop, double *obj, int nPop, int nVar, int toursize);
void generationalselect(int **parent_pop, double *parent_obj, int ***parent_sched, int **child_pop, double *child_obj, int ***child_sched, int nPop, int nVar);
void generationalselect2(int **parent_pop, double *parent_obj, double **parent_makespan, int ***parent_sched, int **child_pop, double *child_obj, double **child_makespan, int ***child_sched, int nPop, int nVar);
void replaceworstselect(int **parent_pop, double *parent_obj, int ***parent_sched, int **child_pop, double *child_obj, int ***child_sched, int nPop, int nVar);
void eval_pop_crisp(double *makespan, int ***sched, int **pop, int nPop, int nVar);
void eval_pop_fuzzy_modal(double *makespan, int ***sched, int **pop, int nPop, int nVar);
void eval_pop_fuzzy_expected(double **makespan, double *obj, int ***sched, int **pop, int nPop, int nVar);
void eval_pop_fuzzy_mo(double **makespan, double **obj, int ***sched, int **pop, int nPop, int nVar);

void check_const_params(void);
void read_instance(const char *command, int arg);
void generate_instance(void);
void removeRow1DInt(int *V, int len, int idx);
void removeRow2DInt(int **M, int nrow, int ncol, int idx);
int min_vectorDbl(double *V, int n);
int max_vectorDbl(double *V, int n);
int min_vectorInt(int *V, int n);
void min_2DarrayDbl(int *idx, double **A, int r, int c, int dir);
void QuickSort_vecInt(int *idx, int *V, int n);
void QuickSort_vecDbl(int *idx, double *V, int n);


#ifdef	__cplusplus
}
#endif

#endif
