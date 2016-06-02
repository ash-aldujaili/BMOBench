#include <stdio.h>
#include <stdbool.h>
#include "globaldeclare.h"
#include "recorder.h"

// Global declarations
size_t nVars;
size_t nObjs;
size_t maxFuncEvals;
bool lastRecord = false;
size_t problem_id;
size_t iRun;
size_t *func_id;
size_t *inst_id;
char *algoName; // name of the running algorithm
struct SolutionsArchive *myarchive;
char outputDir[LENFULLFNAME] = "../EXP_RESULTS"; // string storing the output directory
size_t nPop =20;
double lobound;
double upbound;
double pcross;
double eta_cross;
double eta_mut;


/**
 * Create a duplicate copy of ${string} and return a pointer to it
 * The caller is responsible for free()ing the memory allocated.
 */
char *my_strdup(const char *string) {
    if (string == NULL)
        return NULL;
    size_t len = strlen(string);
    char *dup = malloc(len + 1);
    memcpy(dup, string, len + 1);
    return dup;
}
