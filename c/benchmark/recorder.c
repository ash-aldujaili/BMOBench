#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "globaldeclare.h"
#include "bmobenchProblems.h"
#include "checkDominance.h"
#include "recorder.h"
#include <sys/types.h>
#include <sys/stat.h>
/*
 Recorder
*/
#define IGNORE_DECISION_VARIABLES
//#define WRITE_WITH_MINIMAL_SPACING



static archive_t *archive = NULL;

void freeArchive(){
  node_t *temp_node;
  if (archive == NULL) return;
  temp_node = archive->head;
  while(temp_node !=NULL){
    // free head
    free(temp_node->obj);
    free(temp_node->var);
    archive->head = temp_node->next;
    free(temp_node);
    temp_node = archive->head;
  }
  // free tail
  //free(archive->tail->obj);
  //free(archive->tail->var);
  //free(archive->tail);
  free(archive);
  archive = NULL;
  lastRecord = false;
  //printf("cleared\n");
  //if (archive == NULL) printf("nullified\n");
}

void updateArchive(double *var, double *obj,size_t timestamp){
   //printf("time stamp %zu\n", timestamp);
   if (archive == NULL){
     //printf("problem %zu\n", problem_id);
     archive = malloc(sizeof(archive_t));
     archive->head = malloc(sizeof(node_t));
     archive->tail = archive->head;
     archive->tail->obj = malloc(sizeof(double) * nObjs);
     archive->tail->var = malloc(sizeof(double) * nVars);
     archive->tail->timestamp = timestamp;
     archive->tail->next = NULL;
   }
   else {
     // check for dominance:
     node_t *temp_node = archive->head;
     while(temp_node != NULL){
       // choose variables
       if(checkDominance(obj,temp_node->obj, nObjs) < 0) return;
       temp_node = temp_node->next;
     }
     // append
     archive->tail->next = malloc(sizeof(node_t));
     archive->tail = archive->tail->next;
     archive->tail->obj = malloc(sizeof(double) * nObjs);
     archive->tail->var = malloc(sizeof(double) * nVars);
     archive->tail->timestamp = timestamp;
     archive->tail->next = NULL;
   }
   for(size_t j=0; j < nObjs; j++) archive->tail->obj[j]=obj[j];
   for(size_t j=0; j < nVars; j++) archive->tail->var[j]=var[j];
   //printf("Update done\n");
}


void printArchive(){
  node_t *temp_node;
  // file stuff
  FILE *logfile = NULL;
  char *buf;
  char *fullLogFilename;
  buf = (char *) malloc(sizeof(char) * LENFNAME);
  fullLogFilename = (char*) malloc(LENFULLFNAME * sizeof(char)); // filename of the output text file
  //printf("Printing ..\n");
  sprintf(buf, "%s_%zuD_%s_nfev%.1e_run%zu.txt", pInfo[problem_id].name, nVars, algoName, (double) BUDGET_MULTIPLIER, iRun);
 //printf("Printing ..\n");
  sprintf(fullLogFilename, "%s/%s", outputDir, buf); // system dependent for "/"
  logfile = fopen(fullLogFilename, "w");
 //printf("Printing ..\n");
  // header print
/*   fprintf(logfile, "# Logging data for %s.\n", algoName);
  fprintf(logfile, "# Problem setting  : problem = %s\n", pInfo[problem_id].name);
  fprintf(logfile, "# Algorithm setting: %s with popsize = %zu\n", algoName, nPop);
  fprintf(logfile, "# Budget setting   : stop at max FEs = %.1e\n", (double)maxFuncEvals);
  fprintf(logfile, "# *****************************************************************************\n"); */
  //printf("Printing ..\n");
  #ifndef IGNORE_DECISION_VARIABLES
    fprintf(logfile, "# timestamp | # %zu variables  |  %zu objectives    \n", nVars, nObjs);
  #else
    fprintf(logfile, "# timestamp |# %zu objectives  \n", nObjs);
  #endif
  // scan and print
  temp_node = archive->head;
  while(temp_node != NULL){
    #ifdef WRITE_WITH_MINIMAL_SPACING
    fprintf(logfile, "%zu ", temp_node->timestamp);  // its timestamp (FEval)
    #else
    fprintf(logfile, "%zu ", temp_node->timestamp);  // its timestamp (FEval)
    #endif
    #ifndef IGNORE_DECISION_VARIABLES
    for (size_t j=0; j < nVars; j++)  // all decision variables of a solution
        #ifdef WRITE_WITH_MINIMAL_SPACING
        fprintf(logfile, "%.9e ", temp_node->var[j]);
        #else
        fprintf(logfile, "%16.9e ", temp_node->var[j]);
        #endif
    #endif
    for (size_t k=0; k < nObjs; k++)  // all objective values of a solution
        #ifdef WRITE_WITH_MINIMAL_SPACING
        fprintf(logfile, "%.10e ", temp_node->obj[k]);
        #else
        fprintf(logfile, "%17.10e ", temp_node->obj[k]);
        #endif
     temp_node = temp_node->next;
     fprintf(logfile, "\n");
  }

  // free STUFF
  fflush(logfile);
  free(logfile);
  free(buf);
  free(fullLogFilename);
}


/* from mobbob-gecco 2016*/
void creatFolder(char *folderName) {
    struct stat st = {0};
    if (stat(folderName, &st) == -1)
        mkdir(folderName, 0775);
}
