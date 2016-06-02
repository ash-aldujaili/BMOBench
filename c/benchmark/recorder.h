

/*
 Recorder
*/
#ifndef RECORD_H
#define RECORD_H

typedef struct node_t node_t;
struct node_t {
  double *obj;
  double *var;
  size_t timestamp;
  node_t *next;
};

typedef struct archive {
    node_t *head;
    node_t *tail;
} archive_t ;

//static archive_t *archive;
/* free allocation */
void freeArchive();
/* add a new sample if it is non-dominated*/
void updateArchive(double *var, double *obj,  size_t timestamp);
/* print to a txt file */
void printArchive();
/* create a folder method */
void creatFolder(char *folderName);
#endif
