
#ifndef RandomFieldsUtilsxport_H
#define RandomFieldsUtilsxport_H 1

extern int CORES;
extern int PL;

bool parallel();


extern utilsoption_type OPTIONS;

#define prefixN 2
extern const char * prefixlist[prefixN], **all[prefixN];
extern int allN[prefixN];

typedef
struct KEY_type KEY_type;
struct KEY_type {
  KEY_type *next;
  utilsoption_type global_utils;
  int pid,  visitingpid;
  bool ok, doshow;
  errorloc_type error_location;

  int *ToIntDummy;
  int ToIntN, ToRealN ;
  double *ToRealDummy;

  double loggamma1old, nu1old,
    loggamma2old, nu2old,
    loggamma_old,nuOld,
    gamma, nuAlt;
};
extern KEY_type *PIDKEY[PIDMODULUS];
KEY_type *KEYT();


#endif
