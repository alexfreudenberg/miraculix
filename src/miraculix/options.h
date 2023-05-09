
/*
 Authors 
 Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2022-2023 Martin Schlather

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/


#ifndef miraculix_options_H
#define miraculix_options_H 1


#define geneticsN 8
struct genetics_options {
  usr_bool haplo_set1;
  bool squared, interprete_coding_as_is, prefer_external_freq;
  normalizing_type normalized;
  centering_type centered;
  // internal:
  coding_type coding;
  double digits;
  double *pcentered;
  int ncentered;
};
extern const char *genetics_names[geneticsN];

#define genetics_START {					\
    Nan,							\
      false, false, false,						\
      AccordingCentering, RowMeans,					\
      UnknownSNPcoding,						\
      3.0,							\
      NULL, 0 /* internal */					\
      }


#define tuningN 14
struct tuning_options {
  bool savegeno, oldversion, smoothSnpGroupSize, missingsFully0,
    meanVsubstract, meanSxIsubstract, gpu, addtransposed;
  int variant, logSnpGroupSize, indivGroupSize, miniRows, miniCols,
    floatLoop // 0 : no; > 0 : number of iterations befor summing up
    ;
};
extern const char *tuning_names[tuningN];


#define tuning_START {				\
    false, false, false, false,				\
      true, true, true, false, 			\
      VARIANT_DEFAULT, 0, 0, 0, 0, 0			\
      }


#define messagesN 1
extern const char * messages[messagesN];
struct messages_options{ 
  bool warn_address;
};
#define messages_START {			\
    true					\
      }

struct option_type {
  genetics_options genetics;
  tuning_options tuning;
  messages_options messages;
};
extern option_type OPTIONS_MIRACULIX;

#define prefixMN 3
extern const char * prefixMlist[prefixMN], **allMoptions[prefixMN];
extern int allMoptionsN[prefixMN];

int main_variant(int variant);

typedef
struct KEY_type KEY_type;
struct KEY_type {
  KEY_type *next;
  option_type global;
  utilsoption_type global_utils;
  int pid,  visitingpid;
  bool ok, doshow;
  errorstring_type error_location;

  int *ToIntDummy;
  int  ToIntN = 0;
};
extern KEY_type *PIDKEY_M[PIDMODULUS];
KEY_type *KEYT_M();

bool exists_tiling(coding_type coding, int variant, bool tiling);
bool exists_variant(coding_type coding, int variant, bool indeed,
		    bool efficient);
bool exists_crossprod(coding_type coding);
bool exists_allelefreq(coding_type coding);


#define NvariantsSHR5 32
const bool variantsSHR5[NvariantsSHR5] = {
  /*<=64*/true, true, true, false, /*128*/true,  true,  true, false,//last:224
  /*256*/ true, true, true, false,       false, false, false, false,//last:480
  /*512*/ true, true, true, false,       false, false, false, false,//last:736
  /*GPU*/ true, false,false,false,  /*R*/true, true, false, false //last:992
};

typedef bool (*exists_coding_t) (coding_type );
bool exists_internal_coding(coding_type m) ;
bool exists_user_coding(coding_type m) ;
bool exists_coding(coding_type m);
void getMoptions(SEXP sublist, int i, bool local);
void setMoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		 bool isList, bool local);

#endif
  
