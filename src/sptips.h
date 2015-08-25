#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "adesub.h"


int intAinB(int a, int *b, int lengthB);
void intANotInB(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize);
void unionInt(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize);
void intersectInt(int *a, int *b, int lengthA, int lengthB, int *res, int *resSize);
void pathTipToRoot(int tip, int *ances, int *desc, int N, int *res, int *resSize);
int mrca2tips(int *ances, int*desc, int a, int b, int N);
void sp2tips(int *ances, int *desc, int N, int tipA, int tipB, int *res, int *resSize);
void spalltips(int *ances, int *desc, int *N, int *nTips, int *res, int *resId, int *resSize);
