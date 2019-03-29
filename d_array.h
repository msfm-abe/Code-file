#ifndef _INCLUDE_dARRAY_
#define _INCLUDE_dARRAY_


/*動的配列を確保する関数
　　1次元の配列を確保するときは
  　　・int型->ikakuho1に(int型で)要素数を渡す。戻り値はint型へのポインタ
	　・double型->dkakuho1に(int型で)要素数を渡す。戻り値はdouble型へのポインタ
	2次元の配列を確保するときは
	　・int型->ikakuho2に(int型)行数,(int型)列数の順で要素数を渡す。戻り値はint型へのポインタへのポインタ
	　・double型->dkakuho2に(int型)行数,(int型)列数の順で要素数を渡す。戻り値はdouble型へのポインタへのポインタ
　動的配列を解放する関数
　　1次元の配列を解放するときは
  　　・int型->ikaiho1に解放する動的配列の先頭のポインタを渡す。
	　・double型->dkaiho1に解放する動的配列の先頭のポインタを渡す。
	2次元の配列を開放するときは
	　・int型->ikaiho2に解放する動的配列の先頭のポインタへのポインタと(int型)行数をこの順で渡す。
	　・double型->dkaiho2に解放する動的配列の先頭のポインタへのポインタと(int型)行数をこの順で渡す。
*/


extern int *ikakuho1(int size);

extern double *dkakuho1(int size);

extern int **ikakuho2(int row,int column);

extern double **dkakuho2(int row,int column);

extern void ikaiho1(int *p);

extern void dkaiho1(double *p);

extern void ikaiho2(int **p,int row);

extern void dkaiho2(double **p,int row);


#endif