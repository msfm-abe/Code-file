#ifndef _INCLUDE_dARRAY_
#define _INCLUDE_dARRAY_


/*���I�z����m�ۂ���֐�
�@�@1�����̔z����m�ۂ���Ƃ���
  �@�@�Eint�^->ikakuho1��(int�^��)�v�f����n���B�߂�l��int�^�ւ̃|�C���^
	�@�Edouble�^->dkakuho1��(int�^��)�v�f����n���B�߂�l��double�^�ւ̃|�C���^
	2�����̔z����m�ۂ���Ƃ���
	�@�Eint�^->ikakuho2��(int�^)�s��,(int�^)�񐔂̏��ŗv�f����n���B�߂�l��int�^�ւ̃|�C���^�ւ̃|�C���^
	�@�Edouble�^->dkakuho2��(int�^)�s��,(int�^)�񐔂̏��ŗv�f����n���B�߂�l��double�^�ւ̃|�C���^�ւ̃|�C���^
�@���I�z����������֐�
�@�@1�����̔z����������Ƃ���
  �@�@�Eint�^->ikaiho1�ɉ�����铮�I�z��̐擪�̃|�C���^��n���B
	�@�Edouble�^->dkaiho1�ɉ�����铮�I�z��̐擪�̃|�C���^��n���B
	2�����̔z����J������Ƃ���
	�@�Eint�^->ikaiho2�ɉ�����铮�I�z��̐擪�̃|�C���^�ւ̃|�C���^��(int�^)�s�������̏��œn���B
	�@�Edouble�^->dkaiho2�ɉ�����铮�I�z��̐擪�̃|�C���^�ւ̃|�C���^��(int�^)�s�������̏��œn���B
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