/* ===================================
** test.c
** started on Wed Jan 24 14:40:45 2007 
** ===================================
*/

 
typedef union ptrvar {
  int i;
  double d;
  char c;
} uptrvar;

int main() {
  uptrvar ptrvar[3];
  
  ptrvar[0].i = 1;
  ptrvar[1].d = 2;
  // ptrvar[1].c = 'a';
  ptrvar[2].d = 3;
  
  printf("\nptrvar[0].i=%d",ptrvar[0].i);
   printf("\nptrvar[0].d=%f",ptrvar[0].d);
  // printf("\nptrvar[1].c=%s",ptrvar[1].c);
  printf("\nptrvar[2].d=%f",ptrvar[2].d);
}
