#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* constants */
#define MAX_LINKS 100
#define N 6
/* const int MAX_LINKS = 100; */
/* const int N = 6; */

/* structures */
typedef struct {
  int id;
  int nb_links;
  int links[MAX_LINKS];
} webpage_s;

void print_vect(double vect[N]){
  for(int i=0;i<N;i++)
    printf("%f\n", vect[i]);
  printf("\n");
  return;
}

void test_print_vect(){
  double vect[N];
  for(int i=0;i<N;i++)
    vect[i]=i;
  print_vect(vect);
  return;
}

void all_ones_vect(double vect[N]){
  for(int i=0;i<N;i++)
      vect[i] = 1;
  return;
}

void test_all_ones_vect(){
  double vect[N];
  all_ones_vect(vect);
  print_vect(vect);
  return;
}

void print_mat(double mat[N][N]){
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++)
      printf("%.3f ", mat[i][j]);
    printf("\n");
  }
  printf("\n");
  return;
}

void test_print_mat(){
  double mat[N][N];
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      mat[i][j]=i;
  print_mat(mat);
  return;
}

void all_ones_mat(double mat[N][N]){
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      mat[i][j] = 1;
  return;
}

void test_all_ones_mat(){
  double mat[N][N];
  all_ones_mat(mat);
  print_mat(mat);
  return;
}

void add_mat(double mat1[N][N], double mat2[N][N], double res[N][N]){
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      res[i][j] = mat1[i][j] + mat2[i][j];
  return;
}

void test_add_mat(){
  double mat1[N][N];
  double mat2[N][N];
  double res[N][N];
  
  all_ones_mat(mat1);
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      mat2[i][j] = i;
  add_mat(mat1,mat2,res);
  print_mat(res);
  return;
}

void scalar_mult_vect(double s, double vect[N], double res[N]){
  for(int i=0;i<N;i++)
      res[i] = s*vect[i];
  return;
}

void test_scalar_mult_vect(){
  double s=4;
  double vect[N];
  double res[N];
  all_ones_vect(vect);
  scalar_mult_vect(s,vect,res);
  print_vect(res);
  return;
}

void scalar_mult_mat(double s, double mat[N][N], double res[N][N]){
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      res[i][j] = s*mat[i][j];
  return;
}

void test_scalar_mult_mat(){
  double s=4;
  double mat[N][N];
  double res[N][N];
  all_ones_mat(mat);
  scalar_mult_mat(s,mat,res);
  print_mat(res);
  return;
}

void mat_vect_lmult(double mat[N][N], double vect[N], double res[N]){
  for(int i=0;i<N;i++){
    double sum = 0;
    for(int j=0;j<N;j++)
      sum += mat[i][j]*vect[j];
    res[i] = sum;
  }
  return;
}

void test_mat_vect_lmult(){
  double mat[N][N];
  double vect[N];
  double res[N];
  all_ones_mat(mat);
  all_ones_vect(vect);
  mat_vect_lmult(mat,vect,res);
  print_vect(res);
  return;
}

webpage_s *init_webpage(int id, int nb_links, int *links){
  // (1) verification of the number of links
  assert(nb_links <= MAX_LINKS && nb_links >=0);
  
  // (2) malloc of the structure
  webpage_s *webpage = (webpage_s *)malloc(sizeof(webpage_s));
  assert(webpage!=NULL);

  // (3) initialization of id
  webpage->id = id;

  // (4) initialization of nb_links
  webpage->nb_links = nb_links;
  
  // (5) initialization of the link array
  for (int i=0;i<nb_links;i++){
    webpage->links[i] = links[i];
  }

  return webpage;
}

void init_trans_mat(double P[N][N], webpage_s *web_graph[N]){
  all_ones_mat(P);
  scalar_mult_mat(0,P,P);
  
  for(int j=0;j<N;j++){
    if (web_graph[j]->nb_links == 0) {
      for(int i=0;i<N;i++)
	P[i][j] = 1.0/(double)N;
    } else {
      for(int k=0;k<web_graph[j]->nb_links;k++){
	P[web_graph[j]->links[k]][j] = 1.0/(double)web_graph[j]->nb_links;
      }
    }
  }
  return;
}

void adj_trans_mat(double alpha, double P[N][N], double M[N][N]){
  double E[N][N];
  double tmp1[N][N];
  double tmp2[N][N];

  all_ones_mat(E);
  scalar_mult_mat((1.0-alpha)/(double)N, E, tmp1);
  scalar_mult_mat(alpha, P, tmp2);
  add_mat(tmp1, tmp2, M);
  return;
}

void init_rank_vect(double R[N]){
  all_ones_vect(R);
  scalar_mult_vect(1.0/(double)N, R, R);
  return;
}

void update_pagerank_it(double M[N][N], double R[N], int nb_iter, double Rf[N]){
  for(int i=0;i<nb_iter;i++){
    mat_vect_lmult(M, R, Rf);
    for(int j=0;j<N;j++)
      R[j] = Rf[j];
  }
  return;
}

double norm_diff(double U[N], double V[N]){
  double n = 0;
  for(int i=0;i<N;i++){
    n += (U[i] - V[i])*(U[i] - V[i]);
  }
  return sqrt(n);
}

int update_pagerank_asy(double M[N][N], double R[N], double eps, double Rf[N]){
  mat_vect_lmult(M, R, Rf);
  int count = 0;
  while(norm_diff(R, Rf)>eps){
    for(int i=0;i<N;i++)
      R[i] = Rf[i];
    mat_vect_lmult(M, R, Rf);
    count += 1;
  }
  return count;
}

int main(int argc, char **argv){
  printf("PageRank Algorithm.\n");

  /* Auxiliary test functions */
  test_print_vect();
  test_all_ones_vect();
  test_print_mat();  
  test_all_ones_mat();
  test_add_mat();
  test_scalar_mult_vect();
  test_scalar_mult_mat();
  test_mat_vect_lmult();

  /* Web graph example */
  /* Declaration of webpages */
  int wp0_links[3] = {1,2,3};
  webpage_s *wp0 = init_webpage(0,3,wp0_links);
  int wp1_links[3] = {0,2,4};
  webpage_s *wp1 = init_webpage(1,3,wp1_links);
  int wp2_links[1] = {3};
  webpage_s *wp2 = init_webpage(2,1,wp2_links);
  int wp3_links[2] = {0,2};
  webpage_s *wp3 = init_webpage(3,2,wp3_links);
  int wp4_links[3] = {2,3,5};
  webpage_s *wp4 = init_webpage(4,3,wp4_links);
  int wp5_links[0] = {};
  webpage_s *wp5 = init_webpage(5,0,wp5_links);

  /* Declaration of the web graph */
  webpage_s *web_graph[N] = {wp0, wp1, wp2, wp3, wp4, wp5};
  
  /* PageRank algorithm */
  /* Variables */
  double P[N][N];
  double M[N][N];
  double alpha = 0.85;
  double R_it[N];
  double pagerank_it[N];
  double R_asy[N];
  double pagerank_asy[N];

  /* Transition matrix */
  printf("Computing transition matrix\n");
  init_trans_mat(P, web_graph);
  print_mat(P);
  adj_trans_mat(alpha, P, M);
  print_mat(M);

  /* Rank vector */
  int nb_iter = 20;

  printf("Computing Rank vector (1)\n");
  init_rank_vect(R_it);
  update_pagerank_it(M, R_it, nb_iter, pagerank_it);
  printf("After %d iterations :\n",nb_iter);
  print_vect(pagerank_it);

  printf("Computing Rank vector (2)\n");
  init_rank_vect(R_asy);
  nb_iter = update_pagerank_asy(M, R_asy, 1e-12, pagerank_asy);
  printf("After %d iterations :\n",nb_iter);
  print_vect(pagerank_asy);
  
  return 0;
}
