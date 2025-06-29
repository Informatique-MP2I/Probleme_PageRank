#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* constants */
#define N 6

/* structures */
typedef struct {
} webpage_s;


void print_vect(double vect[N]){
  /* body of the function */
  return;
}

void test_print_vect(){
  /* write your tests here */
  return;
}

void all_ones_vect(double vect[N]){
  /* body of the function */
  return;
}

void test_all_ones_vect(){
  /* write your tests here */
  return;
}

void print_mat(double mat[N][N]){
  /* body of the function */
  return;
}

void test_print_mat(){
  /* write your tests here */
  return;
}

void all_ones_mat(double mat[N][N]){
  /* body of the function */
  return;
}

void test_all_ones_mat(){
  /* write your tests here */
  return;
}

void add_mat(double mat1[N][N], double mat2[N][N], double res[N][N]){
  /* body of the function */
  return;
}

void test_add_mat(){
  /* write your tests here */
  return;
}

void scalar_mult_vect(double s, double vect[N], double res[N]){
  /* body of the function */
  return;
}

void test_scalar_mult_vect(){
  /* write your tests here */
  return;
}

void scalar_mult_mat(double s, double mat[N][N], double res[N][N]){
  /* body of the function */
  return;
}

void test_scalar_mult_mat(){
  /* write your tests here */
  return;
}

void mat_vect_lmult(double mat[N][N], double vect[N], double res[N]){
  /* body of the function */
  return;
}

void test_mat_vect_lmult(){
  /* write your tests here */
  return;
}

webpage_s *init_webpage(int id, int nb_links, int *links){
  /* body of the function */
  return NULL; /* to change */
}

void init_trans_mat(double P[N][N], webpage_s *web_graph[N]){
  /* body of the function */
  return;
}

void adj_trans_mat(double alpha, double P[N][N], double M[N][N]){
  /* body of the function */
  return;
}

void init_rank_vect(double R[N]){
  /* body of the function */
  return;
}

void update_pagerank_it(double M[N][N], double R[N], int nb_iter, double Rf[N]){
  /* body of the function */
  return;
}

int update_pagerank_asy(double M[N][N], double R[N], double eps, double Rf[N]){
  /* body of the function */
  return -1; /* to change */
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

  /* Declaration of the web graph */
  webpage_s *web_graph[N] = {}; /* to change */
  
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
