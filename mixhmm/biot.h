#define EPS_R 0.1
#define EPS_B 0.1 //eq. 16
#define NORMAL 3 //state of Normal DNA
#define MAGIC 4.5 //used in eq. 18,19
#define NA 999 //used for n/a mean, sd, ..
#define MAX_GTYPE 5 //the max number of genotypes of a state
#define NORMAL_GTYPE 3 //number of normal genotypes
#define N_STATE 6 //number of hidden states
#define UF 0.01 //uf

double b1iot_new(int state, double mean[], double sd[], double uf, double r, double p);
