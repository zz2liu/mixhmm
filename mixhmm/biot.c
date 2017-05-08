/*Modified density function of conditional probability for mixture model*/
/*Version: 0.90, Start Date: Oct/6/2008, Author: Ao Li,                 */
/*Comments: Detailed info about this function can be found in the method*/
/*of the mixture model, "Conditional distribution of mixed LRR values"  */

#include <stdio.h>
#include <math.h>
#include "kc.h"
#include "util.h"
#include "biot.h"

double b1iot_new(int state, double mean[], double sd[],
    double uf, double r, double p)
//return the log likihood of observing mix LRR value r, given the tumor state.
// - state: 1..6
// - mean, sd:
// - uf:
// - p: proportion of Normal DNA
//todo: move the state calculation out for performance.
{
    double muN = mean[NORMAL]; //mean of normal state  
    double muT = mean[state]; //mean of tumor state
    double varN = sd[NORMAL]*sd[NORMAL]; //variance of normal state
    double varT = sd[state]*sd[state]; //variance of tumor state

    double muN_p = mean[NORMAL]+log(p); //+c; // eq.(8)
    double muT_p = mean[state]+log(1-p); //+c; //eq.(9)

    //following formulas are from eq.(6, 7)
    double part1 = exp(2*muN_p+varN)*(exp(varN)-1)+exp(2*muT_p+varT)*(exp(varT)-1);
    double part2 = (exp(muN_p+varN/2)+exp(muT_p+varT/2))*(exp(muN_p+varN/2)+exp(muT_p+varT/2));
    double varM = log(part1/part2+1);
    double muM = log(exp(muN_p+varN/2)+exp(muT_p+varT/2))-varM/2;//eq.(7)
    
    //to keep in b1iot(varM, muM)
    r *= log(2); //log2 -> ln
    //following formulas are from eq.(12)
    double cond_p = uf*EPS_R;
    //cond_p += (1-uf) * pdf_normal (r, muM-c, varM);
    cond_p += (1-uf) * pdf_normal (r, muM, varM);
    if (cond_p==0) cond_p=FLOAT_MINIMUM;
    printf("%f\n", cond_p);
    return log(cond_p);
}

/*void b2_mix_distr(double mean[], double sd[], double p,
        double ***meanM, double ***sdM, double ***gtype_state_p) //return
//return mix_mean, mix_sd for each state, normal genotype and tumor genotype.
// - mean, sd: array of {void, B0, B33, B50, B50D}
// - meanM, sdM: 3darray of (tumor state, normal genotype, tumor genotype)
{
    // modify the original BAF SD to make the model work??
    double sd_n[] = {0, sd[1]/sqrt(p), sd[2]/sqrt(p), sd[3]/sqrt(p),
            sd[4]/sqrt(p), sd[5]/sqrt(p)}; //in a hard way
    double sd_t[] = {0, sd[1]/sqrt(1-p), sd[2]/sqrt(1-p), sd[3]/sqrt(1-p),
            sd[4]/sqrt(1-p), sd[5]/sqrt(1-p)}; //in a hard way

    enum INDEX {B0=1, B25=2, B33=3, B50=4, B50D=5};
    double normal_mean[] = {mean[B0],mean[B50],1-mean[B0]};
    double normal_sd[] = {sd_n[B0],sd_n[B50],sd_n[B0]};

    //parameters of the normal distribution of mix BAF values
    int n_gtype[]= {1,2,3,2,4,5}; //number of genotypes for each state,Null genotype is 1
    int n_copy[] = {0,1,2,2,3,4}; //number of total copies
    double A = (p*(n_copy[NORMAL]))/(p*(n_copy[NORMAL])+(1-p)*(n_copy[state])); //eq.(15)
    double B = 1-A; //eq.(15)


    //initiate
    double tumor_mean[N_STATE+1][MAX_GTYPE];
    double tumor_sd[N_STATE+1][MAX_GTYPE];
    double gtype_state_p[N_STATE+1][MAX_GTYPE][n_gtype[NORMAL]];
            //for p(g'|g,z,zn) - prob of tumor genotype given tumer state and
            //normal genotype, zn unecessary??
    tumor_mean[1]    = {mean[B50D], NA, NA, NA, NA};
    tumor_sd[1]      = {sd_t[B50D], NA, NA, NA, NA};
    gtype_state_p[1] =  {{1,1,1},   NA, NA, NA, NA};

    tumor_mean[2]    = {mean[B0],1-mean[B0],     NA, NA, NA };
    tumor_sd[2]      = {sd_t[B0],sd_t[B0],       NA, NA, NA };
    gtype_state_p[2] = {{1,0.5,0},{0,0.5,1},     NA, NA, NA };

    tumor_mean   [3] = {mean[B0], mean[B50], 1-mean[B0], NA, NA};
    tumor_sd     [3] = {sd_t[B0], sd_t[B50], sd_t[B0],   NA, NA};
    gtype_state_p[3] = {{1,0,0},{0,1,0},{0,0,1},         NA, NA};

    tumor_mean   [4] = {mean[B0],1-mean[B0],  NA, NA, NA};
    tumor_sd     [4] = {sd_t[B0],sd_t[B0],    NA, NA, NA};
    gtype_state_p[4] = {{1,0.5,0},{0,0.5,1},  NA, NA, NA};

    tumor_mean   [5] = {mean[B0],mean[B33],1-mean[B33],1-mean[B0],NA};
    tumor_sd     [5] = {sd_t[B0],sd_t[B33],sd_t[B33],sd_t[B0],    NA};
    gtype_state_p[5] = {{1,0,0},{0,0.5,0},{0,0.5,0},{0,0,1},      NA};

    tumor_mean   [6] = {mean[B0],mean[B25],mean[B50],1-mean[B25],1-mean[B0]};
    tumor_sd     [6] = {sd_t[B0],sd_t[B25],sd_t[B50],sd_t[B25],sd_t[B0]};
    gtype_state_p[6] = {{1,0,0},{0,0.25,0},{0,0.5,0},{0,0.25,0},{0,0,1}};

    double muM, varM; //to be moved out
    int i,j;
    for (i=0; i<n_gtype[NORMAL]; i++) { //genotype of normal
        for (j=0; j<n_gtype[state]; j++) { //genotype of tumor
            muM =  A*normal_mean[i]+B*tumor_mean[j];
            varM = (A*normal_sd[i])*(A*normal_sd[i])+(B*tumor_sd[j])*(B*tumor_sd[j]);

        }
    }
}

double b2iot_new (int state, int n_gtype[],
        double **gtype_mean, double **gtype_sd, double **gtype_p,
        double pfb, double b)
//return the likelihood of observing mix BAF value of b, given the hidden state
//of tumor.
//- n_gtype: 1darray of number of gtypes of each state
//- gtype_p: 2darray of [MAX_GTYPE, NORMAL_GTYPE] of probs of observing tumor
//genotype giving normal genotype.
{
    double cond_p, last;
    if (b==0) {
            cond_p = UF*EPS_B*MAGIC;//eq (18)
            last_term = cdf_normal(0,muM, sqrt(varM))
    } else if (b==1) {
            cond_p = UF*EPS_B*MAGIC;//eq (19)
            last_term = 1-cdf_normal(1,muM, sqrt(varM))
    } else {
            cond_p = UF*EPS_B;//eq (16)
            last_term = pdf_normal(b,muM, sqrt(varM));
    }

    double cond_p;
    //for p(g|zn) - prob of a normal genotype given its state.
    double gtype_normal_p[] = {(1-pfb)*(1-pfb),2*pfb*(1-pfb),pfb*pfb};//AA,AB,BB
    int i,j;
    for (i=0; i<n_gtype[NORMAL]; i++) {//normal genotype
        for (j=0; j<n_gtype[state]; j++) {//tumor genotype
            cond_p += (1-UF)*gtype_normal_p[i]*gtype_state_p[j][i] * last_term
        }
    }
    if (p==0) p=FLOAT_MINIMUM;
    return log(p);
}
*/
