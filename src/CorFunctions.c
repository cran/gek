#include <R.h>
#include <Rmath.h>
#include <string.h>
#define square(x) ((x)*(x))


/* correlation functions */

// Matern-5/2
double Matern5_2(const double *x1, const int *n1, const double *x2, const int *n2, const int *d, const int *i1, const int *i2, const double *param) {
	double res = 0.;
	double tmp = 0.;
	for (int k = 0; k < *d; k++) {
		tmp = fabs(x1[*i1 + *n1 * k] - x2[*i2 + *n2 * k]) / param[k] * sqrt(5.);
		res += tmp - log(1 + tmp + square(tmp) / 3);
	}
	return(exp(-res));
	
}


// Matern-3/2
double Matern3_2(const double *x1, const int *n1, const double *x2, const int *n2, const int *d, const int *i1, const int *i2, const double *param) {
	double res = 0.;
	double tmp = 0.;
	for (int k = 0; k < *d; k++) {
		tmp = fabs(x1[*i1 + *n1 * k] - x2[*i2 + *n2 * k]) / param[k] * sqrt(3.);
		res += tmp - log(1 + tmp);
	}
	return(exp(-res));
}

// Gaussian
double Gaussian(const double *x1, const int *n1, const double *x2, const int *n2, const int *d, const int *i1, const int *i2, const double *param) {
	double res = 0.;
	for (int k = 0; k < *d; k++) {
		res += 0.5 * square((x1[*i1 + *n1 * k] - x2[*i2 + *n2 * k]) / param[k]);
	}
	return(exp(-res));
}


/* correlation matrices  */

void corMat(const double *x, const int *n, const int *d, const double *param, const char **type, double *ans) {
	
	double (*corFun)(const double *, const int *, const double *, const int *, const int *, const int *, const int *, const double *);

	if (strcmp(*type, "gaussian") == 0) corFun = Gaussian;
	else if (strcmp(*type, "matern3_2") == 0) corFun = Matern3_2;
	else if (strcmp(*type, "matern5_2") == 0) corFun = Matern5_2;
	
	for (int i = 0; i < *n; i++) {
		for (int j = 0; j < i; j++) {	
			ans[j + *n * i] = ans[i + *n * j] = (*corFun)(x, n, x, n, d, &i, &j, param);
		}
	ans[i + *n *i] = 1;
	}
} 

void corMat2(const double *x1, const int *n1, const double *x2, const int *n2, const int *d, const double *param, const char **type, double *ans) {
	
	double (*corFun)(const double *, const int *, const double *, const int *, const int *, const int *, const int *, const double *);
	
	if (strcmp(*type, "gaussian") == 0) corFun = Gaussian;
	else if (strcmp(*type, "matern3_2") == 0) corFun = Matern3_2;
	else if (strcmp(*type, "matern5_2") == 0) corFun = Matern5_2;
	
	for (int i1 = 0; i1 < *n1; i1++) {
		for (int i2 = 0; i2 < *n2; i2++) {				
			ans[i1 + *n1 * i2] = (*corFun)(x1, n1, x2, n2, d, &i1, &i2, param);
		}
	}
	
} 



/* first partial derivatives of correlation function with respect to x1 */

// Matern-5/2
double Matern5_2_dx(const double *x1, const int *n1, const double *x2, const int *n2, const int *i1, const int *i2, const double *param, const int *k, const double *C) {
	double res = 0.;
	double diff = x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k];
	double absdiff = sqrt(5.) * fabs(diff);
	res = -5 * diff * (absdiff + param[*k]) / (param[*k] * (3 * square(param[*k]) + 3 * param[*k] * absdiff + square(absdiff)));
	return(res * C[*i1 + *n1 * *i2]);
}

// Matern-3/2
double Matern3_2_dx(const double *x1, const int *n1, const double *x2, const int *n2, const int *i1, const int *i2, const double *param, const int *k, const double *C) {
	double res = 0.;
	double diff = x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k]; 
	res = -3 * diff / (param[*k] * (param[*k] + sqrt(3.) * fabs(diff)));
	return(res * C[*i1 + *n1 * *i2]);
}

// Gaussian
double Gaussian_dx(const double *x1, const int *n1, const double *x2, const int *n2, const int *i1, const int *i2, const double *param, const int *k, const double *C) {
	double res = 0.;
	res = -(x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k]) / square(param[*k]);
	return(res * C[*i1 + *n1 * *i2]);
}


/* correlation matrix with first derivatives with respect to x1 */

void corMat_dx(const double *x, const int *n, const double *param, int *k, const char **type, double *C, double *ans) {
	
	double (*corFun_dx)(const double *, const int *, const double *, const int *, const int *, const int *, const double *, const int *, const double *);

	if (strcmp(*type, "gaussian") == 0) corFun_dx = Gaussian_dx;
	else if (strcmp(*type, "matern3_2") == 0) corFun_dx = Matern3_2_dx;
	else if (strcmp(*type, "matern5_2") == 0) corFun_dx = Matern5_2_dx;

	for (int i = 0; i < *n; i++) {
		for (int j = 0; j < i; j++) {
			ans[i + *n * j] = (*corFun_dx)(x, n, x, n, &i, &j, param, k, C);
			ans[j + *n * i] = -ans[i + *n * j];
		}
		ans[i + *n * i] = 0;
	}
}

void corMat_dx_all(const double *X, const int *n, const int *d, const double *param, const char **type, double *C, double *tmp, double *ans) {
	for(int k = 0; k < *d; k++){
		corMat_dx(X, n, param, &k, type, C, tmp);
		for(int i = 0; i < *n * *n; i++){
			ans[k + *d * i] = tmp[i];
		}
	}
}


void corMat2_dx(const double *x1, const int *n1, const double *x2, const int *n2, const double *param, int *k, const char **type, double *C, double *ans) {

	double (*corFun_dx)(const double *, const int *, const double *, const int *, const int *, const int *, const double *, const int *, const double *);

	if (strcmp(*type, "gaussian") == 0) corFun_dx = Gaussian_dx;
	else if (strcmp(*type, "matern3_2") == 0) corFun_dx = Matern3_2_dx;
	else if (strcmp(*type, "matern5_2") == 0) corFun_dx = Matern5_2_dx;

	for(int i1 = 0; i1 < *n1; i1++) {
		for(int i2 = 0; i2 < *n2; i2++) {				
			ans[i1 + *n1 * i2] = (*corFun_dx)(x1, n1, x2, n2, &i1, &i2, param, k, C);
		}
	}
}

void corMat2_dx_all(const double *x1, const int *n1, const double *x2, const int *n2, const int *d, const double *param, const char **type, double *C, double *tmp, double *ans) {
	for(int k = 0; k < *d; k++){
		corMat2_dx(x1, n1, x2, n2, param, &k, type, C, tmp);
		for(int i = 0; i < *n1 * *n2; i++){
			ans[k + *d * i] = tmp[i];
		}
	}
}



/* second partial derivatives of correlation function with respect to x1 and x2 */

// Matern-5/2
double Matern5_2_dxdy(const double *x1, const int *n1, const double *x2, const int *n2, const int *i1, const int *i2, const double *param, const int *k, const int *l, const double *C) {
	double res = 0.;
	if(*k == *l){
		double absdiff = sqrt(5.) * fabs(x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k]);
		double num = 5 * (square(param[*k]) + param[*k] * absdiff - square(absdiff));
		double den = square(param[*k]) * (3 * square(param[*k]) + 3 * param[*k] * absdiff + square(absdiff));
		res = num / den;
	}else{
		double diffk = sqrt(5.) * (x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k]);
		double diffl = sqrt(5.) * (x1[*i1 + *n1 * *l] - x2[*i2 + *n2 * *l]);
		double num = 5 * diffk * diffl * (fabs(diffk) + param[*k]) * (fabs(diffl) + param[*l]);
		double den = param[*k] * param[*l] * (3 * square(param[*k]) + 3 * param[*k] * fabs(diffk) + square(diffk)) * (3 * square(param[*l]) + 3 * param[*l] * fabs(diffl) + square(diffl));
		res = -num / den;
	}
	return(res * C[*i1 + *n1 * *i2]);
}


// Matern-3/2
double Matern3_2_dxdy(const double *x1, const int *n1, const double *x2, const int *n2, const int *i1, const int *i2, const double *param, const int *k, const int *l, const double *C) {
	double res = 0.;
	if(*k == *l){
		double diff = sqrt(3.) * fabs(x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k]);
		res = 3 / (square(param[*k]) * (param[*k] + diff)) * (param[*k] - diff);
	}else{
		double diffk = x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k];
		double diffl = x1[*i1 + *n1 * *l] - x2[*i2 + *n2 * *l];
		res = -9 * diffk * diffl / (param[*k] * param[*l] * (param[*k] + sqrt(3.) * fabs(diffk)) * (param[*l] + sqrt(3.) * fabs(diffl)));
	}
	return(res * C[*i1 + *n1 * *i2]);
}

// Gaussian
double Gaussian_dxdy(const double *x1, const int *n1, const double *x2, const int *n2, const int *i1, const int *i2, const double *param, const int *k, const int *l, const double *C) {
	double res = 0.;
	if(*k == *l){
		res = (1 - square((x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k]) / param[*k])) / square(param[*k]);
	}else{
		res = -(x1[*i1 + *n1 * *k] - x2[*i2 + *n2 * *k]) * (x1[*i1 + *n1 * *l] - x2[*i2 + *n2 * *l]) / square(param[*k]) / square(param[*l]);
	}
	return(res * C[*i1 + *n1 * *i2]);
}


/* correlation matrix with second derivatives with respect to x1 and x2 */

void corMat_dxdy_all(const double *x, const int *n, const int *d, const double *param, const char **type, double *C, double *ans) {
	
	double (*corFun_dxdy)(const double *, const int *, const double *, const int *, const int *, const int *, const double *, const int *, const int *, const double *);

	if (strcmp(*type, "gaussian") == 0) corFun_dxdy = Gaussian_dxdy;
	else if (strcmp(*type, "matern3_2") == 0) corFun_dxdy = Matern3_2_dxdy;
	else if (strcmp(*type, "matern5_2") == 0) corFun_dxdy = Matern5_2_dxdy;

	for(int i = 0; i < *n; i++){
		for(int k = 0; k < *d; k++){
			for(int j = 0; j < i; j++){
				for(int l = 0; l < k; l++){	
					ans[(k + *d * i) * *n * *d + j * *d + l] = ans[(l + *d * j) * *n * *d + i * *d + k] = ans[(l + *d * i) * *n * *d + j * *d + k] = ans[(k + *d * j) * *n * *d + i * *d + l] = 
					(*corFun_dxdy)(x, n, x, n, &i, &j, param, &k, &l, C);
				}
				ans[(k + *d * i) * *n * *d + j * *d + k] = ans[(k + *d * j) * *n * *d + i * *d + k] = (*corFun_dxdy)(x, n, x, n, &i, &j, param, &k, &k, C);
			}
			ans[(k + *d * i) * *n * *d + i * *d + k] = (*corFun_dxdy)(x, n, x, n, &i, &i, param, &k, &k, C);
		}
	}
}


void corMat2_dxdy_all(const double *x1, const int *n1, const double *x2, const int *n2, const int *d, const double *param, const char **type, double *C, double *ans) {
	
	double (*corFun_dxdy)(const double *, const int *, const double *, const int *, const int *, const int *, const double *, const int *, const int *, const double *);
	
	if (strcmp(*type, "gaussian") == 0) corFun_dxdy = Gaussian_dxdy;
	else if (strcmp(*type, "matern3_2") == 0) corFun_dxdy = Matern3_2_dxdy;
	else if (strcmp(*type, "matern5_2") == 0) corFun_dxdy = Matern5_2_dxdy;
	
	for(int i1 = 0; i1 < *n1; i1++){
		for(int k = 0; k < *d; k++){
			for(int i2 = 0; i2 < *n2; i2++){
				for(int l = 0; l < *d; l++){	
					ans[(k + *d * i1) * *n1 * *d + i2 * *d + l] = ans[(l + *d * i2) * *n1 * *d + i1 * *d + k] = ans[(l + *d * i1) * *n1 * *d + i2 * *d + k] = ans[(k + *d * i2) * *n1 * *d + i1 * *d + l] = 
					(*corFun_dxdy)(x1, n1, x2, n2, &i1, &i2, param, &k, &l, C);
				}
				ans[(k + *d * i1) * *n1 * *d + i2 * *d + k] = ans[(k + *d * i2) * *n1 * *d + i1 * *d + k] = (*corFun_dxdy)(x1, n1, x2, n2, &i1, &i2, param, &k, &k, C);
			}
			ans[(k + *d * i1) * *n1 * *d + i1 * *d + k] = (*corFun_dxdy)(x1, n1, x2, n2, &i1, &i1, param, &k, &k, C);
		}
	}
	
} 


/* frist partial derivative of correlation functions with respect to theta */

// Matern-5/2
double Matern5_2_dp(const double *X, const int *n, const int *i, const int *j, const double *param, const int *k, const double *C) {
	double res = 0.;
	double diff = sqrt(5.) * fabs(X[*i + *n * *k] - X[*j + *n * *k]);
	res = square(diff) * (diff + param[*k]) / (square(param[*k]) * (3 * square(param[*k]) + 3 * param[*k] * diff + square(diff)));
	return(res * C[*i + *n * *j]);
}

//  Matern-3/2
double Matern3_2_dp(const double *X, const int *n, const int *i, const int *j, const double *param, const int *k, const double *C) {
	double res = 0.;
	double diff = sqrt(3.) * fabs(X[*i + *n * *k] - X[*j + *n * *k]);
	res = square(diff) / (square(param[*k]) * (diff + param[*k]));
	return(res * C[*i + *n * *j]);
}

// Gaussian
double Gaussian_dp(const double *X, const int *n, const int *i, const int *j, const double *param, const int *k, const double *C) {
	double res = 0.;
	res = square((X[*i + *n * *k] - X[*j + *n * *k]) / param[*k]) / param[*k];
	return(res * C[*i + *n * *j]);
}


/* correlation matrix with first derivatives with respect to theta */


void corMat_dp(const double *X, const int *n, const double *param, const char **type, int *k, double *C, double *ans) {
	(*k)--;  
	double (*corFun_dp)(const double *, const int *, const int *, const int *, const double *, const int *, const double *);

	if (strcmp(*type, "gaussian") == 0)  corFun_dp = Gaussian_dp;
	else if (strcmp(*type, "matern3_2") == 0) corFun_dp = Matern3_2_dp;
	else if (strcmp(*type, "matern5_2") == 0) corFun_dp = Matern5_2_dp;
	
	for (int i = 0; i < *n; i++) {
		for (int j = 0; j < i; j++) {
			ans[j + *n * i] = ans[i + *n * j] = (*corFun_dp)(X, n, &i, &j, param, k, C);
		}
		ans[i + *n * i] = 0;
	}
}


/* second partial derivatives of correlation function with respect to x1 and theta */


//  Matern-5/2
double Matern5_2_dxdp(const double *X, const int *n, const int *d, const int *i, const int *j, const double *param, const int *k, const int *m, const double *D) {
	double res = 0.;
	double diff = sqrt(5) * fabs(X[*i + *n * *m] - X[*j + *n * *m]);
	if(*k == *m){
		res = (2 * square(param[*m]) + 2 * param[*m] * diff - square(diff)) / (square(param[*m]) * (param[*m] + diff));
	}else{
		res = -square(diff) * (diff + param[*m]) / (square(param[*m]) * (3 * square(param[*m]) + 3 * param[*m] * diff + square(diff))); 
	}
	return(res * D[*k + *d * (*i + *n * *j)]);
}

// Matern-3/2
double Matern3_2_dxdp(const double *X, const int *n, const int *d, const int *i, const int *j, const double *param, const int *k, const int *m, const double *D) {
	double res = 0.;
	double diff = sqrt(3.) * fabs(X[*i + *n * *m] - X[*j + *n * *m]);
	if(*k == *m){
		res = -(square(diff) - diff * param[*m] - 2 * square(param[*m])) / (square(param[*m]) * (diff + param[*m]));
	}else{
		res = -square(diff) / (square(param[*m]) * (diff + param[*m]));
	}
	return(res * D[*k + *d * (*i + *n * *j)]);
}

// Gaussian
double Gaussian_dxdp(const double *X, const int *n, const int *d, const int *i, const int *j, const double *param, const int *k, const int *m, const double *D) {
	double res = 0.;
	if(*k == *m){
		res = (2 - square((X[*i + *n * *m] - X[*j + *n * *m]) / param[*m])) / param[*m];
	}else{
		res = -square((X[*i + *n * *m] - X[*j + *n * *m]) / param[*m]) / param[*m];
	}
	return(res * D[*k + *d * (*i + *n * *j)]);
}


/* correlation matrix with second derivatives with respect to x1 and theta */


void corMat_dxdp(const double *X, const int *n, const int *d, const double *param, const char **type, int *m, const double *D, double *ans) {
	(*m)--;  
	
	double (*corFun_dxdp)(const double *, const int *, const int *, const int *, const int *, const double *, const int *, const int *, const double *);

	if (strcmp(*type, "gaussian") == 0)  corFun_dxdp = Gaussian_dxdp;
	else if (strcmp(*type, "matern3_2") == 0) corFun_dxdp = Matern3_2_dxdp;
	else if (strcmp(*type, "matern5_2") == 0) corFun_dxdp = Matern5_2_dxdp;
	
	for (int k = 0; k < *d; k++) {
		for(int i = 0; i < *n; i++){
			for(int j = 0; j < i; j++){
				ans[k + *d * (j + *n * i)] = (*corFun_dxdp)(X, n, d, &i, &j, param, &k, m, D);
				ans[k + *d * (i + *n * j)] = -ans[k + *d * (j + *n * i)];
			}
			ans[k + *d * (i + *n * i)] = 0;
		}
	}
}


/* third partial derivatives of correlation function with respect to x1, x2 and theta */


// Matern-5/2
double Matern5_2_dxdydp(const double *X, const int *n, const int *d, const int *i, const int *j, const double *param, const int *k, const int *l, const int *m, const double *D, const double *E) {
	double res = 0.;
	double diff = sqrt(5.) * fabs(X[*i + *n * *m] - X[*j + *n * *m]);
	if(*k == *m && *m == *l){
		res = (-2 * square(param[*m]) - 2 * param[*m] * diff + square(diff)) / (square(param[*m]) * (param[*m] + diff)) * E[(*m + *d * *i) + *n * *d * (*m + *d * *j)] +
			5 * (X[*i + *n * *m] - X[*j + *n * *m]) * (2 * param[*m] + diff) / square(param[*m] * (param[*m] + diff)) * D[*m + *d * (*j + *n * *i)];		
	}
	if(*k != *m && *m == *l){
		res = (-2 * square(param[*m]) - 2 * param[*m] * diff + square(diff)) / (square(param[*m]) * (param[*m] + diff)) * E[(*m + *d * *i) + *n * *d * (*k + *d * *j)];
	}
	if(*k == *m && *m != *l){
		res = (-2 * square(param[*m]) - 2 * param[*m] * diff + square(diff)) / (square(param[*m]) * (param[*m] + diff)) * E[(*l + *d * *i) + *n * *d * (*m + *d * *j)];		
	}
	if(*k != *m && *m != *l){
		res = square(diff) * (diff + param[*m]) / (square(param[*m]) * (3 * square(param[*m]) + 3 * param[*m] * diff + square(diff))) * E[(*l + *d * *i) + *n * *d * (*k + *d * *j)];
	}
	return(res);
}

// Matern-3/2
double Matern3_2_dxdydp(const double *X, const int *n, const int *d, const int *i, const int *j, const double *param, const int *k, const int *l, const int *m, const double *D, const double *E) {
	double res = 0.;
	double diff = sqrt(3.) * fabs(X[*i + *n * *m] - X[*j + *n * *m]);
	if(*k == *m && *m == *l){
		res = (square(diff) - param[*m] * (diff + 2 * param[*m])) / (square(param[*m]) * (diff + param[*m])) * E[(*m + *d * *i) + *n * *d * (*m + *d * *j)] +
			sqrt(3.) / square(param[*m]) * D[*m + *d * (*j + *n * *i)];			
	}
	if(*k != *m && *m == *l){
		res = (square(diff) - param[*m] * (diff + 2 * param[*m])) / (square(param[*m]) * (diff + param[*m])) * E[(*m + *d * *i) + *n * *d * (*k + *d * *j)];
	}
	if(*k == *m && *m != *l){
		res = (square(diff) - param[*m] * (diff + 2 * param[*m])) / (square(param[*m]) * (diff + param[*m])) * E[(*l + *d * *i) + *n * *d * (*m + *d * *j)];
	}
	if(*k != *m && *m != *l){
		res = square(diff) / (square(param[*m]) * (diff + param[*m])) * E[(*l + *d * *i) + *n * *d * (*k + *d * *j)];
	}
	return(res);
}


// Gaussian
double Gaussian_dxdydp(const double *X, const int *n, const int *d, const int *i, const int *j, const double *param, const int *k, const int *l, const int *m, const double *D, const double *E) {
	double res = 0.;
	if(*k == *m && *m == *l){
		res = (square((X[*i + *n * *m] - X[*j + *n * *m]) / param[*m]) - 2) / param[*m] * E[(*m + *d * *i) + *n * *d * (*m + *d * *j)] +
			2 * (X[*i + *n * *m] - X[*j + *n * *m]) / square(param[*m]) / param[*m] * D[*m + *d * (*j + *n * *i)];	
	}
	if(*k != *m && *m == *l){
		res = (square((X[*i + *n * *m] - X[*j + *n * *m]) / param[*m]) - 2) / param[*m] * E[(*k + *d * *i) + *n * *d * (*m + *d * *j)];
	}
	if(*k == *m && *m != *l){
		res = (square((X[*i + *n * *m] - X[*j + *n * *m]) / param[*m]) - 2) / param[*m] * E[(*m + *d * *i) + *n * *d * (*l + *d * *j)];
	}
	if(*k != *m && *m != *l){
		res = square((X[*i + *n * *m] - X[*j + *n * *m]) / param[*m]) / param[*m] * E[(*l + *d * *i) + *n * *d * (*k + *d * *j)];	
	}
	return(res);
}


/* correlation matrix with third derivatives with respect to x1, x2 and theta */


void corMat_dxdydp(const double *X, const int *n, const int *d, const double *param, const char **type, int *m, const double *D, const double *E, double *ans) {
	(*m)--;
	double (*corFun_dxdydp)(const double *, const int *, const int *, const int *, const int *, const double *, const int *, const int *, const int *, const double *, const double *);

	if (strcmp(*type, "gaussian") == 0)  corFun_dxdydp = Gaussian_dxdydp;
	else if (strcmp(*type, "matern3_2") == 0) corFun_dxdydp = Matern3_2_dxdydp;
	else if (strcmp(*type, "matern5_2") == 0) corFun_dxdydp = Matern5_2_dxdydp;

	for(int i = 0; i < *n; i++){
		for(int k = 0; k < *d; k++){
			for(int j = 0; j < i; j++){
				for(int l = 0; l < k; l++){	
					ans[(k + *d * i) * *n * *d + j * *d + l] = ans[(l + *d * j) * *n * *d + i * *d + k] = ans[(l + *d * i) * *n * *d + j * *d + k] = ans[(k + *d * j) * *n * *d + i * *d + l] = 
					(*corFun_dxdydp)(X, n, d, &i, &j, param, &k, &l, m, D, E);
				}
				ans[(k + *d * i) * *n * *d + j * *d + k] = ans[(k + *d * j) * *n * *d + i * *d + k] = (*corFun_dxdydp)(X, n, d, &i, &j, param, &k, &k, m, D, E);
			}
			ans[(k + *d * i) * *n * *d + i * *d + k] = (*corFun_dxdydp)(X, n, d, &i, &i, param, &k, &k, m, D, E);
		}
	}
}

