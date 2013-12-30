// Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
// <klambauer@bioinf.jku.at>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>


double mydnbinom (double x, double r, double mu){
	double res;
	res = exp(lgamma(x+r)-lgamma(x+1.0)-lgamma(r)+x*log(mu/(mu+r))+r*log(r/(r+mu)));
	return(res);
}


double mydpois (double x, double lambda){
	double res;
	res = exp(x*log(lambda)-lgamma(x+1)-lambda);
	return(res);
}


double rEquation (double* x, double rr, double* A, double* alphaRowMean,
		int N, int n, int componentIdx, double wm, int gg, double eta){
	double sumDigammaXR=0.0;
	int i;
	double resEquation;
	for (i=0;i<N;i++){
		sumDigammaXR = sumDigammaXR+digamma(x[i]+rr)*A[componentIdx+n*i+gg*n*N];
	}
	resEquation = sumDigammaXR-N*alphaRowMean[componentIdx]*digamma(rr)+
			N*alphaRowMean[componentIdx]*log(rr/(rr+wm)) - eta;

	return(resEquation);
}

/* New find_r_bisection */
double find_r_bisection(double* x, double* A, double* alphaRowMean,
		int N, int n, int j, double* wm, int gg, int maxInnerCyc,
		double eta, double rmax) {
	double retval = 0;

	double EPS = 1e-12,EPS2=1e-6;
	double leftPoint=0.01;
	double rightPoint=leftPoint+EPS;
	double middlePoint = (leftPoint + rightPoint) / 2;

	double resEquLeft = rEquation(x, leftPoint, A, alphaRowMean, N, n, j, wm[j], gg, eta);
	double resEquRight = 0;

	double signLeft=sign(resEquLeft);
	double signRight=signLeft;

	int innerCyc = 0;

	while((signRight==signLeft) && (innerCyc<maxInnerCyc)){
		rightPoint=2*rightPoint;
		resEquRight = rEquation(x, rightPoint, A, alphaRowMean, N,n, j, wm[j], gg, eta);
		signRight=sign(resEquRight);
		innerCyc++;
	}

	// cannot solve equation, return upper bound
	if ((innerCyc==maxInnerCyc) && rmax > 0)
		return rmax;
	// cannot solve equation, without upper bound, return high value
	if ((innerCyc==maxInnerCyc) && rmax < 0)
		return 1e8;

	//start bisection
	innerCyc=0;
	middlePoint=0.5*(rightPoint+leftPoint);
	double resEquMiddle = rEquation(x, middlePoint, A, alphaRowMean, N, n, j, wm[j], gg, eta);
	double signMiddle = sign(resEquMiddle);

	while ((innerCyc<maxInnerCyc) && (fabs(resEquMiddle)>EPS2)){

		if (signLeft!=signMiddle){
			rightPoint=middlePoint;
			resEquRight = resEquMiddle;
			signRight = signMiddle;
		} else {
			leftPoint=middlePoint;
			resEquLeft = resEquMiddle;
			signLeft = signMiddle;
		}
		middlePoint = 0.5*(rightPoint+leftPoint);
		resEquMiddle = rEquation(x, middlePoint, A, alphaRowMean, N,n,  j, wm[j], gg, eta);
		signMiddle = sign(resEquMiddle);

		innerCyc = innerCyc+1;
	}
	//retval=middlePoint+eps3;
	retval = middlePoint;
	if (rmax >0 && retval>rmax)
		retval=rmax;

	return(retval);
}



/*  find r by regular falsi */

double find_r_regulafalsi(double* x, double* A, double* alphaRowMean,
		int N, int n, int j, double* wm, int gg, int maxInnerCyc,
		double eta, double rmax) {
	double retval = 0;

	double EPS = 1e-12,EPS2=1e-6;
	double leftPoint=0.01;
	double rightPoint=leftPoint+EPS;
	double middlePoint = (leftPoint + rightPoint) / 2;

	double resEquLeft = rEquation(x, leftPoint, A, alphaRowMean, N, n, j, wm[j], gg, eta);
	double resEquRight = 0;

	double signLeft=sign(resEquLeft);
	double signRight=signLeft;

	int innerCyc = 0;

	while((signRight==signLeft) && (innerCyc<maxInnerCyc)){
		rightPoint=2*rightPoint;
		resEquRight = rEquation(x, rightPoint, A, alphaRowMean, N,n, j, wm[j], gg, eta);
		signRight=sign(resEquRight);
		innerCyc++;
	}

	// can't solve equation, bail out
	if ((innerCyc==maxInnerCyc))
		return rmax;

	//start bisection
	innerCyc=0;
	middlePoint=0.5*(rightPoint+leftPoint);
	double resEquMiddle = rEquation(x, middlePoint, A, alphaRowMean, N, n, j, wm[j], gg, eta);
	double signMiddle = sign(resEquMiddle);

	int direction = 0; // direction we used in the last step, needed for Illinois algorithm.
	while ((innerCyc<maxInnerCyc) && (fabs(resEquMiddle)>EPS2)){

		if (signRight == signMiddle){
			rightPoint=middlePoint;
			resEquRight = resEquMiddle;
			signRight = signMiddle;
			++direction;
		} else {
			leftPoint=middlePoint;
			resEquLeft = resEquMiddle;
			signLeft = signMiddle;
			--direction;
		}

		if (abs(direction) >= 2) { // if we take the same direction twice, do a bisection step instead
			middlePoint = 0.5*(rightPoint+leftPoint);
			direction = 0;
		} else
			middlePoint = (resEquRight*leftPoint - resEquLeft*rightPoint) / (resEquRight - resEquLeft);
		resEquMiddle = rEquation(x, middlePoint, A, alphaRowMean, N,n,  j, wm[j], gg, eta);
		signMiddle = sign(resEquMiddle);

		innerCyc = innerCyc+1;
	}
	//retval=middlePoint+eps3;
	retval = middlePoint;


	if (rmax >0 && retval>rmax)
		retval=rmax;

	return(retval);

}


// if classes==NULL, then we run in unsupervised mode!
SEXP dexus_impl(SEXP XS, SEXP nrowS, SEXP ncolS,
		SEXP alphaINITS, SEXP rINITS, SEXP meansINITS,
		SEXP gammaS, SEXP nS, SEXP cycS, SEXP varToMeanTS, int* classes,
		SEXP etaS, SEXP minMuS, SEXP rmaxS) {
	// calculating the log-likelihood with the basic algorithm
	double* X = REAL(XS);
	int g = INTEGER(nrowS)[0];
	int N =INTEGER(ncolS)[0];

	const double EPS = 1e-100;
	const double DEFAULT_R = 20; // R value used when re-setting clusters

	double* alphaINIT=REAL(alphaINITS);
	double* rINIT=REAL(rINITS);
	double* meansINIT=REAL(meansINITS);

	double varToMeanT=REAL(varToMeanTS)[0];
	double* eta = REAL(etaS);
	double minMu = REAL(minMuS)[0];
	double rmax = REAL(rmaxS)[0];

	double* gamma=REAL(gammaS);

	int n = INTEGER(nS)[0];
	int cyc = INTEGER(cycS)[0];

	double gammaSum;
	int iter,i,j,gg,maxInnerCyc=100;

	double alphaColSum;
	double *alphaRowMean;
	alphaRowMean = (double*) R_alloc(n, sizeof(double));

	double *V2; // squared sum of weights
	V2 = (double*) R_alloc(n, sizeof(double));

	double *x;
	x = (double*) R_alloc(N, sizeof(double));

	double *wm;
	wm = (double*) R_alloc(n, sizeof(double));

	double *wvar;
	wvar = (double*) R_alloc(n, sizeof(double));

	//SEXP p_RET;
	//PROTECT(p_RET = allocVector(REALSXP, g*n));
	//double* p=REAL(p_RET);

	SEXP alpha_RET;
	PROTECT(alpha_RET = allocVector(REALSXP, g*n));
	double* alpha=REAL(alpha_RET);

	SEXP r_RET;
	PROTECT(r_RET = allocVector(REALSXP, g*n));
	double* r=REAL(r_RET);

	SEXP A_RET;
	PROTECT(A_RET = allocVector(REALSXP, g*n*N));
	double* A=REAL(A_RET);

	SEXP means_RET;
	PROTECT(means_RET = allocVector(REALSXP, g*n));
	double* means=REAL(means_RET);

	gammaSum=0.0;
	for (j=0;j<n;j++){
		gammaSum = gammaSum+gamma[j];
	}

	GetRNGstate();

	//loop along rows of the data matrix, genes, ROI
	for (gg=0;gg<g;gg++){
		//Rprintf("New Gene..\n");

		//initializing alpha, p, r
		for (j=0;j<n;j++){
			alpha[j+gg*n]=alphaINIT[j];
			r[j+gg*n]=rINIT[j+gg*n];
			means[j+gg*n] = meansINIT[j+gg*n];
		}

		x = X + gg*N; // current gene

		for (iter=0;iter<=cyc;iter++){
			// E-step
			for (j=0;j<n;j++){
				alphaRowMean[j]=0.0;
				V2[j]=0.0;
			}

			for (i=0;i<N;i++){
				alphaColSum=0.0;
				for (j=0;j<n;j++){

					if (classes != NULL) {  // supervised!
						if (classes[i]==j){
							A[j+n*i+gg*n*N] = 1.0;
						} else {
							A[j+n*i+gg*n*N] = 0.0;
						}
					} else { // unsupervised!
						//if (x[i] < minMu) x[i]=minMu;
						if (r[j+gg*n] > 0)
							A[j+n*i+gg*n*N] = alpha[j+gg*n]*mydnbinom(x[i],r[j+gg*n],means[j+gg*n]);
						else
							A[j+n*i+gg*n*N] = alpha[j+gg*n]*mydpois(x[i],means[j+gg*n]);
					}

					if (A[j+n*i+gg*n*N] < EPS){
						A[j+n*i+gg*n*N] = EPS;
					}
					alphaColSum = alphaColSum+A[j+n*i+gg*n*N];
				}
				for (j=0;j<n;j++){
					A[j+n*i+gg*n*N] = A[j+n*i+gg*n*N]/alphaColSum;
					alphaRowMean[j]=alphaRowMean[j]+A[j+n*i+gg*n*N]/N;
					V2[j]=V2[j]+A[j+n*i+gg*n*N]*A[j+n*i+gg*n*N];
				}
			}

			// Calculating weighted mean and weighted variance for each component
			// Substitute with online update for var to make efficient
			for (j=0;j<n;j++){
				//Rprintf("Component: %lf\n", (double) j);

				wm[j]=EPS;
				for (i=0;i<N;i++){
					wm[j] = wm[j]+x[i]*A[j+n*i+gg*n*N];
				}
				wm[j]=wm[j]/(N*alphaRowMean[j]);
				//Rprintf("weightedMean: %lf\n", wm[j]);

				wvar[j]=0.0;
				for (i=0;i<N;i++){
					wvar[j] = wvar[j]+A[j+n*i+gg*n*N]*(x[i]-wm[j])*(x[i]-wm[j]);
				}
				//unbiased??
				//wvar[j]=wvar[j]*(N*alphaRowMean[j]/(pow(N*alphaRowMean[j],2)+V2[j]) );
				//biased:
				wvar[j]=wvar[j]*N*alphaRowMean[j]/(N*N*alphaRowMean[j]*alphaRowMean[j]-V2[j]);
				//(wvar <- sum(a*(x-wm)^2)*sum(a)/(sum(a)^2-a2))

				//Rprintf("weightedVar: %lf\n", wvar[j]);
			}

			// updating components r, p, alpha
			for (j=0;j<n;j++) {

				// updating components r, p, alpha
				means[j+gg*n] = wm[j];
				if (iter<cyc){
					alpha[j+gg*n] = (alphaRowMean[j]+1.0/n*(gamma[j]-1.0))/(1.0+1.0/n*(gammaSum-n));
				} else {
					alpha[j+gg*n] = alphaRowMean[j];
				}


				if (wvar[j]/wm[j] < varToMeanT) { // Poisson component
					r[j+gg*n] = -1.0;
					//r[j+gg*n] = rmax;
				} else {
					r[j+gg*n] = find_r_bisection(x, A, alphaRowMean, N, n, j, wm,
							gg, maxInnerCyc, eta[gg], -1.0);
					//r[j+gg*n] = find_r_regulafalsi(x, A, alphaRowMean, N, n, j, wm,
					//						gg, maxInnerCyc, eta[gg], -1.0);
					//r[j+gg*n] = find_r_bisection(x, A, alphaRowMean, N, n, j, wm,
					//							gg, maxInnerCyc, eta[gg], rmax);

					//p[j+gg*n] = wm[j]/(wm[j]+r[j+gg*n]+EPS);
					if (r[j+gg*n] > 10000) {
						r[j+gg*n] = -1.0; //Approx with Poisson
					}
				}


				if (means[j+gg*n]< minMu){
					means[j+gg*n] = minMu;
				}
			}
			/*
			// look for overlapping clusters
			for (j = 0; j < n; j++) {
				for (int k = j+1; k < n; k++) {
					double d = means[gg*n + j] - means[gg*n + k];
					if (d*d < 1.0) { // replace cluster center with random datapoint
						int idx = int(unif_rand() * N);
						means[gg*n + k] = std::max(x[idx], minMu);
						r[gg*n + k] = DEFAULT_R;
					}
				}
			}
			 */
		}

		if (rmax > 0){
			for (j=0;j<n;j++) {
				if ((r[j+gg*n] > rmax) || (r[j+gg*n] < 0)) {
					r[j+gg*n] = rmax;
				}
			}
		}



		R_CheckUserInterrupt();
	} // end over data rows

	PutRNGstate();

	SEXP namesRET;
	PROTECT(namesRET = allocVector(STRSXP, 4));
	//SET_STRING_ELT(namesRET, 0, mkChar("p"));
	SET_STRING_ELT(namesRET, 0, mkChar("r"));
	SET_STRING_ELT(namesRET, 1, mkChar("alpha"));
	SET_STRING_ELT(namesRET, 2, mkChar("A"));
	SET_STRING_ELT(namesRET, 3, mkChar("means"));

	SEXP RET;
	PROTECT(RET = allocVector(VECSXP, 4));
	//SET_VECTOR_ELT(RET, 0, p_RET);
	SET_VECTOR_ELT(RET, 0, r_RET);
	SET_VECTOR_ELT(RET, 1, alpha_RET);
	SET_VECTOR_ELT(RET, 2, A_RET);
	SET_VECTOR_ELT(RET, 3, means_RET);
	setAttrib(RET, R_NamesSymbol, namesRET);

	UNPROTECT(6);
	return(RET);
}


extern "C" SEXP dexus(SEXP XS, SEXP nrowS, SEXP ncolS,
		SEXP alphaINITS, SEXP rINITS, SEXP meansINITS,
		SEXP gammaS, SEXP nS, SEXP cycS, SEXP varToMeanTS,
		SEXP etaS, SEXP minMuS, SEXP rmaxS) {
	return dexus_impl(XS, nrowS, ncolS, alphaINITS, rINITS, meansINITS, gammaS,
			nS, cycS, varToMeanTS, NULL, etaS, minMuS, rmaxS);
}


extern "C" SEXP dexusS(SEXP XS, SEXP nrowS, SEXP ncolS,
		SEXP alphaINITS, SEXP rINITS, SEXP meansINITS,
		SEXP gammaS, SEXP nS, SEXP cycS, SEXP varToMeanTS, SEXP classesS,
		SEXP etaS, SEXP minMuS, SEXP rmaxS) {
	int* classes=INTEGER(classesS);
	return dexus_impl(XS, nrowS, ncolS, alphaINITS, rINITS, meansINITS, gammaS, nS,
			cycS, varToMeanTS, classes, etaS, minMuS, rmaxS);
}


/*
system("R CMD SHLIB /home/klambaue/Dropbox/WorkspaceSequencing/projects/mRNASeq/mixnb/src/mixnb.cpp")
dyn.load("/home/klambaue/Dropbox/WorkspaceSequencing/projects/mRNASeq/mixnb/src/mixnb.so")

#system("R CMD SHLIB /files/work/bioinf/mixnb/code/mixnb/src/mixnb.cpp")
#dyn.load("/files/work/bioinf/mixnb/code/mixnb/src/mixnb.so")

x <- rnbinom(mu=10,size=1/2, n=100)
res <- .C("find_r_regulafalsi_ext",retval=as.numeric(0), as.numeric(x),
				as.numeric(rep(1,length(x))),
				as.numeric(1),
				as.integer(length(x)),
				as.integer(1),
				as.integer(0),
				as.numeric(mean(x)),
				as.integer(0),
				as.integer(10000),
				as.numeric(0))$retval
(res)
getRNBrobust(x)

 */

//find_r_regulafalsi(double* x, double* A, double* alphaRowMean, int N, int n, int j, double* wm, int gg, int maxInnerCyc, double eta)

extern "C" void find_r_bisection_ext(double* returnvalue, double* x,
		double* A, double* alphaRowMean,
		int* N, int* n, int* j, double* wm, int* gg, int* maxInnerCyc,
		double* eta, double* rmax) {
	*returnvalue = find_r_bisection(x, A, alphaRowMean, *N, *n, *j, wm,
			*gg, *maxInnerCyc, *eta, *rmax);
}

extern "C" void find_r_regulafalsi_ext(double* returnvalue, double* x,
		double* A, double* alphaRowMean,
		int* N, int* n, int* j, double* wm, int* gg, int* maxInnerCyc,
		double* eta, double* rmax) {
	*returnvalue = find_r_regulafalsi(x, A, alphaRowMean, *N, *n, *j, wm,
			*gg, *maxInnerCyc, *eta, *rmax);
}

/**
 * This is a C re-implementation of nbinomTestForMatrices 
 * from DESeq (Anders & Huber, 2010), based on DESeq version 1.8.3
 * 
 * This runs ~ 10-30% faster than the R version.
 *
 * Going from the R code to the C code:
 * aS => kAs
 * muA => mu * sum(sizeFactorsA)
 * sizeAS => 1 / sumDispsA[i]
 * sumSizeFactorsAS = sum(sizeFactorsA)
 */
extern "C" SEXP dexus_pval_calculation(SEXP aS, SEXP bS, SEXP muAS,
		SEXP muBS, SEXP sizeAS, SEXP sizeBS, SEXP sumSizeFactorsAS,
		SEXP sumSizeFactorsBS) {
	const int* a = INTEGER(aS);
	const int* b = INTEGER(bS);
	const double* muA = REAL(muAS);
	const double* muB = REAL(muBS);
	const double* sizeA = REAL(sizeAS);
	const double* sizeB = REAL(sizeBS);
	const double sumSizeFactorsA = REAL(sumSizeFactorsAS)[0];
	const double sumSizeFactorsB = REAL(sumSizeFactorsBS)[0];
	int n = LENGTH(aS);

	SEXP pvalS;
	PROTECT(pvalS = allocVector(REALSXP, n));
	double* pval = REAL(pvalS);

	// find maximum value of a[i] + b[i], so we only have to alloc ps once
	int maxSum = 0;
	for (int i = 0; i < n; ++i)
		if (a[i] + b[i] > maxSum) maxSum = a[i] + b[i];
	double* psums = (double*) R_alloc(maxSum+1, sizeof(double));


	for (int i = 0; i < n; ++i) {
		if (a == 0 && b == 0) {
			pval[i] = NA_REAL;
			continue;
		}

		int m = a[i] + b[i];
		double sum = 0;
		for (int ks = 0; ks <= m; ++ks) {
			psums[ks] = dnbinom(ks, sizeA[i], muA[i], 0) * dnbinom(m - ks, sizeB[i], muB[i], 0);
			sum += psums[ks];			
		}

		//double pobs = dnbinom_mu(a[i], sizeA[i], muA[i], 0) * dnbinom_mu(b[i], sizeB[i], muB[i], 0);
		//assert(pobs == ps[a[i]]);

		double numer = 0;
		if (a[i] * sumSizeFactorsB < b[i] * sumSizeFactorsA){
			for (int j = 0; j < a[i] + 1; ++j)
				numer += psums[j];
		} else {
			for (int j = a[i]; j <= m ; ++j)
				numer += psums[j];	
		}
		pval[i] = 2 * (numer / sum);
		pval[i] = pval[i] <  1.0 ? pval[i] : 1.0;
	}

	UNPROTECT(1);
	return pvalS;
}

/*
int main() {
	//double x[] = {220, 232, 236, 238, 242, 209, 217, 254, 267, 221, 237, 237};
	double x[] = {511, 528, 526, 465, 512, 468, 500, 478, 531, 495, 495, 530};
	double A[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	double alphaRowMean[] = {1};
	int N = 12;
	int n = 1;
	int j = 0;
	double wm = 0;
	for (int i = 0; i < N; ++i) wm += x[i];
	wm /= N;
	int maxInnerCyc = 10000;
	double r1 = find_r_bisection(x, A, alphaRowMean, N, n, j, &wm, 0, maxInnerCyc, 0);
	double r2 = find_r_regulafalsi(x, A, alphaRowMean, N, n, j, &wm, 0, maxInnerCyc, 0);
	printf("bs: %f\nrf: %f\n", r1, r2);
	return 0;	
}
 */
// DEXUS SEMI-SUPERVISED:: DEXSS ###############################################

SEXP dexss_impl(SEXP XS, SEXP nrowS, SEXP ncolS,
		SEXP alphaINITS, SEXP rINITS, SEXP meansINITS,
		SEXP gammaS, SEXP nS, SEXP cycS, SEXP varToMeanTS, int* classes,
		SEXP etaS, SEXP minMuS, SEXP rmaxS) {
	// calculating the log-likelihood with the basic algorithm
	double* X = REAL(XS);
	int g = INTEGER(nrowS)[0];
	int N =INTEGER(ncolS)[0];

	const double EPS = 1e-100;
	const double DEFAULT_R = 20; // R value used when re-setting clusters

	double* alphaINIT=REAL(alphaINITS);
	double* rINIT=REAL(rINITS);
	double* meansINIT=REAL(meansINITS);

	double varToMeanT=REAL(varToMeanTS)[0];
	double* eta = REAL(etaS);
	double minMu = REAL(minMuS)[0];
	double rmax = REAL(rmaxS)[0];

	double* gamma=REAL(gammaS);

	int n = INTEGER(nS)[0];
	int cyc = INTEGER(cycS)[0];

	double gammaSum;
	int iter,i,j,gg,maxInnerCyc=100;

	double alphaColSum;
	double *alphaRowMean;
	alphaRowMean = (double*) R_alloc(n, sizeof(double));

	double *V2; // squared sum of weights
	V2 = (double*) R_alloc(n, sizeof(double));

	double *x;
	x = (double*) R_alloc(N, sizeof(double));

	double *wm;
	wm = (double*) R_alloc(n, sizeof(double));

	double *wvar;
	wvar = (double*) R_alloc(n, sizeof(double));

	//SEXP p_RET;
	//PROTECT(p_RET = allocVector(REALSXP, g*n));
	//double* p=REAL(p_RET);

	SEXP alpha_RET;
	PROTECT(alpha_RET = allocVector(REALSXP, g*n));
	double* alpha=REAL(alpha_RET);

	SEXP r_RET;
	PROTECT(r_RET = allocVector(REALSXP, g*n));
	double* r=REAL(r_RET);

	SEXP A_RET;
	PROTECT(A_RET = allocVector(REALSXP, g*n*N));
	double* A=REAL(A_RET);

	SEXP means_RET;
	PROTECT(means_RET = allocVector(REALSXP, g*n));
	double* means=REAL(means_RET);

	gammaSum=0.0;
	for (j=0;j<n;j++){
		gammaSum = gammaSum+gamma[j];
	}

	GetRNGstate();

	//loop along rows of the data matrix, genes, ROI
	for (gg=0;gg<g;gg++){
		//Rprintf("New Gene..\n");

		//initializing alpha, p, r
		for (j=0;j<n;j++){
			alpha[j+gg*n]=alphaINIT[j];
			r[j+gg*n]=rINIT[j+gg*n];
			means[j+gg*n] = meansINIT[j+gg*n];
		}

		x = X + gg*N; // current gene

		for (iter=0;iter<=cyc;iter++){
			// E-step
			for (j=0;j<n;j++){
				alphaRowMean[j]=0.0;
				V2[j]=0.0;
			}

			for (i=0;i<N;i++){
				alphaColSum=0.0;
				for (j=0;j<n;j++){

					// semi-supervised!
					if (classes[i] >= 0) {
						if (classes[i]==j){
							A[j+n*i+gg*n*N] = 1.0;
						} else {
							A[j+n*i+gg*n*N] = 0.0;
						}
					} else {
						//if (x[i] < minMu) x[i]=minMu;
						if (r[j+gg*n] > 0)
							A[j+n*i+gg*n*N] = alpha[j+gg*n]*mydnbinom(x[i],r[j+gg*n],means[j+gg*n]);
						else
							A[j+n*i+gg*n*N] = alpha[j+gg*n]*mydpois(x[i],means[j+gg*n]);
					}

					if (A[j+n*i+gg*n*N] < EPS){
						A[j+n*i+gg*n*N] = EPS;
					}
					alphaColSum = alphaColSum+A[j+n*i+gg*n*N];
				}
				for (j=0;j<n;j++){
					A[j+n*i+gg*n*N] = A[j+n*i+gg*n*N]/alphaColSum;
					alphaRowMean[j]=alphaRowMean[j]+A[j+n*i+gg*n*N]/N;
					V2[j]=V2[j]+A[j+n*i+gg*n*N]*A[j+n*i+gg*n*N];
				}
			}

			// Calculating weighted mean and weighted variance for each component
			// Substitute with online update for var to make efficient
			for (j=0;j<n;j++){
				//Rprintf("Component: %lf\n", (double) j);

				wm[j]=EPS;
				for (i=0;i<N;i++){
					wm[j] = wm[j]+x[i]*A[j+n*i+gg*n*N];
				}
				wm[j]=wm[j]/(N*alphaRowMean[j]);
				//Rprintf("weightedMean: %lf\n", wm[j]);

				wvar[j]=0.0;
				for (i=0;i<N;i++){
					wvar[j] = wvar[j]+A[j+n*i+gg*n*N]*(x[i]-wm[j])*(x[i]-wm[j]);
				}
				//unbiased??
				//wvar[j]=wvar[j]*(N*alphaRowMean[j]/(pow(N*alphaRowMean[j],2)+V2[j]) );
				//biased:
				wvar[j]=wvar[j]*N*alphaRowMean[j]/(N*N*alphaRowMean[j]*alphaRowMean[j]-V2[j]);
				//(wvar <- sum(a*(x-wm)^2)*sum(a)/(sum(a)^2-a2))

				//Rprintf("weightedVar: %lf\n", wvar[j]);
			}

			// updating components r, p, alpha
			for (j=0;j<n;j++) {

				// updating components r, p, alpha
				means[j+gg*n] = wm[j];
				if (iter<cyc){
					alpha[j+gg*n] = (alphaRowMean[j]+1.0/n*(gamma[j]-1.0))/(1.0+1.0/n*(gammaSum-n));
				} else {
					alpha[j+gg*n] = alphaRowMean[j];
				}

				if (wvar[j]/wm[j] < varToMeanT) { // Poisson component
					r[j+gg*n] = -1.0;
					//r[j+gg*n] = rmax;
				} else {
					r[j+gg*n] = find_r_bisection(x, A, alphaRowMean, N, n, j, wm,
							gg, maxInnerCyc, eta[gg], -1.0);
					//r[j+gg*n] = find_r_regulafalsi(x, A, alphaRowMean, N, n, j, wm,
					//						gg, maxInnerCyc, eta[gg], -1.0);
					//r[j+gg*n] = find_r_bisection(x, A, alphaRowMean, N, n, j, wm,
					//							gg, maxInnerCyc, eta[gg], rmax);

					//p[j+gg*n] = wm[j]/(wm[j]+r[j+gg*n]+EPS);
					if (r[j+gg*n] > 10000) {
						r[j+gg*n] = -1.0; //Approx with Poisson
					}
				}


				if (means[j+gg*n]< minMu){
					means[j+gg*n] = minMu;
				}
			}
			/*
			// look for overlapping clusters
			for (j = 0; j < n; j++) {
				for (int k = j+1; k < n; k++) {
					double d = means[gg*n + j] - means[gg*n + k];
					if (d*d < 1.0) { // replace cluster center with random datapoint
						int idx = int(unif_rand() * N);
						means[gg*n + k] = std::max(x[idx], minMu);
						r[gg*n + k] = DEFAULT_R;
					}
				}
			}
			 */
		}

		if (rmax > 0){
			for (j=0;j<n;j++) {
				if ((r[j+gg*n] > rmax) || (r[j+gg*n] < 0)) {
					r[j+gg*n] = rmax;
				}
			}
		}



		R_CheckUserInterrupt();
	} // end over data rows

	PutRNGstate();

	SEXP namesRET;
	PROTECT(namesRET = allocVector(STRSXP, 4));
	//SET_STRING_ELT(namesRET, 0, mkChar("p"));
	SET_STRING_ELT(namesRET, 0, mkChar("r"));
	SET_STRING_ELT(namesRET, 1, mkChar("alpha"));
	SET_STRING_ELT(namesRET, 2, mkChar("A"));
	SET_STRING_ELT(namesRET, 3, mkChar("means"));

	SEXP RET;
	PROTECT(RET = allocVector(VECSXP, 4));
	//SET_VECTOR_ELT(RET, 0, p_RET);
	SET_VECTOR_ELT(RET, 0, r_RET);
	SET_VECTOR_ELT(RET, 1, alpha_RET);
	SET_VECTOR_ELT(RET, 2, A_RET);
	SET_VECTOR_ELT(RET, 3, means_RET);
	setAttrib(RET, R_NamesSymbol, namesRET);

	UNPROTECT(6);
	return(RET);
}


extern "C" SEXP dexss(SEXP XS, SEXP nrowS, SEXP ncolS,
		SEXP alphaINITS, SEXP rINITS, SEXP meansINITS,
		SEXP gammaS, SEXP nS, SEXP cycS, SEXP varToMeanTS, SEXP classesS,
		SEXP etaS, SEXP minMuS, SEXP rmaxS) {

	int* classes=INTEGER(classesS);
	return dexss_impl(XS, nrowS, ncolS, alphaINITS, rINITS, meansINITS, gammaS, nS,
			cycS, varToMeanTS, classes, etaS, minMuS, rmaxS);

}

