#include <math.h>
#include <R.h>
#include <Rinternals.h>
//#include <Rinterface.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>
//SEXP La_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v, SEXP method); 
SEXP La_svd(SEXP jobu, SEXP x, SEXP s, SEXP u, SEXP v); // 3.0.3 http://fossies.org/dox/R-3.0.3/Lapack_8c_source.html


#define MAXCLASSES	1024
int knn_class[MAXCLASSES];

#define MAXSAMPLES	1024*10
int knn_samples[MAXSAMPLES];
double knn_distances[MAXSAMPLES];

#define distanceForIJ(d, n, i, j)	d[n*(i-1) - i*(i-1)/2 + j-i     - 1]
// i, j are one-based indexes, d is 0-based index array

#ifndef R_STATS_H
#define R_STATS_H

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif


#endif


#define both_FINITE(a,b) (R_FINITE(a) && R_FINITE(b))
#ifdef R_160_and_older
#define both_non_NA both_FINITE
#else
#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))
#endif


SEXP galgo_error(char *msg) {
	SEXP r;
	Rprintf(msg);
	PROTECT(r = allocVector(INTSXP,1));
	INTEGER(r)[0] = 0;
	UNPROTECT(1);
	return r;
}


// this does only one evaluation
SEXP galgoKNN(SEXP distance, SEXP samplesClasses, SEXP trainingSamplesVector, SEXP testSamplesVector, SEXP neighbours, SEXP min_neighbours, SEXP resultType) {
	// distance is a distance object between all samples
	// training & test vectors are 1 based indexes (R-like usage)
	// result are the same length than test vector


	// the dissimilarity between (row) i and j is do[n*(i-1) - i*(i-1)/2 + j-i]. IN 1-based index array

	if (!isReal(distance)) {
		return (galgo_error("distance should be double. Use as.double().\n"));
	}
	if (!isInteger(samplesClasses)) {
		return (galgo_error("samplesClasses should be integer. Use as.integer().\n"));
	}
	if (!isInteger(trainingSamplesVector)) {
		return (galgo_error("trainingSamplesVector should be integer. Use as.integer().\n"));
	}
	if (!isInteger(testSamplesVector)) {
		return (galgo_error("testSamplesVector should be integer. Use as.integer().\n"));
	}
	if (!isInteger(neighbours)) {
		return (galgo_error("neighbours should be integer. Use as.integer().\n"));
	}
	if (!isInteger(min_neighbours)) {
		return (galgo_error("min_neighbours should be integer. Use as.integer().\n"));
	}
	if (!isInteger(resultType)) {
		return (galgo_error("resultType should be integer. Use as.integer().\n"));
	}

	int knn = INTEGER(neighbours)[0];
	int mnn = INTEGER(min_neighbours)[0];

	int i, j, k, l, maxclass, minclass, nclasses;
	int nl = LENGTH(samplesClasses);
	int ntr = LENGTH(trainingSamplesVector);
	int nte = LENGTH(testSamplesVector);
	int samI, samJ;
	double d;
	SEXP r;

	minclass = maxclass = INTEGER(samplesClasses)[0];
	for (i=1; i < nl; i++) {
		j = INTEGER(samplesClasses)[i];
		if (j > maxclass) maxclass = j;
		if (j < minclass) minclass = j;
	}
	nclasses = maxclass - minclass + 1;

	// resultType != 0, return the number of matched classes
	// otherwise return samplesInClass Vectorz
	PROTECT(r = allocVector(INTSXP,(INTEGER(resultType)[0] ? 1 : nte)));
	INTEGER(r)[0] = 0; // initialize first when resultyType != 0
	if (INTEGER(resultType)[0] == 0) for (i=0; i < nte; i++) INTEGER(r)[i] = -2;
	//  -2 means CLASS WERE NOT PREDICTED    ::: DEPRACATED
	//  -1 means CLASS COULD NOT BE DECIDED
	// >=0 means CLASS PREDICTED


	for (i=0; i < nte; i++) {
		samI = INTEGER(testSamplesVector)[i];
		//Rprintf("Sample [%d], original class:[%d]",samI,INTEGER(samplesClasses)[samI-1]);

		// initialize nearest neighbours
		for (j=0; j < knn; j++)	{ knn_samples[j]=0; knn_distances[j]=9e99; }

		// find nearest neighbours
		for (j=0; j < ntr; j++)	{
			samJ = INTEGER(trainingSamplesVector)[j]; // always in training
			if (samI != samJ) {
				d = (samI > samJ ? distanceForIJ(REAL(distance), nl, samJ, samI) : distanceForIJ(REAL(distance), nl, samI, samJ));
				//Rprintf("Verifying distance [%f] sample [%d]\n",d,samJ);
				if (d < knn_distances[0]) {
					//Rprintf("Adding distance [%f] sample [%d]\n",d,samJ);
					// there has to be a place for sample samJ
					// find out position: the place for samJ is k-1
					for (k=1; (k < knn) && (d < knn_distances[k]); k++);

					// shift all list to mantain it ordered
					for (l=1; l < k; l++) {
						knn_distances[l-1] = knn_distances[l];
						knn_samples[l-1] = knn_samples[l];
					}
					knn_distances[k-1] = d;
					knn_samples[k-1] = samJ;
				}
			}
		}

		//knn_samples has the nearest neighbours

		// initialize counter for classes
		for (j=0; j < nclasses; j++) knn_class[j] = 0;

		// count classes from nearest samples 
		for (j=0; (j < knn); j++) if (knn_samples[j] > 0) knn_class[INTEGER(samplesClasses)[knn_samples[j]-1] - minclass]++;

		// find out the greatest
		// j is the "ok" flag, if it is 0 then class is not trustable
		for (j=1,l=0,k=1; k < nclasses; k++) 
			if (knn_class[k] > knn_class[l]) j = l = k; 
		else if (knn_class[k] == knn_class[l]) j = 0;
		
		k = ((knn_class[l] >= mnn) && (j != 0) ? l+minclass : -1);
		if (INTEGER(resultType)[0]) {
			INTEGER(r)[0] += (k == INTEGER(samplesClasses)[samI-1] ? 1 : 0);
		} else {
			INTEGER(r)[i] = k;
			//Rprintf(", predicted class:[%d]\n",INTEGER(r)[samI-1]);
		}
	}

	UNPROTECT(1);
	return (r);
}







#define MAXCLASSES_MLHD		128
#define MAXVARS_MLHD		1024*4
double means[MAXCLASSES_MLHD*MAXVARS_MLHD];
int samples_class[MAXCLASSES_MLHD];
int samplesInClass[MAXSAMPLES];
double maxMLHD[MAXSAMPLES];
double fac1[MAXVARS_MLHD], fac2[MAXVARS_MLHD];


 /* get the list element named str, or return NULL */
 
 SEXP getListElement(SEXP list, char *str)
 {
   SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
   int i;
 
   for (i = 0; i < length(list); i++)
	 if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	   elmt = VECTOR_ELT(list, i);
	   break;
	 }
   return elmt;
 }


SEXP galgoMLHD(SEXP data, SEXP samplesClasses, SEXP trainingSamplesVector, SEXP testSamplesVector, SEXP jobu, SEXP jobv, SEXP u, SEXP v, SEXP dummy, SEXP method, SEXP cov, SEXP inv, SEXP resultType) {
// data[rows=samples, cols=genes=variables]

	if (!isInteger(samplesClasses)) {
		return (galgo_error("samplesClasses should be integer. Use as.integer().\n"));
	}
	if (!isInteger(trainingSamplesVector)) {
		return (galgo_error("trainingSamplesVector should be integer. Use as.integer().\n"));
	}
	if (!isInteger(testSamplesVector)) {
		return (galgo_error("testSamplesVector should be integer. Use as.integer().\n"));
	}
	if (!isInteger(resultType)) {
		return (galgo_error("resultType should be integer. Use as.integer().\n"));
	}

	int i, j, k, l, maxclass, minclass, nclasses, c, n, ncovclasses;
	int nl = length(samplesClasses);
	int ntr = LENGTH(trainingSamplesVector);
	int nte = LENGTH(testSamplesVector);
	int rows, nvars;
	double d;
	SEXP r;
	rows = nrows(data);
	nvars = ncols(data); // number of variables
	double s = 0.0, xa=0.0, xb=0.0, b=0.0, a=0.0;


	// Total number of classes
	minclass = maxclass = INTEGER(samplesClasses)[0];
	for (i=1; i < nl; i++) {
		j = INTEGER(samplesClasses)[i];
		if (j > maxclass) maxclass = j;
		if (j < minclass) minclass = j;
	}
	nclasses = maxclass - minclass + 1;

	// number of samples per class (in training)
	for (i=0; i < nclasses; i++) samples_class[i]=0;
	for (i=0; i < ntr; i++) samples_class[INTEGER(samplesClasses)[INTEGER(trainingSamplesVector)[i]-1]-minclass]++;

	// mean vectors per class: rows in means are classes, columns are variables
	for (c=0; c < nclasses; c++) for (j=0; j < nvars; j++) means[j*nclasses + c]=0; // initialize to 0
	for (i=0; i < ntr; i++) {
		l = INTEGER(trainingSamplesVector)[i]-1;
		c = (INTEGER(samplesClasses)[l]-minclass);
		for (j=0; j < nvars; j++) means[j*nclasses + c] += REAL(data)[j*rows + l];
	}
	for (j=0; j < nvars; j++) for (c=0; c < nclasses; c++) means[j*nclasses + c] /= samples_class[c];
	//for (c=0; c < nclasses; c++) {
	//	Rprintf("class=%d ",c);
	//	for (j=0; j < nvars; j++) Rprintf("%f ",means[j*nclasses + c]);
	//	Rprintf("\n");
	//}

	// compute covariance matrix
	for (i=0; i < nvars; i++)	for (j=i; j < nvars; j++) REAL(cov)[i*nvars + j] = 0;
	ncovclasses=0;
	for (c=0; c < nclasses; c++) {
		// n=# samples in class c, samplesInClass are indexes for samples in class c
		n=0;
		for (k=0; k < ntr; k++)	{
			l = INTEGER(trainingSamplesVector)[k]-1;
			if ((INTEGER(samplesClasses)[l]-minclass) == c) samplesInClass[n++] = l;
		}
		if (n > 1)	{
			ncovclasses++;
			for (i=0; i < nvars; i++)	{
				for (j=i; j < nvars; j++) {
					s = b = a = 0.0;
					for (k=0; k < n; k++)	{
						l = samplesInClass[k];
						xa = REAL(data)[i*rows + l];
						xb = REAL(data)[j*rows + l];
						s += xa*xb;
						b += xb;
						a += xa;
					}
					s /= (double) (n-1);
					a /= (double) (n-1);
					b /= (double) (n);
					REAL(cov)[i*nvars + j] += s - a*b;
					//if (i != j)	REAL(cov)[j*nvars + i] += s-a*b;
				}
			}
		}
	}
	s = (double) (ntr - ncovclasses);
	for (i=0; i < nvars; i++)	for (j=i; j < nvars; j++) REAL(cov)[j*nvars + i] = (REAL(cov)[i*nvars + j] /= s);
	//Rprintf("cov={\n");
	//for (i=0; i < nvars; i++) { for (j=0; j < nvars; j++) Rprintf("%f ",REAL(cov)[i*nvars + j]); Rprintf("\n"); }
	//Rprintf("}\n");

	//SVD
	//r = La_svd(jobu, jobv, cov, dummy, u, v, method);
	r = La_svd(jobu, cov, dummy, u, v); // R 3.0.3
	//r$d, r$u, r$vt
	//svd: v=t(r$vt)
	//ginv: v %*% (1/d * t(u)) [only "d" positives ( > max(tol * d[0], 0) # THIS IS IMPORTANT, di[i] > 0.00001 <- VERIFY THIS

	//inv=t(r$vt) %*% (1/r$d * t(r$u))
	
	SEXP divar = getListElement(r, "d");
	double *di = REAL(divar);
	SEXP iu = getListElement(r, "u");
	SEXP vi = getListElement(r, "vt");

	
	//compute cov = 1/di * t(iu)
	n = LENGTH(divar);
	for (i=0; i < n; i++) if(di[i] > 0.00001) di[i] = 1.0/di[i];  // by the moment consider all d positives
	for (j=0; j < n; j++) {
		s = 0;
		for (i=0; i < n; i++) REAL(cov)[j*n + i] = di[i] * REAL(iu)[i*n + j];// by the moment consider all d positives
	}

	//compute inv = t(vi) %*% cov
	for (i=0; i < n; i++) {
		for (j=0; j < n; j++) {
			s = 0;
			for (c=0; c < n; c++) s += REAL(vi)[i*n + c] * REAL(cov)[j*n + c];
			REAL(inv)[j*n + i] = s;
		}
	}

	// predict class
	for (i=0; i < nte; i++) { maxMLHD[i] = -9e99; samplesInClass[i] = -1; }
	for (c=0; c < nclasses; c++) {
		// compute fac1 = means[c,] %*% inv
		//Rprintf("t(m) * inv =");
		for (k=0; k < nvars; k++) {
			for (s=0, j = 0; j < nvars; j++) s += means[j*nclasses + c] * REAL(inv)[j*nvars + k];
			fac1[k] = s;
			//Rprintf("%f ",fac1[k]);
		}
		//Rprintf("\n");
		// compute d = .5 * fac1 %*% means[c,]
		for (d=0,k=0; k < nvars; k++) d += fac1[k] * means[k*nclasses + c];
		d *= 0.5;

		// compute s = fac1 %*% datos - d
		for (k=0; k < nte; k++) {
			l = INTEGER(testSamplesVector)[k]-1;
			for (s=0, j=0; j < nvars; j++) s += fac1[j] * REAL(data)[j*rows + l];
			s -= d;
			if (s > maxMLHD[k]) {
				maxMLHD[k] = s;
				samplesInClass[k] = c;
			}
		}
	}

	// samplesInClass is the bloody result

	if (INTEGER(resultType)[0]) {
		// return the number of matched classes
		c = 0;
		for (k=0; k < nte; k++) if ((INTEGER(samplesClasses)[INTEGER(testSamplesVector)[k]-1]-minclass) == samplesInClass[k]) c++;
		PROTECT(r = allocVector(INTSXP, 1));
		INTEGER(r)[0] = c;
	} else {
		// return samplesInClass Vector
		PROTECT(r = allocVector(INTSXP, nte));
		for (k=0; k < nte; k++) INTEGER(r)[k] = samplesInClass[k] + minclass;
	}
	UNPROTECT(1);
	return r;
}

/*
n=length(chr)
.method="dgesdd"
.jobu="S"
.jobv=""
.u <- matrix(0, n, n)
.v <- matrix(0, n, n)
.inv <- matrix(0, n, n)
.d <- double(n)


.d[] <- 0
.v[] <- 0
.u[] <- 0
*/




#define MAXCLASSES_NEARCENT			128
#define MAXVARIABLES_NEARCENT		1024*10
double centroids[MAXCLASSES_NEARCENT * MAXVARIABLES_NEARCENT];
int    samplesInNearClass[MAXCLASSES_NEARCENT];

SEXP galgoNearCent(SEXP data, SEXP samplesClasses, SEXP trainingSamplesVector, SEXP testSamplesVector, SEXP nearMethod, SEXP resultType) {

	if (!isReal(data)) {
		return (galgo_error("data should be double. Use as.double().\n"));
	}
	if (!isInteger(samplesClasses)) {
		return (galgo_error("samplesClasses should be integer. Use as.integer().\n"));
	}
	if (!isInteger(trainingSamplesVector)) {
		return (galgo_error("trainingSamplesVector should be integer. Use as.integer().\n"));
	}
	if (!isInteger(testSamplesVector)  &&  !isNull(testSamplesVector)) {
		return (galgo_error("testSamplesVector should be integer. Use as.integer().\n"));
	}
	if (!isInteger(nearMethod)) {
		return (galgo_error("nearMethod should be integer. Use as.integer().\n"));
	}
	if (!isInteger(resultType)) {
		return (galgo_error("resultType should be integer. Use as.integer().\n"));
	}

	
	int i, j, k, l, maxclass, minclass, nclasses, c, p;
	int nl = LENGTH(samplesClasses);
	int ntr = LENGTH(trainingSamplesVector);
	int nte = LENGTH(testSamplesVector);
	int nvar = ncols(data);
	int nrow = nrows(data);
	int minC;
	double minD;
	double d;
	SEXP r;

	minclass = maxclass = INTEGER(samplesClasses)[0];
	for (i=1; i < nl; i++) {
		j = INTEGER(samplesClasses)[i];
		if (j > maxclass) maxclass = j;
		if (j < minclass) minclass = j;
	}
	nclasses = maxclass - minclass + 1;

	// set centroids to zero
	for (k=nvar*nclasses-1; k >= 0; k--) centroids[k]=0; // centroids rows = variables, cols = classes
	for (k=nclasses-1; k >= 0; k--) samplesInNearClass[k]=0;

	if (INTEGER(nearMethod)) {
		// mean: a) sum to centroids
		for (i=0; i < ntr; i++) {
			l = INTEGER(trainingSamplesVector)[i] - 1;
			k = (INTEGER(samplesClasses)[l] - minclass);
			samplesInNearClass[k]++;
			c = k * nvar;
			for (j=0; j < nvar; j++) centroids[c+j] += REAL(data)[j*nrow + l];
		}
		// mean: b) divide sum by classes
		for (k=nclasses-1; k >= 0; k--) {
			c = k * nvar;
			l = samplesInNearClass[k];
			for (j=0; j < nvar; j++) centroids[c+j] /= l;
		}
	} else {
		// median : medioids
		for (j=0; j < nvar; j++) {
			p = 0;
			for (k=nclasses-1; k >= 0; k--) {
				for (i=0; i < ntr; i++) {
					l = INTEGER(trainingSamplesVector)[i] - 1;
					if ((INTEGER(samplesClasses)[l] - minclass) == k) knn_distances[p++] = REAL(data)[j*nrow + l];
				}
				R_rsort(knn_distances, p);
				if (p % 2 == 1) {
					centroids[k*nvar + j] = knn_distances[(p+1)/2];
				} else {
					centroids[k*nvar + j] = (knn_distances[p/2] + knn_distances[p/2+1])/2;
				}
			}
		}
	}

	if (INTEGER(resultType)[0]) {
		PROTECT(r = allocVector(INTSXP, 1));
		INTEGER(r)[0] = 0;
	} else {
		PROTECT(r = allocVector(INTSXP, nte));
	}
	for (i=0; i < nte; i++) {
		l = INTEGER(testSamplesVector)[i] - 1;
		minC = -1;
		minD = 9e99;
		for (k=nclasses-1; k >= 0; k--) {
			c = k*nvar;
			d = 0;
			for (j=0; j < nvar; j++) d += pow((centroids[c+j] - REAL(data)[j*nrow + l]),2);
			if (d <= minD) {
				minD = d;
				minC = k;
			}
		}
		if (INTEGER(resultType)[0]) {
			if (minC == (INTEGER(samplesClasses)[l] - minclass)) INTEGER(r)[0]++;
		} else {
			INTEGER(r)[i] = minC + minclass;
		}
	}
	UNPROTECT(1);
	return r;
}



static double G_euclidean(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist;
    int count, j;

    count= 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	//if(both_non_NA(x[i1], x[i2])) {
	    dev = (x[i1] - x[i2]);
	//    if(!ISNAN(dev)) {
		dist += dev * dev;
		count++;
	//    }
	//}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return sqrt(dist);
}

static double G_maximum(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist = -DBL_MAX;
    for(j = 0 ; j < nc ; j++) {
	//if(both_non_NA(x[i1], x[i2])) {
	    dev = fabs(x[i1] - x[i2]);
	//    if(!ISNAN(dev)) {
		if(dev > dist)
		    dist = dev;
		count++;
	//    }
	//}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    return dist;
}

static double G_manhattan(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	//if(both_non_NA(x[i1], x[i2])) {
	    dev = fabs(x[i1] - x[i2]);
	//    if(!ISNAN(dev)) {
		dist += dev;
		count++;
	//    }
	//}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return dist;
}

static double G_canberra(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist, sum, diff;
    int count, j;

    count = 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	//if(both_non_NA(x[i1], x[i2])) {
	    sum = fabs(x[i1] + x[i2]);
	    diff = fabs(x[i1] - x[i2]);
	    if (sum > DBL_MIN || diff > DBL_MIN) {
		dev = diff/sum;
//		if(!ISNAN(dev) ||
//		   (!R_FINITE(diff) && diff == sum &&
//		    /* use Inf = lim x -> oo */ (dev = 1.))) {
			if (diff==sum) dev=1.0;
		    dist += dev;
		    count++;
//		}
	    }
	//}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return dist;
}

static double G_dist_binary(double *x, int nr, int nc, int i1, int i2)
{
    int total, count, dist;
    int j;

    total = 0;
    count = 0;
    dist = 0;

    for(j = 0 ; j < nc ; j++) {
	//if(both_non_NA(x[i1], x[i2])) {
	//    if(!both_FINITE(x[i1], x[i2])) {
	//	warning(_("dist(.,\"binary\"): treating non-finite values as NA"));
	//    }
	//    else {
		if(x[i1] || x[i2]) {
		    count++;
		    if( ! (x[i1] && x[i2]) ) dist++;
		}
		total++;
	//    }
	//}
	i1 += nr;
	i2 += nr;
    }

    if(total == 0) return NA_REAL;
    if(count == 0) return 0;
    return (double) dist / count;
}

static double G_minkowski(double *x, int nr, int nc, int i1, int i2, double p)
{
    double dev, dist;
    int count, j;

    count= 0;
    dist = 0;
    for(j = 0 ; j < nc ; j++) {
	//if(both_non_NA(x[i1], x[i2])) {
	    dev = (x[i1] - x[i2]);
	//    if(!ISNAN(dev)) {
		dist += R_pow(fabs(dev), p);
		count++;
	//    }
	//}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    if(count != nc) dist /= ((double)count/nc);
    return R_pow(dist, 1.0/p);
}

enum { EUCLIDEAN=1, MAXIMUM, MANHATTAN, CANBERRA, BINARY, MINKOWSKI };
/* == 1,2,..., defined by order in the R function dist */

void G_distance(double *x, int *nr, int *nc, double *d, int *diag, 
		int *method, double *p)
{
    int dc, i, j, ij;
    double (*distfun)(double*, int, int, int, int) = NULL;

    switch(*method) {
    case EUCLIDEAN:
	distfun = G_euclidean;
	break;
    case MAXIMUM:
	distfun = G_maximum;
	break;
    case MANHATTAN:
	distfun = G_manhattan;
	break;
    case CANBERRA:
	distfun = G_canberra;
	break;
    case BINARY:
	distfun = G_dist_binary;
	break;
    case MINKOWSKI:
	if(!R_FINITE(*p) || *p <= 0)
	    error(_("distance(): invalid p"));
	break;
    default:
	error(_("distance(): invalid distance"));
    }
    dc = (*diag) ? 0 : 1; /* diag=1:  we do the diagonal */
    ij = 0;
	if (*method != MINKOWSKI) {
		for(j = 0 ; j <= *nr ; j++)
		for(i = j+dc ; i < *nr ; i++)
			d[ij++] = distfun(x, *nr, *nc, i, j);
	} else {
		for(j = 0 ; j <= *nr ; j++)
		for(i = j+dc ; i < *nr ; i++)
			d[ij++] = G_minkowski(x, *nr, *nc, i, j, *p);
	}
}

/* La.svd, called from svd */
// I was unable to call La_svd from my C code. Therefore I´m using the La_svd routine from R source.
// NOTE: This code is copyright of R. It was used here to avoid compilation errors.
// We used the R source code as shown in the below address accessed on 27-March-2014
// https://svn.r-project.org/R/trunk/src/modules/lapack/Lapack.c
// We appreaciate the help and comments from R team.

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001--2013  The R Core Team.
 *  Copyright (C) 2003--2010  The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

// static was removed for compilation
/* La.svd, called from svd */
SEXP La_svd(SEXP jobu, SEXP x, SEXP s, SEXP u, SEXP vt)
{
    int n, p, info = 0;

    if (!isString(jobu))
	error("'jobu' must be a character string");
    int *xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0]; p = xdims[1];

    /* work on a copy of x  */
    double *xvals;
    if (!isReal(x)) {
       x = coerceVector(x, REALSXP);
       xvals = REAL(x);
    } else {
	xvals = (double *) R_alloc(n * (size_t) p, sizeof(double));
	Memcpy(xvals, REAL(x), n * (size_t) p);
    }
    PROTECT(x);

    SEXP dims = getAttrib(u, R_DimSymbol);
    if (TYPEOF(dims) != INTSXP) error("non-integer dims");
    int ldu = INTEGER(dims)[0];
    dims = getAttrib(vt, R_DimSymbol);
    if (TYPEOF(dims) != INTSXP) error("non-integer dims");
    int ldvt = INTEGER(dims)[0];
    double tmp;
    /* min(n,p) large is implausible, but cast to be sure */
    int *iwork= (int *) R_alloc(8*(size_t)(n < p ? n : p), sizeof(int));

    /* ask for optimal size of work array */
    const char *ju = CHAR(STRING_ELT(jobu, 0));
    int lwork = -1;
    F77_CALL(dgesdd)(ju, &n, &p, xvals, &n, REAL(s),
		     REAL(u), &ldu, REAL(vt), &ldvt,
		     &tmp, &lwork, iwork, &info);
    if (info != 0)
	error(_("error code %d from Lapack routine '%s'"), info, "dgesdd");
    lwork = (int) tmp;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgesdd)(ju, &n, &p, xvals, &n, REAL(s),
		     REAL(u), &ldu, REAL(vt), &ldvt,
		     work, &lwork, iwork, &info);
    if (info != 0)
	error(_("error code %d from Lapack routine '%s'"), info, "dgesdd");

    SEXP val = PROTECT(allocVector(VECSXP, 3));
    SEXP nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("d"));
    SET_STRING_ELT(nm, 1, mkChar("u"));
    SET_STRING_ELT(nm, 2, mkChar("vt"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, s);
    SET_VECTOR_ELT(val, 1, u);
    SET_VECTOR_ELT(val, 2, vt);
    UNPROTECT(3);
    return val;
}

