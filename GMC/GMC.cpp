/* GMC: input a matrix (X,Y), output a vector of (GMC(X|Y),GMC(Y|X)).
 * Yida Yin
 * yyin44@wisc.edu
 */

// README
// For fast calculation,
// I comment those lines calculates GMC(Y|X).
// In order to calculate this you need to uncomment those lines.

#include <Rcpp.h>
#include <string.h>
#include <math.h> 
#include <algorithm>
#include <assert.h>

#include <iostream>
using namespace std;

void sortByOrder(double a[], const double b[],const int n){
// Sort a[] by the order of b[]
// Bubble Sort!
// TODO if there are two elements in b[] which have the same value.
    assert(n>0);
    double *temp_b = new double[n];
    memcpy(temp_b,b,n*sizeof(double));
    int i, j;
    double temp, temp2;
    for (j = 0; j < n - 1; j++){
        for (i = 0; i < n - 1 - j; i++){
            if(temp_b[i] > temp_b[i + 1]){
                temp = temp_b[i];
                temp_b[i] =temp_b[i + 1];
                temp_b[i + 1] = temp;
                temp2 = a[i];
                a[i] = a[i + 1];
                a[i + 1] = temp2;                
            }
        }
    }
    delete[] temp_b;
}


double CalculateMean(const double value[], int n){
    assert(n>0);
    double sum = 0;
    for(int i = 0; i < n; i++)
        sum += value[i];
    return (sum / n);
}

double CalculateSampleVariane(const double value[], int n){
    assert(n>0);
    double mean = CalculateMean(value, n);
    double temp = 0;
    for(int i = 0; i < n; i++){
        temp += (value[i] - mean) * (value[i] - mean) ;
    }
    return temp / (n - 1);
}
/*
void logitLikeTransformation(double a[], const int n){
// Calsulate 1/(1+abs(a - t(a))),   a = n*n matrix
// Only works on square matrix! 
    assert(n>0);
    double * ta = new double[n*n];
    // Construct t(a)
    // IT IS THE BOTTLENECK OF THIS PROGRAM!
    memcpy(ta, a, n*n*sizeof(double));
    for (int i = 0; i<n; i++) {
        for (int j = i+1; j<n; j++) {
            swap(ta[i*n + j], ta[j*n + i]);
        }
    }
    // Calsulate 1/(1+abs(a - t(a)))
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            a[i*n+j] = 1/(1+abs(a[i*n+j] - ta[i*n+j]));
        }
    }    
    delete[] ta;
}
*/
void logitLikeTransformation(double a[], const int n) {
    // Calsulate 1/(1+abs(a - t(a))),   a = n*n matrix
    // Only works on square matrix! 
    assert(n>0);
    // make a copy of a
    double * ta = new double[n*n];
    memcpy(ta, a, n*n*sizeof(double));
    // Calsulate 1/(1+abs(a - t(a)))
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            a[i*n + j] = 1 / (1 + abs(a[i*n + j] - ta[j*n + i]));
        }
    }
    delete[] ta;
}

double sumOfRowOfMatrix(const double a[],const int n, const int rol, const int colbegin, const int colend){
// Calculate the sum of rol(th) row , from colbegin(th) colunm to colend(th) colunm.
    assert(colbegin>=0 && colbegin<colend && rol >=0 && n>=0 && colend<=n);
    double sum = 0;
    for(int i=colbegin; i<colend; i++){
        sum += a[rol*n+i];
    }
    return sum;
}

double sumOfProductOfMatrix(const double a[],const double b[],const int n, const int rol, const int colbegin, const int colend){
// Calculate the sum of matrix(a[i, colbegin:colend]) %*% vector(b[colbegin:colend])
    assert(colbegin>=0 && colbegin<colend && rol >=0 && n>0 && colend<=n);
    double sum = 0;
    for(int i=colbegin; i<colend; i++){
        sum += a[rol*n+i] * b[i];
    }
    return sum;
}

void vectorABS(double a[], const int n){
    for(int i=0; i<n; i++){
        a[i] = abs(a[i]);
    }
}

double vectorSum(const double a[], const int n){
    double sum=0;
    for(int i=0; i<n; i++){
        sum += a[i];
    }
    return sum;
}

double cov(const double a[], const double b[], const int n){
    double prod = 0;
    for(int i=0; i<n; i++){
        prod += a[i]*b[i];
    }
    return prod/n - (vectorSum(a,n)/n) * (vectorSum(b,n)/n);
}



// [[Rcpp::export]]
double GMCcpp(const Rcpp::NumericMatrix data, const int nlocal=25) {
    const int n = data.nrow(), ln = nlocal;
    assert(data.ncol() == 2 && nlocal>0 && nlocal<n);
    // TODO: check if the mem allocation failed
    double *x = new double[n];
    double *y = new double[n];
    double GMC1 = 0, GMC2 = 0;
    for(int i=0; i<n; i++){
        x[i] = data(i,0);
        y[i] = data(i,1);
    }
    const double varx = CalculateSampleVariane(x,n);
    const double vary = CalculateSampleVariane(y,n);
    double *E_xy = new double[n];
    double *E_yx = new double[n];

    //--------- Create xdata & ydata in another way.  Actually, xdata = [sortByOrder(xdata_y), sorted(y)]
    // xdata_y is the second colunm of xdata; To save memory.  
    double *xdata_y = new double[n];
    double *ydata_x = new double[n];
    memcpy(xdata_y, y, n*sizeof(double));
    memcpy(ydata_x, x ,n*sizeof(double));
    sortByOrder(xdata_y, x, n);
    sortByOrder(ydata_x, y, n);
    sort(x,x+n);  sort(y,y+n);

    // -------- Create matrix X & Y, they are both n*n matrices
    double *XX = new double [n*n];
    double *YY = new double [n*n];
    for(int i=0; i<n; i++){
        memcpy(XX + i*n, x, n*sizeof(double));
        memcpy(YY + i*n, y, n*sizeof(double));
    }
    // Notice: We do not need to do the transpose

    //logitLikeTransformation(XX, n);
    logitLikeTransformation(YY, n);
    // -------- for loop
    for(int i=0; i<n; i++){
        int li = max(0, i-ln);
        int ui = min(n-1, i+ln);

        // Notice: In R, [1:3] means 1,2,3. It includes the last element!
        //E_yx[i] = sumOfProductOfMatrix(XX, xdata_y, n, i, li, ui+1) / sumOfRowOfMatrix(XX, n, i, li, ui+1);
        E_xy[i] = sumOfProductOfMatrix(YY, ydata_x, n, i, li, ui+1) / sumOfRowOfMatrix(YY, n, i, li, ui+1);
    }
    GMC1 = CalculateSampleVariane(E_xy,n) / varx;
    //GMC2 = CalculateSampleVariane(E_yx,n) / vary;
    delete[] x;
    delete[] y;
    delete[] E_xy;
    delete[] E_yx;
    delete[] XX;
    delete[] YY;

    return GMC1;
}

// [[Rcpp::export]]
double GMCLinearcpp(const Rcpp::NumericVector Xb, const Rcpp::NumericVector e, const Rcpp::NumericVector beta, const double lambda1, const double lambda2){
    double GMC1=0;
    if(Xb.size() != e.size()) {cerr<<"ERROR!"<<endl; return 0;}
    const int n = Xb.size();
    double *cXb = new double[n];
    double *ce = new double[n];
    double *cbeta = new double[n];
    for(int i=0; i<n; i++){
        cXb[i] = Xb(i);
        ce[i] = e(i);
    }
    for( int i=0; i<beta.size(); i++){
        cbeta[i] = beta(i);
    }
    double var_cXb = CalculateSampleVariane(cXb,n);
    vectorABS(cbeta,n);



    GMC1 = (-1) * var_cXb / (var_cXb + CalculateSampleVariane(ce,n)) + lambda1 * abs(cov(cXb,ce,n)) + lambda2 * vectorSum(cbeta,beta.size());

    delete[] cXb;
    delete[] ce;
    delete[] cbeta;

    return GMC1;
}
