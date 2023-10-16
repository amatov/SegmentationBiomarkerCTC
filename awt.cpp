#include "mex.h"

% A Trou Wavelet Transform. J.-L. Starck, F. Murtagh, A. Bijaoui, "Image Processing and Data Analysis: The Multiscale Approach (2000)"
	
# Nik Mihaylov (2014)

#include <Eigen>
#include <vector>
#include <cmath>
#include <algorithm>


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;

double log2( double n ){
	return log(n) / log(double(2));
}

Matrix convolve(Matrix I, int k) {

    const int N = I.rows(), M = I.cols();
    const int k1 = 1 << k;
    const int k2 = 1 << (k+1);

    // IMPLEMENTS: tmp = padarray(I, [k2 0], 'replicate');
    Matrix tmp(N + 2*k2, M);
    for( int i = 0; i < k2; ++i ) {
        tmp.row(i) = I.row(0);
        tmp.row(k2 + N + i) = I.row(N-1);
    }
    tmp.middleRows(k2,N) = I;

    // convolve the columns
    for( int i = k2; i < k2 + N; ++i ) {
        I.row(i - k2) = 6*tmp.row(i) + 4*tmp.row(i + k1) + 4*tmp.row(i - k1) + tmp.row(i + k2) + tmp.row(i - k2);
    }

    // IMPLEMENTS: tmp = padarray(I * .0625, [0 k2], 'replicate');
    Matrix J = I*0.0625;
    tmp.resize(N, M + 2*k2);
    for( int i = 0; i < k2; ++i ) {
        tmp.col(i) = J.col(0);
        tmp.col(k2 + M + i) = J.col(M - 1);
    }
    tmp.middleCols(k2, M) = J;

    // convolve the rows
    for( int i = k2; i < k2 + M; ++i ) {
        I.col(i - k2) = 6*tmp.col(i) + 4*tmp.col(i + k1) + 4*tmp.col(i - k1) + tmp.col(i + k2) + tmp.col(i - k2);
    }

    return I * 0.0625;
}


std::vector<Matrix> awt(const Matrix& m, int nBands) {

    const int N = m.rows(), M = m.cols();

    const int K = ceil(std::max(log2(N), log2(M)));

    if( nBands < 1 || nBands > K ) {
        nBands = K;
    }

    std::vector<Matrix> res(nBands + 1);

    Matrix lastA = m;
    for( int k = 0; k < nBands; ++k ) {
        const Matrix newA = convolve(lastA, k);
        res[k] = lastA - newA;
        lastA = newA;
    }
    res[nBands] = lastA;

    return res;
}



#ifdef __cplusplus
extern "C"
#endif
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* Check that the number of input and output parameters is valid */
    if( nrhs != 1 && nrhs != 2) {
        mexErrMsgTxt("One or two input parameters required.");
        return;
    }

    if( nlhs > 1 ) {
        mexErrMsgTxt("One output parameter required.");
        return;
    }

    /* Read input parameter dimensions */
    const int mrows = mxGetM(prhs[0]);
    const int mcols = mxGetN(prhs[0]);

    /* Get a pointer to the input */
    double * const M = mxGetPr(prhs[0]);

    /* Get the number of bands */
    int nBands = 0;
    if( nrhs == 2 ) {
        double * p = mxGetPr(prhs[1]);
        nBands = static_cast<int>(*p);
    }


    // initialize the matrix with the input data
    Matrix I(mrows, mcols);
    const int n = mrows*mcols;
    for( int i = 0; i < n; ++i ) {
        I(i) = M[i];
    }

    std::vector<Matrix> ress = awt(I, nBands);

    // create the 3D array for the result
    const mwSize dims[] = { mrows, mcols, static_cast<mwSize>(ress.size()) };
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);

    // copy the data into the result array
    double * const res = mxGetPr(plhs[0]);
    for( int k = 0; k < ress.size(); ++k ) {
        const Matrix& rm = ress[k];
        const int offset = k * n;
        for( int i = 0; i < n; ++i ) {
            res[offset + i] = rm(i);
        }
    }

}
