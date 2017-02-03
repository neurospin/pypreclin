#include <stdlib.h>
#include <string.h>
/*#include <cmath>*/

#include <math.h>
#include "mex.h"


template <typename T> void interp1_lin(
	mwSize mY, mwSize mXI, mwSize nN, T *YI,
	const T *Y, const T *XI) {
	
	mwSize ix,n,m;
	T xi, dx, yi;
	
	T *YI_local = YI;
	const T *Y_local = Y;
	const T *XI_local = XI;
	
	for(n = 0; n < nN; n++) {
		for(m = 0; m < mXI; m++) {
			xi = XI_local[m];
			
			if(xi >= 1) {			
				if(xi < mY) {
					ix = XI_local[m];
					dx = XI_local[m] - (T)ix;
					yi = Y_local[ix-1];
					YI_local[m] = yi + dx*(Y_local[ix] - yi);
				}
				else {
					if(xi > mY)
						YI_local[m] == NAN;
					else
						YI_local[m] = Y_local[mY-1];
				}
			}
			else
				YI_local[m] = NAN;
		}
		
		Y_local += mY;
		XI_local += mXI;
		YI_local += mXI;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
		const mxArray *prhs[])  {
		
	int isComplex = 0;
	int isDouble = 0;
	mxClassID classID;
	const mxArray *mxY, *mxXI;
	mxArray *mxYI;
	mwSize mY, mXI, nN;
	void *XI;
	
	/* Check for proper number of arguments */
	if (nrhs != 2)
		mexErrMsgTxt("Function requires 2 input arguments:\n\tYI = n_interp1_mc(Y, XI);");
	else if (nlhs > 1)
		mexErrMsgTxt("Function requires one output argument.");
			
	mxY = prhs[0];
	mxXI = prhs[1];
	
	/* check for argument type */
	classID = mxGetClassID(mxY);
	if(classID != mxGetClassID(mxXI))
		mexErrMsgTxt("All parameters have to be of the same type (''double'' or ''single'')");
	if(classID != mxDOUBLE_CLASS & classID != mxSINGLE_CLASS)
		mexErrMsgTxt("''n_interp1_lin'' only supports arguments of type ''double'' or ''single''.");

	/* Input arguments */
	/* Y */
	if(mxGetNumberOfElements(mxY) == 0)
		mexErrMsgTxt("First parameter is empty.");
	isComplex = mxIsComplex(mxY);
	mY = mxGetM(mxY);
	nN = mxGetN(mxY);

	/* XI */
	if(mxIsComplex(mxXI))
		mexErrMsgTxt("Third parameter has to be non-complex.");
	mXI = mxGetM(mxXI);
	if(nN != mxGetN(mxXI))
		mexErrMsgTxt("Both Y & XI have to have the same number of columns.");
	XI = mxGetPr(mxXI);
	
	/* Output argument */
	if(isComplex)
		plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(mxXI),
				mxGetDimensions(mxXI), classID, mxCOMPLEX);
	else
		plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(mxXI),
				mxGetDimensions(mxXI), classID, mxREAL);
	mxYI = plhs[0];

	switch(classID) {
		case mxDOUBLE_CLASS:
			interp1_lin <double> (mY, mXI, nN, (double*)mxGetPr(mxYI),
					(const double*)mxGetPr(mxY), (const double*)XI);
			if(isComplex)
				interp1_lin <double> (mY, mXI, nN, (double*)mxGetPi(mxYI),
						(const double*)mxGetPi(mxY), (const double*)XI);
			break;
		case mxSINGLE_CLASS:
			interp1_lin <float> (mY, mXI, nN, (float*)mxGetData(mxYI),
					(const float*)mxGetData(mxY), (const float*)XI);
			if(isComplex)
				interp1_lin <float> (mY, mXI, nN, (float*)mxGetImagData(mxYI),
						(const float*)mxGetImagData(mxY), (const float*)XI);
			break;
	}

	return;
}
