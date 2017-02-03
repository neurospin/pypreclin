#include <stdlib.h>
#include <string.h>
/*#include <cmath>*/

#include <math.h>
#include "mex.h"


template <typename T> void rescale_XI(mwSize nX, mwSize nXI,
	mwSize *XIi, T *XId, const T *X, const T *XI) {
	
	mwSize ix, ixi;
	T x1, x2, xi;
	
	for(ixi = 0; ixi < nXI; ixi++) {
		xi = XI[ixi];
		XIi[ixi] = nX;
		XId[ixi] = 0;
		for (ix = 0; ix < nX-1; ix++) {
			x1 = X[ix];
			x2 = X[ix+1];
		
			if(xi >= x1 && xi <= x2) {
				XIi[ixi] = ix;
				XId[ixi] = (xi - x1)/(x2 - x1);
				break;
			}
		}
	}
}

template <typename T> void interp1_lin(
	mwSize nX, mwSize nXI, const mwSize *XIi, const T *XId,
	mwSize nN, const T *Y, T *YI, T fill) {
	
	mwSize i,j;
	unsigned int buff_size = nX*sizeof(T);
	T *Y_buff = (T*)mxMalloc(buff_size + 2*sizeof(T));
	Y_buff[nX] = fill;
	Y_buff[nX+1] = fill;
	
	T *YI_local = YI;
	const T *Y_local = Y;
	
	for(j = 0; j < nN; j++) {
		memcpy((void*)Y_buff, (void*)Y_local, buff_size);
		
		for(i = 0; i < nXI; i++)
			YI_local[i] = Y_buff[XIi[i]] + (Y_buff[XIi[i]+1]-Y_buff[XIi[i]])*XId[i];
		
		YI_local += nXI;
		Y_local += nX;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
		const mxArray *prhs[])  {
		
	double fill = 0;
	double fill_i = 0;
	int isComplex = 0;
	int isDouble = 0;
	int emptyX;
	mxClassID classID;
	const mxArray *mxX, *mxY, *mxXI, *mxFill;
	mxArray *mxYI;
	mwSize ndims, ndims_Y, nX, nXI, i, *dims, *dims_Y, *XIi;
	void *XId;
	
	/* Check for proper number of arguments */
	if (nrhs != 3 && nrhs != 4)
		mexErrMsgTxt("Function requires 3 or 4 input arguments:\n\tYI = n_interp1_lin(X,Y, XI[,fill]);");
	else if (nlhs > 1)
		mexErrMsgTxt("Function requires one output argument.");
			
	mxX = prhs[0];
	mxY = prhs[1];
	mxXI = prhs[2];
	
	/* check for argument type */
	classID = mxGetClassID(mxY);
	if(classID == mxDOUBLE_CLASS)
		isDouble = 1;
	else
		if(classID == mxSINGLE_CLASS)
			isDouble = 0;
		else
			mexErrMsgTxt("''n_interp1_lin'' only supports arguments of type ''double'' or ''single''.");

	if(mxGetClassID(mxX) != mxGetClassID(mxY) ||
		mxGetClassID(mxY) != mxGetClassID(mxXI))
		mexErrMsgTxt("All parameters have to be of the same type (''double'' or ''single'')");

	/* Input arguments */
	/* X */
	if(mxGetNumberOfElements(mxX) == 0)
		mexErrMsgTxt("First parameter is empty.");
	ndims = mxGetNumberOfDimensions(mxX);
	if((ndims > 2 || mxIsComplex(mxX)))
		mexErrMsgTxt("First parameter has to be a non-complex vector.");
	else  {
		dims = (mwSize*)mxGetDimensions(mxX); /* d_psi_1 dims */
		if(dims[0] != 1 && dims[1] != 1) {
			printf("dims = [%d %d]\n",dims[0],dims[1]);
			mexErrMsgTxt("First parameter is not a vector.");
		}
		else
			nX = dims[0]*dims[1];
	}

	/* Y */
	isComplex = mxIsComplex(mxY);
	ndims_Y = mxGetNumberOfDimensions(mxY);
	dims_Y = (mwSize*)mxGetDimensions(mxY); /* d_psi_1 dims */
	if(dims_Y[0] != nX)
		mexErrMsgTxt("Number of lines of Y has to be the same as the number of elements in X for ''n_interp1_lin(X,Y,XI,fill)''");

	/* XI */
	ndims = mxGetNumberOfDimensions(mxXI);
	if((ndims > 2 || mxIsComplex(mxXI)))
		mexErrMsgTxt("Third parameter has to be a non-complex vector.");
	else {
		dims = (mwSize*)mxGetDimensions(mxXI); /* d_psi_1 dims */
		if(dims[0] != 1 && dims[1] != 1)
			mexErrMsgTxt("Third parameter is not a vector.");
		else
			nXI = dims[0]*dims[1];
	}

	/* fill value */
	if(nrhs == 4) {
		mxFill = prhs[3];
		if((mxGetClassID(mxFill) != mxDOUBLE_CLASS) && (mxGetClassID(mxFill) != mxSINGLE_CLASS))
			mexErrMsgTxt("Fill value has to be a ''double'' or  a ''single''");

		if(mxGetClassID(mxFill) == mxDOUBLE_CLASS) {
			fill = ((double*)mxGetData(mxFill))[0];
			if(mxIsComplex(mxFill))
				fill_i = ((double*)mxGetImagData(prhs[3]))[0];
		}
		else {
			fill = (double)((float*)mxGetData(mxFill))[0];
			if(mxIsComplex(mxFill))
				fill_i = (double)((float*)mxGetImagData(mxFill))[0];
		}
	}

	/* Output argument */
	mwSize *dims_YI = (mwSize*)mxCalloc(ndims_Y, sizeof(mwSize));
	dims_YI[0] = nXI;
	for(i = 1; i < ndims_Y; i++) dims_YI[i] = dims_Y[i];
	if(isComplex)
		plhs[0] = mxCreateNumericArray(ndims_Y, dims_YI, classID, mxCOMPLEX);
	else
		plhs[0] = mxCreateNumericArray(ndims_Y, dims_YI, classID, mxREAL);
	mxYI = plhs[0];

	/* interpolation arrays */
	XIi = (mwSize*)mxCalloc(nXI, sizeof(mwSize));
	if(isDouble)
		XId = mxCalloc(nXI, sizeof(double));
	else
		XId = mxCalloc(nXI, sizeof(float));

	if(isDouble)
		rescale_XI <double> (nX, nXI, XIi, (double*)XId,
			(double*)mxGetData(mxX), (double*)mxGetData(mxXI));
	else
		rescale_XI <float> (nX, nXI, XIi, (float*)XId,
			(float*)mxGetData(mxX), (float*)mxGetData(mxXI));

	if(isDouble) {
		interp1_lin <double> (nX, nXI, XIi, (double*)XId, mxGetN(mxY),
			(double*)mxGetPr(mxY), (double*)mxGetPr(mxYI), fill);
		if(isComplex)
			interp1_lin <double> (nX, nXI, XIi, (double*)XId, mxGetN(mxY),
				(double*)mxGetPi(mxY), (double*)mxGetPi(mxYI), fill_i);
	}
	else {
		interp1_lin <float> (nX, nXI, XIi, (float*)XId, mxGetN(mxY),
			(float*)mxGetData(mxY), (float*)mxGetData(mxYI), (float)fill);
		if(isComplex)
			interp1_lin <float> (nX, nXI, XIi, (float*)XId, mxGetN(mxY),
				(float*)mxGetImagData(mxY), (float*)mxGetImagData(mxYI), (float)fill_i);
	}

	mxFree(XIi);
	mxFree(XId);

	return;
}
