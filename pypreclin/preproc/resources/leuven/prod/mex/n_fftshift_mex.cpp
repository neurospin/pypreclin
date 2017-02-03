#include <stdlib.h>
#include <string.h>
/*#include <cmath>*/

#include <math.h>
#include "mex.h"

// // template <typename T>
// // void fftshift1(const mxArray *mxX, mxArray *mxY, bool isComplex) {
// // 	mwSize M = mxGetM(mxX);
// // 	mwSize N = mxGetN(mxX);
// // 	
// // 	if(M == 1) {
// // 		memcpy(mxGetData(mxY), mxGetData(mxX), M*N*sizeof(T));
// // 		if(isComplex)
// // 			memcpy(mxGetImagData(mxY), mxGetImagData(mxX), M*N*sizeof(T));
// // 	}
// // 	else {
// // 		mwSize shift = (M+1)/2;
// // 		mwSize i;
// // 		
// // 		size_t shiftSize1 = (M-shift)*sizeof(T);
// // 		size_t shiftSize2 = shift*sizeof(T);
// // 		
// // 		T *x2 = (T*)mxGetData(mxX);
// // 		T *x1 = x2 + shift;
// // 		T *y1 = (T*)mxGetData(mxY);
// // 		T *y2 = y1 + M-shift;
// // 		
// // 		for(i = 0; i < N; i++) {
// // 			memcpy(y1, x1, shiftSize1);
// // 			memcpy(y2, x2, shiftSize2);
// // 			
// // 			x1 += M; x2 += M;
// // 			y1 += M; y2 += M;
// // 		}
// // 		if(isComplex) {
// // 			x2 = (T*)mxGetImagData(mxX);
// // 			x1 = x2 + shift;
// // 			y1 = (T*)mxGetImagData(mxY);
// // 			y2 = y1 + M-shift;
// // 			
// // 			for(i = 0; i < N; i++) {
// // 				memcpy(y1, x1, shiftSize1);
// // 				memcpy(y2, x2, shiftSize2);
// // 				x1 += M; x2 += M;
// // 				y1 += M; y2 += M;
// // 			}
// // 		}
// // 	}
// // }
// // 
// // template <typename T>
// // void fftshift2(const mxArray *mxX, mxArray *mxY, bool isComplex) {
// // 	
// // 	mwSize ndims = mxGetNumberOfDimensions(mxX);
// // 	const mwSize *dims = mxGetDimensions(mxX);
// // 	
// // 	mwSize M = dims[1];
// // 	
// // 	if(M == 1) {
// // 		memcpy(mxGetData(mxY), mxGetData(mxX),
// // 				mxGetM(mxX)*mxGetN(mxX)*sizeof(T));
// // 		if(isComplex)
// // 			memcpy(mxGetImagData(mxY), mxGetImagData(mxX),
// // 					mxGetM(mxX)*mxGetN(mxX)*sizeof(T));
// // 	}
// // 	else {
// // 		mwSize i;
// // 		mwSize P = dims[0];
// // 		mwSize Q = 1; for(i = 2; i < ndims; i++) Q *= dims[i];
// // 		
// // 		mwSize shift = (M+1)/2;
// // 		
// // 		size_t shiftSize1 = P*(M-shift)*sizeof(T);
// // 		size_t shiftSize2 = P*shift*sizeof(T);
// // 		
// // 		T *x2 = (T*)mxGetData(mxX);
// // 		T *x1 = x2 + P*shift;
// // 		T *y1 = (T*)mxGetData(mxY);
// // 		T *y2 = y1 + P*(M-shift);
// // 		
// // 		for(i = 0; i < Q; i++) {
// // 			memcpy(y1, x1, shiftSize1);
// // 			memcpy(y2, x2, shiftSize2);
// // 			
// // 			x1 += P*M; x2 += P*M;
// // 			y1 += P*M; y2 += P*M;
// // 		}
// // 		if(isComplex) {
// // 			x2 = (T*)mxGetImagData(mxX);
// // 			x1 = x2 + P*shift;
// // 			y1 = (T*)mxGetImagData(mxY);
// // 			y2 = y1 + P*(M-shift);
// // 			
// // 			for(i = 0; i < Q; i++) {
// // 				memcpy(y1, x1, shiftSize1);
// // 				memcpy(y2, x2, shiftSize2);
// // 				x1 += P*M; x2 += P*M;
// // 				y1 += P*M; y2 += P*M;
// // 			}
// // 		}
// // 	}
// // }

template <typename T>
void fftshiftn(const mxArray *mxX, mxArray *mxY, bool isComplex, mwSize dim) {
	
	mwSize ndims = mxGetNumberOfDimensions(mxX);
	const mwSize *dims = mxGetDimensions(mxX);
	
	if((dim > ndims) || dim < 1 || (dims[dim-1] == 1)) {
		memcpy(mxGetData(mxY), mxGetData(mxX),
				mxGetM(mxX)*mxGetN(mxX)*sizeof(T));
		if(isComplex)
			memcpy(mxGetImagData(mxY), mxGetImagData(mxX),
					mxGetM(mxX)*mxGetN(mxX)*sizeof(T));
	}
	else {
		mwSize i;
		
		mwSize M = dims[dim-1];
		mwSize P = 1; for(i = 0; i < dim-1; i++) P *= dims[i];
		mwSize Q = 1; for(i = dim; i < ndims; i++) Q *= dims[i];
		
		mwSize shift = (M+1)/2;
		
		size_t shiftSize1 = P*(M-shift)*sizeof(T);
		size_t shiftSize2 = P*shift*sizeof(T);
		
		T *x2 = (T*)mxGetData(mxX);
		T *x1 = x2 + P*shift;
		T *y1 = (T*)mxGetData(mxY);
		T *y2 = y1 + P*(M-shift);
		
		for(i = 0; i < Q; i++) {
			memcpy(y1, x1, shiftSize1);
			memcpy(y2, x2, shiftSize2);
			
			x1 += P*M; x2 += P*M;
			y1 += P*M; y2 += P*M;
		}
		if(isComplex) {
			x2 = (T*)mxGetImagData(mxX);
			x1 = x2 + P*shift;
			y1 = (T*)mxGetImagData(mxY);
			y2 = y1 + P*(M-shift);
			
			for(i = 0; i < Q; i++) {
				memcpy(y1, x1, shiftSize1);
				memcpy(y2, x2, shiftSize2);
				x1 += P*M; x2 += P*M;
				y1 += P*M; y2 += P*M;
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
		const mxArray *prhs[])  {
	
	const mxArray *mxX = prhs[0];
	mxClassID classID = mxGetClassID(mxX);
	bool isComplex = mxIsComplex(mxX);
	mwSize dim = 0;
	mwSize i;
	
	if(classID != mxDOUBLE_CLASS && classID != mxSINGLE_CLASS)
		mexErrMsgTxt("'n_fftshift1_mex(x)' only accepts 'double' or 'single' data type.");
	
	
	plhs[0] = mxDuplicateArray(mxX);
	mxArray *mxY  = plhs[0];
	mxArray *mxTmp, *mxT1, *mxT2;
	
	if(nrhs >= 2) {
		if(mxIsNumeric(prhs[1]))
			switch(mxGetClassID(prhs[1])) {
				case mxDOUBLE_CLASS:
					dim = (mwSize)(((double*)mxGetData(prhs[1]))[0]);
					break;
				case mxSINGLE_CLASS:
					dim = (mwSize)(((float*)mxGetData(prhs[1]))[0]);
					break;
				case mxINT8_CLASS:
					dim = (mwSize)(((char*)mxGetData(prhs[1]))[0]);
					break;
				case mxUINT8_CLASS:
					dim = (mwSize)(((unsigned char*)mxGetData(prhs[1]))[0]);
					break;
				case mxINT16_CLASS:
					dim = (mwSize)(((short*)mxGetData(prhs[1]))[0]);
					break;
				case mxUINT16_CLASS:
					dim = (mwSize)(((unsigned short*)mxGetData(prhs[1]))[0]);
					break;
				case mxINT32_CLASS:
					dim = (mwSize)(((int*)mxGetData(prhs[1]))[0]);
					break;
				case mxUINT32_CLASS:
					dim = (mwSize)(((unsigned int*)mxGetData(prhs[1]))[0]);
					break;
				case mxINT64_CLASS:
					dim = (mwSize)(((long*)mxGetData(prhs[1]))[0]);
					break;
				case mxUINT64_CLASS:
					dim = (mwSize)(((unsigned long*)mxGetData(prhs[1]))[0]);
					break;
				default:
					dim = -1;
			}
			dim = dim < -1 ? -1 : dim;
	}

	const mwSize ndims = mxGetNumberOfDimensions(mxX);
	switch(dim) {
		case -1:
			break;
		case 0:
			mxTmp = mxDuplicateArray(mxX);
			if(ndims%2 == 0) {
				mxT1 = mxTmp;
				mxT2 = mxY;
			}
			else {
				mxT1 = mxY;
				mxT2 = mxTmp;				
			}
			if(classID == mxDOUBLE_CLASS) {
				fftshiftn <double> (mxX, mxT1, isComplex, 1);
				for(dim = 2; dim <= ndims; dim++)
					if(dim%2 == 1)
						fftshiftn <double> (mxT2, mxT1, isComplex, dim);
					else
						fftshiftn <double> (mxT1, mxT2, isComplex, dim);
			}
			else {
				fftshiftn <float> (mxX, mxT1, isComplex, 1);
				for(dim = 2; dim <= ndims; dim++)
					if(dim%2 == 1)
						fftshiftn <float> (mxT2, mxT1, isComplex, dim);
					else
						fftshiftn <float> (mxT1, mxT2, isComplex, dim);
			}
			mxDestroyArray(mxTmp);
			break;
// 		case 1:
// 			if(classID == mxDOUBLE_CLASS)
// 				fftshift1 <double> (mxX, mxY, isComplex);
// 			else
// 				fftshift1 <float> (mxX, mxY, isComplex);
// 			break;
// 		case 2:
// 			if(classID == mxDOUBLE_CLASS)
// 				fftshift2 <double> (mxX, mxY, isComplex);
// 			else
// 				fftshift2 <float> (mxX, mxY, isComplex);
// 			break;
		default:
			if(classID == mxDOUBLE_CLASS)
				fftshiftn <double> (mxX, mxY, isComplex, dim);
			else
				fftshiftn <float> (mxX, mxY, isComplex, dim);
			break;
	}
	
	return;
}
