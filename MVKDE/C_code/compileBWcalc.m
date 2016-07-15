function compileBWcalc()
% compile mex functions for bandwidth

disp('Start compiling...')

%mex -outdir ../ -v  -largeArrayDims mex_getIntSquaredHessian.cpp
mex -outdir ../ -v  -largeArrayDims G:\Anna\matlabFunctions\MVKDE\C_code\mex_getIntSquaredHessian.cpp


disp('Done!')