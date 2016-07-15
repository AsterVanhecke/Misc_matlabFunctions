function checkcompiledBWestimator() 

init_path = pwd ;
code_path = [init_path,'/C_code' ] ;
try
   tmp = mex_getIntSquaredHessian(1, 1, 1, 1) ;
catch
    disp('It appears that the code was not yet compiled for your system. Let us compile it...') ;
    try
        cd(code_path) ; compileBWcalc() ; 
        cd(init_path) ; tmp = mex_getIntSquaredHessian(1, 1, 1, 1) ;
        disp('Your code is properly compiled!') ;        
    catch
        try
            disp('Okay, something is still wrong. Have you set your compiler?') ;
            disp('We will try to set it up. Just chose a Visual studio compliler if you have one, or chose an lcc.') ;
            mex -v ;
            cd(code_path) ; compileBWcalc() ; 
            cd(init_path) ; tmp = mex_getIntSquaredHessian(1, 1, 1, 1) ;
            disp('Yes! Worked!') ;
            disp('Your code is properly compiled!') ;
        catch
            disp('Okay, still not compiled... Is your system 32 bit?') ;  
            try
                cd(code_path) ;                
                mex -outdir ../ -v  mex_getIntSquaredHessian.cpp
                cd(init_path) ; tmp = mex_getIntSquaredHessian(1, 1, 1, 1) ;                
                disp('Your code is properly compiled!') ;
            catch
                disp('I have no clue what is wrong. Exiting demo.') ;
                return ;
            end
        end
    end
end
cd(init_path) ;