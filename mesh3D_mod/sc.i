%module sc

%{
#include "sc.h"
%}

%include "std_vector.i"
// Templates
namespace std {
    %template(IntVector)  vector<int>;
    %template(DoubleVector) vector<double>;
}  

%include "sc.h"

    


