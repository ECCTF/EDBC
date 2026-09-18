This code is an implementation of Algorithm 1 and the Bernstein's algorithm for many coefficient set such as C^{\ast}={\pm 1, \pm 5\}$ and $C^{\ast}={\pm 1, \pm 3,\ldots,\pm C_{max}}$ with $1\le C_{max}\le 31$. It is compiled by by Microsoft visual studio 2015 on Win 10. One can generate curve.exe to show a canonic EDBC with $C^{\ast}={\pm 1, \pm 3,\ldots,\pm C_{7}}$. 
The code is shown in main.cpp, AdaptiveDBC.cpp，and AdaptiveDBC.h.

Miracl lib  is used to implement  big number arithmetic. And one should use the correct path of the miral.lib. Our experiments  are compiled  using C++. TheTimeRatioofProduceEDBC.exe shows the ratio of time cost of our Dynamic Programming to Bernstein's algorithm is 1/6 for random 256-bit integers.
curve.exe shows  a canonic EDBC for a positive integer less than 760-bit.   For other coefficient set, one should choose different functions in AdaptiveDBC.h.
For producing different optimal EDBC, one can use the cooresponding value of A,mA,D, and T.
   