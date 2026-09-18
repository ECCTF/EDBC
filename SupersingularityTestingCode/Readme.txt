This code is an implementation of our new product tree supersingularity testing in Sect. 6.2. It is compiled by by Microsoft visual studio 2015 on Win 10. One can generate Supersingularity.exe to show whether the curve is supersingular using our double-base product tree supersingularity testing by input A value of curve. 
The code is shown in main.cpp, doublebasechain.cpp,  doublebasechain.h, preComputation.cpp, preComputation.h, pointOperation.cpp, pointOperation.h, scalarMultiplication.cpp,  scalarMultiplication.h, Supersingularity.cpp, and Supersingularity.h.

Miracl lib  is used to implement  big number arithmetic. One should use the correct path of the miral.lib. Our experiments  are compiled  using C++.
SupersingularityTesting.exe shows  the cost of scalar multiplication  on Edwards curve with $a=-1$ using a combination of projective coordinates and extended coordinates.  
   