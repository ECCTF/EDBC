This code is an implementation of our new product tree supersingularity testing in Sect. 6.2 and Pope, Reijnders, Robert, and Sferlazza's result where they presented a state-of-the-art
method to compute supersingularity testing by Robert's method [56] that computes
pairing from theMontgomery ladder.. It is compiled by by Microsoft visual studio 2015 on Win 10.  
The code is shown in main.cpp, Supersingularity.cpp, doublebasechain.cpp,  doublebasechain.h, preComputation.cpp, preComputation.h, pointOperation.cpp, pointOperation.h, scalarMultiplication.cpp, and scalarMultiplication.h.

Miracl lib  is used to implement  big number arithmetic. One should use the correct path of the miral.lib. Our experiments  are compiled  using C++.
SupersingularityTesting.exe shows  our supersingularity testing that obtains
17.6% of gain compared with Pope, Reijnders, Robert, and Sferlazza's result.   
   