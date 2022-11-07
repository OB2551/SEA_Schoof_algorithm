# Elliptic Curves

This code is an implementation of the Schoof and SEA algorithm for counting points on elliptic curves over finite fields. Python was chosen for ease of implementation at the cost of performance. SEA starts slowing significantly down past 128-bit primes. Anything larger than 32 bits for Schoof and the Baby-step giant-step algorithm is very slow and generally infeasible.

EllipticCurve.py  contains an EllipticCurve class which requires a prime p, and integers A and B to initialise an instance of the class. 
An EllipticCurve object has the methods schoof() and sea() to count the number of points on the elliptic curve ![Untitled3](https://user-images.githubusercontent.com/67613774/200216559-d157c584-0b8c-45fa-8f9f-dea281c6dd05.png)


EllipticCurve.py  has dependencies on the other files in this repository, except for the GUI.py and .png file.

Other packages used are in the implementation are standard python packages: itertools, math, random, functools and operator.
The gmpy2 package is also used to help speed up arithmetic and will need to be installed.
Example.py contains a few example computations.

The GUI.py file provides a basic gui to use the SEA or Schoof algorithm for point counting, where the user can specify the curve parameters p,A,B.
customtkinter and tkinter need to be installed for GUI.py to run.
Alternatively it can be downloaded as a standalone application here:

https://file.io/CRsqJYctm0vU

The large file size is largely due to text files for precomputed modular polynomials, which have enormous coefficients.
The application opens with the terminal window, as it is possible there are bugs that I am unaware of. The SEA algorithm is intended for large primes - for very small primes or very rare edge cases it can potentially fail.

While using large curves, the app window may change to "Not repsonding" while the computation is still runnning.

![Untitled2](https://user-images.githubusercontent.com/67613774/200215308-df3f9062-c35b-41a2-a2cf-74ab13810e4e.png)

![Untitled](https://user-images.githubusercontent.com/67613774/200215314-c0a96679-cd48-4821-910c-4c162d25ce79.png)

There are several options available that alter the SEA algorithm, by default it is set to the options that generally give the best performance. 

![output1](https://user-images.githubusercontent.com/67613774/200216582-bdd250a8-a30e-48bc-a0f9-df46d303a83f.jpg)

Point counting is of interest in elliptic curve cryptography.
For an introduction to point counting on elliptic curves, see for example:
https://en.wikipedia.org/wiki/Counting_points_on_elliptic_curves
