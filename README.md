# installation:
1. For building the project just run $make$.

2. For running the tests $make tests && runtest$

3. For running the programm $./svd -method=<METHOD> -type=<TYPE> file$

Where METHOD is one of the following:

- gr (Golub Reisnch Method)



# References:
1. William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery - 
Numerical recipes_ the art of scientific computing-Cambridge University Press (2007).pdf

2. https://www.cs.utexas.edu/~inderjit/public_papers/DesignMRRR_toms06.pdf

3. Matrix generator for testing: https://math.nist.gov/MatrixMarket/generators.html

4. Input data (matrices): https://math.nist.gov/MatrixMarket/formats.html
This is the same data format we used in assigment 7 and 8, I suggest to use the same format. Read 3. for more info on how to generate matrices (LAPACK Test Software?). I would ask later on isis what we should use anyway.

5. https://www.cs.utexas.edu/users/inderjit/public_papers/HLA_SVD.pdf

6. READ THE PART ABOUT NUMERICAL APPROACH https://en.wikipedia.org/wiki/Singular_value_decomposition.

# Files:
- matrix.h: Matrix object implementation.
- vector.h: Vector object implementation.

# Methods

1. Golub & Kahan (Alfonso)


