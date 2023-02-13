#include "matrix.h"
#include "vector.h"
#include "methods_solver3a.h"
#include <algorithm>

using namespace std;

void thresholdcheck(VectorD &supdiagonalvector)
{
    #if DEBUG
    cout << "starting 'thresholdcheck' on sup diagvector: ";
    supdiagonalvector.display();
    #endif
    // zero out any value on the super-diagonal that meets the convergence criteria
    double threshold = 0.00001;
    for (int i = 0; i < supdiagonalvector.size(); i++)
    {
        if(supdiagonalvector[i] <= threshold)
        {
            supdiagonalvector[i] = 0.0;
        }
    }
    #if DEBUG
    cout << "finished 'thresholdcheck'. output vector:";
    supdiagonalvector.display();
    #endif
}


int reversecountzeros(VectorD &countingvector, int start, int end)
{
    #if DEBUG
    cout << "starting 'reversecountzeros' on supdiag vector:";
    countingvector.display();
    #endif
    int zerocount = 0;
    for (int i = start; i >= end; --i)
    {
        //cout << "reverse count zeros start = " << start << " end = " << end << " iteration = " << i << " value = " << countingvector[i + countingvector.size()] << endl;
        if (countingvector[i + countingvector.size()] == 0.0)
        {
            zerocount++;
        }
        else
        {
            break;
        }
    }
    #if DEBUG
    cout << "finished 'reversecountzeros'. zero count = " << zerocount << endl;
    #endif
    return zerocount;
}

int reversecountnonzeros(VectorD &countingvector, int start, int end)
{
    #if DEBUG
    cout << "starting 'reversecountnonzeros' on supdiag vector:";
    countingvector.display();
    #endif
    int nonzerocount = 0;
    for (int i = start; i >= end; --i)
    {
        if (countingvector[i + countingvector.size()] == 0.0)
        {
            break;
        }
        else
        {
            nonzerocount++;
        }
    }
    #if DEBUG
    cout << "finished 'reversecount non zeros'. non zero count = " << nonzerocount << endl;
    #endif
    return nonzerocount;
}

void calculatesingularvalues(VectorD &diagonalvector, int &B_size)
{
    #if DEBUG
    cout << "begining Singular Value calculation! on vector: ";
    diagonalvector.display();
    #endif
    for (int i = 0; i < B_size; i++)
    {
        diagonalvector[i] = sqrt(diagonalvector[i]);
    }
    #if DEBUG
    cout << "Finished Calculating Singular Values!";
    diagonalvector.display();
    #endif
}

void segmentandprocess(VectorD &diagonal_vector, VectorD &supdiagonal_vector, int &B_size, int &B2_size, int &B3_size)
{

    #if DEBUG
    cout << "starting 'segmentandprocess' iteration\n"
        "B_size = " << B_size << "\n"
        "B2_size = " << B2_size << "\n"
        "B3_size = " << B3_size << "\n"
        "sup diag vector = ";
        supdiagonal_vector.display();
    cout << "diag vector = ";
        diagonal_vector.display();
    #endif

    //VectorD B2_diag(B2_size);
    VectorD B2_supdiag(max(B2_size - 1, 0));
    VectorD B3_diag(B3_size);
    VectorD diagonal_vector_copy(B_size);
    VectorD supdiagonal_vector_copy(B_size);
    bool shiftsuccess = false;
    VectorD B2_diag_new(B2_size);

    #if DEBUG
    cout << ">>>>> segmenting B2 vector" << endl;
    #endif
    // segment the B vectors into 2 vectors, B2/3
    for (int i = 0, j = B_size - B3_size; j < B_size; i++, j++)
    {
        B3_diag[i] = diagonal_vector[j];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size; i++, j++)
    {
        B2_diag_new[i] = diagonal_vector[j];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size - 1; i++, j++)
    {
        B2_supdiag[i] = supdiagonal_vector[j];
    }
    #if DEBUG
    cout << ">>>>> finished segmenting B2 vectors from B vectors" << endl;
    #endif
    //reverse loop, starting from the end of B2_diag
    for(int i = -1; i >= -1 * B2_size; --i)
    {
        if (B2_diag_new[i + B2_size] == 0)
        {
            #if DEBUG
            cout << ">>>>> applying givens roation to: \n"
                    ">>>>> B2_supdiag ";
                    B2_supdiag.display();
            cout << ">>>>> B2_diag ";
                    B2_diag_new.display();
            #endif
            //apply givens roation so that B2_supdiag[i-1] == 0, and update the original diag/supdiag vectors (note, i is backwards indexing)
            cout << "*******************************************************\n"
                    "*******************************************************\n"
                    "*******************************************************\n"
                    "*     givens roation indexing is wrong!!!!!           *\n"
                    "*******************************************************\n"
                    "*******************************************************\n"
                    "*******************************************************\n" << endl;
            int x = B2_supdiag[(i + B2_size) - 1];
            int y = B2_diag_new[i + B2_size];
            int r = sqrt(pow(x, 2) + pow(y, 2));
            int c = y / r;
            int s = x / r;
            B2_diag_new[i + B2_size] = s * B2_supdiag[(i + B2_size) - 1] + c * B2_diag_new[i + B2_size];
            B2_diag_new[(i + B2_size) - 1] = c * B2_diag_new[(i + B2_size) - 1];
            B2_supdiag[(i + B2_size) - 1] = c * B2_supdiag[(i + B2_size) - 1] - s * B2_diag_new[i + B2_size]; // should be 0 after rotation
            #if DEBUG
            cout << ">>>>> finished givens roation. output vectors: \n"
                    ">>>>> B2_supdiag ";
                    B2_supdiag.display();
            cout << ">>>>> B2_diag ";
                    B2_diag_new.display();
            #endif

            //copy all changes to B2_vectors to the original diag & supdiag vectors
            for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size; i++, j++)
            {
            diagonal_vector[j] = B2_diag_new[i];
            }
            for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size - 1; i++, j++)
            {
            supdiagonal_vector[j] = B2_supdiag[i];
            }
        }
        else

        {
            continue;
        }
    }
    #if DEBUG
    cout << ">>>>> applying shifting strategy \n"
            ">>>>> supdiagonal_vector = ";
            supdiagonal_vector.display();
    cout << ">>>>> diagonal_vector = ";
            diagonal_vector.display();
    #endif
    //apply shifting strategy logic to the B2 vectors. so B2_diag[-1] decreases below tollerence

    //set parameter mu = to the minimum singular value (i.e. minimum of B3_diag. but zero if B3_diag is empty)
    //double mu = min_element(*diagonal_vector.begin(), *diagonal_vector.end());
    //double minimum = (0 == B3_size) ? 0 : B3_diag[0];
    double minimum = diagonal_vector[0];
    for (int i = 0; i < B_size; i++)
    {
        double a = diagonal_vector[i];
        if (a < minimum)
        {
            minimum = a;
        }
    }
    double mu = minimum;

    //create vector copies to revert changes if shift fails.
    for (int i = 0; i < B_size; i++)
    {
        diagonal_vector_copy[i] = diagonal_vector[i];
    }

    for (int i = 0; i < B_size - 1; i++)
    {
        supdiagonal_vector_copy[i] = supdiagonal_vector[i];
    }
    
    while ((shiftsuccess == false))
    {
        double d = diagonal_vector[0]; // - mu;
        for (int k = 0; k < B_size - 1; k++)
        {
            double t = 0.0;
            #if DEBUG
            cout << ">>>>> shifting strategy starting iteration " << k << " with values mu = " << mu << " d = " << d << " t = " << t << "\n"
                    ">>>>>      supdiagonal_vector = ";
                    supdiagonal_vector.display();
            cout << ">>>>>      diagonal_vector = ";
                    diagonal_vector.display();
            #endif
            diagonal_vector[k] = d + supdiagonal_vector[k];
            t = diagonal_vector[k+1] / diagonal_vector[k];
            supdiagonal_vector[k] = supdiagonal_vector[k] * t;
            d = d * t; // - mu;
            #if DEBUG
            cout << ">>>>> shifting strategy finishing values: mu = " << mu << " d = " << d << " t = " << t << "\n"
                    ">>>>>      supdiagonal_vector = ";
                    supdiagonal_vector.display();
            cout << ">>>>>      diagonal_vector = ";
                    diagonal_vector.display();
            #endif
            if (d < 0)
            {
                #if DEBUG
                cout << ">>>>> Shift failed! d was < 0 before finishing the iterations. decrease mu and iterate again." << endl;
                #endif
                for (int i = 0; i < B_size; i++)
                {
                    diagonal_vector[i] = diagonal_vector_copy[i];
                }
                for (int i = 0; i < B_size - 1; i++)
                {
                    supdiagonal_vector[i] = supdiagonal_vector_copy[i];
                }
                //mu = mu / 2;
                shiftsuccess = false;
                break;
            }
            shiftsuccess = true;
        }
        // if the shift was successful, set the last diagonal vector element = d
        if(shiftsuccess == true)
        {
            diagonal_vector[diagonal_vector.size() - 1] = d;
        }
    }



    #if DEBUG
    cout << ">>>>> completed 'segmentandprocess' iteration. \n"
            "sup diagonal = ";
            supdiagonal_vector.display();
    cout << "diag vector = ";
            diagonal_vector.display();
    #endif
}


void Householder2(MatrixD &A)
{
    #if DEBUG
    cout << "starting householder transformation with matrix:" << endl;
    A.display();
    #endif

    int m = A.num_rows(); // 4
    int n = A.num_cols(); // 3

    double scale;
    double f, g, s, h;
    int k, l, j;
    VectorD hv(n);

    g = 0.0;

    int i;
    for (i = 0; i < n; ++i) // iterate through all the rows
    {
        hv[i] = g; // hv[0] = 0
        scale = 0.0;
        //Transformations for the cols
        if (i < m) // if we haven't iterated past the number of columns (so the HH is only done on a sq matrix)
        {

            //loop through rows of column i starting from the top and sum(square(elements in col i))
            for (l = i; l < m; l++)
            {
                scale += A(l, i) * A(l, i); //1
            }
            //if the sum(square(element)) = 0, so every element was 0, then g=0
            if (scale == 0.0)
                g = 0.0;
            //if the sum(square(element)) was not 0
            else
            {
                f = A(i, i); //1
                //g with the value of sqrt(scale), but the inverse sign of the diagonal element A(i,i)
                g = -copysign(sqrt(scale), f); //-1
                //h is the diagonal element A(i,i) * norm of col i - sum(square(elements in col i))
                h = f * g - scale; //-2

                //diagonal element A(i,i) = diagonal element A(i,i) - norm with the same sign
                A(i, i) = f - g; //2
                
                //sorta looks like matrix multiplication?
                for (j = i + 1; j < n; j++)//j=1 iterating through columns
                {
                    s = 0.0;
                    for (k = i; k < m; k++)//k=0 iterating through rows
                    {
                        //increment s with each element in column i * each element in column i+1
                        s += A(k, i) * A(k, j);//6+0+0+0
                    }
                    //diagonal element A(i,i) = s / (the diagonal element A(i,i) * norm of col i - sum(square(elements in col i)))
                    f = s / h;//-3
                    //for each element in the column i+1, add f * the element in column i
                    for (k = i; k < m; k++)
                    {
                      A(k, j) += f * A(k, i);//-3*6
                    }
                }
            }
        }

        g = 0.0;
        scale = 0.0;

        // Transformations for the rows
        if (i + 1 != n && i + 1 <= m) // we haven't run out of rows
        {
            for (l = i + 1; l < n; l++)
            {
                //sum the sq of every element in row i starting to the right of the diagonal
                scale += A(i, l) * A(i, l);//36
            }

            if (scale == 0.0)
                g = 0.0;

            else
            {
                f = A(i, i + 1);//6

                g = -copysign(sqrt(scale), f);//-6
                h = f * g - scale;//-72
                //set the superdiagonal value
                A(i, i + 1) = f - g;//0

                for (j = i + 1; j < n; j++)
                {
                    //set remaining hv elements = row element / h
                    hv[j] = A(i, j) / h;//0,0
                }

                // Update the matrix
                for (l = i + 1; l < m; l++)
                {
                    s = 0.0;
                    for (k = i + 1; k < n; k++)
                    {
                        s += A(l, k) * A(i, k);
                    }
                    for (k = i + 1; k < n; k++)
                    {
                        A(l, k) += s * hv[k];
                    }
                }
            }
        }
    }
    #if DEBUG
    cout << "completed householder transformation. Output matrix:" << endl;
    //VectorD householderdiagvv(n);
    //for (int i = 0; i < n; i++) {
    //    householderdiagvv[i] = A(i, i);
    //}
    //householderdiagvv.display_h();
    A.display();
    #endif
}





void Solver3_methods(VectorD &diagvector, VectorD &supdiagvector, int &B_size)
{
    #if DEBUG
    cout << "starting 'solver3_methods'" << endl;
    #endif
    int B3_size;
    int B2_size;
    int loop = 1;
    int max_loops = 10000;

    thresholdcheck(supdiagvector);
    B3_size = (B_size - 1 == reversecountzeros(supdiagvector, -1, -1 * (B_size - 1))) ? B_size : reversecountzeros(supdiagvector, -1, -1 * (B_size - 1));

    while (B3_size < B_size && loop <= max_loops)
    {
        B2_size = (0 == reversecountnonzeros(supdiagvector, -1 * B3_size - 1, -1 * (B_size - 1))) ? 0 : 1 + reversecountnonzeros(supdiagvector, -1 * B3_size - 1, -1 * (B_size - 1));
        segmentandprocess(diagvector, supdiagvector, B_size, B2_size, B3_size);
        thresholdcheck(supdiagvector);
        B3_size = (B_size - 1 == reversecountzeros(supdiagvector, -1, -1 * (B_size - 1))) ? B_size : reversecountzeros(supdiagvector, -1, -1 * (B_size - 1));

        #if DEBUG        
                cout << "************* loop " << loop << " completed. new B3_size = " << B3_size << endl;
                cout << "*************" << endl;
                cout << "*************" << endl;
        #endif

        loop++;
    }

    calculatesingularvalues(diagvector, B_size);
    #if DEBUG
    cout << "finished 'solver3_methods'" << endl;
    #endif
}

void Solver3_main(MatrixD &matrixA, VectorD &diagV)
{
    #if DEBUG
    cout << "starting method ph, 'Solver3_main'" << endl;
    #endif
    int n = min(matrixA.num_cols(), matrixA.num_rows());
    VectorD supdiagV(n - 1);
    
    Householder2(matrixA);

    // convert and square bidiagonal matrix into two vectors: diag, and supdiag
    for (int i = 0; i < n; i++) {
      diagV[i] = (matrixA(i, i) * matrixA(i, i)) / 4;
    }
    for (int i = 0; i < n - 1; i++) {
      supdiagV[i] = (matrixA(i, i + 1) * matrixA(i, i + 1)) / 4;
    }
    #if DEBUG
    cout << "initial diagonal and supdiagonal vectors:" << endl;
    cout << "diagV = \n";
    diagV.display_h();
    cout << "supdiagV = \n";
    supdiagV.display_h();
    #endif

    Solver3_methods(diagV, supdiagV, n);

    #if DEBUG
    cout << "completed 'Sovler3_main'" << endl;
    #endif

}