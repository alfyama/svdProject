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
    double threshold = 0.01;
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
    VectorD B2_diag_copy(B_size);
    VectorD B2_supdiag_copy(B_size);
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
        if (B2_diag_new[i + B2_diag_new.size()] == 0)
        {
            #if DEBUG
            cout << ">>>>> applying givens roation";
            #endif
            //apply givens roation so that B2_supdiag[i-1] == 0, and update the original diag/supdiag vectors (note, i is backwards indexing)
            int x = B2_supdiag[i-1];
            int y = B2_diag_new[i];
            int r = sqrt(pow(x, 2) + pow(y, 2));
            int c = y / r;
            int s = x / r;
            B2_diag_new[i] = s * B2_supdiag[i-1] + c * B2_diag_new[i];
            B2_diag_new[i-1] = c * B2_diag_new[i-1];
            B2_supdiag[i-1] = c * B2_supdiag[i-1] - s * B2_diag_new[i]; // should be 0 after rotation

        }
        else

        {
            break;
        }
    }
    #if DEBUG
    cout << ">>>>> applying shifting strategy \n"
            ">>>>> B2_supdiag = ";
            B2_supdiag.display();
    cout << ">>>>> B2_diag = ";
            B2_diag_new.display();
    #endif
    //apply shifting strategy logic to the B2 vectors. so B2_diag[-1] decreases below tollerence

    //set parameter mu = to the minimum singular value (i.e. minimum of B3_diag. but zero if B3_diag is empty)
    //double mu = min_element(*diagonal_vector.begin(), *diagonal_vector.end());
    double minimum = (0 == B3_size) ? 0 : B3_diag[0];
    for (int i = 0; i < B3_size; i++)
    {
        double a = B3_diag[i];
        if (a < minimum)
        {
            minimum = a;
        }
    }
    double mu = minimum;

    //create vector copies to revert changes if shift fails.
    for (int i = 0; i < B2_size; i++)
    {
        B2_diag_copy[i] = B2_diag_new[i];
    }

    for (int i = 0; i < B2_size - 1; i++)
    {
        B2_supdiag_copy[i] = B2_supdiag[i];
    }
    
    while ((shiftsuccess == false))
    {
        double d = B2_diag_new[0] - mu;
        for (int k = 0; k < B2_size - 1; k++)
        {
            double t = 0.0;
            #if DEBUG
            cout << ">>>>> shifting strategy starting iteration " << k << " with values mu = " << mu << " d = " << d << " t = " << t << "\n"
                    ">>>>>      B2_supdiag = ";
                    B2_supdiag.display();
            cout << ">>>>>      B2_diag_new = ";
                    B2_diag_new.display();
            #endif
            B2_diag_new[k] = d + B2_supdiag[k];
            t = B2_diag_new[k+1] / B2_diag_new[k];
            B2_supdiag[k] = B2_supdiag[k] * t;
            d = d * t - mu;
            #if DEBUG
            cout << ">>>>> shifting strategy finishing values: mu = " << mu << " d = " << d << " t = " << t << "\n"
                    ">>>>>      B2_supdiag = ";
                    B2_supdiag.display();
            cout << ">>>>>      B2_diag_new = ";
                    B2_diag_new.display();
            #endif
            if (d < 0)
            {
                #if DEBUG
                cout << ">>>>> Shift failed! d was < 0 before finishing the iterations. decrease mu and iterate again." << endl;
                #endif
                for (int i = 0; i < B2_size; i++)
                {
                    B2_diag_new[i] = B2_diag_copy[i];
                }
                for (int i = 0; i < B2_size - 1; i++)
                {
                    B2_supdiag[i] = B2_supdiag_copy[i];
                }
                mu = mu / 2;
                shiftsuccess = false;
                break;
            }
            shiftsuccess = true;
        }
        // if the shift was successful, set the last diagonal vector element = d
        if(shiftsuccess == true)
        {
            B2_diag_new[B2_diag_new.size() - 1] = d;
        }
    }

    //copy all changes to B2_vectors to the original diag & supdiag vectors
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size; i++, j++)
    {
        diagonal_vector[j] = B2_diag_new[i];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size - 1; i++, j++)
    {
        supdiagonal_vector[j] = B2_supdiag[i];
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
    cout << "starting householder transformation1" << endl;
    #endif

    int m = A.num_rows();
    int n = A.num_cols();

    double scale;
    double f, g, s, h;
    int k, l, j;
    VectorD hv(n);

    g = 0.0;

    int i;
    for (i = 0; i < n; ++i)
    {
        hv[i] = g;
        scale = 0.0;
        //Transformations for the cols
        if (i < m)
        {

            //loop through rows of column i starting from the top and sum(square(elements in col i))
            for (l = i; l < m; l++)
            {
                scale += A(l, i) * A(l, i);
            }
            //if the sum(square(element)) = 0, so every element was 0, then g=0
            if (scale == 0.0)
                g = 0.0;
            //if the sum(square(element)) was not 0
            else
            {
                f = A(i, i);
                //g with the value of sqry(scale), but the sign of the diagonal element A(i,i)
                g = -copysign(sqrt(scale), f);
                //h is the diagonal element A(i,i) * norm of col i - sum(square(elements in col i))
                h = f * g - scale;

                //diagonal element A(i,i) = diagonal element A(i,i) - norm with the same sign
                A(i, i) = f - g;

                for (j = i + 1; j < n; j++)
                {
                    s = 0.0;
                    for (k = i; k < m; k++)
                    {
                        //increment s with each element in column i * each element in column i+1
                        s += A(k, i) * A(k, j);
                    }
                    //diagonal element A(i,i) = s / (the diagonal element A(i,i) * norm of col i - sum(square(elements in col i)))
                    f = s / h;
                    //for each element in the column i+1, add f * the element in column i
                    for (k = i; k < m; k++)
                    {
                      A(k, j) += f * A(k, i);
                    }
                }
            }
        }

        g = 0.0;
        scale = 0.0;

        // Transformations for the rows
        if (i + 1 != n && i + 1 <= m)
        {
            for (l = i + 1; l < n; l++)
            {
              scale += A(i, l) * A(i, l);
            }

            if (scale == 0.0)
                g = 0.0;

            else
            {
                f = A(i, i + 1);

                g = -copysign(sqrt(scale), f);
                h = f * g - scale;
                A(i, i + 1) = f - g;

                for (j = i + 1; j < n; j++)
                {
                    hv[j] = A(i, j) / h;
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
    //cout << "completed householder transformation. Output matrix:" << endl;
    //VectorD householderdiagvv(n);
    //for (int i = 0; i < n; i++) {
    //    householderdiagvv[i] = A(i, i);
    //}
    //householderdiagvv.display_h();
    //A.display();
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
      diagV[i] = matrixA(i, i) * matrixA(i, i);
    }
    for (int i = 0; i < n - 1; i++) {
      supdiagV[i] = matrixA(i, i + 1) * matrixA(i, i + 1);
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