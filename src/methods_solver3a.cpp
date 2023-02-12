#include "matrix.h"
#include "vector.h"
#include "methods_solver3a.h"
#include <algorithm>

using namespace std;

void thresholdcheck(VectorF &supdiagonalvector)
{
#if DEBUG
cout << "starting 'thresholdcheck' on sup diagvector: " << endl;
supdiagonalvector.display();
#endif
    // zero out any value on the super-diagonal that meets the convergence criteria
    float threshold = 0.01;
    for (int i = 0; i < supdiagonalvector.size(); i++)
    {
        if(supdiagonalvector[i] <= threshold)
        {
            supdiagonalvector[i] = 0.0;
        }
    }
#if DEBUG
cout << "finished 'thresholdcheck'. output vector:" << endl;
supdiagonalvector.display();
#endif
}


int reversecountzeros(VectorF &countingvector, int start, int end)
{
#if DEBUG
cout << "starting 'reversecountzeros' on supdiag vector:" << endl;
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

int reversecountnonzeros(VectorF &countingvector, int start, int end)
{
#if DEBUG
cout << "starting 'reversecountnonzeros' on supdiag vector:" << endl;
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

void calculatesingularvalues(VectorF &diagonalvector, int &B_size)
{
#if DEBUG
cout << "begining Singular Value calculation! on vector:" << endl;
diagonalvector.display();
#endif
    for (int i = 0; i < B_size; i++)
    {
        diagonalvector[i] = sqrt(diagonalvector[i]);
    }
#if DEBUG
cout << "Finished Calculating Singular Values!" << endl;
diagonalvector.display();
#endif
}

void segmentandprocess(VectorF &diagonal_vector, VectorF &supdiagonal_vector, int &B_size, int &B2_size, int &B3_size)
{

#if DEBUG
cout << "starting 'segmentandprocess' iteration" << endl;
    //"B_size = " << B_size << "\n"
    //"B2_size = " << B2_size << "\n"
    //"B3_size = " << B3_size << endl;
    //supdiagonal_vector.display();
    //diagonal_vector.display();
#endif

    //VectorF B2_diag(B2_size);
    VectorF B2_supdiag(max(B2_size - 1, 0));
    //VectorF B3_diag(B3_size);
    //VectorF B3_supdiag(max(B3_size - 1, 0));
    VectorF diag_copy(B_size);
    VectorF supdiag_copy(B_size);
    bool shiftsuccess = false;
    VectorF B2_diag_new(B2_size);

#if DEBUG
cout << "segmenting B2 and B3 vectors" << endl;
//B3_diag.display();
//B3_supdiag.display();
//B2_diag_new.display();
//B2_supdiag.display();
//
//B2_diag_new[0] = 33;
//B2_diag_new.display();
//B2_supdiag.display();

#endif
    // segment the B vectors into 2 vectors, B2/3
    //for (int i = 0, j = B_size - B3_size; j < B_size; i++, j++)
    //{
    //    B3_diag[i] = diagonal_vector[j];
    //}
    //for (int i = 0, j = B_size - B3_size; j < B_size - 1; i++, j++)
    //{
    //    B3_supdiag[i] = supdiagonal_vector[j];
    //}
    //diagonal_vector.display();
    //supdiagonal_vector.display();
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size; i++, j++)
    {
        B2_diag_new[i] = diagonal_vector[j];
    }
    //diagonal_vector.display();
    //supdiagonal_vector.display();
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size - 1; i++, j++)
    {
        B2_supdiag[i] = supdiagonal_vector[j];
    }
#if DEBUG
cout << "finished segmenting B2 vectors from B vectors" << endl;
//B2_supdiag.display();
//B2_diag_new.display();
//supdiagonal_vector.display();
#endif
    //reverse loop, starting from the end of B2_diag
    for(int i = -1; i >= -1 * B2_size; --i)
    {
        if (B2_diag_new[i + B2_diag_new.size()] == 0)
        {
#if DEBUG
cout << "applying givens roation";
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

            //// inverse the segmenting B logic so all changes are transfered to diagvector and supdiagvector
            //for (int i = 0, j = B_size - B3_size; j < B_size; i++, j++)
            //{
            //    diagonal_vector[j] = B3_diag[i];
            //}
            //for (int i = 0, j = B_size - B3_size; j < B_size - 1; i++, j++)
            //{
            //    supdiagonal_vector[j] = B3_supdiag[i];
            //}
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
            break;
        }
    }
#if DEBUG
cout << "applying shifting strategy. Starting vectors:" << endl;
supdiagonal_vector.display();
diagonal_vector.display();
#endif
    //apply shifting strategy logic to diag & supdiag vectors so diag[-1] decreases below tollerence
    //QUESTION: should the shift be on B1/2, B2/3 or the entire B?
    //float mu = min_element(*diagonal_vector.begin(), *diagonal_vector.end());
    float minimum = diagonal_vector[0];
    for (int i = 0; i < diagonal_vector.size(); i++)
    {
        float a = diagonal_vector[i];
        if (a < minimum)
        {
            minimum = a;
        }
    }
    float mu = minimum;

    //create vector copies to revert changes if shift fails.
    for (int i = 0; i < B_size; i++)
    {
        diag_copy[i] = diagonal_vector[i];
    }

    for (int i = 0; i < B_size - 1; i++)
    {
        supdiag_copy[i] = supdiagonal_vector[i];
    }
    
    while ((shiftsuccess == false))
    {
        float d = diagonal_vector[0] - mu;
        cout << "shifting strategy parameters: mu = " << mu << " d = " << d << endl;
        for (int k = 0; k < B_size - 1; k++)
        {
            double t = 0.0;

            cout << "shifting strategy iteration k = " << k << " , starting values:\n"
                "diagonal_vector = " << endl;
                diagonal_vector.display();
            cout << "t = " << t << "\n"
                "supdiagonal_vector[k] = " << endl;
                supdiagonal_vector.display();

            diagonal_vector[k] = d + supdiagonal_vector[k];
            t = diagonal_vector[k+1] / diagonal_vector[k];
            //cout << "sup diagonal vector [k] before update: " << supdiagonal_vector[k] << endl;
            supdiagonal_vector[k] = supdiagonal_vector[k] * t;
            //cout << "sup diagonal vector [k] after update: " << supdiagonal_vector[k] << endl;
            d = d * t - mu;
            cout << "shifting strategy finishing values: mu = " << mu << " d = " << d << "\n"
                "diagonal_vector = " << endl;
                diagonal_vector.display();
            cout << "t = " << t << "\n"
                "supdiagonal_vector[k] = " << endl;
                supdiagonal_vector.display();
            if (d < 0)
            {
                cout << "Shift failed! d was < 0 before finishing the iterations. decrease mu and iterate again." << endl;
                for (int i = 0; i < B_size; i++)
                {
                    diagonal_vector[i] = diag_copy[i];
                }
                for (int i = 0; i < B_size - 1; i++)
                {
                    supdiagonal_vector[i] = supdiag_copy[i];
                }
                mu = mu / 2;
                shiftsuccess = false;
                break;
            }
            shiftsuccess = true;
            //cout << "vectors after shift " << k << " :" << endl;
            //diagonal_vector.display();
            //supdiagonal_vector.display();
        }
        // if the shift was successful, set the last diagonal vector element = d
        if(shiftsuccess == true)
        {
            diagonal_vector[diagonal_vector.size() - 1] = d;
        }
    }

#if DEBUG
cout << "completed 'segmentandprocess' iteration. output vectors:" << endl;
supdiagonal_vector.display();
diagonal_vector.display();
#endif
}


void Householder2(MatrixF &A)
{
#if DEBUG
cout << "starting householder transformation1" << endl;
#endif

    int m = A.num_rows();
    int n = A.num_cols();

    float scale;
    float f, g, s, h;
    int k, l, j;
    VectorF hv(n);

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
cout << "completed householder transformation. Output matrix:" << endl;
//A.display();
//cout << endl;
#endif
}





void Solver3_methods(VectorF &diagvector, VectorF &supdiagvector, int &B_size)
{
#if DEBUG
cout << "starting 'solver3_methods'" << endl;
#endif
    int B3_size;
    int B2_size;
    int loop = 0;
    int max_loops = 20;

    thresholdcheck(supdiagvector);
    B3_size = reversecountzeros(supdiagvector, -1, -1 * (B_size - 1));

    while (B3_size < B_size && loop < max_loops)
    {
        B2_size = (0 == reversecountnonzeros(supdiagvector, -1 * B3_size - 1, -1 * (B_size - 1))) ? 0 : reversecountnonzeros(supdiagvector, -1 * B3_size - 1, -1 * (B_size - 1));
        cout << "reverse count non zeros complete. B2_size = " << B2_size << endl;
        segmentandprocess(diagvector, supdiagvector, B_size, B2_size, B3_size);
        thresholdcheck(supdiagvector);
        B3_size = (B_size - 1 == reversecountzeros(supdiagvector, -1, -1 * (B_size - 1))) ? B_size : reversecountzeros(supdiagvector, -1, -1 * (B_size - 1));
        
        cout << "*************loop " << loop << " completed. new B3_size = " << B3_size << endl;
        cout << "*************" << endl;
        cout << "*************" << endl;

        loop ++;
    }

    calculatesingularvalues(diagvector, B_size);
#if DEBUG
cout << "finished 'solver3_methods'" << endl;
#endif
}

void Solver3_main(MatrixF &matrixA, VectorF &diagV)
{
#if DEBUG
cout << "starting method ph, 'Solver3_main'" << endl;
#endif
    int n = min(matrixA.num_cols(), matrixA.num_rows());
    VectorF supdiagV(n - 1);
    
    Householder2(matrixA);
    // assert that matrix is bi-diagonal

    // convert and square bidiagonal matrix into two vectors: diag, and supdiag
    for (int i = 0; i < n; i++) {
      diagV[i] = matrixA(i, i) * matrixA(i, i);
    }
    for (int i = 0; i < n - 1; i++) {
      supdiagV[i] = matrixA(i, i + 1) * matrixA(i, i + 1);
    }
#if DEBUG
cout << "initial diagonal and supdiagonal vectors:" << endl;
diagV.display();
supdiagV.display();
#endif

    Solver3_methods(diagV, supdiagV, n);

#if DEBUG
cout << "completed 'Sovler3_main'" << endl;
#endif

}