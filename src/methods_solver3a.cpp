#include "matrix.h"
#include "vector.h"
#include "methods_solver3a.h"

using namespace std;

void thresholdcheck(VectorF &supdiagonalvector)
{

    // zero out any value on the super-diagonal that meets the convergence criteria
    float threshold = 0.001;
    auto set_to_zero = [threshold](float& x) { x = (x < threshold) ? 0 : x; };
    for (int i = 1; i < supdiagonalvector.size(); i++)
    {
        if(supdiagonalvector[i] < threshold)
        {
            set_to_zero(supdiagonalvector[i]);
        }
    }
}


int reversecountzeros(VectorF &countingvector, int start, int end)
{
    int zerocount = 0;
    for (int i = start; i <= end; --i)
    {
        if (countingvector[i] == 0)
        {
            zerocount++;
        }
        else
        {
            break;
        }
    }
    return zerocount;
}

int reversecountnonzeros(VectorF countingvector, int start, int end)
{
    int nonzerocount = 0;
    for (int i = start; i <= end; --i)
    {
        if (countingvector[i] == 0)
        {
            break;
        }
        else
        {
            nonzerocount++;
        }
    }
    return nonzerocount;
}

void calculatesingularvalues(VectorF &diagonalvector, int &B_size)
{
    for (int i = 0; i <= B_size; i++)
    {
        diagonalvector[i] = sqrt(diagonalvector[i]);
    }
}

void segmentandprocess(VectorF &diagonal_vector, VectorF &supdiagonal_vector, int &B_size, int &B2_size, int &B3_size)
{
    VectorF B2_diag(B2_size);
    VectorF B2_supdiag(B2_size - 1);
    VectorF B3_diag(B3_size);
    VectorF B3_supdiag(B3_size - 1);
    VectorF diag_copy(B_size);
    VectorF supdiag_copy(B_size);
    bool shiftsuccess = false;


      // segment the B vectors into 2 vectors, B2/3
    for (int i = 0, j = B_size - B3_size; j <= B_size; i++, j++)
    {
        B3_diag[i] = diagonal_vector[j];
    }
    for (int i = 0, j = B_size - B3_size; j <= B_size - 1; i++, j++)
    {
        B3_supdiag[i] = supdiagonal_vector[j];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size; i++, j++)
    {
        B2_diag[i] = diagonal_vector[j];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size - 1; i++, j++)
    {
        B2_supdiag[i] = supdiagonal_vector[j];
    }
    
    //reverse loop, starting from the end of B2_diag
    for(int i = -1; i >= -1 * B2_size; --i)
    {
        if (B2_diag[i] == 0)
        {
            //apply givens roation so that B2_supdiag[i-1] == 0
            int x = B2_supdiag[i-1];
            int y = B2_diag[i];
            int r = sqrt(pow(x, 2) + pow(y, 2));
            int c = y / r;
            int s = x / r;
            B2_diag[i] = s * B2_supdiag[i-1] + c * B2_diag[i];
            B2_diag[i-1] = c * B2_diag[i-1];
            B2_supdiag[i-1] = c * B2_supdiag[i-1] - s * B2_diag[i]; // should be 0 after rotation
        }
        else

        {
            break;
        }
    }
    //apply shifting strategy logic to diag & supdiag vectors so diag[-1] decreases below tollerence

    //QUESTION: should the shift be on B1/2, B2/3 or the entire B?
    
    double mu = min(*diagonal_vector.begin(), *diagonal_vector.end());
    double d = diagonal_vector[0] - mu;

    for (int i = 0; i <= B_size; i++)
    {
        diag_copy[i] = diagonal_vector[i];
    }

    for (int i = 0; i <= B_size - 1; i++)
    {
        supdiag_copy[i] = supdiagonal_vector[i];
    }

    while ((shiftsuccess = false))
    {
        for (int k = 0; k < B_size - 1; k++)
        {
            double t;

            diagonal_vector[k] = d + supdiagonal_vector[k];
            t = diagonal_vector[k+1] / diagonal_vector[k];
            supdiagonal_vector[k] = supdiagonal_vector[k] * t;
            d = d * t - mu;

            if (d < 0)
            {
                for (int i = 0; i <= B_size; i++)
                {
                    diagonal_vector[i] = diag_copy[i];
                }
                for (int i = 0; i <= B_size - 1; i++)
                {
                    supdiagonal_vector[i] = supdiag_copy[i];
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
            diagonal_vector[-1] = d;
        }
    }


    // inverse the segmenting B logic so all changes are transfered to diagvector and supdiagvector
    for (int i = 0, j = B_size - B3_size; j <= B_size; i++, j++)
    {
        diagonal_vector[j] = B3_diag[i];
    }
    for (int i = 0, j = B_size - B3_size; j <= B_size - 1; i++, j++)
    {
        supdiagonal_vector[j] = B3_supdiag[i];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size; i++, j++)
    {
        diagonal_vector[j] = B2_diag[i];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size - 1; i++, j++)
    {
        supdiagonal_vector[j] = B2_supdiag[i];
    }
}


void householder(MatrixF &A)
{
    int n = A.num_rows();
    int m = A.num_cols();

    for (int k = 0; k < m; ++k)
    {
        int s = 0;
        for (int i = k; i < n; ++i)
        {
            s += A(i, k) * A(i, k);
        }
        int sign = (A(k, k) > 0) ? 1 : -1;
        int mu = sqrt(s) * sign;
        A(k, k) -= mu;
        int norm = 0;
        for (int i = k; i < n; ++i)
        {
            norm += A(i, k) * A(i, k);
        }
        norm = sqrt(norm);
        if (norm == 0)
        {
            continue;
        }
        for (int i = k; i < n; ++i)
        {
            A(i, k) /= norm;
        }
        for (int j = k + 1; j < m; ++j)
        {
            s = 0;
            for (int i = k; i < n; ++i)
            {
                s += A(i, k) * A(i, j);
            }
            for (int i = k; i < n; ++i)
            {
                A(i, j) -= 2 * s * A(i, k);
            }
        }
    }
}


void Solver3_methods(VectorF &diagvector, VectorF &supdiagvector, int &B_size)
{
    int B3_size;
    int B2_size;
    VectorD singularvalues(B_size);

    thresholdcheck(supdiagvector);
    B3_size = reversecountzeros(supdiagvector, -1, -1 * (B_size - 1));

    while (B3_size < B_size)
    {
        B2_size = reversecountnonzeros(supdiagvector, -1 * B3_size, -1 * (B_size - 1));
        segmentandprocess(diagvector, supdiagvector, B_size, B2_size, B3_size);
        thresholdcheck(supdiagvector);
        B3_size = reversecountzeros(supdiagvector, -1, -1 * (B_size - 1));
    }

    //Assert that (B3_size == B_size), otherwise something isn't right...
    calculatesingularvalues(diagvector, B_size);
}

void Solver3_main(MatrixF &matrixA, VectorF &diagV) {
    
    int n = matrixA.num_cols();
    VectorF supdiagV(n);
    VectorD solution(n);
    
    householder(matrixA);
    // assert that matrix is bi-diagonal

    // convert and square bidiagonal matrix into two vectors: diag, and supdiag
    for (int i = 0; i < n; i++) {
      diagV[i] = matrixA(i, i) * matrixA(i, i);
    }
    for (int i = 0; i < n - 1; i++) {
      supdiagV[i] = matrixA(i, i + 1) * matrixA(i, i + 1);
    }

    Solver3_methods(diagV, supdiagV, n);
}