#include "matrix.h"
#include "vector.h"

using namespace std;

void thresholdcheck(VectorD supdiagonalvector)
{

    // zero out any value on the super-diagonal that meets the convergence criteria
    double threshold = 0.001;
    auto set_to_zero = [threshold](double& x) { x = (x < threshold) ? 0 : x; };
    for (int i = 1; i < supdiagonalvector.size(); i++)
    {
        if(supdiagonalvector[i] < threshold)
        {
            set_to_zero(supdiagonalvector[i]);
        }
    }
}


int reversecountzeros(VectorD countingvector, int start, int end)
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

int reversecountnonzeros(VectorD countingvector, int start, int end)
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

VectorD calculatesingularvalues(VectorD diagonalvector, int B_size)
{
    VectorD sigma(B_size);

    for (int i = 0; i <= B_size; i++)
    {
        sigma[i] = sqrt(diagonalvector[i]);
    }

    return sigma;
}

void segmentandprocess(VectorD diagvector, VectorD supdiagvector, int B_size, int B2_size, int B3_size)
{
    VectorD B2_diag(B2_size);
    VectorD B2_supdiag(B2_size - 1);
    VectorD B3_diag(B3_size);
    VectorD B3_supdiag(B3_size - 1);
    VectorD diag_copy(B_size);
    VectorD supdiag_copy(B_size);
    bool shiftsuccess = false;


      // segment the B vectors into 2 vectors, B2/3
    for (int i = 0, j = B_size - B3_size; j <= B_size; i++, j++)
    {
        B3_diag[i] = diagvector[j];
    }
    for (int i = 0, j = B_size - B3_size; j <= B_size - 1; i++, j++)
    {
        B3_supdiag[i] = supdiagvector[j];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size; i++, j++)
    {
        B2_diag[i] = diagvector[j];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size - 1; i++, j++)
    {
        B2_supdiag[i] = supdiagvector[j];
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
    // TODO:
    // add .begin() and .end()
    
    double mu = min(diagvector.begin(), diagvector.end());
    double d = diagvector[0] - mu;

    for (int i = 0; i <= B_size; i++)
    {
        diag_copy[i] = diagvector[i];
    }

    for (int i = 0; i <= B_size - 1; i++)
    {
        supdiag_copy[i] = supdiagvector[i];
    }

    while (shiftsuccess = false)
    {
        for (int k = 0; k < B_size - 1; k++)
        {
            double t;

            diagvector[k] = d + supdiagvector[k];
            t = diagvector[k+1] / diagvector[k];
            supdiagvector[k] = supdiagvector[k] * t;
            d = d * t - mu;

            if (d < 0)
            {
                for (int i = 0; i <= B_size; i++)
                {
                    diagvector[i] = diag_copy[i];
                }
                for (int i = 0; i <= B_size - 1; i++)
                {
                    supdiagvector[i] = supdiag_copy[i];
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
            diagvector[-1] = d;
        }
    }


    // inverse the segmenting B logic so all changes are transfered to diagvector and supdiagvector
    for (int i = 0, j = B_size - B3_size; j <= B_size; i++, j++)
    {
        diagvector[j] = B3_diag[i];
    }
    for (int i = 0, j = B_size - B3_size; j <= B_size - 1; i++, j++)
    {
        supdiagvector[j] = B3_supdiag[i];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size; i++, j++)
    {
        diagvector[j] = B2_diag[i];
    }
    for (int i = 0, j = B_size - B3_size - B2_size; j < B_size - B3_size - 1; i++, j++)
    {
        supdiagvector[j] = B2_supdiag[i];
    }
}


VectorD Solver3_methods(VectorD diagvector, VectorD supdiagvector, int B_size)
{
    int B3_size;
    int B2_size;
    VectorD sigma;

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
    sigma = calculatesingularvalues(diagvector, B_size);
}


VectorD Solver3_main(MatrixD bidiagmatrix)
{
    //assert that matrix is bi-diagonal

    int n = bidiagmatrix.num_cols();
    VectorD diag(n);
    VectorD supdiag(n);
    VectorD solution(n);
    VectorD B3_diag;

    int q;
    int p;
    double mu;
    double d;
    bool shiftsuccess;


    // convert and square bidiagonal matrix into two vectors: diag, and supdiag
    for (int i = 0; i < n; i++)
    {
        diag[i] = bidiagmatrix(i, i) * bidiagmatrix(i, i);
    }
    for (int i = 0; i < n - 1; i++)
    {
        supdiag[i] = bidiagmatrix(i, i + 1) * bidiagmatrix(i, i + 1);
    }

    solution = Solver3_methods(diag, supdiag, n);

    return solution;
}