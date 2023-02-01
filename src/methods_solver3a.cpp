#include "matrix.h"
#include "vector.h"



//input = two vectors s (diagonal of B) and 3 (superdiagonal of B)
void svdsolver1(VectorD diag, VectorD supdiag)
{
    VectorD B2_diag;
    VectorD B2_supdiag;
    VectorD B3_diag;
    VectorD B3_supdiag;
    int q;
    int p;

    // zero out any value on the super-diagonal that meets the convergence criteria
    double threshold = 0.001;
    auto set_to_zero = [threshold](double& x) { x = (x < threshold) ? 0 : x; };
    for (int i = 1; i < supdiag.size(); i++)
    {
        if(supdiag[i] < threshold)
        {
            set_to_zero(supdiag[i]);
        }
    }

    //TODO: initialize B3_diag so we can enter the loop and count 

    // begin loop to segment the B matrix into B1, B2, B3 and perform Given's rotations until B3 is nxn
    while (B3_diag.size() < diag.size())
    {
        // determine p & q to segment B so that B3 is diagonal, B2 has no zero superdiag
        int zero_count = 0;
        // count zeros from end of vector
        for (int i = -1; i >= -1 * supdiag.size(); --i)
        {
            if (supdiag[i] == 0)
            {
                ++zero_count;
            }
            else
            {
                break;
            }
        }

        int nonzero_count = 0;
        // count nonzeros starting from the last zero found in previous loop
        for (int i = -1 - zero_count; i >= -1 * supdiag.size(); --i)
        {
            if (supdiag[i] == 0)
            {
                break;
            }
            else
            {
                ++nonzero_count;
            }
        }

        // segment the B vectors into 2 vectors, B2/3
        for (int i = diag.size() - zero_count + 1; i <= diag.size(); i++)
        {
            B3_diag.push_back(diag[i]);
        }
        for (int i = diag.size() - zero_count + 1; i <= diag.size() - 1; i++)
        {
            B3_supdiag.push_back(diag[i]);
        }
        for (int i = diag.size() - zero_count - nonzero_count + 1; i < diag.size() - zero_count; i++)
        {
            B2_diag.push_back(diag[i]);
        }
        for (int i = diag.size() - zero_count - nonzero_count + 1; i < diag.size() - zero_count - 1; i++)
        {
            B2_supdiag.push_back(diag[i]);
        }

        //
        if (B3_diag.size() == diag.size())
        {
            //take sq root of B3_diag and create a matrix 'S' with this new vector on the diagonal
            //matrix 'S' is the Singlualr Value matrix
            break;
        }
        else
        {
            //reverse loop, starting from B2_diag[-1]
            for(int i = -1; i >= -1 * B2_diag.size(); --i)
            {
                if (B2_diag[i] == 0)
                {
                    //TODO:
                    //apply givens roation so that B2_supdiag[i-1] == 0
                }
                else
                {
                    break;
                }
            }

            //TODO:
            //apply shifting strategy logic
            //so that at least the last element in B2_diag == 0
            
            
            //then 'break;' so that we can re-segment the original matrix B
            //after setting some of the supdiag values to 0 (which means B3 will grow)
            break;
        }
    }
}

