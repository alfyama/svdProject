#include "matrix.h"
#include "vector.h"



//input = two vectors s (diagonal of B) and 3 (superdiagonal of B)
void svdsolver1(VectorD diag, VectorD supdiag)
{

    int q;
    int p;

    // square entries
    for (int i = 0; i <= diag.size(); i++)
    {
        diag[i] = diag[i] * diag[i];
    }
    for (int i = 0; i <= supdiag.size(); i++)
    {
        supdiag[i] = supdiag[i] * supdiag[i];
    }

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

        // create B vectors of the right size
        VectorD B2_diag(nonzero_count);
        VectorD B2_supdiag(nonzero_count-1);
        VectorD B3_diag(zero_count);
        VectorD B3_supdiag(zero_count-1);

        // segment the B vectors into 2 vectors, B2/3
        for (int i = 0, j = diag.size() - zero_count + 1; j <= diag.size(); i++, j++)
        {
            B3_diag[i] = diag[j];
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
            VectorD S(diag.size());

            // TODO
            //take sq root of B3_diag and create a matrix 'S' with this new vector on the diagonal
            //matrix 'S' is the Singlualr Value matrix
            for (int i = 0; i <= B3_diag.size(); i++)
            {
                S[i] = sqroot(B3_diag[i]);
            }
            
            break;
        }
        else
        {
            //reverse loop, starting from B2_diag[-1]
            for(int i = -1; i >= -1 * B2_diag.size(); --i)
            {
                if (B2_diag[i] == 0)
                //QUESTION: why does the s need to be zero? for numerical fesibility?
                {
                    //apply givens roation so that B2_supdiag[i-1] == 0
                    int x = B2_supdiag[i-1];
                    int y = B2_diag[i];
                    int r = sqrt(pow(x, 2) + pow(y, 2));
                    int c = y / r;
                    int s = x / r;

                    B2_diag[i] = s * B2_supdiag[i-1] + c * B2_diag[i];
                    B2_diag[i-1] = c * B2_diag[i-1];
                    B2_supdiag[i-1] = c * B2_supdiag[i-1] - s * B2_diag[i]; // should be 0
                }
                else
                {
                    break;
                }
            }

            //TODO:
            //apply shifting strategy logic to diag & supdiag vectors
            // - choose mu as minimum(singular value of B) ** 2
            // - set d = s1 - mu
            // - loop from s1 to sn-1 but if d<0, then choose a new mu
            //      - logic...
            //      - logic...
            // - finally, sn = d
            
            double mu = min(diag);
            double d = diag[0] - mu;
            for (int i = 0; i <= diag.size(); i++)
            {
                VectorD diag_copy(diag.size());
                diag_copy[i] = diag[i];
            }
            for (int i = 0; i <= supdiag.size(); i++)
            {
                VectorD supdiag_copy(supdiag.size());
                supdiag_copy[i] = supdiag[i];
            }

            for (int k = 0; k < diag.size(); k++)
            {
                double t;
                diag[k] = d + supdiag[k];
                t = diag[k+1] / diag[k];
                supdiag[k] = supdiag[k] * t;
                d = d * t - mu;
                if (d < 0)
                {
                    for (int i = 0; i <= diag.size(); i++)
                    {
                        diag[i] = diag_copy[i];
                    }
                    for (int i = 0; i <= supdiag.size(); i++)
                    {
                        supdiag[i] = supdiag_copy[i];
                    }
                    mu = mu / 2;
                    break;
                }
            }

            if (int k == diag.size() - 1)
            {
                diag[k + 1] = d;
                break;
            }
            
        }
    }
}

