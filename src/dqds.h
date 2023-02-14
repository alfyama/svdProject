#include "matrix.h"
#include "vector.h"
#include <algorithm>

using namespace std;

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void thresholdcheck(Vector<T> &supdiagonalvector)
{
#if DEBUG
    cout << "starting 'thresholdcheck' on sup diagvector: ";
    supdiagonalvector.display();
#endif
    // zero out any value on the super-diagonal that meets the convergence criteria
    double threshold = 0.0001;
    for (int i = 0; i < supdiagonalvector.size(); i++)
    {
        if (supdiagonalvector[i] <= threshold)
        {
            supdiagonalvector[i] = 0.0;
        }
    }
#if DEBUG
    cout << "finished 'thresholdcheck'. output vector:";
//    supdiagonalvector.display();
#endif
}

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
int reversecountzeros(Vector<T> &countingvector, int start, int end)
{
#if DEBUG
    cout << "starting 'reversecountzeros' on supdiag vector:";
//    countingvector.display();
#endif
    int zerocount = 0;
    for (int i = start; i >= end; --i)
    {
        // cout << "reverse count zeros start = " << start << " end = " << end << " iteration = " << i << " value = " << countingvector[i + countingvector.size()] << endl;
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

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
int reversecountnonzeros(Vector<T> &countingvector, int start, int end)
{
#if DEBUG
    cout << "starting 'reversecountnonzeros' on supdiag vector:";
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

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void calculatesingularvalues(Vector<T> &diagonalvector, int &B_size)
{
#if DEBUG
    cout << "begining Singular Value calculation! on vector: ";
#endif
    for (int i = 0; i < B_size; i++)
    {
        diagonalvector[i] = sqrt(diagonalvector[i]);
    }
#if DEBUG
    cout << "Finished Calculating Singular Values!";
#endif
}

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void segmentandprocess(Vector<T> &diagonal_vector, Vector<T> &supdiagonal_vector, int &B_size, int &B2_size, int &B3_size)
{

#if DEBUG
    cout << "starting 'segmentandprocess' iteration\n"
#endif

    // Vector<T> B2_diag(B2_size);
    Vector<T> B2_supdiag(max(B2_size - 1, 0));
    Vector<T> B3_diag(B3_size);
    Vector<T> diagonal_vector_copy(B_size);
    Vector<T> supdiagonal_vector_copy(B_size);
    bool shiftsuccess = false;
    Vector<T> B2_diag_new(B2_size);

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


    // reverse loop, starting from the end of B2_diag
    for (int i = -1; i >= -1 * B2_size; --i)
    {
        if (B2_diag_new[i + B2_size] == 0)
        {
#if DEBUG
            cout << ">>>>> applying givens roation to: \n"
#endif
            //cout << "x = " << x << " y = " << y << " r = " << r << " c = " << c << " s = " << s << endl;
            B2_diag_new[i + B2_size] = B2_supdiag[(i + B2_size) - 1];
            B2_diag_new[(i + B2_size) - 1] = -1 * B2_diag_new[(i + B2_size) - 1];
            B2_supdiag[(i + B2_size) - 1] = 0; // should be 0 after rotation
#if DEBUG
            cout << ">>>>> finished givens roation." << endl;
#endif

            // copy all changes to B2_vectors to the original diag & supdiag vectors
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
#endif
    // apply shifting strategy logic to the B2 vectors. so B2_diag[-1] decreases below tollerence

    // create vector copies to revert changes if shift fails.
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
        T d = diagonal_vector[0];
        for (int k = 0; k < B_size - 1; k++)
        {
            double t = 0.0;
            diagonal_vector[k] = d + supdiagonal_vector[k];
            t = diagonal_vector[k + 1] / diagonal_vector[k];
            supdiagonal_vector[k] = supdiagonal_vector[k] * t;
            d = d * t;
#if DEBUG
            cout << ">>>>> shifting strategy finishing" << endl;
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
                shiftsuccess = false;
                break;
            }
            shiftsuccess = true;
        }
        // if the shift was successful, set the last diagonal vector element = d
        if (shiftsuccess == true)
        {
            diagonal_vector[diagonal_vector.size() - 1] = d;
        }
    }

#if DEBUG
    cout << ">>>>> completed 'segmentandprocess' iteration." endl;
#endif
}

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void Householder(Matrix<T> &A, Vector<T> &hv, Vector<T> &w)
{
#if DEBUG
    cout << "starting householder transformation with matrix:" << endl;
#endif

    int m = A.num_rows();
    int n = A.num_cols();

    double scale;
    double f, g, s, h;
    int k, l, j;

    g = 0.0;

    int i;
    for (i = 0; i < n; ++i)
    {
        hv[i] = g;
        scale = 0.0;
        // Transformations for the cols
        if (i < m)
        {
            // loop through rows of column i starting from the top and sum(square(elements in col i))
            for (l = i; l < m; l++)
            {
                scale += A(l, i) * A(l, i);
            }
            // if the sum(square(element)) = 0, so every element was 0, then g=0
            if (scale == 0.0)
                g = 0.0;
            // if the sum(square(element)) was not 0
            else
            {
                f = A(i, i);
                // g with the value of sqrt(scale), but the inverse sign of the diagonal element A(i,i)
                g = -copysign(sqrt(scale), f);
                // h is the diagonal element A(i,i) * norm of col i - sum(square(elements in col i))
                h = f * g - scale;

                // diagonal element A(i,i) = diagonal element A(i,i) - norm with the same sign
                A(i, i) = f - g;

                // sorta looks like matrix multiplication?
                for (j = i + 1; j < n; j++)
                {
                    s = 0.0;
                    for (k = i; k < m; k++)
                    {
                        // increment s with each element in column i * each element in column i+1
                        s += A(k, i) * A(k, j);
                    }
                    // diagonal element A(i,i) = s / (the diagonal element A(i,i) * norm of col i - sum(square(elements in col i)))
                    f = s / h;
                    // for each element in the column i+1, add f * the element in column i
                    for (k = i; k < m; k++)
                    {
                        A(k, j) += f * A(k, i);
                    }
                }
            }
        }
        w[i] = g;

        g = 0.0;
        scale = 0.0;

        // Transformations for the rows
        if (i + 1 != n && i + 1 <= m)
        {
            for (l = i + 1; l < n; l++)
            {
                // sum the sq of every element in row i starting to the right of the diagonal
                scale += A(i, l) * A(i, l);
            }

            if (scale == 0.0)
                g = 0.0;

            else
            {
                f = A(i, i + 1);
                g = -copysign(sqrt(scale), f);
                h = f * g - scale;
                // set the superdiagonal value
                A(i, i + 1) = f - g;

                for (j = i + 1; j < n; j++)
                {
                    // set remaining hv elements = row element / h
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
    cout << "completed householder transformation." << endl;
#endif
}

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void dqds_methods(Vector<T> &diagvector, Vector<T> &supdiagvector, int &B_size)
{
#if DEBUG
    cout << "starting 'dqds_methods'" << endl;
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
    cout << "finished 'dqds_methods'" << endl;
#endif
}

template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
void dqds_main(Matrix<T> &matrixA, Vector<T> &diagV)
{
#if DEBUG
    cout << "starting method ph, 'dqds_main'" << endl;
#endif
    int n = min(matrixA.num_cols(), matrixA.num_rows());
    Vector<T> supdiagV(n - 1);
    Vector<T> householdersupdiag(n);

    Householder(matrixA, householdersupdiag, diagV);
    // square diag and supdiag vectors
    for (int i = 0; i < n; i++)
    {
        diagV[i] = diagV[i] * diagV[i];
    }
    for (int i = 0, j = 1; i < n - 1; i++, j++)
    {
        supdiagV[i] = householdersupdiag[j] * householdersupdiag[j];
    }

#if DEBUG
    cout << "initial diagonal and supdiagonal vectors:" << endl;
    cout << "diagV = \n";
    diagV.display_h();
    cout << "supdiagV = \n";
    supdiagV.display_h();
#endif

    dqds_methods(diagV, supdiagV, n);
    diagV.reorder();

#if DEBUG
    cout << "completed 'dqds_main'" << endl;
#endif
}