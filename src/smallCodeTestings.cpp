#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

class AMatrix
{
	private:
	double *inputData;
    int n_Rows, m_Cols, nm_Elements;
    
    public:
		//  Throughout this code when dealing with a  n*m matrix, n and m are  n:= number of rows and m:= number of columns.
    AMatrix();
    AMatrix(int n, int m);
    AMatrix(int n, int m, const double Data[]);

    //Objects
    AMatrix& set();
	AMatrix& Identity();
    AMatrix& toIdentity(); //converting rectangular matrix to rectangular identity
	AMatrix& zero();
	AMatrix transpose();
    AMatrix norm(); //norm of a given matrix column

    AMatrix operator*(const AMatrix &B);
	AMatrix operator*(const double s);
	AMatrix operator+(const AMatrix &B);
	AMatrix operator-(const AMatrix &B);
	AMatrix &operator=(const AMatrix &B);
    
    double GetElement(int row, int col) const;
    bool SetElement(int row, int col, double elementValue);
    int Sub2Ind(int row, int col) const;
    int GetNumRows() const;
    int GetNumCols() const;

    double norm (int col_index) const;
	double paddingOut (int row_index, int col_index) const;
    void PrintMatrix(int precision);
    void QRdecomposition(AMatrix &Q, AMatrix &R);

    // Destructor
    ~AMatrix();
};


AMatrix::AMatrix()
{
  n_Rows = 1;
  m_Cols = 1;
  nm_Elements = 1;
  inputData = nullptr;
}

// Construct zero matrix (all elements 0)
AMatrix::AMatrix(int n, int m)
{
  n_Rows = n;
  m_Cols = m;
  nm_Elements = n_Rows * m_Cols;
  inputData = new double [nm_Elements];
  for (int i=0; i<nm_Elements; i++)
	  inputData[i] = 0.0;
}

AMatrix::AMatrix(int n, int m, const double Data[])
{
	n_Rows = n;
    m_Cols = m;
    nm_Elements = n_Rows * m_Cols;
	inputData = new double[nm_Elements];
	for (int i=0; i<nm_Elements; i++)
	  inputData[i] = Data[i];
}

AMatrix::~AMatrix() 
{
	// Destructor.
	if (inputData != nullptr)
		delete[] inputData;
	
	inputData = nullptr;
}


AMatrix &AMatrix::zero() 
{
	for (int i = 0; i < n_Rows; i++)
		for (int j = 0; j < m_Cols; j++)
			inputData[j * m_Cols] = 0.0;
	return *this;
}

/*AMatrix AMatrix::norm() 
{
    int j = 0;  // index of column that we want to have its norm
    double sum = 0.0;      
		for (int i = 0; i < n_Rows; i++)
        {
			sum =+ inputData[j + (i * m_Cols)]; //vertical sweeping of the desire column to get the indices
        }  
	return *this;
}*/

AMatrix &AMatrix::Identity() //
{
    int n = this->GetNumRows();
	int m = this->GetNumCols();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        inputData[i * m + j] = (i == j) ? 1.0 : 0.0;
        }
      }
    return *this;
 }

AMatrix &AMatrix::toIdentity() //
{
    int n = GetNumRows();
	int m = GetNumCols();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        inputData[i * m + j] = (i == j) ? 1.0 : 0.0;
        }
      }
    return *this;
 }


double AMatrix::GetElement(int row, int col) const
{
	int linearIndex = Sub2Ind(row, col);
	return inputData[linearIndex];

}

bool AMatrix::SetElement(int row, int col, double elementValue)
{
	int linearIndex = Sub2Ind(row, col);
	return inputData[linearIndex] = elementValue;
}

int AMatrix::Sub2Ind(int row, int col) const
{
		return (row * m_Cols) + col;
}

int AMatrix::GetNumRows() const
{
	return n_Rows;
}

int AMatrix::GetNumCols() const
{
	return m_Cols;
}

// Operators
AMatrix AMatrix::operator*(const AMatrix &B) 
{
	AMatrix R(n_Rows, B.m_Cols);
	for (int j = 0; j < n_Rows; j++)
		for (int i = 0; i < B.m_Cols; i++) 
		{
			R.inputData[j * B.m_Cols + i] = 0.0;
			for (int k = 0; k < m_Cols; k++)
				R.inputData[j * B.m_Cols + i] += inputData[j * m_Cols + k] * B.inputData[k * B.m_Cols + i];
		}
	return R;
}

AMatrix AMatrix::operator*(const double s) 
{
	AMatrix R(n_Rows, m_Cols);
	for (int i = 0; i < n_Rows; i++)
		for (int j = 0; j < m_Cols; j++)
			R.inputData[i * m_Cols + j] = inputData[i * m_Cols + j] * s;
	return R;
}

AMatrix AMatrix::operator+(const AMatrix &B) 
{
	AMatrix R(n_Rows, m_Cols);
	for (int j = 0; j < n_Rows; j++)
		for (int i = 0; i < m_Cols; i++)
			R.inputData[j * m_Cols + i] = inputData[j * m_Cols + i] + B.inputData[j * m_Cols + i];
	return R;
}

AMatrix AMatrix::operator-(const AMatrix &B) 
{
	AMatrix R(n_Rows, m_Cols);
	for (int j = 0; j < n_Rows; j++)
		for (int i = 0; i < m_Cols; i++)
			R.inputData[j * m_Cols + i] = inputData[j * m_Cols + i] - B.inputData[j * m_Cols + i];
	return R;
}

AMatrix &AMatrix::operator=(const AMatrix &B) 
{
	if (n_Rows * m_Cols != B.n_Rows * B.m_Cols) 
	{
		delete [] inputData;
		inputData = new double[B.n_Rows * B.m_Cols];
	}
	n_Rows = B.n_Rows;
	m_Cols = B.m_Cols;

	for (int i = 0; i < n_Rows * m_Cols; i++) inputData[i] = B.inputData[i];
	return *this;
}


double AMatrix::norm (int col_index) const
{
    double sum = 0.0;      
		for (int i = 0; i < n_Rows; i++)
        {
			sum += (inputData[col_index + (i * m_Cols)] * inputData[col_index + (i * m_Cols)]); //vertical sweeping of the desire column to get the indices
        }  
	return sqrt(sum);
}

double AMatrix::paddingOut (int row_index, int col_index) const
{     
	AMatrix paddedA();
	for (int j = 0; j < row_index; j++)
	{
		for (int i = 0; i < col_index; i++)
        {
			paddedA = inputData[j + (i * m_Cols)] ; //vertical sweeping of the desire column to get the indices
        }  
	}
	return sqrt(sum);
}

void AMatrix::PrintMatrix(int precision)
{
	int nRows = GetNumRows();
	int mCols = GetNumCols();
	for (int i = 0; i < nRows; ++i)
    {
	  for (int j = 0; j < mCols; ++j)
        {
            std::cout << std::setw(10)<<std::fixed << std::setprecision(precision) << this->GetElement(i, j) << "  ";
        }
	    std::cout << std::endl;
	}    
}


void AMatrix::QRdecomposition(AMatrix &Q, AMatrix &R) 
{
	double Norm; // For initializing the matrix column Norm
	double alpha; // alpha coefficient
	AMatrix u(m, 1);
    AMatrix v(m, 1);
	AMatrix P(m, m); 
    AMatrix I(m, m);

	Q = AMatrix(m, m);
	R = *this; //making a copy of the input matrix (the current object). It is somehow the modified A through the algorithm to build final R. 
				//R = Hn-1...H3H2H1A       R1=H1A (in the first loop)

    int nRows = R.GetNumRows(); // Getting the number of rows of the input matrix 
	int mCols = R.GetNumCols(); //          ''          columns          ''

	for (int i = 0; i < mCols ; i++) 
	{
		u.zero(); v.zero(); // Initializing u and v vectors
		
		double Norm = 0.0;
		Norm = R.norm(0); // computing the norm of the specified matrix column
		
		alpha = u.A[i] < 0 ? mag : -mag;

		for (int j = i; j < m; j++) v.A[j] /= Norm;

		P = I - (v * v.transpose()) * 2.0;

		R = P * R;
		Q = Q * P;
	}
}

    int main(int argc, char* argv[]) 
{
    int n = 5;
    int m = 4;
    double test[] = {5, 3, 12, 2,    4, 2, 17, 10,    5, 3, 2, 6,   1, 7, 6, 6,   2, -6, 4, -3};
    AMatrix A(n, m, test);
    AMatrix B(n, n);
    B = B.Identity();
    double Norm = A.norm(0);
    A.PrintMatrix(3);
    std::cout << "---------------------------------------------------"<< std::endl;
    B.PrintMatrix(3);
    std::cout << "---------------------------------------------------"<< std::endl;
    std::cout << Norm<< std::endl;
    std::cout << temp[6]<< std::endl;
    //std:: cout<<*temp;
    

	return 0;
}