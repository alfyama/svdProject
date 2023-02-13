#include <iostream> 
#include <iomanip>
#include <math.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)<(b)?(b):(a))

class Matrix 
{
  private:
  protected:
  public:
	double *data_;
	int num_rows_;
	int num_cols_;

	Matrix(const int num_rows_, const int num_cols_);
	Matrix(const int num_rows_, const int num_cols_, const double data_[]);

	~Matrix();

	Matrix& set(const double data_[]);
	Matrix& Identity();
	Matrix& zero();
	Matrix transpose();
	Matrix operator*(const Matrix &B);
	Matrix operator*(const double s);
	Matrix operator+(const Matrix &B);
	Matrix operator-(const Matrix &B);
	Matrix &operator=(const Matrix &B);
	void jacobiRotationSVD(Matrix &Q, Matrix &B);
	void output();

    // Dimensions methods
    inline int num_rows() const { return num_rows_; }
    inline int num_cols() const { return num_cols_; }
};


void Matrix::output() 
{
	for (int j = 0; j < num_rows_; j++) 
	{
		for (int i = 0; i < num_cols_; i++) 
		{
			std::cout << std::setw(16) << data_[j * num_cols_ + i] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

Matrix::Matrix(const int m, const int n) : num_rows_(m), num_cols_(n) 
{
	data_ = new double[m * n];
	this->Identity();
}

Matrix::Matrix(const int m, const int n, const double data_[]) : num_rows_(m), num_cols_(n) 
{
	this->data_ = new double[m * n];
	for (int i = 0; i < m * n; i++) this->data_[i] = data_[i];
}

Matrix::~Matrix() 
{
	delete [] data_;
}

Matrix &Matrix::set(const double data_[]) 
{
	for (int i = 0; i < num_rows_* num_cols_; i++) this->data_[i] = data_[i];
	return *this;
}

Matrix &Matrix::Identity() 
{
	for (int j = 0; j < num_rows_; j++)
		for (int i = 0; i < num_cols_; i++)
			data_[j * num_cols_ + i] = i == j ? 1.0 : 0.0;
	return *this;
}

Matrix &Matrix::zero() 
{
	for (int j = 0; j < num_rows_; j++)
		for (int i = 0; i < num_cols_; i++)
			data_[j * num_cols_ + i] = 0.0;
	return *this;
}

Matrix Matrix::transpose() 
{
	Matrix R(num_cols_, num_rows_);
	for (int j = 0; j < num_rows_; j++)
		for (int i = 0; i < num_cols_; i++)
			R.data_[i * num_rows_ + j] = data_[j * num_cols_ + i];
	return R;
}


Matrix Matrix::operator*(const Matrix &B) 
{
	Matrix R(num_rows_, B.num_cols_);
	for (int j = 0; j < num_rows_; j++)
		for (int i = 0; i < B.num_cols_; i++) 
		{
			R.data_[j * B.num_cols_ + i] = 0.0;
			for (int k = 0; k < num_cols_; k++)
				R.data_[j * B.num_cols_ + i] += data_[j * num_cols_ + k] * B.data_[k * B.num_cols_ + i];
		}
	return R;
}

Matrix Matrix::operator*(const double s) 
{
	Matrix R(num_rows_, num_cols_);
	for (int j = 0; j < num_rows_; j++)
		for (int i = 0; i < num_cols_; i++)
			R.data_[j * num_cols_ + i] = data_[j * num_cols_ + i] * s;
	return R;
}

Matrix Matrix::operator+(const Matrix &B) 
{
	Matrix R(num_rows_, num_cols_);
	for (int j = 0; j < num_rows_; j++)
		for (int i = 0; i < num_cols_; i++)
			R.data_[j * num_cols_ + i] = data_[j * num_cols_ + i] + B.data_[j * num_cols_ + i];
	return R;
}

Matrix Matrix::operator-(const Matrix &B) 
{
	Matrix R(num_rows_, num_cols_);
	for (int j = 0; j < num_cols_; j++)
		for (int i = 0; i < num_rows_; i++)
			R.data_[j * num_cols_ + i] = data_[j * num_cols_ + i] - B.data_[j * num_cols_ + i];
	return R;
}

Matrix &Matrix::operator=(const Matrix &B) 
{
	if (num_rows_ * num_cols_ != B.num_rows_ * B.num_cols_) 
	{
		delete [] data_;
		data_ = new double[B.num_rows_ * B.num_cols_];
	}
	num_rows_ = B.num_rows_;
	num_cols_ = B.num_cols_;

	for (int i = 0; i < num_rows_ * B.num_cols_; i++) data_[i] = B.data_[i];
	return *this;
}

void Matrix::jacobiRotationSVD(Matrix &Q, Matrix &B) 
{
	
    B = *this;
    Matrix U(num_rows_, num_rows_);
    U = U.Identity();
    Matrix V(num_rows_, num_rows_);
    V = V.Identity();
    
    double epsilon = 1e-12;
	double Norm;
	double alpha;
	Matrix u(num_rows_, 1); 
    Matrix v(num_rows_, 1);
	Matrix H(num_rows_, num_rows_);
    Matrix I(num_rows_, num_rows_);
	Q = Matrix(num_rows_, num_rows_);

	for (int i = 0; i < num_cols_; i++) 
	{
		u.zero(); v.zero();
		
		Norm = 0.0;
		for (int j = i; j < num_rows_; j++) 
		{
			u.data_[j] = B.data_[j * num_cols_ + i];
			Norm += u.data_[j] * u.data_[j];
		}

		Norm = sqrt(Norm);
		
		alpha = u.data_[i] < 0 ? Norm : -Norm;

		Norm = 0.0;
		for (int j = i; j < num_rows_; j++) 
		{
			v.data_[j] = j == i ? u.data_[j] + alpha : u.data_[j];
			Norm += v.data_[j] * v.data_[j];
		}
		Norm = sqrt(Norm);

		if (Norm < epsilon) continue;

		for (int j = i; j < num_rows_; j++) v.data_[j] /= Norm;

		H = I - (v * v.transpose()) * 2.0;

		B = H * B;
		Q = Q * H;
	}
	
    double N;
    double sin;
	double cos;
    double r;
    double B_[4];
    double R_[4];

    for (int i = 0; i < num_cols_; i++) 
        for (int j = 0; j < num_cols_; j++)
           N += B.data_[i * num_cols_ + j] * B.data_[i * num_cols_ + j];
    N = sqrt(N); 

    double s = 0.0;
    while (sqrt(s) <= (epsilon * epsilon) * (N * N)) 
    {
        for (int i = 0; i < num_cols_ - 1; i++) 
        {
            for (int j = i + 1; j < num_cols_; j++)
			{
                s += B.data_[i * num_cols_ + j] * B.data_[j * num_cols_ + i];


                r = sqrt (((B.data_[i * num_cols_ + i]) * (B.data_[i * num_cols_ + i])) + (B.data_[j * num_cols_ + i] * B.data_[j * num_cols_ + i]));
                cos = B.data_[i * num_cols_ + i] / r;
			    sin = B.data_[j * num_cols_ + i] / r;
    // X = B.data_[i * num_cols_ + i]
    // Y = B.data_[j * num_cols_ + i]
			    B_[0] = cos * B.data_[i * num_cols_ + i] + sin * B.data_[j * num_cols_ + i];
			    B_[1] = cos * B.data_[i * num_cols_ + j] + sin * B.data_[j * num_cols_ + j];
                B_[2] = -sin * B.data_[i * num_cols_ + i] + cos * B.data_[j * num_cols_ + i];
                B_[3] = -sin * B.data_[i * num_cols_ + j] + cos * B.data_[j * num_cols_ + j];
    // X = B.data_[i * num_cols_ + i]
    // Y = B.data_[j * num_cols_ + i]   

                cos = B_[0]/sqrt((B_[0] * B_[0]) + (B_[1] * B_[1]));
			    sin = B_[1]/sqrt((B_[0] * B_[0]) + (B_[1] * B_[1]));

			    R_[0] = cos * B_[0] + sin * B_[1];
                R_[1] = -sin * B_[0] + cos * B_[1];
			    R_[2] = cos * B_[2] + sin * B_[3];
                R_[3] = -sin * B_[2] + cos * B_[3];

                B.data_[i * num_cols_ + i] =  R_[0];
                B.data_[j * num_cols_ + j] =  R_[3];

			//U = U * R;
			//V = V * R;
			}
			
        }
        
    }
}

int main(int argc, char* argv[]) 
{
	double test[] = {5, 1, 3, 5,    10, 2, 7, 1,    12, 4, 4, 3,   15, 7, 4, 3,   22, 6, 4, 5};
	Matrix A__(5,4,test);
	Matrix Q__(2,2);
    Matrix B__(2,2);
	A__.jacobiRotationSVD(Q__, B__);

	A__.output();
	Q__.output();
	B__.output();
    (Q__*B__).output();
	

	return 0;
}