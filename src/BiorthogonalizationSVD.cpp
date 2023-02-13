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
	void BiorthogonalizationSVD(Matrix &U, Matrix &V);
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

void Matrix::BiorthogonalizationSVD(Matrix &U, Matrix &V) 
{
	double N;
    double sin;
	double cos;
    double epsilon = 0.001;
	double sum[4];
	//double G[4];
	double U_[4];
    Matrix R(num_rows_, num_rows_); 
    Matrix I(num_rows_, num_rows_);
    //Matrix V(num_rows_, num_rows_);
    Matrix sigma(num_rows_, 1);
    I = I.Identity();
    V = V.Identity();
	U = *this;

    for (int i = 0; i < num_cols_; i++) 
        for (int j = 0; j < num_cols_; j++)
            N += U.data_[i * num_cols_ + j] * U.data_[i * num_cols_ + j];
    N = sqrt(N);        
        
    double s = 0.0;
    while (sqrt(s) <= (epsilon * epsilon) * (N * N)) 
    {
        for (int i = 0; i < num_cols_ - 1; i++) 
        {
            for (int j = i + 1; j < num_cols_; j++)
			{
                for (int k = 0; k < num_rows_; k++)
                {
                    s += U.data_[k * num_cols_ + i] * U.data_[k * num_cols_ + j];
                    sum[0] += U.data_[k * num_cols_ + i] * U.data_[k * num_cols_ + i];
					sum[1] += U.data_[k * num_cols_ + i] * U.data_[k * num_cols_ + j];
					sum[2] += U.data_[k * num_cols_ + j] * U.data_[k * num_cols_ + i];
					sum[3] += U.data_[k * num_cols_ + j] * U.data_[k * num_cols_ + j];
                }
			cos = sum[0]/sqrt((sum[0] * sum[0]) + (sum[2] * sum[2]));
			sin = sum[2]/sqrt((sum[0] * sum[0]) + (sum[2] * sum[2]));
			
			U_[0] = cos * sum[0] + sin * sum[2];
			U_[1] = cos * sum[1] + sin * sum[3];
			//U_[2] =	-sin * sum[0] + cos * sum[2];
			//U_[3] = -sin * sum[1] + cos * sum[3];
			cos = U_[0]/sqrt((U_[0] * U_[0]) + (U_[1] * U_[1]));
			sin = U_[1]/sqrt((U_[0] * U_[0]) + (U_[1] * U_[1]));
			//G[0] = cos * U_[0] + sin * U_[1];
			//G[1] = -sin * U_[0] + cos * U_[1];

			U.data_[i * num_cols_ + i] = cos * U.data_[i * num_cols_ + i] + sin * U.data_[i * num_cols_ + j];

            U.data_[i * num_cols_ + j] = -sin * U.data_[i * num_cols_ + i] + cos * U.data_[i * num_cols_ + j];

			U = U * R;
			//V = V * R;
			}
			
        }
    }

    for (int i = 0; i < num_cols_; i++) 
    {
        for (int k = 0; k < num_rows_; k++) 
        {
        sigma.data_[i] += U.data_[k * num_cols_ + i] * U.data_[k * num_cols_ + i];
        }
    sigma.data_[i] = sqrt(sigma.data_[i]);
    }
}


int main(int argc, char* argv[]) 
{
	double test[] = {2, 1, 3, 5,   1, 0, 7, 1,    0, 3, 4, 3,   3, 7, 4, 3,  1, 6, 4, 3};
	Matrix A_(5,4,test);
	Matrix U_(2,2), V_(2,2);
	A_.BiorthogonalizationSVD(U_, V_);

	A_.output();
	U_.output();
	//sigma.output();
	//V_.output();
	//(U_*R__).output();

	return 0;
}
