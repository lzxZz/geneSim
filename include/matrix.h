#include <iostream>
using namespace std;

namespace Matrix
{
//矩阵计算类
class Matrix
{
  public:
    Matrix(size_t rows, size_t cols) : row_count(rows), col_count(cols)
    {
        value = new double[row_count * col_count];
    }
    Matrix(size_t rows, size_t cols, double init_value) : row_count(rows), col_count(cols)
    {
        value = new double[row_count * col_count];
        for (size_t i = 0; i < row_count; i++)
        {
            for (size_t j = 0; j < col_count; j++)
            {
                set_value(i, j, init_value);
            }
        }
    }

    ~Matrix()
    {
        delete[] value;
    }

    double get_value(size_t x, size_t y)
    {
        if (x > row_count || y > col_count)
        {
            throw out_of_range("greater than matrix dimension!");
        }
        return *(value + x * col_count + y);
    }
    void set_value(size_t x, size_t y, double value)
    {
        if (x > row_count || y > col_count)
        {
            throw out_of_range("greater than matrix dimension!");
        }
        *(this->value + x * col_count + y) = value;
    }

    void multi(double number)
    {
        for (size_t i = 0; i < row_count; i++)
        {
            for (size_t j = 0; j < col_count; j++)
            {
                set_value(i,j,get_value(i,j) * number);
            }
            
        }
    }

    static Matrix getE(int dimension){
        Matrix m(dimension,dimension,0);
        for (int i = 0; i < dimension ; i ++)
        {
            m.set_value(i,i,1);
        }
        return m;
    }

    void print()
    {
        for (size_t i = 0; i < row_count; i++)
        {
            for (size_t j = 0; j < col_count; j++)
            {
                cout << get_value(i, j) << "\t";
            }
            cout << endl;
        }
        cout << endl;
        cout << endl;
    }

  private:
    size_t row_count;
    size_t col_count;
    double *value;
};


using Vector = Matrix;

} // namespace Matrix