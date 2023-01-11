#include "mkl.hpp"

void vector()
{
    Casino::MKL::Vector<float> a(10),b(10),c(10);
    std::cout << a.size() << std::endl;
    a.random();
    b.random();
    c = a + b;
    c.print();
    c = a * b;
    c.print();
    c = a - b;
    c.print();
    c = 10.0f * a;
    c.print();
    c = a * 10.0f;
    c.print();
    c = 10.0f + a;
    c.print();
    c = a + 10.0f;
    c.print();
    c = 10.0f - a;
    c.print();
    c = a - 10.0f;
    c.print();
    c = sin(a*M_PI);
    c.print();
    c = cos((float)M_PI*a);
    c.print();
}
int main()
{
    Casino::MKL::Matrix<float> a(3,3),b(1,3),c(3,3);
    a.fill(4);
    c = a * 10.0f;
    c.print();
    c = 1.0f / (1.0f - a;
    c.print();
}