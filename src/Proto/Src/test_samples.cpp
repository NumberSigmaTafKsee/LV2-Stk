#include "Samples.hpp"
#include "mkl.hpp"

using namespace Casino;


using std::chrono::nanoseconds;
using std::chrono::duration_cast;
typedef std::chrono::high_resolution_clock Clock;

void sampleVector()
{
    sample_vector<float> v(10),a(100000),b(100000);
    fill(v,100.0f);
    fill(a,10.0f);
    fill(b,20.0f);
    // mkl needs to warmup like a gpu when it is warmed up it is fast
    v = (a+b)*b;
    auto start_time = Clock::now();
    //v = Casino::Samples::add(a,b);
    for(int i = 0; i < 100; i++) v = (a+b)*b;
    auto end_time = Clock::now();
    std::cout << "Time difference:"
      << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/100 << " nanoseconds" << std::endl;
    //std::cout << v;
}

void mklVector()
{
    Casino::MKL::Vector<float> v(1000000),a(100000),b(100000);
    v.fill(100.0f);
    a.fill(10.0f);
    b.fill(20.0f);
    // mkl needs to warmup like a gpu when it is warmed up it is fast
    v = (a+b)*b;
    auto start_time = Clock::now();
    //v = Casino::Samples::add(a,b);
    for(int i = 0; i < 100; i++) v = (a+b)*b;
    auto end_time = Clock::now();
    std::cout << "Time difference:"
      << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/100 << " nanoseconds" << std::endl;
    //std::cout << v;
}
int main()
{
    
    mkl_set_env_mode(1);    
    mkl_set_threading_layer(MKL_THREADING_GNU);
    mkl_set_interface_layer(MKL_INTERFACE_LP64+MKL_INTERFACE_GNU);
    mkl_enable_instructions(MKL_ENABLE_AVX2);    
    mkl_cbwr_set(MKL_CBWR_AVX2);    
    mkl_set_num_threads(8);
    mkl_set_num_threads_local(8);
    mkl_set_dynamic(0);
    
    
    
    //mklVector();
    //sampleVector();
    complex_matrix<float> v(3,3);
    v.fill(std::complex<float>(1,0.26));
    std::cout << v*v << std::endl;

    Casino::eigen::complex_matrix<float> v2(3,3);
    v2.fill(std::complex<float>(1,0.26));
    std::cout << v2*v2 << std::endl;
}
