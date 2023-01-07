%{
#include<armadillo>
#include<vector>
#include<complex>
#include<string>
%}

// this is for script language indexing matrix like [i][j]
// if you can think of a better way to do it
%inline %{
template<typename T>
struct MatrixView {
    arma::Mat<T> * m;
    size_t row;

    MatrixView(arma::Mat<T> * x, size_t r) : row(r),m(x) {}

    T& operator [](size_t i) { return (*m)(row,i); }

    T    __getitem__(size_t i) { return (*m)(row,i); }
    void __setitem__(size_t i, const T val) { (*m)(row,i) = val; }
};

// will return to the armadillo way soon
enum ApproxType
{
    APPROX_ABSDIFF=0,
    APPROX_RELDIFF,
    APPROX_BOTH
};

%}

namespace arma 
{
    template<typename T>
    class Col {        
    public:
        Col();
        Col(size_t n);
        Col(const Col<T> & v);        
        ~Col();

        // this pretty useless for python but seems to work in Lua
        Col<T>& operator = (const Col<T> & v);
        Col<T>& operator = (const arma::Col<T> & v);
        
        T& at(size_t r);

        T& operator()(size_t i);
        T& operator[](size_t i);

        %extend {
            T    __getitem__(size_t i) { return (*$self)(i); }
            void __setitem__(size_t i, T v) { (*$self)(i) = v; }

        Col<T> operator - () { return -*$self;}
        Col<T> operator + (const Col<T> & b) { return *$self + b; }
        Col<T> operator - (const Col<T> & b) { return *$self - b; }
        Col<T> operator * (const Col<T> & b) { return *$self * b; }
        Col<T> operator / (const Col<T> & b) { return *$self / b; }
/*
        Col<T> operator % (const Col<T> & b) {
            Col<T> r(*this);
            r = vector % b;
            return r;
        }
*/

        Col<T> operator + (const T & b) { return *$self + b; }
        Col<T> operator - (const T & b) { return *$self - b; }
        Col<T> operator * (const T & b) { return *$self * b; }
        Col<T> operator / (const T & b) { return *$self / b; }

            // print is a python keyword
            void println(const std::string & header = "") { $self->print(header); }

            T*   data() { return $self->memptr(); }
        }

        size_t size() const;
        void zeros();
        void zeros(size_t r);
        void ones();
        void ones(size_t r);
        //void clean(const T v);
        void randu();
        void randu(size_t r);
        void randn();
        void randn(size_t r);
        void fill(const T val);
        void replace(const T oldv, const T newv);
        void clamp(const T min, const T max);
        void set_size(size_t r);
        void resize(size_t n);
        void copy_size(const Col<T> & m);
        void reset();      
        void raw_print(const std::string & header = "");
        void brief_print(const std::string & header = "");
        T * memptr();
    };

    template<typename T>
    class Row {
    public:
        Row();
        Row(size_t n);
        Row(const Row<T> & v);
        ~Row();

        Row<T>& operator = (const Row<T> & v);
        Row<T>& operator = (const arma::Row<T> & v);
        
        T& at(size_t r);
        size_t size() const;
        T& operator()(size_t i);
        T& operator[](size_t i);

        %extend {
            T    __getitem__(size_t i) { return (*$self)(i); }            
            void __setitem__(size_t i, T v) { (*$self)(i) = v; }
            void println(const std::string & header = "") { $self->print(header); }
            T*   data() { return $self->memptr(); }

            Row<T> operator - () { return -*$self; }
            Row<T> operator + (const Row<T> & b) { return *$self + b; }
            Row<T> operator - (const Row<T> & b) { return *$self - b; }
            Row<T> operator * (const Row<T> & b) { return *$self * b; }
            Row<T> operator / (const Row<T> & b) { return *$self / b; }
    /*
            Row<T> operator % (const Row<T> & b) {
                Row<T> r(*this);
                r = vector % b;
                return r;
            }
    */
            
            Row<T> operator + (const T & b) { return *$self + b; }
            Row<T> operator - (const T & b) { return *$self - b; }
            Row<T> operator * (const T & b) { return *$self * b; }
            Row<T> operator / (const T & b) { return *$self / b; }

        }

        void zeros();
        void zeros(size_t r);
        void ones();
        void ones(size_t r);
        void randu();
        void randu(size_t r);
        void randn();
        void randn(size_t r);
        void fill(const T val);
        //void clean(const T v);
        void replace(const T oldv, const T newv);
        void clamp(const T min, const T max);
        void set_size(size_t r);
        void resize(size_t n);
        void copy_size(const Row<T> & m);
        void reset();
        T*   memptr();
        void raw_print(const std::string & header = "");
        void brief_print(const std::string & header = "");
    };

    template<typename T>
    class Mat {
    public:
        Mat();
        Mat(size_t rows, size_t cols);
        Mat(const Mat<T> & m);        
        //Mat(const std::initializer_list<T> & l) { matrix = arma::Mat<T>(l); }
        ~Mat();

        Mat<T>& operator = (const Mat<T> & m);
                

        T& operator()(size_t i, size_t j);
        
        
        T& at(size_t r, size_t c);

        void zeros();
        void zeros(size_t r, size_t c);
        void ones();
        void ones(size_t r, size_t c);
        void eye();
        void eye(size_t r, size_t c);
        void randu();
        void randu(size_t r, size_t c);
        void randn();
        void randn(size_t r, size_t c);
        void fill(const T val);
        //void clean( T v);
        void replace(const T oldv, const T newv);
        void clamp(const T min, const T max);
        void set_size(size_t r, size_t c);
        void reshape(size_t r, size_t c);
        void resize(size_t r, size_t c);
        void copy_size(const Mat<T> & m);
        void reset();
        T*   memptr();
        T*   colptr(size_t n);

        Col<T> as_col();
        Row<T> as_row();

        Mat<T> t();
        //Mat<T> i() { return Mat<T>(matrix.i()); }

        T min();
        T max();

        size_t index_min();
        size_t index_max();
        
        Mat<T> eval();

        bool is_empty();
        bool is_trimatu();
        bool is_trimatl();
        bool is_diagmat();
        bool is_square();
        bool is_symmetric();
        bool is_zero();
        bool is_finite();
        bool has_inf();
        bool has_nan();

        void raw_print(const std::string & header = "");
        void brief_print(const std::string & header = "");

        Col<T> diag(size_t k=0);
        //void diag(const Col<T> & v) { matrix.diag(v); }
        
        
        void insert_rows(size_t r, size_t n);
        void insert_cols(size_t c, size_t n);
        
        void shed_row(size_t n);
        void shed_rows(size_t f, size_t l);

        void shed_col(size_t n);
        void shed_cols(size_t f, size_t l);
        
        void swap_rows(size_t r1, size_t r2);
        void swap_cols(size_t c1, size_t c2);

        void swap(Mat<T> & x);

        Mat<T> submat(size_t r1, size_t c1, size_t r2, size_t c2);
        

        void save(const std::string& filename);
        void load(const std::string& filename);
        
        //void save_csv(const std::string& filename, const std::string& header = "") { matrix.save(arma::csv_name(filename,header)); }
        //void load_csv(const std::string& filename, const std::string& header = "") { matrix.load(arma::csv_name(filename,header)); }

        %extend {
            MatrixView<T> __getitem__(size_t row) { return MatrixView<T>($self,row); }

            // not sure why it is not in the class?
            Mat<T> operator - () { return -*$self; }    
            Mat<T> operator + (const Mat<T> & b) { return *$self + b; }
            Mat<T> operator - (const Mat<T> & b) { return *$self - b; }
            Mat<T> operator * (const Mat<T> & b) { return *$self * b; }
            Mat<T> operator / (const Mat<T> & b) { return *$self / b; }
    /*
            Mat<T> operator % (const Mat<T> & b) {
                Mat<T> r(*this);
                r = matrix % b;
                return r;
            }
    */
            Mat<T> operator + (const T & b) { return *$self + b; }
            Mat<T> operator - (const T & b) { return *$self - b; }
            Mat<T> operator * (const T & b) { return *$self * b; }
            Mat<T> operator / (const T & b) { return *$self / b; }

            void __setitem__(size_t col, const Col<T> & val) { $self->col(col) = val; }
            void __setitem__(size_t row, const Row<T> & val) { $self->row(row) = val; }
        
            void copy(const Mat<T> & m) { *$self = m; }

            size_t rows() const { return $self->n_rows; }
            size_t cols() const { return $self->n_cols; }
            
            void set_row(size_t r, const Row<T> & v) {
                $self->row(r) = v;
            }
            void set_col(size_t c, const Col<T> & v) {
                $self->col(c) = v;
            }
            Col<T>    get_col(size_t c) { return Col<T>($self->col(c)); }
            Row<T>    get_row(size_t r) { return Row<T>($self->row(r)); }

            void println(const std::string & header = "") { $self->print(header); }

            void addToEachRow(const Row<T> & v) { $self->each_row() += v; }
            void addToEachCol(const Col<T> & v) { $self->each_col() += v; }

            T*   data() { return $self->memptr(); }
        }
    };

    class SizeMat
    {
    public:
    
    const uword n_rows;
    const uword n_cols;
    
    inline explicit SizeMat(const uword in_n_rows, const uword in_n_cols);
    
    inline uword operator[](const uword dim) const;
    inline uword operator()(const uword dim) const;
    
    inline bool operator==(const SizeMat& s) const;
    inline bool operator!=(const SizeMat& s) const;
    
    inline SizeMat operator+(const SizeMat& s) const;
    inline SizeMat operator-(const SizeMat& s) const;
    
    inline SizeMat operator+(const uword val) const;
    inline SizeMat operator-(const uword val) const;
    
    inline SizeMat operator*(const uword val) const;
    inline SizeMat operator/(const uword val) const;
    };



    class arma_rng
    {
    public:
        
        typedef arma_rng_cxx11::seed_type seed_type;
        static void set_seed(const seed_type val);
        static void set_seed_random();    
    };
}

%inline %{
    namespace ops {
    template<typename T>    
    Col<T> operator + (const T a, const Col<T> & b) {
        Col<T> r(b.size());
        r = a + b;
        return r;
    }
    template<typename T>    
    Col<T> operator - (const T a, const Col<T> & b) {
        Col<T> r(b.size());
        r = a - b;
        return r;
    }
    template<typename T>    
    Col<T> operator * (const T a, const Col<T> & b) {
        Col<T> r(b.size());
        r = a * b;
        return r;
    }
    template<typename T>    
    Col<T> operator / (const T a, const Col<T> & b) {
        Col<T> r(b.size());
        r = a * b;
        return r;
    }

    template<typename T>    
    Mat<T> operator + (const T a, const Mat<T> & b) {
        Mat<T> r(b.size());
        r = a + b;
        return r;
    }
    template<typename T>    
    Mat<T> operator - (const T a, const Mat<T> & b) {
        Mat<T> r(b.size());
        r = a - b;
        return r;
    }
    template<typename T>    
    Mat<T> operator * (const T a, const Mat<T> & b) {
        Mat<T> r(b.size());
        r = a * b;
        return r;
    }

    template<typename T>    
    Col<T> linspace(T start, T end, size_t N=100) { return arma::linspace(start,end,N); }

    template<typename T>    
    Col<T> logspace(T start, T end, size_t N=100) { return arma::logspace(start,end,N); }

    template<typename T>    
    Col<T> regspace(T start, T end, size_t N=100) { return arma::regspace(start,end,N); }


    template<typename T>
    Col<T> randperm(size_t n) { return Col<T>(arma::randperm(n)); }

    template<typename T>
    Mat<T> randperm(size_t r, size_t c) { return Mat<T>(arma::randperm(r,c)); }

    template<typename T>
    Mat<T> eye(size_t r, size_t c) { return Mat<T>(arma::eye(r,c)); }

    template<typename T>
    T zero(const T & x) {
        T r(x);
        r.zeros();
        return r;
    }
    template<typename T>
    T ones(const T & x) {
        T r(x);
        r.ones();
        return r;
    }
    template<typename T>
    T random(const T & x) {
        T r(x);
        r.randu();
        return r;
    }
    template<typename T>
    T randu(const T & x) {
        T r(x);
        r.randu();
        return r;
    }


    template<typename T>
    Col<T> random_vector(size_t n, T min=(T)0, T max = (T)0) {
        Col<T> r(n);
        r.randu();
        r = min + (max-min)*r;
        return r;
    }

    template<typename T>
    Mat<T> toeplitz(const Col<T> & a) { return Mat<T>(arma::toeplitz(a)); }

    template<typename T>
    Mat<T> circ_toeplitz(const Col<T> & a) { return Mat<T>(arma::circ_toeplitz(a)); }

    template<typename T>
    Mat<T> toeplitz(const Col<T> & a, const Col<T> & b) { return Mat<T>(arma::toeplitz(a,b)); }


    template<typename T>
    T accu(const Col<T> & a) { return arma::accu(a); }
    template<typename T>
    T accu(const Row<T> & a) { return arma::accu(a); }
    template<typename T>
    T accu(const Mat<T> & a) { return arma::accu(a); }

    template<typename T>
    Col<T> arg(const Col<T> & a) { return arma::arg(a); }
    template<typename T>
    Row<T> arg(const Row<T> & a) { return arma::arg(a); }
    template<typename T>
    Mat<T> arg(const Mat<T> & a) { return arma::arg(a); }

    template<typename T>
    T sum(const Col<T> & a) { return arma::sum(a); }
    template<typename T>
    T sum(const Row<T> & a) { return arma::sum(a); }
    template<typename T>
    Col<T> row_sum(const Mat<T> & a) { return arma::sum(a,0); }
    template<typename T>
    Row<T> col_sum(const Mat<T> & a) { return arma::sum(a,1); }

    template<typename T>
    Mat<T> affmul(const Mat<T> & a, const Mat<T> & b) { return Mat<T>(arma::affmul(a,b)); }
    template<typename T>
    Mat<T> affmul(const Mat<T> & a, const Col<T> & b) { return Mat<T>(arma::affmul(a,b)); }

    template<typename T>
    std::complex<T> cdot(const Col<std::complex<T>> & a, const Col<std::complex<T>> & b) { return arma::cdot(a,b); }
    template<typename T>
    std::complex<T> cdot(const Row<std::complex<T>> & a, const Col<std::complex<T>> & b) { return arma::cdot(a,b); }

    template<typename T>
    bool approx_equal(const Mat<T> & a, const Mat<T> & b, const std::string& method, T tol) {
        std::string foo;        
        return arma::approx_equal(a,b, method.c_str(), tol);
    }

    template<typename T>
    double as_scalar(Col<T> expr) { return arma::as_scalar(expr); }
    template<typename T>
    double as_scalar(Row<T> expr) { return arma::as_scalar(expr); }

    template<typename T>
    Mat<T> clamp(const Mat<T> & a, T min, T max) { return Mat<T>(arma::clamp(a,min,max)); }

    template<typename T>
    T cond(const Mat<T> & a) { return arma::cond(a); }

    template<typename T>
    Col<T> cumsum(const Col<T> & v) { return Col<T>(arma::cumsum(v)); }
    template<typename T>
    Row<T> cumsum(const Row<T> & v) { return Col<T>(arma::cumsum(v)); }
    template<typename T>
    Mat<T> cumsum(const Mat<T> & m) { return Mat<T>(arma::cumsum(m)); }

    template<typename T>
    Col<T> cumprod(const Col<T> & v) { return Col<T>(arma::cumprod(v)); }
    template<typename T>
    Row<T> cumprod(const Row<T> & v) { return Col<T>(arma::cumprod(v)); }
    template<typename T>
    Mat<T> cumprod(const Mat<T> & m) { return Mat<T>(arma::cumprod(m)); }

    template<typename T>
    T det(const Mat<T> & m) { return arma::det(m); }

    template<typename T>
    Mat<T> diagmat(const Col<T> & v) { return Mat<T>(arma::diagmat(v)); }
    template<typename T>
    Mat<T> diagmat(const Mat<T> & m) { return Mat<T>(arma::diagmat(m)); }

    template<typename T>
    Col<T> diagvec(const Mat<T> & a, size_t k = 0) { return Col<T>(arma::diagvec(a)); }

    template<typename T>
    Col<T> diff(const Col<T> & v, size_t k=1) { return Col<T>(arma::diff(v,k)); }

    template<typename T>
    Mat<T> diff(const Mat<T> & m, size_t k=1, size_t dim = 0) { return Mat<T>(arma::diff(m,k,dim)); }

    template<typename T>
    T dot(const Col<T> & a, const Col<T> & b) { return arma::dot(a,b); }
    template<typename T>
    T dot(const Row<T> & a, const Row<T> & b) { return arma::dot(a,b); }

    template<typename T>
    T norm_dot(const Col<T> & a, const Col<T> & b) { return arma::norm_dot(a,b); }
    template<typename T>
    T norm_dot(const Row<T> & a, const Row<T> & b) { return arma::norm_dot(a,b); }

    template<typename T>
    void inplace_trans(Mat<T> & m) { arma::inplace_trans(m); }

    template<typename T>
    T min(const Col<T> & x) { return arma::min(x); }
    template<typename T>
    T min(const Row<T> & x) { return arma::min(x); }
    
    template<typename T>
    Col<T> col_min(const Mat<T> & x,size_t dim=0) { return arma::min(x,dim); }
    template<typename T>
    Row<T> row_min(const Mat<T> & x,size_t dim=0) { return arma::min(x,dim); }

    template<typename T>
    T max(const Col<T> & x) { return arma::max(x); }
    template<typename T>
    T max(const Row<T> & x) { return arma::max(x); }

    template<typename T>
    Row<T> row_max(const Mat<T> & x, size_t dim=0) { return arma::max(x,dim); }

    template<typename T>
    Col<T> col_max(const Mat<T> & x, size_t dim=0) { return arma::max(x,dim); }

    template<typename T>
    Col<T> normalize(const Col<T> & x) { return Col<T>(arma::normalise(x)); }
    template<typename T>
    Row<T> normalize(const Row<T> & x) { return Col<T>(arma::normalise(x)); }
    template<typename T>
    Mat<T> normalize(const Mat<T> & x) { return Mat<T>(arma::normalise(x)); }

    template<typename T>
    T prod(const Col<T> & x) { return arma::prod(x); }
    template<typename T>
    T prod(const Row<T> & x) { return arma::prod(x); }
    
    template<typename T>
    Col<T> col_prod(const Mat<T> & x, size_t dim=0) { return arma::prod(x,dim); }
    template<typename T>
    Row<T> row_prod(const Mat<T> & x, size_t dim=0) { return arma::prod(x,dim); }

    template<typename T>
    Mat<std::complex<T>> powmat(const Mat<T> & x, T n=0) { return Mat<std::complex<T>>(arma::powmat(x,n)); }

    template<typename T>
    size_t rank(const Mat<T> & x) { return arma::rank(x); }
    template<typename T>
    size_t rank(const Mat<T> & x, T tol) { return arma::rank(x,tol); }

    template<typename T>
    Mat<T> repelem(const Mat<T> & a, size_t r, size_t c) { return Mat<T>(arma::repelem(a,r,c)); }

    template<typename T>
    Mat<T> repmat(const Mat<T> & a, size_t r, size_t c) { return Mat<T>(arma::repmat(a,r,c)); }

    template<typename T>
    Mat<T> reshape(const Mat<T> & a, size_t r, size_t c) { return Mat<T>(arma::reshape(a,r,c)); }

    template<typename T>
    Col<T> reverse(const Col<T> & a) { return Col<T>(arma::reverse(a)); }
    template<typename T>
    Mat<T> reverse(const Mat<T> & a, size_t dim=0) { return Mat<T>(arma::reverse(a,dim)); }

    template<typename T>
    Col<T> shift(const Col<T> & a, size_t n) { return Col<T>(arma::shift(a,n)); }
    template<typename T>
    Mat<T> shift(const Mat<T> & a, size_t n, size_t dim = 0) { return Mat<T>(arma::shift(a,n,dim)); }

    template<typename T>
    Col<T> shuffle(const Col<T> & a, size_t n) { return Col<T>(arma::shuffle(a,n)); }
    template<typename T>
    Mat<T> shuffle(const Mat<T> & a, size_t dim = 0) { return Mat<T>(arma::shuffle(a,dim)); }

    template<typename T>
    const arma::SizeMat size(const Col<T> & v) { return arma::size(v); }
    template<typename T>
    const arma::SizeMat size(const Row<T> & v) { return arma::size(v); }
    template<typename T>
    const arma::SizeMat size(const Mat<T> & m) { return arma::size(m); }

    template<typename T>
    void sort(const Col<T> & v, const std::string& dir = "ascend") { arma::sort(v, dir.c_str()); }
    template<typename T>
    void sort(const Row<T> & v, const std::string& dir = "ascend") { arma::sort(v, dir.c_str()); }
    //template<typename T>
    //void sort(const Mat<T> & m, const std::string& dir = "ascend", int dim=0) { arma::sort(m, dir,dim); }

    template<typename T>
    Mat<T> symmatu(const Mat<T> & m) { return Mat<T>(arma::symmatu(m)); }

    template<typename T>
    T trace(const Mat<T> & m) { return arma::trace(m); }

    template<typename T>
    Mat<T> trans(const Mat<T> & m) { return Mat<T>(arma::trans(m)); }

    template<typename T>
    Mat<T> trimatu(const Mat<T> & m, size_t k = 0) { return Mat<T>(arma::symmatu(m,k)); }

    template<typename T>
    Col<T> vectorize(const Mat<T> & m) { return Mat<T>(arma::vectorise(m)); }


    template<typename T>
    Col<T> abs(const Col<T> & a) { return Col<T>(arma::abs(a)); }
    template<typename T>
    Row<T> abs(const Row<T> & a) { return Row<T>(arma::abs(a)); }
    template<typename T>
    Mat<T> abs(const Mat<T> & a) { return Mat<T>(arma::abs(a)); }

    template<typename T>
    Col<T> exp(const Col<T> & a) { return Col<T>(arma::exp(a)); }
    template<typename T>
    Row<T> exp(const Row<T> & a) { return Row<T>(arma::exp(a)); }
    template<typename T>
    Mat<T> exp(const Mat<T> & a) { return Mat<T>(arma::expmat(a)); }

    template<typename T>
    Col<T> exp2(const Col<T> & a) { return Col<T>(arma::exp2(a)); }
    template<typename T>
    Row<T> exp2(const Row<T> & a) { return Row<T>(arma::exp2(a)); }
    template<typename T>
    Mat<T> exp2(const Mat<T> & a) { return Mat<T>(arma::exp2(a)); }

    template<typename T>
    Col<T> exp10(const Col<T> & a) { return Col<T>(arma::exp10(a)); }
    template<typename T>
    Row<T> exp10(const Row<T> & a) { return Row<T>(arma::exp10(a)); }
    template<typename T>
    Mat<T> exp10(const Mat<T> & a) { return Mat<T>(arma::exp10(a)); }

    template<typename T>
    Col<T> expm1(const Col<T> & a) { return Col<T>(arma::expm1(a)); }
    template<typename T>
    Row<T> expm1(const Row<T> & a) { return Row<T>(arma::expm1(a)); }
    template<typename T>
    Mat<T> expm1(const Mat<T> & a) { return Mat<T>(arma::expm1(a)); }

    template<typename T>
    Col<T> trunc_exp(const Col<T> & a) { return Col<T>(arma::trunc_exp(a)); }
    template<typename T>
    Row<T> trunc_exp(const Row<T> & a) { return Row<T>(arma::trunc_exp(a)); }
    template<typename T>
    Mat<T> trunc_exp(const Mat<T> & a) { return Mat<T>(arma::trunc_exp(a)); }

    template<typename T>
    Col<T> log(const Col<T> & a) { return Col<T>(arma::log(a)); }
    template<typename T>
    Row<T> log(const Row<T> & a) { return Row<T>(arma::log(a)); }
    template<typename T>
    Mat<T> log(const Mat<T> & a) { return Mat<T>(arma::log(a)); }

    template<typename T>
    Col<T> log2(const Col<T> & a) { return Col<T>(arma::log2(a)); }
    template<typename T>
    Row<T> log2(const Row<T> & a) { return Row<T>(arma::log2(a)); }
    template<typename T>
    Mat<T> log2(const Mat<T> & a) { return Mat<T>(arma::log2(a)); }

    template<typename T>
    Col<T> log1p(const Col<T> & a) { return Col<T>(arma::log1p(a)); }
    template<typename T>
    Row<T> log1p(const Row<T> & a) { return Row<T>(arma::log1p(a)); }
    template<typename T>
    Mat<T> log1p(const Mat<T> & a) { return Mat<T>(arma::log1p(a)); }

    template<typename T>
    Col<T> log10(const Col<T> & a) { return Col<T>(arma::log10(a)); }
    template<typename T>
    Row<T> log10(const Row<T> & a) { return Row<T>(arma::log10(a)); }
    template<typename T>
    Mat<T> log10(const Mat<T> & a) { return Mat<T>(arma::log10(a)); }

    template<typename T>
    Col<T> trunc_log(const Col<T> & a) { return Col<T>(arma::trunc_log(a)); }
    template<typename T>
    Row<T> trunc_log(const Row<T> & a) { return Row<T>(arma::trunc_log(a)); }
    template<typename T>
    Mat<T> trunc_log(const Mat<T> & a) { return Mat<T>(arma::trunc_log(a)); }

    template<typename T>
    Col<T> pow(const Col<T> & a, T x) { return Col<T>(arma::pow(a,x)); }
    template<typename T>
    Row<T> pow(const Row<T> & a, T x) { return Row<T>(arma::pow(a,x)); }
    template<typename T>
    Mat<T> pow(const Mat<T> & a, T x) { return Mat<T>(arma::pow(a,x)); }

    template<typename T>
    Col<T> square(const Col<T> & a) { return Col<T>(arma::square(a)); }
    template<typename T>
    Row<T> square(const Row<T> & a) { return Row<T>(arma::square(a)); }
    template<typename T>
    Mat<T> square(const Mat<T> & a) { return Mat<T>(arma::square(a)); }

    template<typename T>
    Col<T> sqrt(const Col<T> & a) { return Col<T>(arma::sqrt(a)); }
    template<typename T>
    Row<T> sqrt(const Row<T> & a) { return Row<T>(arma::sqrt(a)); }
    template<typename T>
    Mat<T> sqrt(const Mat<T> & a) { return Mat<T>(arma::sqrt(a)); }

    template<typename T>
    Col<T> floor(const Col<T> & a) { return Col<T>(arma::floor(a)); }
    template<typename T>
    Row<T> floor(const Row<T> & a) { return Row<T>(arma::floor(a)); }
    template<typename T>
    Mat<T> floor(const Mat<T> & a) { return Mat<T>(arma::floor(a)); }

    template<typename T>
    Col<T> ceil(const Col<T> & a) { return Col<T>(arma::ceil(a)); }
    template<typename T>
    Row<T> ceil(const Row<T> & a) { return Row<T>(arma::ceil(a)); }
    template<typename T>
    Mat<T> ceil(const Mat<T> & a) { return Mat<T>(arma::ceil(a)); }

    template<typename T>
    Col<T> round(const Col<T> & a) { return Col<T>(arma::round(a)); }
    template<typename T>
    Row<T> round(const Row<T> & a) { return Row<T>(arma::round(a)); }
    template<typename T>
    Mat<T> round(const Mat<T> & a) { return Mat<T>(arma::round(a)); }

    template<typename T>
    Col<T> trunc(const Col<T> & a) { return Col<T>(arma::trunc(a)); }
    template<typename T>
    Row<T> trunc(const Row<T> & a) { return Row<T>(arma::trunc(a)); }
    template<typename T>
    Mat<T> trunc(const Mat<T> & a) { return Mat<T>(arma::trunc(a)); }

    template<typename T>
    Col<T> erf(const Col<T> & a) { return Col<T>(arma::erf(a)); }
    template<typename T>
    Row<T> erf(const Row<T> & a) { return Row<T>(arma::erf(a)); }
    template<typename T>
    Mat<T> erf(const Mat<T> & a) { return Mat<T>(arma::erf(a)); }

    template<typename T>
    Col<T> erfc(const Col<T> & a) { return Col<T>(arma::erfc(a)); }
    template<typename T>
    Row<T> erfc(const Row<T> & a) { return Row<T>(arma::erfc(a)); }
    template<typename T>
    Mat<T> erfc(const Mat<T> & a) { return Mat<T>(arma::erfc(a)); }

    template<typename T>
    Col<T> tgamma(const Col<T> & a) { return Col<T>(arma::tgamma(a)); }
    template<typename T>
    Row<T> tgamma(const Row<T> & a) { return Row<T>(arma::tgamma(a)); }
    template<typename T>
    Mat<T> tgamma(const Mat<T> & a) { return Mat<T>(arma::tgamma(a)); }

    template<typename T>
    Col<T> lgamma(const Col<T> & a) { return Col<T>(arma::lgamma(a)); }
    template<typename T>
    Row<T> lgamma(const Row<T> & a) { return Row<T>(arma::lgamma(a)); }
    template<typename T>
    Mat<T> lgamma(const Mat<T> & a) { return Mat<T>(arma::lgamma(a)); }

    template<typename T>
    Col<T> sign(const Col<T> & a) { return Col<T>(arma::sign(a)); }
    template<typename T>
    Row<T> sign(const Row<T> & a) { return Row<T>(arma::sign(a)); }
    template<typename T>
    Mat<T> sign(const Mat<T> & a) { return Mat<T>(arma::sign(a)); }

    template<typename T>
    Col<T> cos(const Col<T> & a) { return Col<T>(arma::cos(a)); }
    template<typename T>
    Row<T> cos(const Row<T> & a) { return Row<T>(arma::cos(a)); }
    template<typename T>
    Mat<T> cos(const Mat<T> & a) { return Mat<T>(arma::cos(a)); }

    template<typename T>
    Col<T> sin(const Col<T> & a) { return Col<T>(arma::sin(a)); }
    template<typename T>
    Row<T> sin(const Row<T> & a) { return Row<T>(arma::sin(a)); }
    template<typename T>
    Mat<T> sin(const Mat<T> & a) { return Mat<T>(arma::sin(a)); }

    template<typename T>
    Col<T> tan(const Col<T> & a) { return Col<T>(arma::tan(a)); }
    template<typename T>
    Row<T> tan(const Row<T> & a) { return Row<T>(arma::tan(a)); }
    template<typename T>
    Mat<T> tan(const Mat<T> & a) { return Mat<T>(arma::tan(a)); }

    template<typename T>
    Col<T> acos(const Col<T> & a) { return Col<T>(arma::acos(a)); }
    template<typename T>
    Row<T> acos(const Row<T> & a) { return Row<T>(arma::acos(a)); }
    template<typename T>
    Mat<T> acos(const Mat<T> & a) { return Mat<T>(arma::acos(a)); }

    template<typename T>
    Col<T> asin(const Col<T> & a) { return Col<T>(arma::asin(a)); }
    template<typename T>
    Row<T> asin(const Row<T> & a) { return Row<T>(arma::asin(a)); }
    template<typename T>
    Mat<T> asin(const Mat<T> & a) { return Mat<T>(arma::asin(a)); }

    template<typename T>
    Col<T> atan(const Col<T> & a) { return Col<T>(arma::atan(a)); }
    template<typename T>
    Row<T> atan(const Row<T> & a) { return Row<T>(arma::atan(a)); }
    template<typename T>
    Mat<T> atan(const Mat<T> & a) { return Mat<T>(arma::atan(a)); }


    template<typename T>
    Col<T> cosh(const Col<T> & a) { return Col<T>(arma::cosh(a)); }
    template<typename T>
    Row<T> cosh(const Row<T> & a) { return Row<T>(arma::cosh(a)); }
    template<typename T>
    Mat<T> cosh(const Mat<T> & a) { return Mat<T>(arma::cosh(a)); }

    template<typename T>
    Col<T> sinh(const Col<T> & a) { return Col<T>(arma::sinh(a)); }
    template<typename T>
    Row<T> sinh(const Row<T> & a) { return Row<T>(arma::sinh(a)); }
    template<typename T>
    Mat<T> sinh(const Mat<T> & a) { return Mat<T>(arma::sinh(a)); }

    template<typename T>
    Col<T> tanh(const Col<T> & a) { return Col<T>(arma::tanh(a)); }
    template<typename T>
    Row<T> tanh(const Row<T> & a) { return Row<T>(arma::tanh(a)); }
    template<typename T>
    Mat<T> tanh(const Mat<T> & a) { return Mat<T>(arma::tanh(a)); }

    template<typename T>
    Col<T> acosh(const Col<T> & a) { return Col<T>(arma::acosh(a)); }
    template<typename T>
    Row<T> acosh(const Row<T> & a) { return Row<T>(arma::acosh(a)); }
    template<typename T>
    Mat<T> acosh(const Mat<T> & a) { return Mat<T>(arma::acosh(a)); }

    template<typename T>
    Col<T> asinh(const Col<T> & a) { return Col<T>(arma::asinh(a)); }
    template<typename T>
    Row<T> asinh(const Row<T> & a) { return Row<T>(arma::asinh(a)); }
    template<typename T>
    Mat<T> asinh(const Mat<T> & a) { return Mat<T>(arma::asinh(a)); }

    template<typename T>
    Col<T> atanh(const Col<T> & a) { return Col<T>(arma::atanh(a)); }
    template<typename T>
    Row<T> atanh(const Row<T> & a) { return Row<T>(arma::atanh(a)); }
    template<typename T>
    Mat<T> atanh(const Mat<T> & a) { return Mat<T>(arma::atanh(a)); }

    template<typename T>
    Col<T> atan2(const Col<T> & a, const Col<T> & b) { return Col<T>(arma::atan2(a,b)); }
    template<typename T>
    Row<T> atan2(const Row<T> & a, const Row<T> & b) { return Row<T>(arma::atan2(a,b)); }
    template<typename T>
    Mat<T> atan2(const Mat<T> & a, const Mat<T> & b) { return Mat<T>(arma::atan2(a,b)); }

    template<typename T>
    Col<T> hypot(const Col<T> & a, const Col<T> & b) { return Col<T>(arma::hypot(a,b)); }
    template<typename T>
    Row<T> hypot(const Row<T> & a, const Row<T> & b) { return Row<T>(arma::hypot(a,b)); }
    template<typename T>
    Mat<T> hypot(const Mat<T> & a, const Mat<T> & b) { return Mat<T>(arma::hypot(a,b)); }

    template<typename T>
    Col<T> sinc(const Col<T> & a) { return Col<T>(arma::sinc(a)); }
    template<typename T>
    Row<T> sinc(const Row<T> & a) { return Row<T>(arma::sinc(a)); }
    template<typename T>
    Mat<T> sinc(const Mat<T> & a) { return Mat<T>(arma::sinc(a)); }

    /*
    template<typename T>
    Col<T> col(const Mat<T> & a) { return Col<T>(a)); }

    template<typename T>
    Row<T> row(const Mat<T> & a) { return Col<T>(arma::row(a)); }
    */

    template<typename T>
    void addToEachRow(Mat<T> & m, const Row<T> & v) {
        for(size_t i = 0; i < m.rows(); i++)
            m.row(i) += v;
    }


    template<typename T>
    void addToEachCol(Mat<T> & m, const Col<T> & v) {
        for(size_t i = 0; i < m.cols(); i++)
            m.col(i) += v;
    }

    // complex

    template<typename T>
    void set_real(Col<std::complex<T>> & v, const Col<T> & value) { v.set_real(value); }
    template<typename T>
    void set_real(Row<std::complex<T>> & v, const Row<T> & value) { v.set_real(value); }
    template<typename T>
    void set_real(Mat<std::complex<T>> & m, const Mat<T> & value) { m.set_real(value); }


    template<typename T>
    void set_imag(Col<std::complex<T>> & v, const Col<T> & value) { v.set_imag(value); }
    template<typename T>
    void set_imag(Row<std::complex<T>> & v, const Row<T> & value) { v.set_imag(value); }
    template<typename T>
    void set_imag(Mat<std::complex<T>> & m, const Mat<T> & value) { m.set_imag(value); }

    template<typename T>
    Mat<std::complex<T>> st(Mat<std::complex<T>> & m) { return (m.st()); }

    template<typename T>
    bool is_hermitian(const Mat<std::complex<T>>& m) { return m.is_hermitian(); }

    template<typename T>
    Mat<T> real(Mat<std::complex<T>> & m) { return Mat<T>(arma::real(m)); }

    template<typename T>
    Mat<T> imag(Mat<std::complex<T>> & m) { return Mat<T>(arma::imag(m)); }

    template<typename T>
    std::complex<double> log_det(Mat<T> & m) { return std::complex<double>(arma::log_det(m)); }

    template<typename T>
    bool log_det(Mat<T> & m, double val, double sign) { return arma::log_det(val,sign,m); }

    template<typename T>
    double log_det_sympd(Mat<T> & m) { return arma::log_det_sympd(m); }

    template<typename T>
    bool   log_det_sympd(Mat<T> & m, double result) { return arma::log_det_sympd(result,m); }

    template<typename T>
    Mat<std::complex<T>> logmat(Mat<T> & m) { return Mat<std::complex<T>>(arma::logmat(m)); }

    template<typename T>
    Mat<T> logmat_sympd(Mat<T> & m) { return Mat<T>(arma::logmat_sympd(m)); }

    template<typename T>
    Mat<std::complex<T>> sqrtmat(Mat<T> & m) { return Mat<std::complex<T>>(arma::sqrtmat(m)); }

    template<typename T>
    Mat<T> sqrtmat_sympd(Mat<T> & m) { return Mat<T>(arma::sqrtmat_sympd(m)); }


    // solvers
    template<typename T>
    Mat<T> chol(const Mat<T> & m, const std::string & layout = "upper") { return Mat<T>(arma::chol(m,layout.c_str())); }

    template<typename T>
    Col<T> eig_sym(Mat<T> & m) { return Col<T>(arma::eig_sym(m)); }

    template<typename T>
    Col<std::complex<T>> eig_gen(Mat<T> & m) { return Col<std::complex<T>>(arma::eig_gen(m)); }

    template<typename T>
    Col<std::complex<T>> eig_pair(Mat<T> & a, Mat<T> & b) { return Col<std::complex<T>>(arma::eig_pair(a,b)); }

    template<typename T>
    Mat<T> hess(Mat<T> & x) { return Mat<T>(arma::hess(x)); }

    template<typename T>
    Mat<T> inv(Mat<T> & x) { return Mat<T>(arma::inv(x)); }

    template<typename T>
    Mat<T> inv_sympd(Mat<T> & x) { return Mat<T>(arma::inv_sympd(x)); }

    template<typename T>
    void lu(Mat<T> & L, Mat<T> & U, Mat<T>& P, Mat<T> & X) { arma::lu(L,U,P,X); }

    template<typename T>
    void lu(Mat<T> & L, Mat<T> & U, Mat<T> & X) { arma::lu(L,U,X); }

    // I don't really know what this does, I don't really use the solvers much.
    template<typename T>
    Mat<T> null(Mat<T> & A) { return Mat<T>(arma::null(A)); }

    template<typename T>
    Mat<T> orth(Mat<T> & A) { return Mat<T>(arma::orth(A)); }

    template<typename T>
    Mat<T> pinv(Mat<T> & A) { return Mat<T>(arma::pinv(A)); }

    template<typename T>
    void qr(Mat<T> & Q, Mat<T> & R, Mat<T> & X) { (arma::qr(Q, R, X)); }

    template<typename T>
    void qr_econ(Mat<T> & Q, Mat<T> & R, Mat<T> & X) { (arma::qr_econ(Q, R, X)); }

    template<typename T>
    Mat<T> schur(Mat<T> & A) { return Mat<T>(arma::schur(A)); }

    template<typename T>
    Col<T> solve(Mat<T> & A, Col<T> & b) { return Mat<T>(arma::solve(A,b)); }

    template<typename T>
    Col<T> svd(Mat<T> & X) { return Mat<T>(arma::svd(X)); }

    // DSP
    template<typename T>
    Col<T> conv(Col<T> & A, Col<T> & b) { return Col<T>(arma::conv(A, b)); }

    template<typename T>
    Mat<T> conv2(Mat<T> & A, Mat<T> & b) { return Mat<T>(arma::conv2(A, b)); }

    template<typename T>
    Col<std::complex<T>> fft(Col<T> & A) { return Col<std::complex<T>>(arma::fft(A)); }

    template<typename T>
    Col<std::complex<T>> fft(Col<std::complex<T>> & A) { return Col<std::complex<T>>(arma::fft(A)); }

    template<typename T>
    Col<std::complex<T>> ifft(Col<std::complex<T>> & A) { return Col<std::complex<T>>(arma::ifft(A)); }

    template<typename T>
    Mat<std::complex<T>> fft2(Mat<T> & A) { return Mat<std::complex<T>>(arma::fft2(A)); }

    template<typename T>
    Mat<std::complex<T>> fft2(Mat<std::complex<T>> & A) { return Mat<std::complex<T>>(arma::fft2(A)); }

    template<typename T>
    Mat<std::complex<T>> ifft2(Mat<std::complex<T>> & A) { return Mat<std::complex<T>>(arma::ifft2(A)); }

    template<typename T>
    void interp1(Col<T> & X, Col<T> & Y, Col<T> & XI, Col<T> & YI, const std::string& method = "linear", double extrap = arma::datum::nan ) { (arma::interp1(X, Y, XI, YI)); }

    template<typename T>
    void interp2(Mat<T> & X, Mat<T> & Y, Mat<T> & Z, Mat<T> & XI, Mat<T> & YI, Mat<T> & ZI, const std::string& method = "linear", double extrap = arma::datum::nan ) { (arma::interp2(X, Y, Z, XI, YI, ZI, method.c_str(),extrap)); }


    template<typename T>
    Col<T> polyfit(Col<T> & X, Col<T> & Y, size_t N) { return Col<T>(arma::polyfit(X, Y, N)); }

    template<typename T>
    Col<T> polyval(Col<T> & P, Col<T> & X) { return Col<T>(arma::polyval(P, X)); }

    // Statistics

    template<typename T>
    double mean(Col<T> & V) { return (arma::mean(V)); }
    template<typename T>
    double median(Col<T> & V) { return (arma::median(V)); }
    template<typename T>
    double stddev(Col<T> & V) { return (arma::stddev(V)); }
    template<typename T>
    double var(Col<T> & V) { return (arma::var(V)); }
    template<typename T>
    double range(Col<T> & V) { return (arma::range(V)); }

    /*
    template<typename T>
    double mean(Mat<T> & V) { return (arma::mean(V)); }
    template<typename T>
    double median(Mat<T> & V) { return (arma::median(V)); }
    template<typename T>
    double stddev(Mat<T> & V) { return (arma::stddev(V)); }
    template<typename T>
    double var(Mat<T> & V) { return (arma::var(V)); }
    template<typename T>
    double range(Mat<T> & V) { return (arma::range(V)); }
    */

    template<typename T>
    Col<T> cov(Col<T> & X, Col<T>& Y) { return Col<T>(arma::cov(X,Y)); }
    template<typename T>
    Col<T> cov(Col<T> & X) { return Col<T>(arma::cov(X)); }
    
    /*
    template<typename T>
    Mat<T> cov(Mat<T> & X, Mat<T>& Y) { return Mat<T>(arma::cov(X,Y)); }
    template<typename T>
    Mat<T> cov(Mat<T> & X) { return Mat<T>(arma::cov(X.Mat)); }
    */

    template<typename T>
    Col<T> cor(Col<T> & X, Col<T>& Y) { return Col<T>(arma::cor(X,Y)); }
    template<typename T>
    Col<T> cor(Col<T> & X) { return Col<T>(arma::cor(X)); }

    /*
    template<typename T>
    Mat<T> cor(Mat<T> & X, Mat<T>& Y) { return Mat<T>(arma::cor(X,Y)); }    
    template<typename T>
    Mat<T> cor(Mat<T> & X) { return Mat<T>(arma::cor(X)); }
    */

    template<typename T>
    uvec hist(Col<T> & X) { return (arma::hist(X)); }    
    template<typename T>
    uvec hist(Col<T> & X, size_t n_bins) { return (arma::hist(X,n_bins)); }
    //template<typename T>
    //umat hist(Mat<T> & X) { return (arma::hist(X)); }
    //template<typename T>
    //umat hist(Mat<T> & X, size_t n_bins) { return (arma::cov(X,n_bins)); }
    /*
    template<typename T>
    uvec histc(Col<T> & X, Col<T> & edges) { return (arma::histc(X,edges)); }
    template<typename T>
    uvec histc(Col<T> & X, Col<T> & edges, size_t n_bins) { return (arma::histc(X,edges,n_bins)); }
    template<typename T>
    umat histc(Mat<T> & X, Col<T> & edges) { return (arma::histc(X,edges)); }
    template<typename T>
    umat histc(Mat<T> & X, Col<T> & edges, size_t n_bins) { return (arma::histc(X,edges,n_bins)); }
    */
    template<typename T>
    Col<T> quantile(Col<T> & V, Col<T> & P) { return Col<T>(arma::quantile(V,P));}
    template<typename T>
    Col<T> quantile(Mat<T> & V, Col<T> & P) { return Mat<T>(arma::quantile(V,P));}

    
    template<typename T>
    void princomp(Mat<T> & MC, Mat<T> & X) { (arma::princomp(MC, X));}
    /*
    template<typename T>
    void princomp(Mat<T> & MC, Mat<T> & mat_score, Mat<T> & X) { (arma::princomp(MC,mat_score, X));}

    template<typename T>
    void princomp(Mat<T> & MC, Mat<T> & mat_score, Col<T> & latent, Mat<T> & X) { (arma::princomp(MC,mat_score, latent, X)); }

    template<typename T>
    void princomp(Mat<T> & MC, Mat<T> & mat_score, Col<T> & latent, Col<T> & tsqr, Mat<T> & X) { (arma::princomp(MC,mat_score, latent, tsqr, X)); }

    template<typename T>
    void princomp(Mat<std::complex<T>> & MC, Mat<std::complex<T>> & mat_score, Col<std::complex<T>> & latent, Col<std::complex<T>> & tsqr, Mat<std::complex<T>> & X) { (arma::princomp(MC,mat_score, latent, tsqr, X)); }
    */
    }
%}
