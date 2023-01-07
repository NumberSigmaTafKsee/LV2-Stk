%module essentia
%{
#include "essentia/essentia.h"
#include "essentia/algorithmfactory.h"
#include "essentia/algorithm.h"
%}
%include "std_vector.i"
%include "std_string.i"

%template (float_vector) std::vector<float>;
%template (double_vector) std::vector<double>;

namespace essentia {
    class Parameter {

    public:

    enum ParamType {
        UNDEFINED,

        REAL,
        STRING,
        BOOL,
        INT,
        STEREOSAMPLE,

        VECTOR_REAL,
        VECTOR_STRING,
        VECTOR_BOOL,
        VECTOR_INT,
        VECTOR_STEREOSAMPLE,

        VECTOR_VECTOR_REAL,
        VECTOR_VECTOR_STRING,
        VECTOR_VECTOR_STEREOSAMPLE,

        VECTOR_MATRIX_REAL,

        MAP_VECTOR_REAL,
        MAP_VECTOR_STRING,
        MAP_VECTOR_INT,
        MAP_REAL,

        MATRIX_REAL
    };


    public:

    // Constructor for just declaring type (not providing a value)
    Parameter(ParamType tp);

/*
  // Constructor for simple parameters
  #define SPECIALIZE_CTOR(valueType, paramType, mName)                         \
  Parameter (const valueType& x) : _type(paramType), _##mName(x), _configured(true) {}

  SPECIALIZE_CTOR(std::string,  STRING,       str);
  SPECIALIZE_CTOR(Real,         REAL,         real);
  SPECIALIZE_CTOR(bool,         BOOL,         boolean);
  SPECIALIZE_CTOR(int,          INT,          real);
  SPECIALIZE_CTOR(double,       REAL,         real);
  SPECIALIZE_CTOR(uint,         INT,          real);
  SPECIALIZE_CTOR(StereoSample, STEREOSAMPLE, ssamp);
*/
    Parameter(const std::string&);
    Parameter(const double&);
    Parameter(const bool&);
    //Parameter(const StereoSample&);
    Parameter(const char* x);
/*
// Constructor for vector parameters
  #define SPECIALIZE_VECTOR_CTOR(valueType, paramType)                         \
  Parameter(const std::vector<valueType>& v) : _type(paramType), _configured(true) {\
    _vec.resize(v.size());                                                     \
    for (int i=0; i<int(v.size()); ++i) { _vec[i] = new Parameter(v[i]); }     \
  }
  SPECIALIZE_VECTOR_CTOR(Real,                      VECTOR_REAL);
  SPECIALIZE_VECTOR_CTOR(std::string,               VECTOR_STRING);
  SPECIALIZE_VECTOR_CTOR(bool,                      VECTOR_BOOL);
  SPECIALIZE_VECTOR_CTOR(int,                       VECTOR_INT);
  SPECIALIZE_VECTOR_CTOR(StereoSample,              VECTOR_STEREOSAMPLE);
  SPECIALIZE_VECTOR_CTOR(std::vector<Real>,         VECTOR_VECTOR_REAL);
  SPECIALIZE_VECTOR_CTOR(std::vector<std::string>,  VECTOR_VECTOR_STRING);
  SPECIALIZE_VECTOR_CTOR(std::vector<StereoSample>, VECTOR_VECTOR_STEREOSAMPLE);
  SPECIALIZE_VECTOR_CTOR(TNT::Array2D<Real>,        VECTOR_MATRIX_REAL);

  // Constructor for map parameters
  #define SPECIALIZE_MAP_CTOR(valueType, paramType)                            \
  Parameter(const std::map<std::string, valueType>& m) : _type(paramType), _configured(true) { \
    for (std::map<std::string, valueType>::const_iterator i = m.begin();       \
         i != m.end();                                                         \
         ++i) { _map[(*i).first] = new Parameter((*i).second); }               \
  }

  SPECIALIZE_MAP_CTOR(std::vector<std::string>, MAP_VECTOR_STRING);
  SPECIALIZE_MAP_CTOR(std::vector<Real>, MAP_VECTOR_REAL);
  SPECIALIZE_MAP_CTOR(std::vector<int>, MAP_VECTOR_INT);
  SPECIALIZE_MAP_CTOR(Real, MAP_REAL);
*/

    Parameter(const std::vector<std::string>&);
    Parameter(const std::vector<double>&);
    Parameter(const std::vector<bool>&);
    Parameter(const std::vector<std::vector<std::string>>&);
    Parameter(const std::vector<std::vector<double>>&);
    Parameter(const std::vector<std::vector<int>>&);
    Parameter(const std::vector<std::vector<bool>>&);
    //Parameter(const std::vector<StereoSample>&);
    //Parameter(const std::vector<TNT::Array2D<Real>>&);
    
    Parameter(const std::map<std::string, std::vector<std::string>> &);
    Parameter(const std::map<std::string, std::vector<Real>> &);
    Parameter(const std::map<std::string, std::vector<int>> &);
    Parameter(const std::map<std::string, Real> &);

    //Parameter(const TNT::Array2D<Real>&);

    Parameter(const Parameter& p);

    // also define ctor with a ptr, which allows a nice trick: we can now construct
    // a Parameter from a map<string, Parameter*> which do not necessarily have the
    // same type
    Parameter(const Parameter* p);
    ~Parameter();

    void clear();

    Parameter& operator=(const Parameter& p);
    bool operator==(const Parameter& p) const;
    bool operator!=(const Parameter& p) const;
    ParamType type() const { return _type; }
    bool isConfigured() const { return _configured; }



    std::string toString(int precision = 12) const;
    std::string toLower() const;

    
    bool toBool();
    double toDouble();
    float toFloat();
    StereoSample toStereoSample();

    int toInt() const;
    Real toReal() const;

    /*
    #define TOVECTOR(fname, valueType, paramType)                               \
    std::vector<valueType > toVector##fname() const {                           \
        if (!_configured)                                                         \
        throw EssentiaException("Parameter: parameter has not been configured yet (ParamType=", _type, ")"); \
        if (_type != paramType)                                                   \
        throw EssentiaException("Parameter: parameter is not of type: ", paramType); \
                                                                                \
        std::vector<valueType > result(_vec.size());                              \
        for (int i=0; i<(int)_vec.size(); ++i) {                                  \
        result[i] = _vec[i]->to##fname();                                       \
        }                                                                         \
        return result;                                                            \
    }

    TOVECTOR(Real, Real, VECTOR_REAL)
    TOVECTOR(String, std::string, VECTOR_STRING)
    TOVECTOR(Int, int, VECTOR_INT)
    TOVECTOR(Bool, bool, VECTOR_BOOL)
    TOVECTOR(StereoSample, StereoSample, VECTOR_STEREOSAMPLE)
    TOVECTOR(VectorReal, std::vector<Real>, VECTOR_VECTOR_REAL)
    TOVECTOR(VectorString, std::vector<std::string>, VECTOR_VECTOR_STRING)
    TOVECTOR(VectorStereoSample, std::vector<StereoSample>, VECTOR_VECTOR_STEREOSAMPLE)
    TOVECTOR(MatrixReal, TNT::Array2D<Real>, VECTOR_MATRIX_REAL)
    //  TOVECTOR(MatrixInt, TNT::Array2D<int>, VECTOR_MATRIX_INT)

    #define TOMAP(fname, valueType, paramType)                                   \
    std::map<std::string, valueType > toMap##fname() const {                     \
        if (!_configured)                                                          \
        throw EssentiaException("Parameter: parameter has not been configured yet (ParamType=", _type, ")"); \
        if (_type != paramType)                                                    \
        throw EssentiaException("Parameter: parameter is not of type: ", paramType); \
                                                                                \
        std::map<std::string, valueType > result;                                  \
                                                                                \
        for (std::map<std::string, Parameter*>::const_iterator i = _map.begin();   \
            i != _map.end();                                                      \
            ++i) {                                                                \
        result[i->first] = i->second->to##fname();                               \
        }                                                                          \
                                                                                \
        return result;                                                             \
    }

    TOMAP(VectorReal, std::vector<Real>, MAP_VECTOR_REAL)
    TOMAP(VectorString, std::vector<std::string>, MAP_VECTOR_STRING)
    TOMAP(VectorInt, std::vector<int>, MAP_VECTOR_INT)
    TOMAP(Real, Real, MAP_REAL)
    //  TOMAP(String, std::string, MAP_STRING)
    //  TOMAP(Int, int, MAP_INT)
    //  TOMAP(Bool, bool, MAP_BOOL)
    //  TOMAP(StereoSample, StereoSample, MAP_STEREOSAMPLE)

    #define TOMATRIX(fname, valueType, paramType)                                \
    TNT::Array2D<valueType> toMatrix##fname() const {                            \
        if (!_configured)                                                          \
        throw EssentiaException("Parameter: parameter has not been configured yet (ParamType=", _type, ")");\
        if (_type != paramType)                                                    \
        throw EssentiaException("Parameter: parameter is not of type: ", paramType);\
        TNT::Array2D<valueType> result(_vec.size(), _vec[0]->_vec.size());         \
                                                                                \
        for (int i=0; i<result.dim1(); ++i) {                                      \
        for (int j=0; j<result.dim2(); ++j) {                                    \
            result[i][j] = _vec[i]->_vec[j]->to##fname();                          \
        }                                                                        \
        }                                                                          \
        return result;                                                             \
    }

    TOMATRIX(Real, Real, MATRIX_REAL)
    //  TOMATRIX(String, std::string, MATRIX_STRING)
    //  TOMATRIX(Int, int, MATRIX_INT)
    //  TOMATRIX(Bool, bool, MATRIX_BOOL)
    */

    };
namespace standard {

class Algorithm 
{

 public:
  static const std::string processingMode;

 public:

  virtual ~Algorithm();

  /**
   * Return the names of all the inputs that have been defined for this object.
   */
  std::vector<std::string> inputNames() const;

  /**
   * Return the names of all the outputs that have been defined for this object.
   */
  std::vector<std::string> outputNames() const;

  /**
   * Do the actual computation once that everything is set and configured.
   * The recommended use for this function is to first get the inputs and
   * outputs into local ref variables (const for the inputs) and then do the
   * processing.
   * This allow you also to write a "classic" function call with parameters
   * which you would just wrap with the parameterless function.
   */
  virtual void compute()=0;

  /**
   * This function will be called when doing batch computations between each
   * file that is processed. That is, if your algorithm is some sort of state
   * machine, it allows you to reset it to its original state to process
   * another file without having to delete and reinstantiate it.
   */
  virtual void reset()=0;


  // methods for having access to the types of the inputs/outputs
  //std::vector<const std::type_info*> inputTypes() const;
  //std::vector<const std::type_info*> outputTypes() const;

};

} // namespace standard
} // namespace essentia

%inline %{
    essentia::standard::AlgorithmFactory &factory = essentia::standard::AlgorithmFactory::instance();
    essentia::standard::Algorithm* create_algorithm(const std::string& name) {
        return factory.create(name);
    }
    
%}