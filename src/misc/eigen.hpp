#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
//#include <Eigen/SPQRSupport>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

namespace pano {
namespace misc {

// Generic functor
template <class InternalFunctorT, class T, int NX = Eigen::Dynamic,
          int NY = Eigen::Dynamic>
struct GenericFunctor {
  typedef T Scalar;
  enum { InputsAtCompileTime = NX, ValuesAtCompileTime = NY };
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime>
      JacobianType;

  const int m_inputs, m_values;
  InternalFunctorT m_fun;
  GenericFunctor(const InternalFunctorT &fun)
      : m_fun(fun), m_inputs(InputsAtCompileTime),
        m_values(ValuesAtCompileTime) {}
  GenericFunctor(const InternalFunctorT &fun, int inputs, int values)
      : m_fun(fun), m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  // you should define that in the subclass :
  inline int operator()(const InputType &x, ValueType &v) const {
    m_fun(x, v);
    return 0;
  }
};

template <class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic,
          class InternalFunctorT>
GenericFunctor<InternalFunctorT, T, NX, NY>
MakeGenericFunctor(const InternalFunctorT &fun) {
  return GenericFunctor<InternalFunctorT, T, NX, NY>(fun);
}

template <class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic,
          class InternalFunctorT>
GenericFunctor<InternalFunctorT, T, NX, NY>
MakeGenericFunctor(const InternalFunctorT &fun, int inputs, int values) {
  return GenericFunctor<InternalFunctorT, T, NX, NY>(fun, inputs, values);
}

template <class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic,
          class InternalFunctorT>
Eigen::NumericalDiff<GenericFunctor<InternalFunctorT, T, NX, NY>>
MakeGenericNumericDiffFunctor(const InternalFunctorT &fun) {
  return Eigen::NumericalDiff<GenericFunctor<InternalFunctorT, T, NX, NY>>(
      GenericFunctor<InternalFunctorT, T, NX, NY>(fun));
}

template <class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic,
          class InternalFunctorT>
Eigen::NumericalDiff<GenericFunctor<InternalFunctorT, T, NX, NY>>
MakeGenericNumericDiffFunctor(const InternalFunctorT &fun, int inputs,
                              int values) {
  return Eigen::NumericalDiff<GenericFunctor<InternalFunctorT, T, NX, NY>>(
      GenericFunctor<InternalFunctorT, T, NX, NY>(fun, inputs, values));
}
}
}
