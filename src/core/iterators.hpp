#pragma once

#include <cassert>
#include <chrono>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "meta.hpp"

namespace pano {
namespace core {

//// ITERATORS
template <class IterT, class T>
struct IsIteratorOfType
    : std::is_same<typename std::iterator_traits<IterT>::value_type, T> {};

template <class ContainerT, class T>
struct IsContainerOfType
    : IsIteratorOfType<decltype(std::begin(std::declval<ContainerT>())), T> {};

// range
template <class IterT> struct Range {
  IterT b, e;
  Range(IterT bb, IterT ee) : b(bb), e(ee) {}
  template <class ContainerT>
  explicit Range(ContainerT &&cont) : b(std::begin(cont)), e(std::end(cont)) {}
  IterT begin() const { return b; }
  IterT end() const { return e; }

  template <class FunT> inline void forEach(FunT &&fun) const {
    IterT i = b;
    while (i != e) {
      fun(*i);
      ++i;
    }
  }
};

template <class IterT> Range<IterT> MakeRange(IterT b, IterT e) {
  return Range<IterT>(b, e);
}

// element of container MUST support PredT(ele) -> bool
// ConditionalIterator will automatically skip elements which DO NOT satisfy
// PredT in iteration
template <class IterT, class PredT>
class ConditionalIterator
    : public std::iterator<
          std::forward_iterator_tag,
          typename std::iterator_traits<IterT>::value_type,
          typename std::iterator_traits<IterT>::difference_type,
          typename std::iterator_traits<IterT>::pointer,
          typename std::iterator_traits<IterT>::reference> {

public:
  using Iterator = IterT;

  inline ConditionalIterator(Iterator it_, Iterator end_, PredT pred_ = PredT())
      : _it(it_), _end(end_), _pred(pred_) {
    if (_it != _end && !_pred(*_it))
      ++(*this);
  }

  inline ConditionalIterator &operator++() {
    assert(_it != _end);
    ++_it;
    while (_it != _end && !_pred(*_it))
      ++_it;
    return *this;
  }

  inline reference operator*() const { return *_it; }

  inline pointer operator->() const { return &(*_it); }

  inline bool operator==(const ConditionalIterator &i) const {
    return _it == i._it;
  }

  inline bool operator!=(const ConditionalIterator &i) const {
    return !(*this == i);
  }

  inline Iterator internalIterator() const { return _it; }

private:
  IterT _it;
  IterT _end;
  PredT _pred;
};

// class ConditionalContainerWrapper
template <class ContainerT, class ElementPredT>
class ConditionalContainerWrapper {
public:
  using OriginalIterator = typename ContainerT::iterator;
  using iterator = ConditionalIterator<OriginalIterator, ElementPredT>;
  using value_type = typename std::iterator_traits<iterator>::value_type;

  inline ConditionalContainerWrapper(ContainerT *cont_,
                                     ElementPredT elePred_ = ElementPredT())
      : _cont(cont_), _elePred(elePred_) {}
  inline iterator begin() {
    return iterator(std::begin(*_cont), std::end(*_cont), _elePred);
  }
  inline iterator end() {
    return iterator(std::end(*_cont), std::end(*_cont), _elePred);
  }
  inline iterator begin() const {
    return iterator(std::begin(*_cont), std::end(*_cont), _elePred);
  }
  inline iterator end() const {
    return iterator(std::end(*_cont), std::end(*_cont), _elePred);
  }

private:
  ContainerT *_cont;
  ElementPredT _elePred;
};

// class ConstConditionalContainerWrapper
template <class ContainerT, class ElementPredT>
class ConstConditionalContainerWrapper {
public:
  using OriginalIterator = typename ContainerT::const_iterator;
  using iterator = ConditionalIterator<OriginalIterator, ElementPredT>;
  using value_type = typename std::iterator_traits<iterator>::value_type;

  inline ConstConditionalContainerWrapper(
      const ContainerT *cont_, ElementPredT elePred_ = ElementPredT())
      : _cont(cont_), _elePred(elePred_) {}
  inline iterator begin() const {
    return iterator(std::begin(*_cont), std::end(*_cont), _elePred);
  }
  inline iterator end() const {
    return iterator(std::end(*_cont), std::end(*_cont), _elePred);
  }

private:
  const ContainerT *_cont;
  ElementPredT _elePred;
};

// make conditional container
template <class ContainerT, class ElementPredT>
inline ConditionalContainerWrapper<ContainerT, ElementPredT>
MakeConditionalContainer(ContainerT *cont_,
                         ElementPredT elePred_ = ElementPredT()) {
  return ConditionalContainerWrapper<ContainerT, ElementPredT>(cont_, elePred_);
}

template <class ContainerT, class ElementPredT>
inline ConstConditionalContainerWrapper<ContainerT, ElementPredT>
MakeConditionalContainer(const ContainerT *cont_,
                         ElementPredT elePred_ = ElementPredT()) {
  return ConstConditionalContainerWrapper<ContainerT, ElementPredT>(cont_,
                                                                    elePred_);
}
}
}
