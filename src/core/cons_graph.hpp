#pragma once

#include "handle.hpp"
#include <numeric>

namespace pano {
namespace core {

// configuration for constraint
template <class ConstraintDataT, class... ComponentTs>
struct ConstraintConfig {};

template <class T> struct IsConstraintConfig : std::false_type {};

template <class ConstraintDataT, class... ComponentTs>
struct IsConstraintConfig<ConstraintConfig<ConstraintDataT, ComponentTs...>>
    : std::true_type {};

namespace {

// get member type of constraint config
template <class ConstraintConfigT> struct MemberTypesOfConstraintConfig {};

template <class ConstraintDataType, class... ComponentTs>
struct MemberTypesOfConstraintConfig<
    ConstraintConfig<ConstraintDataType, ComponentTs...>> {
  using ConstraintData = ConstraintDataType;
  using ComponentDataTuple = std::tuple<ComponentTs...>;
};

// constraint data tuple from constraint config tuple
template <class ConstraintConfigTupleT>
struct ConstraintDataTupleFromConstraintConfigTuple {};

template <class... ConstraintConfigTs>
struct ConstraintDataTupleFromConstraintConfigTuple<
    std::tuple<ConstraintConfigTs...>> {
  using type = std::tuple<typename MemberTypesOfConstraintConfig<
      ConstraintConfigTs>::ConstraintData...>;
};
}

// make a layer content tuple with 2 levels
template <class ComponentTagT, class ConstraintTagT, class ComponentDataT,
          class... ConstraintConfigTs>
struct ComponentTopo {
  using ConstraintDataTuple = std::tuple<typename MemberTypesOfConstraintConfig<
      ConstraintConfigTs>::ConstraintData...>;
  std::tuple<std::set<
      HandleOfType<ConstraintTagT, typename MemberTypesOfConstraintConfig<
                                       ConstraintConfigTs>::ConstraintData>>...>
      allConstraints;

  HandleOfType<ComponentTagT, ComponentDataT> hd;
  explicit inline ComponentTopo(int id = -1) : hd(id) {}

  template <class ConstraintDataT>
  inline std::set<HandleOfType<ConstraintTagT, ConstraintDataT>> &
  constraints() {
    enum {
      _idx =
          TypeFirstLocationInTuple<ConstraintDataT, ConstraintDataTuple>::value
    };
    static_assert(_idx >= 0, "Invalid ConstraintDataT");
    return std::get<_idx>(allConstraints);
  }

  template <class ConstraintDataT>
  inline const std::set<HandleOfType<ConstraintTagT, ConstraintDataT>> &
  constraints() const {
    enum {
      _idx =
          TypeFirstLocationInTuple<ConstraintDataT, ConstraintDataTuple>::value
    };
    static_assert(_idx >= 0, "Invalid ConstraintDataT");
    return std::get<_idx>(allConstraints);
  }

  template <class Archive> inline void serialize(Archive &ar) {
    ar(allConstraints, hd);
  }
};

template <class ComponentTagT, class ConstraintTagT, class ConstraintDataT,
          class... ComponentDataTs>
struct ConstraintTopo {
  static_assert(sizeof...(ComponentDataTs) > 0,
                "ComponentDataTs must be more than zero");
  using ComponentDataTuple = std::tuple<ComponentDataTs...>;
  using ComponentHandlesTuple =
      std::tuple<HandleOfType<ComponentTagT, ComponentDataTs>...>;
  ComponentHandlesTuple allComponents;

  HandleOfType<ConstraintTagT, ConstraintDataT> hd;
  explicit inline ConstraintTopo() : hd(-1) {}
  explicit inline ConstraintTopo(
      int id, HandleOfType<ComponentTagT, ComponentDataTs>... componentHandles)
      : hd(id), allComponents(std::forward_as_tuple(componentHandles...)) {}

  template <class ComponentDataT,
            int Idx = TypeFirstLocationInTuple<ComponentDataT,
                                               ComponentDataTuple>::value>
  inline typename std::tuple_element<Idx, ComponentHandlesTuple>::type
  component() {
    static_assert(Idx >= 0, "Invalid ComponentDataT");
    return std::get<Idx>(allComponents);
  }

  template <class ComponentDataT,
            int Idx = TypeFirstLocationInTuple<ComponentDataT,
                                               ComponentDataTuple>::value>
  inline const typename std::tuple_element<Idx, ComponentHandlesTuple>::type &
  component() const {
    static_assert(Idx >= 0, "Invalid ComponentDataT");
    return std::get<Idx>(allComponents);
  }

  template <int Idx>
  inline typename std::tuple_element<Idx, ComponentHandlesTuple>::type
  component() {
    return std::get<Idx>(allComponents);
  }

  template <int Idx>
  inline const typename std::tuple_element<Idx, ComponentHandlesTuple>::type
  component() const {
    return std::get<Idx>(allComponents);
  }

  template <class Archive> inline void serialize(Archive &ar) {
    ar(allComponents, hd);
  }
};

namespace {

template <class ComponentTagT, class ConstraintTagT, class ComponentDataT,
          class ConstraintConfigTupleT>
struct ComponentTripletArrayFromConstraintConfigTuple {};

template <class ComponentTagT, class ConstraintTagT, class ComponentDataT,
          class... ConstraintConfigTs>
struct ComponentTripletArrayFromConstraintConfigTuple<
    ComponentTagT, ConstraintTagT, ComponentDataT,
    std::tuple<ConstraintConfigTs...>> {
  using type =
      TripletArray<ComponentTopo<ComponentTagT, ConstraintTagT, ComponentDataT,
                                 ConstraintConfigTs...>,
                   ComponentDataT>;
};

template <class ComponentTagT, class ConstraintTagT, class ComponentDataTupleT,
          class ConstraintConfigTupleT>
struct ComponentTripletArrayTupleFromConstraintConfigTuple {};

template <class ComponentTagT, class ConstraintTagT,
          class ConstraintConfigTupleT, class... ComponentDataTs>
struct ComponentTripletArrayTupleFromConstraintConfigTuple<
    ComponentTagT, ConstraintTagT, std::tuple<ComponentDataTs...>,
    ConstraintConfigTupleT> {
  using type =
      std::tuple<typename ComponentTripletArrayFromConstraintConfigTuple<
          ComponentTagT, ConstraintTagT, ComponentDataTs,
          ConstraintConfigTupleT>::type...>;
};

template <class ComponentTagT, class ConstraintTagT, class ConstraintConfigT>
struct ConstraintTripletArrayFromConstraintConfig {};

template <class ComponentTagT, class ConstraintTagT, class ConstraintDataT,
          class... ComponentDataTs>
struct ConstraintTripletArrayFromConstraintConfig<
    ComponentTagT, ConstraintTagT,
    ConstraintConfig<ConstraintDataT, ComponentDataTs...>> {
  using type = TripletArray<ConstraintTopo<ComponentTagT, ConstraintTagT,
                                           ConstraintDataT, ComponentDataTs...>,
                            ConstraintDataT>;
};

template <class ComponentTagT, class ConstraintTagT,
          class ConstraintConfigTupleT>
struct ConstraintTripetArrayTupleFromConstraintConfigTuple {};

template <class ComponentTagT, class ConstraintTagT,
          class... ConstraintConfigTs>
struct ConstraintTripetArrayTupleFromConstraintConfigTuple<
    ComponentTagT, ConstraintTagT, std::tuple<ConstraintConfigTs...>> {
  using type = std::tuple<typename ConstraintTripletArrayFromConstraintConfig<
      ComponentTagT, ConstraintTagT, ConstraintConfigTs>::type...>;
};
}

struct ComponentTag {};
struct ConstraintTag {};

template <class DataT>
using ComponentHandle = HandleOfType<ComponentTag, DataT>;

template <class DataT>
using ConstraintHandle = HandleOfType<ConstraintTag, DataT>;

template <class HandleOfTypeT> struct TagOfHandleOfType {};

template <class TagT, class DataT>
struct TagOfHandleOfType<HandleOfType<TagT, DataT>> {
  using type = TagT;
};

template <class HandleOfTypeT> struct DataTypeOfHandleOfType {};

template <class TagT, class DataT>
struct DataTypeOfHandleOfType<HandleOfType<TagT, DataT>> {
  using type = DataT;
};

template <class ComponentDataTupleT, class ConstraintConfigTupleT>
class ConstraintGraph {
  static_assert(IsTuple<ComponentDataTupleT>::value,
                "ComponentDataTupleT must be std::tuple");
  static_assert(IsTuple<ConstraintConfigTupleT>::value,
                "ConstraintConfigTupleT must be std::tuple");
  using ComponentsTripletArrayTuple =
      typename ComponentTripletArrayTupleFromConstraintConfigTuple<
          ComponentTag, ConstraintTag, ComponentDataTupleT,
          ConstraintConfigTupleT>::type;
  using ConstraintsTripletArrayTuple =
      typename ConstraintTripetArrayTupleFromConstraintConfigTuple<
          ComponentTag, ConstraintTag, ConstraintConfigTupleT>::type;

public:
  using ComponentDataTupleType = ComponentDataTupleT;
  static const int ComponentDataTypeNum =
      std::tuple_size<ComponentDataTupleType>::value;
  using ConstraintDataTupleType =
      typename ConstraintDataTupleFromConstraintConfigTuple<
          ConstraintConfigTupleT>::type;
  static const int ConstraintDataTypeNum =
      std::tuple_size<ConstraintDataTupleType>::value;

  ConstraintGraph() {}
  ConstraintGraph(const ConstraintGraph &) = default;
  ConstraintGraph(ConstraintGraph &&cg)
      : _components(std::move(cg._components)),
        _constraints(std::move(cg._constraints)) {}
  ConstraintGraph &operator=(ConstraintGraph &&cg) {
    _components = std::move(cg._components);
    _constraints = std::move(cg._constraints);
    return *this;
  }

public:
  template <class DataT>
  inline ComponentHandle<std::decay_t<DataT>> addComponent(DataT &&data) {
    using _DataType = std::decay_t<DataT>;
    enum {
      _idx = TypeFirstLocationInTuple<_DataType, ComponentDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    using _TripletType = typename std::tuple_element<
        _idx, ComponentsTripletArrayTuple>::type::value_type;
    using _TopoType = typename _TripletType::TopoType;
    auto &correspondingComponents = std::get<_idx>(_components);
    _TopoType topo;
    topo.hd = ComponentHandle<_DataType>(correspondingComponents.size());
    correspondingComponents.emplace_back(std::move(topo),
                                         std::forward<DataT>(data), true);
    return topo.hd;
  }

  template <class DataT, class... ComponentHandleTs>
  ConstraintHandle<std::decay_t<DataT>>
  addConstraint(DataT &&data, ComponentHandleTs... comps) {
    using _DataType = std::decay_t<DataT>;
    enum {
      _idx = TypeFirstLocationInTuple<_DataType, ConstraintDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    using _TripletType = typename std::tuple_element<
        _idx, ConstraintsTripletArrayTuple>::type::value_type;
    using _TopoType = typename _TripletType::TopoType;

    auto &correspondingConstraints = std::get<_idx>(_constraints);
    _TopoType topo(correspondingConstraints.size(),
                   std::forward<ComponentHandleTs>(comps)...);
    correspondingConstraints.emplace_back(std::move(topo),
                                          std::forward<DataT>(data), true);

    // update constraint handles in each component
    bool dummy[] = {insertConstraintHandleInComponent(topo.hd, comps)...};
    return topo.hd;
  }

  // removed ?
  template <class DataT> inline bool removed(ComponentHandle<DataT> h) const {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    using _TripletType = typename std::tuple_element<
        _idx, ComponentsTripletArrayTuple>::type::value_type;
    using _TopoType = typename _TripletType::TopoType;
    auto &correspondingComponents = std::get<_idx>(_components);
    if (h.invalid())
      return true;
    return !correspondingComponents.at(h.id).exists;
  }

  template <class DataT> inline bool removed(ConstraintHandle<DataT> h) const {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    using _TripletType = typename std::tuple_element<
        _idx, ConstraintsTripletArrayTuple>::type::value_type;
    using _TopoType = typename _TripletType::TopoType;
    auto &correspondingConstraints = std::get<_idx>(_constraints);
    if (h.invalid())
      return true;
    return !correspondingConstraints.at(h.id).exists;
  }

  template <class DataT> inline void remove(ConstraintHandle<DataT> h) {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    using _TripletType = typename std::tuple_element<
        _idx, ConstraintsTripletArrayTuple>::type::value_type;
    using _TopoType = typename _TripletType::TopoType;

    if (removed(h))
      return;

    auto &correspondingConstraints = std::get<_idx>(_constraints);
    // remove constraint handles in each component
    removeConstraintHandleInRelatedComponentsUsingSequence(
        h, typename SequenceGenerator<std::tuple_size<
               typename _TopoType::ComponentDataTuple>::value>::type());

    correspondingConstraints[h.id].exists = false;
  }

  template <class DataT> void remove(ComponentHandle<DataT> h) {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    using _TripletType = typename std::tuple_element<
        _idx, ComponentsTripletArrayTuple>::type::value_type;
    using _TopoType = typename _TripletType::TopoType;

    if (removed(h))
      return;

    auto &correspondingComponents = std::get<_idx>(_components);
    correspondingComponents[h.id].exists = false;

    // remove related constraints
    removeHandlesInContainerTupleUsingSequence(
        topo(h).allConstraints,
        typename SequenceGenerator<std::tuple_size<
            typename _TopoType::ConstraintDataTuple>::value>::type());
  }

  // component info
  template <class DataT>
  inline const DataT &data(ComponentHandle<DataT> hd) const {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    return std::get<_idx>(_components).at(hd.id).data;
  }

  template <class DataT> inline DataT &data(ComponentHandle<DataT> hd) {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    return std::get<_idx>(_components)[hd.id].data;
  }

  template <class DataT>
  inline const typename std::tuple_element<
      TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value,
      ComponentsTripletArrayTuple>::type::value_type::TopoType &
  topo(ComponentHandle<DataT> hd) const {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    return std::get<_idx>(_components).at(hd.id).topo;
  }

  template <class DataT>
  inline typename std::tuple_element<
      TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value,
      ComponentsTripletArrayTuple>::type::value_type::TopoType &
  topo(ComponentHandle<DataT> hd) {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    return std::get<_idx>(_components).at(hd.id).topo;
  }

  // constraint info
  template <class DataT>
  inline const DataT &data(ConstraintHandle<DataT> hd) const {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    return std::get<_idx>(_constraints).at(hd.id).data;
  }

  template <class DataT> inline DataT &data(ConstraintHandle<DataT> hd) {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    return std::get<_idx>(_constraints)[hd.id].data;
  }

  template <class DataT>
  inline const typename std::tuple_element<
      TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value,
      ConstraintsTripletArrayTuple>::type::value_type::TopoType &
  topo(ConstraintHandle<DataT> hd) const {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    return std::get<_idx>(_constraints).at(hd.id).topo;
  }

  template <class DataT>
  inline typename std::tuple_element<
      TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value,
      ConstraintsTripletArrayTuple>::type::value_type::TopoType &
  topo(ConstraintHandle<DataT> hd) {
    enum {
      _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value
    };
    static_assert(_idx >= 0, "DataT is not supported here!");
    return std::get<_idx>(_constraints)[hd.id].topo;
  }

  // all info

  template <class DataT>
  inline typename ConditionalContainerTypeFromTripletArrayType<
      typename std::tuple_element<
          TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value,
          ComponentsTripletArrayTuple>::type>::const_type
  components() const {
    return MakeConditionalContainer(
        std::get<
            TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value>(
            _components));
  }

  template <class DataT>
  inline typename ConditionalContainerTypeFromTripletArrayType<
      typename std::tuple_element<
          TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value,
          ComponentsTripletArrayTuple>::type>::type
  components() {
    return MakeConditionalContainer(
        std::get<
            TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value>(
            _components));
  }

  template <class DataT>
  inline const typename std::tuple_element<
      TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value,
      ComponentsTripletArrayTuple>::type &
  internalComponents() const {
    return std::get<
        TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value>(
        _components);
  }

  template <class DataT>
  inline typename std::tuple_element<
      TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value,
      ComponentsTripletArrayTuple>::type &
  internalComponents() {
    return std::get<
        TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value>(
        _components);
  }

  inline const ComponentsTripletArrayTuple &allComponents() const {
    return _components;
  }
  inline ComponentsTripletArrayTuple &allComponents() { return _components; }

  template <class DataT>
  inline typename ConditionalContainerTypeFromTripletArrayType<
      typename std::tuple_element<
          TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value,
          ConstraintsTripletArrayTuple>::type>::const_type
  constraints() const {
    return MakeConditionalContainer(
        std::get<
            TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value>(
            _constraints));
  }

  template <class DataT>
  inline typename ConditionalContainerTypeFromTripletArrayType<
      typename std::tuple_element<
          TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value,
          ConstraintsTripletArrayTuple>::type>::type
  constraints() {
    return MakeConditionalContainer(
        std::get<
            TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value>(
            _constraints));
  }

  template <class DataT>
  inline const typename std::tuple_element<
      TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value,
      ConstraintsTripletArrayTuple>::type &
  internalConstraints() const {
    return std::get<
        TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value>(
        _constraints);
  }

  template <class DataT>
  inline typename std::tuple_element<
      TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value,
      ConstraintsTripletArrayTuple>::type &
  internalConstraints() {
    return std::get<
        TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value>(
        _constraints);
  }

  inline const ConstraintsTripletArrayTuple &allConstraints() const {
    return _constraints;
  }
  inline ConstraintsTripletArrayTuple &allConstraints() { return _constraints; }

  template <class DataT, class ValueT>
  inline HandledTable<ComponentHandle<DataT>, ValueT>
  createComponentTable(const ValueT &defaultValue = ValueT()) const {
    return HandledTable<ComponentHandle<DataT>, ValueT>(
        internalComponents<DataT>().size(), defaultValue);
  }

  template <class DataT, class ValueT>
  inline HandledTable<ConstraintHandle<DataT>, ValueT>
  createConstraintTable(const ValueT &defaultValue = ValueT()) const {
    return HandledTable<ConstraintHandle<DataT>, ValueT>(
        internalConstraints<DataT>().size(), defaultValue);
  }

  template <class ValueT>
  inline auto createTableForAllComponents() const
      -> decltype(MakeHandledTableForAllComponents<ValueT>(*this)) {
    return MakeHandledTableForAllComponents<ValueT>(*this);
  }

  template <class ValueT>
  inline auto createTableForAllConstraints() const
      -> decltype(MakeHandledTableForAllConstraints<ValueT>(*this)) {
    return MakeHandledTableForAllConstraints<ValueT>(*this);
  }

  // gc stuff
  // void gc() {
  //    auto componentNewLocs = removeAndMapComponentsUsingSequence(
  //        typename SequenceGenerator<ComponentDataTypeNum>::type());
  //    auto constraintNewLocs = removeAndMapConstraintsUsingSequence(
  //        typename SequenceGenerator<ConstraintDataTypeNum>::type());

  //}

  // merge
  void merge(const ConstraintGraph &cg) {
    auto thisComponentSizes = TupleMap(SizeFunctor(), _components);
    auto thisConstraintSizes = TupleMap(SizeFunctor(), _constraints);
    mergeComponentsUsingSequence(
        cg, thisComponentSizes, thisConstraintSizes,
        SequenceGenerator<ComponentDataTypeNum>::type());
    mergeConstraintsUsingSequence(
        cg, thisComponentSizes, thisConstraintSizes,
        SequenceGenerator<ConstraintDataTypeNum>::type());
  }

  inline size_t allComponentsNum() const {
    return allComponentsNumUsingSequence(
        SequenceGenerator<ComponentDataTypeNum>::type());
  }
  inline size_t allConstraintsNum() const {
    return allConstraintsNumUsingSequence(
        SequenceGenerator<ConstraintDataTypeNum>::type());
  }

  template <class Archiver> inline void serialize(Archiver &ar) {
    ar(_components, _constraints);
  }

private:
  // add helpers
  template <class DataT, class ComponentHandleT>
  inline bool insertConstraintHandleInComponent(ConstraintHandle<DataT> consH,
                                                ComponentHandleT compH) {
    topo(compH).constraints<DataT>().insert(consH);
    return true;
  }

  // remove helpers
  template <class DataT, class ComponentHandleT>
  inline bool removeConstraintHandleInComponent(ConstraintHandle<DataT> consH,
                                                ComponentHandleT &compH) {
    if (compH.invalid() || removed(compH))
      return false;
    topo(compH).constraints<DataT>().erase(consH);
    return true;
  }

  template <class DataT, int... Is>
  inline bool removeConstraintHandleInRelatedComponentsUsingSequence(
      ConstraintHandle<DataT> consH, Sequence<Is...>) {
    bool d[] = {removeConstraintHandleInComponent(
        consH, std::get<Is>(topo(consH).allComponents))...};
    return true;
  }

  template <class ContainerT>
  inline bool removeHandlesInContainer(ContainerT &&hs) {
    for (auto &&h : hs) {
      remove(h);
    }
    return true;
  }

  template <class ContainerTupleT, int... Is>
  inline bool removeHandlesInContainerTupleUsingSequence(ContainerTupleT &hs,
                                                         Sequence<Is...>) {
    bool d[] = {removeHandlesInContainer(std::get<Is>(hs))...};
    return true;
  }

  //// gc helpers
  // template <int ... Is>
  // inline std::tuple<std::vector<ComponentHandle<typename
  // std::tuple_element<Is, ComponentDataTupleType>::type>>...>
  //    removeAndMapComponentsUsingSequence(Sequence<Is...>) {
  //    std::tuple<std::vector<ComponentHandle<typename std::tuple_element<Is,
  //    ComponentDataTupleType>::type>>...> hs;
  //    int dummy[] = { RemoveAndMap(std::get<Is>(_components),
  //    std::get<Is>(hs)) ... };
  //    return hs;
  //}

  // template <int ... Is>
  // inline std::tuple<std::vector<ConstraintHandle<typename
  // std::tuple_element<Is, ConstraintDataTupleType>::type>>...>
  //    removeAndMapConstraintsUsingSequence(Sequence<Is...>) {
  //    std::tuple<std::vector<ConstraintHandle<typename std::tuple_element<Is,
  //    ConstraintDataTupleType>::type>>...> hs;
  //    int dummy[] = { RemoveAndMap(std::get<Is>(_constraints),
  //    std::get<Is>(hs)) ... };
  //    return hs;
  //}

  // inline bool updateComponents

  // template <int ...Is>
  // inline void updateComponentsUsingSequence(
  //    const std::tuple<std::vector<ComponentHandle<typename
  //    std::tuple_element<Is, ComponentDataTupleType>::type>>...> &
  //    componentNewLocs,
  //    Sequence<Is...>){
  //    for (auto & c : std::get<Is>(_components)){
  //        c.topo.hd = UpdateOldHandle(std::get<Is>(componentNewLocs),
  //        c.topo.hd);
  //    }
  //}

  // merge helpers
  // udpate component triplet
  template <int ConsId, class ConstraintHandlesSetTupleT,
            class ConstraintSizesTupleT>
  static inline int offsetConstriantHandlesInSetTupleWithIndex(
      ConstraintHandlesSetTupleT &chs,
      const ConstraintSizesTupleT &thisConsSizes) {
    auto &hs = std::get<ConsId>(chs);
    auto oldHs = hs;
    hs.clear();
    for (auto h : oldHs) {
      h.id += std::get<ConsId>(thisConsSizes);
      hs.insert(h);
    }
    return 0;
  }

  template <class ConstraintHandlesSetTupleT, class ConstraintSizesTupleT,
            int... ConsIds>
  static inline int offsetConstriantHandlesInSetTupleUsingSequence(
      ConstraintHandlesSetTupleT &chs,
      const ConstraintSizesTupleT &thisConsSizes, Sequence<ConsIds...>) {
    int dummy[] = {offsetConstriantHandlesInSetTupleWithIndex<ConsIds>(
        chs, thisConsSizes)...};
    return 0;
  }

  template <int CompId, class ComponentSizesTupleT, class ConstraintSizesTupleT>
  static inline int
  offsetComponentTripletsWithIndex(ComponentsTripletArrayTuple &allComps,
                                   const ComponentSizesTupleT &thisCompSizes,
                                   const ConstraintSizesTupleT &thisConsSizes) {
    auto &comps = std::get<CompId>(allComps);
    size_t thisCompSize = std::get<CompId>(thisCompSizes);
    for (size_t i = thisCompSize; i < comps.size();
         i++) { // update new inserted triplet topos
      comps[i].topo.hd.id += thisCompSize;
      offsetConstriantHandlesInSetTupleUsingSequence(
          comps[i].topo.allConstraints, thisConsSizes,
          SequenceGenerator<ConstraintDataTypeNum>::type());
    }
    return 0;
  }

  template <int CompId, class ComponentSizesTupleT, class ConstraintSizesTupleT>
  inline int
  mergeComponentsWithIndex(const ConstraintGraph &cg,
                           const ComponentSizesTupleT &thisCompSizes,
                           const ConstraintSizesTupleT &thisConsSizes) {
    auto &comps = std::get<CompId>(_components);
    const auto &anotherComps = std::get<CompId>(cg._components);
    comps.insert(comps.end(), anotherComps.begin(), anotherComps.end());
    offsetComponentTripletsWithIndex<CompId>(_components, thisCompSizes,
                                             thisConsSizes);
    return 0;
  }

  template <class ComponentSizesTupleT, class ConstraintSizesTupleT,
            int... CompIds>
  inline int mergeComponentsUsingSequence(
      const ConstraintGraph &cg, const ComponentSizesTupleT &thisCompSizes,
      const ConstraintSizesTupleT &thisConsSizes, Sequence<CompIds...>) {
    int dummy[] = {
        mergeComponentsWithIndex<CompIds>(cg, thisCompSizes, thisConsSizes)...};
    return 0;
  }

  // update constraint triplet
  template <int CompPositionInTuple, class ComponentHandlesTupleT,
            class ComponentSizesTupleT>
  static inline int offsetComponentHandlesInTupleWithIndex(
      ComponentHandlesTupleT &chs, const ComponentSizesTupleT &thisCompSizes) {
    using CompType =
        typename DataTypeOfHandleOfType<typename std::tuple_element<
            CompPositionInTuple, ComponentHandlesTupleT>::type>::type;
    enum {
      Idx = TypeFirstLocationInTuple<CompType, ComponentDataTupleType>::value
    };
    std::get<CompPositionInTuple>(chs).id += std::get<Idx>(thisCompSizes);
    return 0;
  }

  template <class ComponentHandlesTupleT, class ComponentSizesTupleT,
            int... CompPositionInTuples>
  static inline int offsetComponentHandlesInTupleUsingSequence(
      ComponentHandlesTupleT &chs, const ComponentSizesTupleT &thisCompSizes,
      Sequence<CompPositionInTuples...>) {
    int dummy[] = {offsetComponentHandlesInTupleWithIndex<CompPositionInTuples>(
        chs, thisCompSizes)...};
    return 0;
  }

  template <int ConsId, class ComponentSizesTupleT, class ConstraintSizesTupleT>
  static inline int offsetConstraintTripletsWithIndex(
      ConstraintsTripletArrayTuple &allCons,
      const ComponentSizesTupleT &thisCompSizes,
      const ConstraintSizesTupleT &thisConsSizes) {
    auto &cons = std::get<ConsId>(allCons);
    size_t thisConsSize = std::get<ConsId>(thisConsSizes);
    for (size_t i = thisConsSize; i < cons.size(); i++) {
      cons[i].topo.hd.id += thisConsSize;
      offsetComponentHandlesInTupleUsingSequence(
          cons[i].topo.allComponents, thisCompSizes,
          SequenceGenerator<std::tuple_size<decltype(
              cons[i].topo.allComponents)>::value>::type());
    }
    return 0;
  }

  template <int ConsId, class ComponentSizesTupleT, class ConstraintSizesTupleT>
  inline int
  mergeConstraintsWithIndex(const ConstraintGraph &cg,
                            const ComponentSizesTupleT &thisCompSizes,
                            const ConstraintSizesTupleT &thisConsSizes) {
    auto &cons = std::get<ConsId>(_constraints);
    const auto &anotherCons = std::get<ConsId>(cg._constraints);
    cons.insert(cons.end(), anotherCons.begin(), anotherCons.end());
    offsetConstraintTripletsWithIndex<ConsId>(_constraints, thisCompSizes,
                                              thisConsSizes);
    return 0;
  }

  template <class ComponentSizesTupleT, class ConstraintSizesTupleT,
            int... ConsIds>
  inline int mergeConstraintsUsingSequence(
      const ConstraintGraph &cg, const ComponentSizesTupleT &thisCompSizes,
      const ConstraintSizesTupleT &thisConsSizes, Sequence<ConsIds...>) {
    int dummy[] = {mergeConstraintsWithIndex<ConsIds>(cg, thisCompSizes,
                                                      thisConsSizes)...};
    return 0;
  }

  template <int... CompIds>
  size_t allComponentsNumUsingSequence(Sequence<CompIds...>) const {
    size_t szs[] = {std::get<CompIds>(_components).size()...};
    return std::accumulate(std::begin(szs), std::end(szs), 0ull);
  }

  template <int... ConsIds>
  size_t allConstraintsNumUsingSequence(Sequence<ConsIds...>) const {
    size_t szs[] = {std::get<ConsIds>(_constraints).size()...};
    return std::accumulate(std::begin(szs), std::end(szs), 0ull);
  }

private:
  ComponentsTripletArrayTuple _components;
  ConstraintsTripletArrayTuple _constraints;
};

template <class DataT, class ConstraintGraphT>
struct ComponentHandledTableFromConstraintGraph {};

template <class DataT, class... ComponentTs, class ConstraintConfigTupleT>
struct ComponentHandledTableFromConstraintGraph<
    DataT,
    ConstraintGraph<std::tuple<ComponentTs...>, ConstraintConfigTupleT>> {
  using type = MixedHandledTable<DataT, ComponentHandle<ComponentTs>...>;
};

template <class DataT, class ConstraintGraphT>
struct ConstraintHandledTableFromConstraintGraph {};

template <class DataT, class ComponentDataTupleT, class... ConstraintConfigTs>
struct ConstraintHandledTableFromConstraintGraph<
    DataT,
    ConstraintGraph<ComponentDataTupleT, std::tuple<ConstraintConfigTs...>>> {
  using type =
      MixedHandledTable<DataT,
                        ConstraintHandle<typename MemberTypesOfConstraintConfig<
                            ConstraintConfigTs>::ConstraintData>...>;
};

template <class DataT, class ComponentT, class... ComponentTs,
          class ConstraintConfigTupleT>
inline HandledTable<ComponentHandle<ComponentT>, DataT>
MakeHandledTableForComponents(const ConstraintGraph<std::tuple<ComponentTs...>,
                                                    ConstraintConfigTupleT> &cg,
                              const DataT &defaultValue = DataT()) {
  return HandledTable<ComponentHandle<ComponentT>, DataT>(
      cg.internalComponents<ComponentT>().size(), defaultValue);
}

template <class DataT, class ConstraintT, class... ComponentTs,
          class ConstraintConfigTupleT>
inline HandledTable<ConstraintHandle<ConstraintT>, DataT>
MakeHandledTableForConstraints(
    const ConstraintGraph<std::tuple<ComponentTs...>, ConstraintConfigTupleT>
        &cg,
    const DataT &defaultValue = DataT()) {
  return HandledTable<ConstraintHandle<ConstraintT>, DataT>(
      cg.internalConstraints<ConstraintT>().size(), defaultValue);
}

template <class DataT, class... ComponentTs, class ConstraintConfigTupleT>
inline MixedHandledTable<DataT, ComponentHandle<ComponentTs>...>
MakeHandledTableForAllComponents(
    const ConstraintGraph<std::tuple<ComponentTs...>, ConstraintConfigTupleT>
        &cg,
    const DataT &d = DataT()) {
  return MixedHandledTable<DataT, ComponentHandle<ComponentTs>...>(
      {cg.internalComponents<ComponentTs>().size()...}, d);
}

template <class DataT, class ComponentDataTupleT, class... ConstraintConfigTs>
inline MixedHandledTable<
    DataT, ConstraintHandle<typename MemberTypesOfConstraintConfig<
               ConstraintConfigTs>::ConstraintData>...>
MakeHandledTableForAllConstraints(
    const ConstraintGraph<ComponentDataTupleT,
                          std::tuple<ConstraintConfigTs...>> &cg,
    const DataT &d = DataT()) {
  return MixedHandledTable<
      DataT, ConstraintHandle<typename MemberTypesOfConstraintConfig<
                 ConstraintConfigTs>::ConstraintData>...>(
      {cg.internalConstraints<typename MemberTypesOfConstraintConfig<
          ConstraintConfigTs>::ConstraintData>()
           .size()...},
      d);
}
}
}
