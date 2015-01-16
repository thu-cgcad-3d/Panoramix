#ifndef PANORAMIX_CORE_CONS_GRAPH_HPP
#define PANORAMIX_CORE_CONS_GRAPH_HPP
 
#include "handle.hpp" 

namespace panoramix {
    namespace core {


        // used to register all variables into one vector or 
        // or register all constraints into one sparse matrix (represented as a vector of nonzero triplets)
        template <class IdentityT, class ValueT>
        class Registration {
        public:
            inline Registration() {}
            template <class IdentityIteratorT, class GetValuesNumFunT>
            Registration(IdentityIteratorT begin, IdentityIteratorT end,
                GetValuesNumFunT getValuesNum, const ValueT & defaultValue = ValueT()){
                size_t sz = 0;
                size_t identityNum = 0;
                for (auto i = begin; i != end; ++i){
                    sz += getValuesNum(*i);
                    identityNum++;
                }
                _values.reserve(sz);
                _valueIdentities.reserve(sz);
                _startPositions.reserve(identityNum);
                while (begin != end){
                    auto & identity = *begin;
                    _startPositions[identity] = _values.size();
                    _values.insert(_values.end(), getValuesNum(*i), defaultValue);
                    _valueIdentities.insert(_valueIdentities.end(), getValuesNum(*i), identity);
                    ++begin;
                }
            }
            template <class IdentityIteratorT, class GetValuesNumFunT, class GetValuesFunT>
            Registration(IdentityIteratorT begin, IdentityIteratorT end,
                GetValuesNumFunT getValuesNum, GetValuesFunT getValues){
                size_t sz = 0;
                size_t identityNum = 0;
                for (auto i = begin; i != end; ++i){
                    sz += getValuesNum(*i);
                    identityNum++;
                }
                _values.reserve(sz);
                _valueIdentities.reserve(sz);
                _startPositions.reserve(identityNum);
                while (begin != end){
                    auto & identity = *begin;
                    _startPositions[identity] = _values.size();
                    for (auto && v : getValues(identity)){
                        _values.push_back(v);
                        _valueIdentities.push_back(identity);
                    }
                    ++begin;
                }
            }

        public:
            size_t size() const { return _values.size(); }

            typename std::vector<ValueT>::iterator begin() { return _values.begin(); }
            typename std::vector<ValueT>::iterator end() { return _values.end(); }

            typename std::vector<ValueT>::const_iterator begin() const { return _values.begin(); }
            typename std::vector<ValueT>::const_iterator end() const { return _values.end(); }

            typename std::vector<ValueT>::const_iterator cbegin() const { return _values.cbegin(); }
            typename std::vector<ValueT>::const_iterator cend() const { return _values.cend(); }

            const ValueT & operator[](size_t index) const { return _values[index]; }
            ValueT & operator[](size_t index) { return _values[index]; }

            const std::vector<ValueT> & values() const { return _values; }

            const IdentityT & identity(size_t index) const { return _valueIdentities.at(index); }
            const size_t startPositionOfIdentity(const IdentityT & id) const { return _startPositions.at(id); }
            size_t orderInIdentity(size_t index) const { return index - startPositionOfIdentity(identity(index)); }

            const ValueT & operator()(const IdentityT & id, size_t nthValue) const {
                return _values[startPositionOfIdentity(id) + nthValue];
            }
            ValueT & operator()(const IdentityT & id, size_t nthValue) {
                return _values[startPositionOfIdentity(id) + nthValue];
            }

            bool empty() const { return _values.empty(); }
            void clear() {
                _values.clear();
                _valueIdentities.clear();
                _startPositions.clear();
            }

            Registration & concat(const Registration & v) {
                size_t oldSize = _values.size();
                _values.insert(_values.end(), v._values.begin(), v._values.end());
                _valueIdentities.insert(_valueIdentities.end(), v._valueIdentities.begin(), v._valueIdentities.end());
                for (auto && idPos : v._startPositions){
                    assert(!Contains(_startPositions, idPos.first));
                    _startPositions[idPos.first] = idPos.second + oldSize;
                }
                return *this;
            }

        public:
            template <class VectorT>
            void convertTo(VectorT & v) const {
                v.resize(_values.size());
                for (size_t i = 0; i < _values.size(); i++){
                    v(i) = _values[i];
                }
            }

        private:
            std::vector<ValueT> _values;
            std::vector<IdentityT> _valueIdentities;
            std::unordered_map<IdentityT, size_t> _startPositions;
        };









        // configuration for constraint
        template <class ComponentDataT, int N>
        struct ComponentOccupation {};

        template <class ComponentOccupationT>
        struct IsComponentOccupation : std::false_type {};

        template <class ComponentDataT, int N>
        struct IsComponentOccupation<ComponentOccupation<ComponentDataT, N>> : std::true_type{};

        template <class ConstraintDataT, class ...ComponentOccupationTs>
        struct ConstraintConfig {};

        template <class T>
        struct IsConstraintConfig : std::false_type {};

        template <class ...ComponentOccupationTs>
        struct IsConstraintConfig<ConstraintConfig<ComponentOccupationTs...>> : std::true_type{};

        namespace {

            // get member type of component occupation
            template <class ComponentOccupationT>
            struct ComponentDataTypeFromComponentOccupation {};

            template <class ComponentDataT, int N>
            struct ComponentDataTypeFromComponentOccupation<ComponentOccupation<ComponentDataT, N>> {
                using type = ComponentDataT;
            };

            // get member type of constraint config
            template <class ConstraintConfigT>
            struct ConstraintDataTypeFromConstraintConfig {};

            template <class ConstraintDataType, class ...ComponentOccupationTs>
            struct ConstraintDataTypeFromConstraintConfig<ConstraintConfig<ConstraintDataType, ComponentOccupationTs...>> {
                using type = ConstraintDataType;
            };


            // constraint data tuple from constraint config tuple
            template <class ConstraintConfigTupleT>
            struct ConstraintDataTupleFromConstraintConfigTuple {};

            template <class ... ConstraintConfigTs>
            struct ConstraintDataTupleFromConstraintConfigTuple<std::tuple<ConstraintConfigTs ...>> {
                using type = std::tuple<typename ConstraintDataTypeFromConstraintConfig<ConstraintConfigTs>::type...>;
            };

            // component handles array
            template <class ComponentTagT, class ComponentOccupationT>
            struct ComponentHandlesArrayFromComponentOccupation {};

            template <class ComponentTagT, class DataT, int N>
            struct ComponentHandlesArrayFromComponentOccupation<ComponentTagT, ComponentOccupation<DataT, N>> {
                using type = std::array<HandleOfType<ComponentTagT, DataT>, N>;
            };

            template <class ComponentTagT, class DataT>
            struct ComponentHandlesArrayFromComponentOccupation<ComponentTagT, ComponentOccupation<DataT, Dynamic>> {
                using type = std::vector<HandleOfType<ComponentTagT, DataT>>;
            };

        }


        template <int N, class HandleT, class ... HandleTs, class = std::enable_if_t<N != Dynamic>>
        inline std::array<HandleT, N> Depends(HandleT h, HandleTs ... hs) {
            return { h, hs ... };
        }

        template <int N, class HandleT, class ... HandleTs, class = std::enable_if_t<N == Dynamic>, class = void>
        inline std::vector<HandleT> Depends(HandleT h, HandleTs ... hs) {
            return { h, hs ... };
        }


        // make a layer content tuple with 2 levels
        template <class ComponentTagT, class ConstraintTagT, class ComponentDataT, class ... ConstraintConfigTs>
        struct ComponentTopo {
            using ConstraintDataTuple = std::tuple<typename ConstraintDataTypeFromConstraintConfig<ConstraintConfigTs>::type ...>;
            std::tuple<std::set<HandleOfType<ConstraintTagT, typename ConstraintDataTypeFromConstraintConfig<ConstraintConfigTs>::type>>...> allConstraints;
            
            HandleOfType<ComponentTagT, ComponentDataT> hd;
            explicit inline ComponentTopo(int id = -1) : hd(id){}
            
            template <class ConstraintDataT>
            inline std::set<HandleOfType<ConstraintTagT, ConstraintDataT>> & constraints() {
                enum { _idx = TypeFirstLocationInTuple<ConstraintDataT, ConstraintDataTuple>::value };
                static_assert(_idx >= 0, "Invalid ConstraintDataT");
                return std::get<_idx>(allConstraints);
            }

            template <class ConstraintDataT>
            inline const std::set<HandleOfType<ConstraintTagT, ConstraintDataT>> & constraints() const {
                enum { _idx = TypeFirstLocationInTuple<ConstraintDataT, ConstraintDataTuple>::value };
                static_assert(_idx >= 0, "Invalid ConstraintDataT");
                return std::get<_idx>(allConstraints);
            }

            template <class Archive> inline void serialize(Archive & ar) { ar(allConstraints, hd); }
        };

        template <class ComponentTagT, class ConstraintTagT, class ConstraintDataT, class ... ComponentOccupationTs>
        struct ConstraintTopo {
            static_assert(sizeof...(ComponentOccupationTs) > 0, "ComponentOccupationTs must be more than zero");
            using ComponentDataTuple = std::tuple<typename ComponentDataTypeFromComponentOccupation<ComponentOccupationTs>::type ...>;
            using ComponentHandlesArrayTuple = std::tuple<typename ComponentHandlesArrayFromComponentOccupation<ComponentTagT, ComponentOccupationTs>::type ...>;
            std::tuple<typename ComponentHandlesArrayFromComponentOccupation<ComponentTagT, ComponentOccupationTs>::type ...> allComponents;

            HandleOfType<ConstraintTagT, ConstraintDataT> hd;
            explicit inline ConstraintTopo() : hd(-1){}
            explicit inline ConstraintTopo(int id, 
                typename ComponentHandlesArrayFromComponentOccupation<ComponentTagT, ComponentOccupationTs>::type && ... handleArrays)
                : hd(id), 
                allComponents(std::forward<typename ComponentHandlesArrayFromComponentOccupation<ComponentTagT, ComponentOccupationTs>::type>(handleArrays) ...) {}

            template <class ComponentDataT, int Idx = TypeFirstLocationInTuple<ComponentDataT, ComponentDataTuple>::value>
            inline typename std::tuple_element<Idx, ComponentHandlesArrayTuple>::type & components() {
                static_assert(Idx >= 0, "Invalid ComponentDataT");
                return std::get<Idx>(allConstraints);
            }

            template <class ComponentDataT, int Idx = TypeFirstLocationInTuple<ComponentDataT, ComponentDataTuple>::value>
            inline const typename std::tuple_element<Idx, ComponentHandlesArrayTuple>::type & components() const {
                static_assert(Idx >= 0, "Invalid ComponentDataT");
                return std::get<Idx>(allConstraints);
            }

            template <class Archive> inline void serialize(Archive & ar) { ar(allComponents, hd); }
        };


        namespace {

            template <class ComponentTagT, class ConstraintTagT, class ComponentDataT, class ConstraintConfigTupleT>
            struct ComponentTripletArrayFromConstraintConfigTuple {};

            template <class ComponentTagT, class ConstraintTagT, class ComponentDataT, class ... ConstraintConfigTs>
            struct ComponentTripletArrayFromConstraintConfigTuple<ComponentTagT, ConstraintTagT, ComponentDataT, std::tuple<ConstraintConfigTs...>> {
                using type = TripletArray<ComponentTopo<ComponentTagT, ConstraintTagT, ComponentDataT, ConstraintConfigTs...>, ComponentDataT>;
            };

            template <class ComponentTagT, class ConstraintTagT, class ComponentDataTupleT, class ConstraintConfigTupleT>
            struct ComponentTripletArrayTupleFromConstraintConfigTuple {};

            template <class ComponentTagT, class ConstraintTagT, class ConstraintConfigTupleT, class ... ComponentDataTs>
            struct ComponentTripletArrayTupleFromConstraintConfigTuple<ComponentTagT, ConstraintTagT, std::tuple<ComponentDataTs...>, ConstraintConfigTupleT> {
                using type = std::tuple<typename ComponentTripletArrayFromConstraintConfigTuple<ComponentTagT, ConstraintTagT, ComponentDataTs, ConstraintConfigTupleT>::type ...>;
            };

            template <class ComponentTagT, class ConstraintTagT, class ConstraintConfigT>
            struct ConstraintTripletArrayFromConstraintConfig {};

            template <class ComponentTagT, class ConstraintTagT, class ConstraintDataT, class ...ComponentOccupationTs>
            struct ConstraintTripletArrayFromConstraintConfig<ComponentTagT, ConstraintTagT, ConstraintConfig<ConstraintDataT, ComponentOccupationTs...>> {
                using type = TripletArray<ConstraintTopo<ComponentTagT, ConstraintTagT, ConstraintDataT, ComponentOccupationTs...>, ConstraintDataT>;
            };

            template <class ComponentTagT, class ConstraintTagT, class ConstraintConfigTupleT>
            struct ConstraintTripetArrayTupleFromConstraintConfigTuple {};

            template <class ComponentTagT, class ConstraintTagT, class ... ConstraintConfigTs>
            struct ConstraintTripetArrayTupleFromConstraintConfigTuple<ComponentTagT, ConstraintTagT, std::tuple<ConstraintConfigTs...>> {
                using type = std::tuple<typename ConstraintTripletArrayFromConstraintConfig<ComponentTagT, ConstraintTagT, ConstraintConfigTs>::type...>;
            };

        }



        struct ComponentTag {};
        struct ConstraintTag {};

        template <class DataT>
        using ComponentHandle = HandleOfType<ComponentTag, DataT>;

        template <class DataT>
        using ConstraintHandle = HandleOfType<ConstraintTag, DataT>;

        template <class HandleOfTypeT>
        struct TagOfHandleOfType {};

        template <class TagT, class DataT>
        struct TagOfHandleOfType<HandleOfType<TagT, DataT>> {
            using type = TagT;
        };


        template <class ComponentDataTupleT, class ConstraintConfigTupleT>
        class ConstraintGraph {

            static_assert(IsTuple<ComponentDataTupleT>::value, "ComponentDataTupleT must be std::tuple");
            static_assert(IsTuple<ConstraintConfigTupleT>::value, "ConstraintConfigTupleT must be std::tuple");
            using ComponentsTripletArrayTuple =
                typename ComponentTripletArrayTupleFromConstraintConfigTuple<ComponentTag, ConstraintTag,
                ComponentDataTupleT, ConstraintConfigTupleT>::type;
            using ConstraintsTripletArrayTuple =
                typename ConstraintTripetArrayTupleFromConstraintConfigTuple<ComponentTag, ConstraintTag,
                ConstraintConfigTupleT>::type;

            using ComponentDataTupleType = ComponentDataTupleT;
            using ConstraintDataTupleType = typename ConstraintDataTupleFromConstraintConfigTuple<ConstraintConfigTupleT>::type;

        public:

            template <class DataT>
            inline ComponentHandle<std::decay_t<DataT>> addComponent(DataT && data) {
                using _DataType = std::decay_t<DataT>;
                enum { _idx = TypeFirstLocationInTuple<_DataType, ComponentDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                using _TripletType = typename std::tuple_element<_idx, ComponentsTripletArrayTuple>::type::value_type;
                using _TopoType = typename _TripletType::TopoType;
                auto & correspondingComponents = std::get<_idx>(_components);
                _TopoType topo;
                topo.hd = ComponentHandle<_DataType>(correspondingComponents.size());
                correspondingComponents.emplace_back(std::move(topo), std::forward<DataT>(data), true);
                return topo.hd;
            }


            template <class DataT, class ... ComponentHandleArrayTs>
            ConstraintHandle<std::decay_t<DataT>> addConstraint(DataT && data, ComponentHandleArrayTs && ... comps) {
                using _DataType = std::decay_t<DataT>;
                enum { _idx = TypeFirstLocationInTuple<_DataType, ConstraintDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                using _TripletType = typename std::tuple_element<_idx, ConstraintsTripletArrayTuple>::type::value_type;
                using _TopoType = typename _TripletType::TopoType;
               
                auto & correspondingConstraints = std::get<_idx>(_constraints);
                _TopoType topo(correspondingConstraints.size(), std::forward<ComponentHandleArrayTs>(comps) ...);
                correspondingConstraints.emplace_back(std::move(topo), std::forward<DataT>(data), true);
                
                // update constraint handles in each component
                bool dummy[] = { updateConstraintHandleInComponent(topo.hd, comps) ... };
                return topo.hd;
            }

            // component info
            template <class DataT>
            inline const DataT & data(ComponentHandle<DataT> hd) const {
                enum { _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                return std::get<_idx>(_components).at(hd.id).data;
            }

            template <class DataT>
            inline DataT & data(ComponentHandle<DataT> hd){
                enum { _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                return std::get<_idx>(_components)[hd.id].data;
            }

            template <class DataT>
            inline const typename std::tuple_element<TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value, 
                ComponentsTripletArrayTuple>::type::value_type::TopoType & topo(ComponentHandle<DataT> hd) const {
                enum {_idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value};
                static_assert(_idx >= 0, "DataT is not supported here!");
                return std::get<_idx>(_components).at(hd.id).topo;
            }

            template <class DataT>
            inline typename std::tuple_element<TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value, 
                ComponentsTripletArrayTuple>::type::value_type::TopoType & topo(ComponentHandle<DataT> hd) {
                enum { _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                return std::get<_idx>(_components).at(hd.id).topo;
            }

            // constraint info
            template <class DataT>
            inline const DataT & data(ConstraintHandle<DataT> hd) const {
                enum { _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                return std::get<_idx>(_constraints).at(hd.id).data;
            }

            template <class DataT>
            inline DataT & data(ConstraintHandle<DataT> hd) {
                enum { _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                return std::get<_idx>(_constraints)[hd.id].data;
            }

            template <class DataT>
            inline const typename std::tuple_element<TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value, 
                ConstraintsTripletArrayTuple>::type::value_type::TopoType & topo(ConstraintHandle<DataT> hd) const {
                enum { _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                return std::get<_idx>(_constraints).at(hd.id).topo;
            }

            template <class DataT>
            inline typename std::tuple_element<TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value, 
                ConstraintsTripletArrayTuple>::type::value_type::TopoType & topo(ConstraintHandle<DataT> hd) {
                enum { _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                return std::get<_idx>(_constraints)[hd.id].topo;
            }


            template <class Archiver> inline void serialize(Archiver & ar) { ar(_components, _constraints); }
        
        private:
            template <class DataT, class ComponentHandleArrayT>
            inline bool updateConstraintHandleInComponent(HandleOfType<ConstraintTag, DataT> consH, ComponentHandleArrayT && compArr){
                for (auto && compH : compArr){
                    topo(compH).constraints<DataT>().insert(consH);
                }
                return true;
            }
        
        
        private:
            ComponentsTripletArrayTuple _components;
            ConstraintsTripletArrayTuple _constraints;
        };
    	
    }
}
 
#endif