#ifndef PANORAMIX_CORE_CONS_GRAPH_HPP
#define PANORAMIX_CORE_CONS_GRAPH_HPP
 
#include <numeric>
#include "handle.hpp" 

namespace panoramix {
    namespace core {


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

            template <class T>
            struct TempContainer {
                template <class ... Ts>
                inline TempContainer(Ts && ... elements) : data({elements ...}) {}
                template <int N>
                inline operator std::array<T, N>() const {
                    std::array<T, N> arr;
                    std::copy(data.begin(), data.end(), arr.begin());
                    return arr;
                }
                template <class AllocT>
                inline operator std::vector<T, AllocT> () {
                    return std::move(data);
                }
                inline typename std::vector<T>::const_iterator begin() const { return data.begin(); }
                inline typename std::vector<T>::const_iterator end() const { return data.end(); }
                std::vector<T> data;
            };

        }

        template <class HandleT, class ... HandleTs>
        inline TempContainer<HandleT> Depends(HandleT h, HandleTs ... hs) {
            return TempContainer<HandleT>(h, hs ... );
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
                const typename ComponentHandlesArrayFromComponentOccupation<ComponentTagT, ComponentOccupationTs>::type & ... handleArrays)
                : hd(id), 
                allComponents(std::forward_as_tuple(handleArrays...)) {}
           /* template <class ... ComponentHandleArrayTs>
            explicit inline ConstraintTopo(int id, ComponentHandleArrayTs && ... handleArrays) 
                : hd(id), allComponents(std::forward_as_tuple(handleArrays...)) {}*/

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

        public:
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
                bool dummy[] = { insertConstraintHandleInComponent(topo.hd, comps) ... };
                return topo.hd;
            }

            // removed ?
            template <class DataT>
            inline bool removed(ComponentHandle<DataT> h) const {
                enum { _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                using _TripletType = typename std::tuple_element<_idx, ComponentsTripletArrayTuple>::type::value_type;
                using _TopoType = typename _TripletType::TopoType;
                auto & correspondingComponents = std::get<_idx>(_components);
                if (h.invalid())
                    return true;
                return !correspondingComponents.at(h.id).exists;
            }

            template <class DataT>
            inline bool removed(ConstraintHandle<DataT> h) const {
                enum { _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                using _TripletType = typename std::tuple_element<_idx, ConstraintsTripletArrayTuple>::type::value_type;
                using _TopoType = typename _TripletType::TopoType;
                auto & correspondingConstraints = std::get<_idx>(_constraints);
                if (h.invalid())
                    return true;
                return !correspondingConstraints.at(h.id).exists;
            }


            template <class DataT>
            inline void remove(ConstraintHandle<DataT> h){
                enum { _idx = TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                using _TripletType = typename std::tuple_element<_idx, ConstraintsTripletArrayTuple>::type::value_type;
                using _TopoType = typename _TripletType::TopoType;

                if (removed(h))
                    return;

                auto & correspondingConstraints = std::get<_idx>(_constraints);
                // remove constraint handles in each component
                removeConstraintHandleInRelatedComponentsUsingSequence(h,
                    typename SequenceGenerator<std::tuple_size<typename _TopoType::ComponentDataTuple>::value>::type());

                correspondingConstraints[h.id].exists = false;
            }

            template <class DataT>
            void remove(ComponentHandle<DataT> h){
                enum { _idx = TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value };
                static_assert(_idx >= 0, "DataT is not supported here!");
                using _TripletType = typename std::tuple_element<_idx, ComponentsTripletArrayTuple>::type::value_type;
                using _TopoType = typename _TripletType::TopoType;

                if (removed(h))
                    return;

                auto & correspondingComponents = std::get<_idx>(_components);
                correspondingComponents[h.id].exists = false;

                // remove related constraints
                removeHandlesInContainerTupleUsingSequence(topo(h).allConstraints,
                    typename SequenceGenerator<std::tuple_size<typename _TopoType::ConstraintDataTuple>::value>::type());
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

            // all info
            
            template <class DataT>
            inline typename ConditionalContainerTypeFromTripletArrayType<typename std::tuple_element<TypeFirstLocationInTuple<DataT, 
                ComponentDataTupleType>::value, ComponentsTripletArrayTuple>::type>::const_type components() const{
                return MakeConditionalContainer(std::get<TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value>(_components));
            }

            template <class DataT>
            inline typename ConditionalContainerTypeFromTripletArrayType<typename std::tuple_element<TypeFirstLocationInTuple<DataT, 
                ComponentDataTupleType>::value, ComponentsTripletArrayTuple>::type>::type components() {
                return MakeConditionalContainer(std::get<TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value>(_components));
            }

            template <class DataT>
            inline const typename std::tuple_element<TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value,
                ComponentsTripletArrayTuple>::type & internalComponents() const {
                return std::get<TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value>(_components);
            }

            template <class DataT>
            inline typename std::tuple_element<TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value,
                ComponentsTripletArrayTuple>::type & internalComponents() {
                return std::get<TypeFirstLocationInTuple<DataT, ComponentDataTupleType>::value>(_components);
            }

            inline const ComponentsTripletArrayTuple & allComponents() const { return _components; }
            inline ComponentsTripletArrayTuple & allComponents() { return _components; }


            template <class DataT>
            inline typename ConditionalContainerTypeFromTripletArrayType<typename std::tuple_element<TypeFirstLocationInTuple<DataT, 
                ConstraintDataTupleType>::value, ConstraintsTripletArrayTuple>::type>::const_type constraints() const{
                return MakeConditionalContainer(std::get<TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value>(_constraints));
            }

            template <class DataT>
            inline typename ConditionalContainerTypeFromTripletArrayType<typename std::tuple_element<TypeFirstLocationInTuple<DataT,
                ConstraintDataTupleType>::value, ConstraintsTripletArrayTuple>::type>::type constraints() {
                return MakeConditionalContainer(std::get<TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value>(_constraints));
            }

            template <class DataT>
            inline const typename std::tuple_element<TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value, 
                ConstraintsTripletArrayTuple>::type & internalConstraints() const {
                return std::get<TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value>(_constraints);
            }

            template <class DataT>
            inline typename std::tuple_element<TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value, 
                ConstraintsTripletArrayTuple>::type & internalConstraints() {
                return std::get<TypeFirstLocationInTuple<DataT, ConstraintDataTupleType>::value>(_constraints);
            }

            inline const ConstraintsTripletArrayTuple & allConstraints() const { return _constraints; }
            inline ConstraintsTripletArrayTuple & allConstraints() { return _constraints; }


            template <class Archiver> inline void serialize(Archiver & ar) { ar(_components, _constraints); }
        
        private:
            template <class DataT, class ComponentHandleArrayT>
            inline bool insertConstraintHandleInComponent(ConstraintHandle<DataT> consH, ComponentHandleArrayT && compArr){
                for (auto && compH : compArr){
                    topo(compH).constraints<DataT>().insert(consH);
                }
                return true;
            }
            
            template <class DataT, class ComponentHandleArrayT>
            inline bool removeConstraintHandleInComponent(ConstraintHandle<DataT> consH, ComponentHandleArrayT & compArr){
                for (auto && compH : compArr){
                    if (compH.invalid() || removed(compH))
                        continue;
                    topo(compH).constraints<DataT>().erase(consH);
                }
                return true;
            }

            template <class DataT, int ... Is>
            inline bool removeConstraintHandleInRelatedComponentsUsingSequence(ConstraintHandle<DataT> consH, Sequence<Is...>){
                bool d[] = { removeConstraintHandleInComponent(consH, std::get<Is>(topo(consH).allComponents))... };
                return true;
            }

            template <class ContainerT>
            inline bool removeHandlesInContainer(ContainerT && hs){
                for (auto && h : hs){
                    remove(h);
                }
                return true;
            }

            template <class ContainerTupleT, int ... Is>
            inline bool removeHandlesInContainerTupleUsingSequence(ContainerTupleT & hs, Sequence<Is...>){
                bool d[] = { removeHandlesInContainer(std::get<Is>(hs))... };
                return true;
            }
        
        private:
            ComponentsTripletArrayTuple _components;
            ConstraintsTripletArrayTuple _constraints;
        };
    	
    }
}
 
#endif