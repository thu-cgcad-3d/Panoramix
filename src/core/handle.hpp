#ifndef PANORAMIX_CORE_HANDLE_HPP
#define PANORAMIX_CORE_HANDLE_HPP

#include <vector>
#include <cstdint>
#include <array>
#include <set>

#include "meta.hpp"
#include "iterators.hpp"

namespace panoramix {
    namespace core {


        static const int Dynamic = -1;


        template <class Tag>
        struct Handle {
            int id;
            inline Handle(int id_ = -1) : id(id_){}
            inline bool operator == (Handle h) const { return id == h.id; }
            inline bool operator != (Handle h) const { return id != h.id; }
            inline void reset() { id = -1; }
            inline bool valid() const { return id >= 0; }
            inline bool invalid() const { return id < 0; }
            template <class Archive> inline void serialize(Archive & ar) { ar(id); }
        };

        template <class ValueT, int N = Dynamic>
        struct SizedContainer {
            using type = std::array<ValueT, N>;
        };

        template <class ValueT>
        struct SizedContainer<ValueT, Dynamic> {
            using type = std::vector<ValueT>;
        };

        template <class ValueT, int N = Dynamic>
        using SizedContainerType = typename SizedContainer<ValueT, N>::type;


        template <class Tag>
        using HandleArray = std::vector<Handle<Tag>>;
        template <class Tag>
        using HandlePtrArray = std::vector<Handle<Tag>*>;

        template <class Tag>
        inline bool operator < (const Handle<Tag> & a, const Handle<Tag> & b){
            return a.id < b.id;
        }


        // is handle ?
        template <class T>
        struct IsHandle : no {};

        template <class Tag>
        struct IsHandle<Handle<Tag>> : yes{};

        // decorated handle types

        template <int L>
        struct AtLevel {
            static const int Level = L;
        };
        template <int L>
        using HandleAtLevel = Handle<AtLevel<L>>;


        template <class TypeTag, class Tag>
        struct OfType {};

        template <class TypeTag, class Tag>
        using HandleOfType = Handle<OfType<TypeTag, Tag>>;

        template <class TypeTag, int L>
        using HandleOfTypeAtLevel = Handle<OfType<TypeTag, AtLevel<L>>>;


        // handled table
        template <class HandleT, class DataT, class ContainerT = std::vector<DataT>>
        struct HandledTable {
            static_assert(IsHandle<HandleT>::value, "HandleT must be a Handle!");
            HandledTable() {}
            explicit HandledTable(size_t maxSize) : data(maxSize) {}
            HandledTable(size_t maxSize, const DataT & d) : data(maxSize, d) {}

            HandledTable(HandledTable && t) : data(std::move(t.data)) {}
            HandledTable & operator = (HandledTable && t) { data = std::move(t.data); return *this; }

            void resize(size_t sz) { data.resize(sz); }
            const DataT & operator[](HandleT h) const { return data[h.id]; }
            const DataT & at(HandleT h) const { return data.at(h.id); }
            DataT & operator[](HandleT h) { return data[h.id]; }

            using value_type = typename ContainerT::value_type;
            struct iterator {
                ContainerT & cont;
                typename ContainerT::iterator i;
                HandleT hd() const { return HandleT(std::distance(std::begin(cont), i)); }
                value_type & operator * () const { return *i; }
                iterator & operator ++ () { ++i; return *this; }
                bool operator == (iterator it) const { return i == it.i; }
                bool operator != (iterator it) const { return i != it.i; }
            };

            struct const_iterator {
                const ContainerT & cont;
                typename ContainerT::const_iterator i;
                HandleT hd() const { return HandleT(std::distance(std::begin(cont), i)); }
                const value_type & operator * () const { return *i; }
                const_iterator & operator ++ () { ++i; return *this; }
                bool operator == (const_iterator it) const { return i == it.i; }
                bool operator != (const_iterator it) const { return i != it.i; }
            };
            
            inline iterator begin() { return iterator{ data, data.begin() }; }
            inline iterator end() { return iterator{ data, data.end() }; }
            inline const_iterator begin() const { return const_iterator{ data, data.begin() }; }
            inline const_iterator end() const { return const_iterator{ data, data.end() }; }

            template <class MappingFunT>
            inline HandledTable<HandleT, std::decay_t<typename FunctionTraits<MappingFunT>::ResultType>> map(MappingFunT && fun) const {
                HandledTable<HandleT, std::decay_t<typename FunctionTraits<MappingFunT>::ResultType>> result(data.size());
                for (int i = 0; i < data.size(); i++){
                    result.data[i] = fun(data[i]);
                }
                return result;
            }

            template <class Archive> inline void serialize(Archive & ar) { ar(data); }
            ContainerT data;
        };


        // mixed handled table
        template <class DataT, class ... HandleTs>
        struct MixedHandledTable {
            static const int HandlesNum = sizeof...(HandleTs);
            using HandleTupleType = std::tuple<HandleTs...>;

            MixedHandledTable() {}

            explicit MixedHandledTable(std::initializer_list<size_t> maxSizes) {
                assert(maxSizes.size() <= HandlesNum);
                int i = 0;
                for (size_t sz : maxSizes){
                    data[i++] = std::vector<DataT>(sz);
                }
            }

            MixedHandledTable(std::initializer_list<size_t> maxSizes, const DataT & d) {
                assert(maxSizes.size() <= HandlesNum);
                int i = 0;
                for (size_t sz : maxSizes){
                    data[i++] = std::vector<DataT>(sz, d);
                }
            }

            template <class ... IntTs>
            explicit MixedHandledTable(const std::tuple<IntTs...> & maxSizes)
                : data({ { std::vector<DataT>(std::get<TypeFirstLocationInTuple<HandleTs, HandleTupleType>::value>(maxSizes)) } }) {}

            MixedHandledTable(MixedHandledTable && t) : data(std::move(t.data)) {}
            MixedHandledTable & operator = (MixedHandledTable && t) { data = std::move(t.data); return *this; }

            template <class HandleT>
            const DataT & operator[](HandleT h) const {
                return data[TypeFirstLocationInTuple<HandleT, HandleTupleType>::value][h.id];
            }
            template <class HandleT>
            const DataT & at(HandleT h) const {
                return data[TypeFirstLocationInTuple<HandleT, HandleTupleType>::value].at(h.id);
            }
            template <class HandleT>
            DataT & operator[](HandleT h) {
                return data[TypeFirstLocationInTuple<HandleT, HandleTupleType>::value][h.id];
            }
            std::array<DataT &, sizeof...(HandleTs)> at(size_t i) {
                return std::array<DataT &, sizeof...(HandleTs)>{{ data[TypeFirstLocationInTuple<HandleTs, HandleTupleType>::value][i]... }};
            }
            std::array<const DataT &, sizeof...(HandleTs)> at(size_t i) const {
                return std::array<DataT &, sizeof...(HandleTs)>{{ data[TypeFirstLocationInTuple<HandleTs, HandleTupleType>::value][i]... }};
            }
            template <class HandleT>
            std::vector<DataT> & dataOfType() {
                return data[TypeFirstLocationInTuple<HandleT, HandleTupleType>::value];
            }


            template <class Archive> inline void serialize(Archive & ar) { ar(data); }
            std::array<std::vector<DataT>, sizeof...(HandleTs)> data;
        };





        // triplet 
        template <class TopoT, class DataT>
        struct Triplet {
            using TopoType = TopoT;
            using DataType = DataT;

            TopoT topo;
            bool exists;
            DataT data;
            inline Triplet(){}
            inline Triplet(const TopoT & t, const DataT & d, bool e = true)
                : topo(t), exists(e), data(d){}
            inline Triplet(const TopoT & t, DataT && d, bool e = true)
                : topo(t), exists(e), data(std::forward<DataT>(d)) {}
            template <class Archive> inline void serialize(Archive & ar) { ar(topo, exists, data); }
        };
        template <class TopoT, class DataT>
        struct TripletExistsPred {
            inline bool operator()(const Triplet<TopoT, DataT> & t) const {
                return t.exists;
            }
        };
        template <class TripletT>
        struct TripletExistsPredFromTripletType {};
        template <class TopoT, class DataT>
        struct TripletExistsPredFromTripletType<Triplet<TopoT, DataT>> {
            using type = TripletExistsPred<TopoT, DataT>;
        };

        template <class TopoT, class DataT>
        using TripletArray = std::vector<Triplet<TopoT, DataT>>;
        template <class TripletArrayT>
        struct ConditionalContainerTypeFromTripletArrayType {};
        template <class TopoT, class DataT>
        struct ConditionalContainerTypeFromTripletArrayType<TripletArray<TopoT, DataT>> {
            using type = ConditionalContainerWrapper<TripletArray<TopoT, DataT>, TripletExistsPred<TopoT, DataT>>;
            using const_type = ConstConditionalContainerWrapper<TripletArray<TopoT, DataT>, TripletExistsPred<TopoT, DataT>>;
        };

        template <class TopoT, class DataT>
        inline ConditionalContainerWrapper<TripletArray<TopoT, DataT>, TripletExistsPred<TopoT, DataT>> 
            MakeConditionalContainer(TripletArray<TopoT, DataT> & arr) {
            return ConditionalContainerWrapper<TripletArray<TopoT, DataT>, TripletExistsPred<TopoT, DataT>>(&arr, TripletExistsPred<TopoT, DataT>());
        }
        template <class TopoT, class DataT>
        inline ConstConditionalContainerWrapper<TripletArray<TopoT, DataT>, TripletExistsPred<TopoT, DataT>> 
            MakeConditionalContainer(const TripletArray<TopoT, DataT> & arr) {
            return ConstConditionalContainerWrapper<TripletArray<TopoT, DataT>, TripletExistsPred<TopoT, DataT>>(&arr, TripletExistsPred<TopoT, DataT>());
        }


        template <class TopoT, class DataT>
        inline auto BoundingBox(const Triplet<TopoT, DataT> & t) -> decltype(BoundingBox(t.data)) {
            return BoundingBox(t.data);
        }



        // helper functions
        template <class ComponentTableT, class UpdateHandleTableT>
        int RemoveAndMap(ComponentTableT & v, UpdateHandleTableT & newlocations) {
            // ComponentTableT : std::vector<Triplet<TopoT, DataT>>
            // UpdateHandleTableT: std::vector<Handle<TopoT>>
            newlocations.resize(v.size());
            int index = 0;
            for (size_t i = 0; i < v.size(); i++){
                newlocations[i] = { v[i].exists == false ? -1 : (index++) };
            }
            v.erase(std::remove_if(v.begin(), v.end(), [](const typename ComponentTableT::value_type & t){
                return !t.exists;
            }), v.end());
            return 0;
        }


        template <class UpdateHandleTableT, class T>
        inline int UpdateOldHandle(const UpdateHandleTableT & newlocationTable, Handle<T> & h) {
            // UpdateHandleTableT: std::vector<Handle<T>>
            if (h.valid())
                h = newlocationTable[h.id];
            return 0;
        }       


        template <class UpdateHandleTableT, class ContainerT, class = std::enable_if_t<IsContainer<ContainerT>::value>>
        inline int UpdateOldHandleContainer(const UpdateHandleTableT& newlocationTable, ContainerT & hs) {
            // UpdateHandleTableT: std::vector<Handle<TopoT>>
            for (auto & h : hs){
                if (h.valid())
                    h = newlocationTable[h.id];
            }
            return 0;
        }

        template <class UpdateHandleTableT, class K, class CompareK, class AllocK>
        inline int UpdateOldHandleContainer(const UpdateHandleTableT& newlocationTable, std::set<Handle<K>, CompareK, AllocK> & hs) {
            // UpdateHandleTableT: std::vector<Handle<TopoT>>
            std::set<Handle<K>, CompareK, AllocK> oldhs = hs;
            hs.clear();
            for (const auto & oldh : oldhs) {
                hs.insert(newlocationTable[oldh.id]);
            }
            return 0;
        }

       /* namespace {
            template <class UpdateHandleTableTupleT, class ContainerT, int ... Is>
            inline int UpdateOldHandleUsingSequence1(const UpdateHandleTableTupleT & newlocationTable, ContainerT & hs, Sequence<Is...>){
                int dummy[] = { UpdateOldHandle(std::get<Is>(newlocationTable), hs)... };
                return 0;
            }
            template <class UpdateHandleTableT, class ContainerTupleT, int ... Is>
            inline int UpdateOldHandleUsingSequence2(const UpdateHandleTableT & newlocationTable, ContainerTupleT & hs, Sequence<Is...>){
                int dummy[] = { UpdateOldHandle(newlocationTable, std::get<Is>(hs))... };
                return 0;
            }
        }

        template <class ... UpdateHandleTableTs, class ContainerT, class = std::enable_if_t<!IsTuple<ContainerT>::value>>
        inline int UpdateOldHandle(const std::tuple<UpdateHandleTableTs...> & newlocationTable, ContainerT & hs){
            return UpdateOldHandleUsingSequence1(newlocationTable, hs, SequenceGenerator<sizeof...(UpdateHandleTableTs)>::type());
        }

        template <class UpdateHandleTableT, class ... ContainerTs>
        inline int UpdateOldHandle(const UpdateHandleTableT & newlocationTable, std::tuple<ContainerTs...> & hs){
            return UpdateOldHandleUsingSequence2(newlocationTable, hs, SequenceGenerator<sizeof...(ContainerTs)>::type());
        }*/

        template <class ContainerT>
        inline void RemoveInValidHandleFromContainer(ContainerT & hs) {
            auto invalid = typename std::iterator_traits<decltype(std::begin(hs))>::value_type();
            invalid.reset();
            hs.erase(std::remove(std::begin(hs), std::end(hs), invalid), std::end(hs));
        }
        template <class T, int N>
        inline void RemoveInValidHandleFromContainer(std::array<T, N> & hs){}
        template <class K, class CompareK, class AllocK>
        inline void RemoveInValidHandleFromContainer(std::set<K, CompareK, AllocK> & hs) {
            auto invalid = typename std::iterator_traits<decltype(std::begin(hs))>::value_type();
            invalid.reset();
            hs.erase(invalid);
        }
        template <class K, class HasherK, class EqualK, class AllocK>
        inline void RemoveInValidHandleFromContainer(std::unordered_set<K, HasherK, EqualK, AllocK> & hs){
            auto invalid = typename std::iterator_traits<decltype(std::begin(hs))>::value_type();
            invalid.reset();
            hs.erase(invalid);
        }


    }
}
 


namespace std {

    template <class Tag>
    struct hash<panoramix::core::Handle<Tag>> {
        inline uint64_t operator()(panoramix::core::Handle<Tag> a) const {
            return static_cast<uint64_t>(a.id);
        }
    };

}


#endif