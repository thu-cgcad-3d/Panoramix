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

        struct Dummy {
            template <class Archive> inline void serialize(Archive & ar) {}
        };


        static const int Dynamic = -1;


        /**
        * @brief Handle struct
        */
        template <class Tag>
        struct Handle {
            int64_t id;
            inline Handle(int64_t id_ = -1) : id(id_){}
            inline bool operator == (Handle h) const { return id == h.id; }
            inline bool operator != (Handle h) const { return id != h.id; }
            inline void reset() { id = -1; }
            inline bool isValid() const { return id >= 0; }
            inline bool isInvalid() const { return id < 0; }
            template <class Archive> inline void serialize(Archive & ar) { ar(id); }
        };
        template <class Tag>
        using HandleArray = std::vector<Handle<Tag>>;
        template <class Tag>
        using HandlePtrArray = std::vector<Handle<Tag>*>;
        template <class Tag>
        inline bool operator < (const Handle<Tag> & a, const Handle<Tag> & b){
            return a.id < b.id;
        }
        template <class Tag>
        struct HandleHasher {
            inline uint64_t operator()(Handle<Tag> a) const {
                return static_cast<uint64_t>(a.id);
            }
        };

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


        // handled table
        template <class HandleT, class DataT>
        struct HandledTable {
            static_assert(IsHandle<HandleT>::value, "HandleT must be a Handle!");
            HandledTable() {}
            explicit HandledTable(size_t maxSize) : data(maxSize) {}
            HandledTable(size_t maxSize, const DataT & d) : data(maxSize, d) {}
            void resize(size_t sz) { data.resize(sz); }
            const DataT & operator[](HandleT h) const { return data[h.id]; }
            const DataT & at(HandleT h) const { return data[h.id]; }
            DataT & operator[](HandleT h) { return data[h.id]; }

            template <class Archive> inline void serialize(Archive & ar) { ar(data); }
            std::vector<DataT> data;
        };





        // triplet 
        template <class TopoT, class DataT>
        struct Triplet {
            using TopoType = TopoT;
            using DataType = DataT;
            TopoT topo;
            uint8_t exists;
            DataT data;
            inline Triplet(){}
            inline Triplet(const TopoT & t, const DataT & d, uint8_t e = true)
                : topo(t), exists(e), data(d){}
            inline Triplet(const TopoT & t, DataT && d, uint8_t e = true)
                : topo(t), exists(e), data(std::forward<DataT>(d)) {}
            template <class Archive> inline void serialize(Archive & ar) { ar(topo, exists, data); }
        };
        template <class TopoT, class DataT>
        struct TripletExistsPred {
            inline uint8_t operator()(const Triplet<TopoT, DataT> & t) const {
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


        template <class TopoT, class DataT>
        inline auto BoundingBox(const Triplet<TopoT, DataT> & t) -> decltype(BoundingBox(t.data)) {
            return BoundingBox(t.data);
        }



        namespace {
            // helper functions
            template <class ComponentTableT, class UpdateHandleTableT>
            int RemoveAndMap(ComponentTableT & v, UpdateHandleTableT & newlocations) {
                // ComponentTableT : std::vector<Triplet<TopoT, DataT>>
                // UpdateHandleTableT: std::vector<Handle<TopoT>>
                newlocations.resize(v.size());
                int64_t index = 0;
                for (size_t i = 0; i < v.size(); i++){
                    newlocations[i] = { v[i].exists == false ? -1 : (index++) };
                }
                //for (int i = int(v.size() - 1); i >= 0; --i){
                //    if (!v[i].exists){
                //        v.erase(v.begin() + i);
                //    }
                //}
                v.erase(std::remove_if(v.begin(), v.end(), [](const typename ComponentTableT::value_type & t){
                    return !t.exists;
                }), v.end());
                return 0;
            }
            template <class UpdateHandleTableT, class TopoT>
            inline void UpdateOldHandle(const UpdateHandleTableT & newlocationTable, Handle<TopoT> & h) {
                // UpdateHandleTableT: std::vector<Handle<TopoT>>
                if (h.isValid())
                    h = newlocationTable[h.id];
            }
            template <class UpdateHandleTableT, class ContainerT>
            inline void UpdateOldHandleContainer(const UpdateHandleTableT& newlocationTable, ContainerT & hs) {
                // UpdateHandleTableT: std::vector<Handle<TopoT>>
                for (auto & h : hs){
                    if (h.isValid())
                        h = newlocationTable[h.id];
                }
            }
            template <class UpdateHandleTableT, class K, class CompareK, class AllocK>
            inline void UpdateOldHandleContainer(const UpdateHandleTableT& newlocationTable, std::set<K, CompareK, AllocK> & hs) {
                // UpdateHandleTableT: std::vector<Handle<TopoT>>
                std::set<K, CompareK, AllocK> oldhs = hs;
                hs.clear();
                for (const auto & oldh : oldhs) {
                    hs.insert(newlocationTable[oldh.id]);
                }
            }
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