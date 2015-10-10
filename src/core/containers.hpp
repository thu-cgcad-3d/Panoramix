#pragma once


#include <nanoflann.hpp>

#include "meta.hpp"
#include "basic_types.hpp"
#include "utility.hpp"
#include "handle.hpp"
#include "iterators.hpp"
 
namespace pano {
    namespace core {



        template <class IteratorT>
        class PointCloudWrapper {
            using PointType = typename std::iterator_traits<IteratorT>::value_type;
            static_assert(IsVecOrPoint<PointType>::value, "value of ContainerT must be core::Point/core::Vec");
            using ValueType = typename PointType::value_type;
            static const int Dimension = PointType::channels;

        public:
            inline explicit PointCloudWrapper(IteratorT b, IteratorT e) : begin(b), end(e) {}

            inline size_t kdtree_get_point_count() const { return std::distance(begin, end); }
            inline ValueType kdtree_distance(const ValueType * p1, const size_t idx_p2, size_t size) const {
                auto iter = begin;
                std::advance(iter, idx_p2);
                ValueType dist = 0.0;
                for (int i = 0; i < size; i++) {
                    dist += Square((*iter)(i)-p1[i]);
                }
                return dist;
            }
            inline ValueType kdtree_get_pt(const size_t idx, int dim) const {
                auto iter = begin;
                std::advance(iter, idx);
                return (*iter)(dim);
            }
            template <class BBOX>
            bool kdtree_get_bbox(BBOX &bb) const { return false; }
        private:
            IteratorT begin, end;
        };




        // RTreeSet
        template <class T, class BoundingBoxFunctorT = DefaultBoundingBoxFunctor>
        class RTreeSet {
        public:
            using BoxType = decltype(std::declval<BoundingBoxFunctorT>()(std::declval<T>()));
            using ValueType = typename BoxType::Type;
            static const int Dimension = BoxType::Dimension;

            inline explicit RTreeSet(const BoundingBoxFunctorT & bboxFun = BoundingBoxFunctorT())
                : _rtree(std::make_unique<third_party::RTree<T, ValueType, Dimension>>()), _bbox(bboxFun) {}

            template <class IteratorT>
            inline RTreeSet(IteratorT begin, IteratorT end, const BoundingBoxFunctorT & bboxFun = BoundingBoxFunctorT())
                : _rtree(std::make_unique<third_party::RTree<T, ValueType, Dimension>>()), _bbox(bboxFun) {
                insert(begin, end);
            }

            inline RTreeSet(RTreeSet && r)
                : _rtree(std::move(r._rtree)), _bbox(std::move(r._bbox)) {}
            inline RTreeSet & operator = (RTreeSet && r) {
                _rtree = std::move(r._rtree);
                _bbox = std::move(r._bbox);
                return *this;
            }

            RTreeSet(const RTreeSet &) = delete;
            RTreeSet & operator = (const RTreeSet &) = delete;

        public:
            inline size_t size() const { return _rtree->Count(); }
            inline bool empty() const { return size() == 0; }

            inline void clear() { return _rtree->RemoveAll(); }

            inline void insert(const T & t) {
                auto box = bbox(t);
                _rtree->Insert(box.minCorner.val, box.maxCorner.val, t);
            }

            template <class IteratorT>
            void insert(IteratorT begin, IteratorT end) {
                while (begin != end) {
                    insert(*begin);
                    ++begin;
                }
            }

            template <class CallbackFunctorT>
            inline int search(const BoxType & b, CallbackFunctorT && callback) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, callback);
            }
            inline int count(const BoxType & b) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, StaticConstantFunctor<bool, true>());
            }

        private:
            std::unique_ptr<third_party::RTree<T, ValueType, Dimension>> _rtree;
            BoundingBoxFunctorT _bbox;
        };


        // RTreeMap
        template <class T, class ValT, class BoundingBoxFunctorT = DefaultBoundingBoxFunctor>
        class RTreeMap {
        public:
            using BoxType = decltype(std::declval<BoundingBoxFunctorT>()(std::declval<T>()));
            using ValueType = typename BoxType::Type;
            static const int Dimension = BoxType::Dimension;

            inline explicit RTreeMap(const BoundingBoxFunctorT & bboxFun = BoundingBoxFunctorT())
                : _rtree(std::make_unique<third_party::RTree<std::pair<T, ValT>, ValueType, Dimension>>()), _bbox(bboxFun) {}

            template <class IteratorT>
            inline RTreeMap(IteratorT begin, IteratorT end, const BoundingBoxFunctorT & bboxFun = BoundingBoxFunctorT())
                : _rtree(std::make_unique<third_party::RTree<std::pair<T, ValT>, ValueType, Dimension>>()), _bbox(bboxFun) {
                insert(begin, end);
            }

            inline RTreeMap(RTreeMap && r)
                : _rtree(std::move(r._rtree)), _bbox(std::move(r._bbox)) {}
            inline RTreeMap & operator = (RTreeMap && r) {
                _rtree = std::move(r._rtree);
                _bbox = std::move(r._bbox);
                return *this;
            }

            RTreeMap(const RTreeMap &) = delete;
            RTreeMap & operator = (const RTreeMap &) = delete;

        public:
            inline size_t size() const { return _rtree->Count(); }
            inline bool empty() const { return size() == 0; }

            inline void clear() { return _rtree->RemoveAll(); }

            inline void insert(const std::pair<T, ValT> & p) {
                auto box = _bbox(p.first);
                _rtree->Insert(box.minCorner.val, box.maxCorner.val, p);
            }
            inline void emplace(const T & key, const ValT & val) {
                auto box = _bbox(key);
                _rtree->Insert(box.minCorner.val, box.maxCorner.val, std::make_pair(key, val));
            }

            template <class IteratorT>
            void insert(IteratorT begin, IteratorT end) {
                while (begin != end) {
                    insert(*begin);
                    ++begin;
                }
            }

            template <class CallbackFunctorT>
            inline int search(const BoxType & b, CallbackFunctorT && callback) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, callback);
            }
            inline int count(const BoxType & b) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, StaticConstantFunctor<bool, true>());
            }

        private:
            std::unique_ptr<third_party::RTree<std::pair<T, ValT>, ValueType, Dimension>> _rtree;
            BoundingBoxFunctorT _bbox;
        };




        // simple RTree
        template <class BoxT, class T>
        class RTree {
        public:
            using BoxType = BoxT;
            using ValueType = typename BoxType::Type;
            static const int Dimension = BoxType::Dimension;

            inline explicit RTree()
                : _rtree(std::make_unique<third_party::RTree<T, ValueType, Dimension>>()) {}

            template <class IteratorT>
            inline RTree(IteratorT begin, IteratorT end)
                : _rtree(std::make_unique<third_party::RTree<T, ValueType, Dimension>>()){
                insert(begin, end);
            }

            inline RTree(RTree && r)
                : _rtree(std::move(r._rtree)) {}
            inline RTree & operator = (RTree && r){
                _rtree = std::move(r._rtree);
                return *this;
            }

            RTree(const RTree &) = delete;
            RTree & operator = (const RTree &) = delete;

        public:
            inline size_t size() const { return _rtree->Count(); }
            inline bool empty() const { return size() == 0; }

            inline void clear() { return _rtree->RemoveAll(); }

            inline void insert(const BoxType & box, const T & t) {
                _rtree->Insert(box.minCorner.val, box.maxCorner.val, t);
            }

            inline void insert(const std::pair<BoxType, T> & p){
                _rtree->Insert(p.first.minCorner.val, p.first.maxCorner.val, p.second);
            }

            template <class IteratorT>
            void insert(IteratorT begin, IteratorT end) {
                while (begin != end){
                    insert(*begin);
                    ++begin;
                }
            }

            template <class CallbackFunctorT>
            inline int search(const BoxType & b, CallbackFunctorT && callback) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, callback);
            }

            inline int count(const BoxType & b) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, StaticConstantFunctor<bool, true>());
            }

        private:
            std::unique_ptr<third_party::RTree<T, ValueType, Dimension>> _rtree;
        };


        
        // RTree Wrapper
        template <class T, class BoundingBoxFunctorT = DefaultBoundingBoxFunctor>
        class RTreeWrapper {
        public:
            using BoxType = decltype(std::declval<BoundingBoxFunctorT>()(std::declval<T>()));
            using ValueType = typename BoxType::Type;
            static const int Dimension = BoxType::Dimension;

            inline explicit RTreeWrapper(BoundingBoxFunctorT getBB = BoundingBoxFunctorT())
                : _rtree(std::make_unique<third_party::RTree<T, ValueType, Dimension>>()),
                _getBoundingBox(getBB) {}

            template <class IteratorT>
            inline RTreeWrapper(IteratorT begin, IteratorT end, BoundingBoxFunctorT getBB = BoundingBoxFunctorT())
                : _rtree(std::make_unique<third_party::RTree<T, ValueType, Dimension>>()),
                _getBoundingBox(getBB){
                insert(begin, end);
            }

            inline RTreeWrapper(RTreeWrapper && r) 
                : _rtree(std::move(r._rtree)), _getBoundingBox(std::move(r._getBoundingBox)) {}
            inline RTreeWrapper & operator = (RTreeWrapper && r){
                _rtree = std::move(r._rtree);
                _getBoundingBox = std::move(r._getBoundingBox);
                return *this;
            }

            RTreeWrapper(const RTreeWrapper &) = delete;
            RTreeWrapper & operator = (const RTreeWrapper &) = delete;

        public:
            inline size_t size() const { return _rtree->Count(); }
            inline bool empty() const { return size() == 0; }

            inline void clear() { return _rtree->RemoveAll(); }
            inline const BoundingBoxFunctorT & getBoundingBox() const { return _getBoundingBox; }

            inline void insert(const BoxType & box, const T & t) {
                _rtree->Insert(box.minCorner.val, box.maxCorner.val, t);
            }

            inline void insert(const T & t) {
                BoxType box = _getBoundingBox(t);
                for (int i = 0; i < Dimension; i++){
                    if (isnan(box.minCorner[i]) || isnan(box.maxCorner[i])){
#ifdef _DEBUG
                        std::cout << "invalid box type (NaN value), ignore this element" << std::endl;
#endif
                        return;
                    }
                    if (!(box.minCorner[i] <= box.maxCorner[i])) {
#ifdef _DEBUG
                        std::cout << "invalid box type (minCorner[i] > maxCorner[i]), ignore this element" << std::endl;
#endif
                        return;
                    }
                }
                _rtree->Insert(box.minCorner.val, box.maxCorner.val, t);
            }

            template <class IteratorT>
            void insert(IteratorT begin, IteratorT end) {
                while (begin != end){
                    insert(*begin);
                    ++begin;
                }
            }

            template <class CallbackFunctorT>
            inline int search(const BoxType & b, CallbackFunctorT && callback) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, callback);
            }

            inline int count(const BoxType & b) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, StaticConstantFunctor<bool, true>());
            }

            template <class CallbackFunctorT>
            inline int searchNear(const T & t, CallbackFunctorT && callback) const {
                auto b = _getBoundingBox(t);
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, callback);
            }

            template <class IsEqualT = std::equal_to<T>>
            bool contains(const T & t, IsEqualT && cmp = IsEqualT()) const {
                bool exists = false;
                search(_getBoundingBox(t), [&exists, &t, &cmp](const T & ele) {
                    if (cmp(ele, t)) {
                        exists = true;
                        return false;
                    }
                    return true;
                });
                return exists;
            }

        private:
            std::unique_ptr<third_party::RTree<T, ValueType, Dimension>> _rtree;
            BoundingBoxFunctorT _getBoundingBox;
        };



        template <class T, int N>
        class VecSet {
        public:
            inline explicit VecSet(const T & ir = Epsilon)
                : _influenceRange(ir), _rtree(std::make_shared<third_party::RTree<int, T, N>>()) {
            }
            template <class IteratorT,
                class = std::enable_if_t<std::is_same<std::iterator_traits<IteratorT>::value_type, Vec<T, N>>::value >>
            inline VecSet(IteratorT begin, IteratorT end, const T & ir = Epsilon)
                : _influenceRange(ir), _rtree(std::make_shared<third_party::RTree<int, T, N>>()){
                for (auto i = begin; i != end; ++i){
                    insert(*i);
                }
            }

            inline const T & influenceRange() const { return _influenceRange; }
            inline int size() const { return _vecData.size(); }
            inline void clear() { _vecData.clear(); _rtree->RemoveAll(); }

            using const_iterator = typename std::vector<Vec<T, N>>::const_iterator;
            inline const_iterator begin() const { return _vecData.begin(); }
            inline const_iterator end() const { return _vecData.end(); }
            inline const_iterator cbegin() const { return _vecData.cbegin(); }
            inline const_iterator cend() const { return _vecData.cend(); }
            
            void insert(const Vec<T, N> & v){
                Box<T, N> box(v, v);
                for (int i = 0; i < N; i++){
                    box.minCorner[i] -= _influenceRange;
                    box.maxCorner[i] += _influenceRange;
                }
                bool existed = false;
                _rtree->Search(box.minCorner.val, box.maxCorner.val, [this, &v, &existed](int storedId){
                    T dist = Distance(v, _vecData[storedId].front());
                    if (dist <= _influenceRange){
                        existed = true;
                        return false;
                    }
                    return true;
                });
                if (!existed){
                    _vecData.push_back(v);
                    _rtree->Insert(v.val, v.val, static_cast<int>(_vecData.size() - 1));
                }
            }

        private:
            T _influenceRange;
            std::vector<Vec<T, N>> _vecData;
            std::shared_ptr<third_party::RTree<int, T, N>> _rtree;
        };




        template <class T, int N, class ValueT>
        class VecMap {
        public:
            inline explicit VecMap(const T & ir = 1e-5)
                : _influenceRange(ir), _rtree(std::make_shared<third_party::RTree<int, T, N>>()) {
            }
            template <class IteratorT,
                class = std::enable_if_t<std::is_same<std::iterator_traits<IteratorT>::value_type, std::pair<Vec<T, N>, ValueT>>::value >>
            inline VecMap(IteratorT begin, IteratorT end, const T & ir = Epsilon)
                : _influenceRange(ir), _rtree(std::make_shared<third_party::RTree<int, T, N>>()){
                for (auto i = begin; i != end; ++i){
                    insert(*i);
                }
            }

            inline const T & influenceRange() const { return _influenceRange; }
            inline int size() const { return _vecData.size(); }
            inline void clear() { _vecData.clear(); _rtree->RemoveAll(); }

            using const_iterator = typename std::vector<std::pair<Vec<T, N>, ValueT>>::const_iterator;
            inline const_iterator begin() const { return _vecData.begin(); }
            inline const_iterator end() const { return _vecData.end(); }
            inline const_iterator cbegin() const { return _vecData.cbegin(); }
            inline const_iterator cend() const { return _vecData.cend(); }

            using iterator = typename std::vector<std::pair<Vec<T, N>, ValueT>>::iterator;
            inline iterator begin() { return _vecData.begin(); }
            inline iterator end() { return _vecData.end(); }

            inline iterator find(const Vec<T, N> & v) const {
                int nid = findNearestStoredIdByRTree(v);
                return nid == -1 ? _vecData.end() : _vecData.begin() + nid;
            }

            inline bool contains(const Vec<T, N> & v) const { return findNearestStoredIdByRTree(v) != -1; }
            inline const ValueT & at(const Vec<T, N> & v) const {
                int nid = findNearestStoredIdByRTree(v);
                assert(nid != -1);
                return _vecData[nid].second;
            }
            inline const ValueT & operator[](const Vec<T, N> & v) const { return at(v); }
            inline ValueT & operator[](const Vec<T, N> & v) { return insert(std::make_pair(v, ValueT()))->second;}

            iterator insert(const std::pair<Vec<T, N>, ValueT> & d){
                int nid = findNearestStoredIdByRTree(d.first);
                if (nid == -1){
                    _vecData.emplace_back(d.first, d.second);
                    return _vecData.end() - 1;
                }
                else{
                    _vecData[nid].second = d.second;
                    return _vecData.begin() + nid;
                }
            }

        private:
            int findNearestStoredIdByRTree(const Vec<T, N> & v) const {
                Box<T, N> box(v, v);
                for (int i = 0; i < N; i++){
                    box.minCorner[i] -= _influenceRange;
                    box.maxCorner[i] += _influenceRange;
                }
                int nearesetStoredId = -1;
                T minDist = _influenceRange;
                _rtree->Search(box.minCorner.val, box.maxCorner.val, [this, &v, &minDist, &nearesetStoredId](int storedId){
                    T dist = Distance(v, _vecData[storedId].first);
                    if (dist <= minDist){
                        minDist = dist;
                        nearesetStoredId = storedId;
                    }
                    return true;
                });
                return nearesetStoredId;
            }


        private:
            T _influenceRange;
            std::vector<std::pair<Vec<T, N>, ValueT>> _vecData;
            std::shared_ptr<third_party::RTree<int, T, N>> _rtree;
        };



        template <class T, int N>
        class VecMultiSet {
        public:
            inline explicit VecMultiSet(const T & ir = 1e-5)
                : _influenceRange(ir), _fullSize(0), _rtree(std::make_shared<third_party::RTree<int, T, N>>()) {
            }
            template <class IteratorT, 
                class = std::enable_if_t<std::is_same<std::iterator_traits<IteratorT>::value_type, Vec<T, N>>::value>>
            inline VecMultiSet(IteratorT begin, IteratorT end, const T & ir = Epsilon)
                : _influenceRange(ir), _fullSize(0), _rtree(std::make_shared<third_party::RTree<int, T, N>>()){
                for (auto i = begin; i != end; ++i){
                    insert(*i);
                }
            }

            inline const T & influenceRange() const { return _influenceRange; }
            inline int size() const { return _vecData.size(); }
            inline int fullSize() const { return _fullSize; }
            inline void clear() { _fullSize = 0; _vecData.clear(); _rtree->RemoveAll(); }

            using const_iterator = typename std::vector<std::vector<Vec<T, N>>>::const_iterator;
            inline const_iterator begin() const { return _vecData.begin(); }
            inline const_iterator end() const { return _vecData.end(); }
            inline const_iterator cbegin() const { return _vecData.cbegin(); }
            inline const_iterator cend() const { return _vecData.cend(); }

            void insert(const Vec<T, N> & v){
                Box<T, N> box(v, v);
                for (int i = 0; i < N; i++){
                    box.minCorner[i] -= _influenceRange;
                    box.maxCorner[i] += _influenceRange;
                }       
                int nearesetStoredId = -1;
                T minDist = _influenceRange;
                _rtree->Search(box.minCorner.val, box.maxCorner.val, [this, &v, &minDist, &nearesetStoredId](int storedId){
                    T dist = Distance(v, _vecData[storedId].front());
                    if (dist <= minDist){
                        minDist = dist;
                        nearesetStoredId = storedId;
                    }
                    return true;
                });
                if (nearesetStoredId == -1){ // not storage at all
                    _vecData.emplace_back(1, v);
                    _rtree->Insert(v.val, v.val, static_cast<int>(_vecData.size() - 1));
                }
                else{
                    _vecData[nearesetStoredId].push_back(v);
                }
                _fullSize++;
            }
        private:
            T _influenceRange;
            int _fullSize;
            std::vector<std::vector<Vec<T, N>>> _vecData;
            std::shared_ptr<third_party::RTree<int, T, N>> _rtree;
        };





        // max heap
        template <class KeyT, class ScoreT = double, class ScoreCompareT = std::less<ScoreT>, class KeyToIdT = std::unordered_map<KeyT, int>>
        class MaxHeap {            
        public:
            using iterator = typename std::vector<Scored<KeyT, ScoreT>>::iterator;
            using const_iterator = typename std::vector<Scored<KeyT, ScoreT>>::const_iterator;
            using value_type = Scored<KeyT, ScoreT>;

            inline MaxHeap(const ScoreCompareT & cmp = ScoreCompareT()) : _scoreCompare(cmp) {}
            
            template <class IteratorT, 
            class = std::enable_if_t<
                std::is_same<std::iterator_traits<IteratorT>::value_type, core::Scored<KeyT, ScoreT>>::value>>
            inline MaxHeap(IteratorT begin, IteratorT end, const ScoreCompareT & cmp = ScoreCompareT()) : _scoreCompare(cmp) {
                _data.reserve(std::distance(begin, end));
                while (begin != end){
                    _data.push_back(*begin);
                    _keyToId[_data.back().component] = _data.size() - 1;
                    *begin++;
                }
                makeMaxHeap();
            }

            template <class IteratorT,
            class = std::enable_if_t<
                std::is_same<std::iterator_traits<IteratorT>::value_type, KeyT>::value>>
            inline MaxHeap(IteratorT vbegin, IteratorT vend, const ScoreT & defaultScore = ScoreT(), 
                const ScoreCompareT & cmp = ScoreCompareT()) : _scoreCompare(cmp) {
                _data.reserve(std::distance(vbegin, vend));
                while (vbegin != vend){
                    _data.push_back(core::ScoreAs(*vbegin, defaultScore));
                    _keyToId[_data.back().component] = _data.size() - 1;
                    ++vbegin;
                }
                //makeMaxHeap(); no need to make heap since all scores are same
            }

            template <class IteratorT, class FuncT,
            class = std::enable_if_t<
                std::is_same<std::iterator_traits<IteratorT>::value_type, KeyT>::value &&
                std::is_same<decltype(std::declval<FuncT>()(*std::declval<IteratorT>())), ScoreT>::value>>
            inline MaxHeap(IteratorT vbegin, IteratorT vend, FuncT && fun) {
                _data.reserve(std::distance(vbegin, vend));
                while (vbegin != vend){
                    _data.push_back(core::ScoreAs(*vbegin, fun(*vbegin)));
                    _keyToId[_data.back().component] = _data.size() - 1;
                    ++vbegin;
                }
                makeMaxHeap(); //need to make heap
            }


            inline const_iterator begin() const { return _data.begin(); }
            inline const_iterator end() const { return _data.end(); }
            inline const_iterator cbegin() const { return _data.cbegin(); }
            inline const_iterator cend() const { return _data.cend(); }

            inline const KeyT & top() const { return _data.front().component; }
            inline const ScoreT & topScore() const { return _data.front().score; }
            inline const ScoreT & operator[](const KeyT & key) const { return _data[_keyToId.at(key)].score; }
            inline const ScoreT & at(const KeyT & key) const { return _data[_keyToId.at(key)].score; }
            inline size_t size() const { return _data.size(); }
            inline bool empty() const { return _data.empty(); }
            inline size_t height() const { return static_cast<size_t>(log2(_data.size())); }
            
            void pop(){
                if (_data.empty())
                    return;
                swapKeys(0, _data.size() - 1);
                _keyToId.erase(_data.back().component);
                _data.erase(_data.end() - 1, _data.end());
                maxHeapify(0);
            }
            
            void setScore(const KeyT & key, const ScoreT & newScore){
                auto & oldScore = at(key);
                if (oldScore == newScore)
                    return;
                else if (newScore > oldScore){ // increase key
                    int id = _keyToId[key];
                    _data[id].score = newScore;
                    while (id > 0 && _scoreCompare(_data[parentId(id)].score, _data[id].score)){
                        swapKeys(id, parentId(id));
                        id = parentId(id);
                    }
                }
                else{ // decrease key
                    int id = _keyToId[key];
                    _data[id].score = newScore;
                    maxHeapify(id);
                }
            }
            
            void push(const Scored<KeyT, ScoreT> & e) {
                _data.push_back(e);
                _keyToId[e.component] = _data.size() - 1;
                int id = _data.size() - 1;
                while (id > 0 && _scoreCompare(_data[parentId(id)].score, _data[id].score)){
                    swapKeys(id, parentId(id));
                    id = parentId(id);
                }
            }
            
            inline void push(const KeyT & t, const ScoreT & s){
                push(core::ScoreAs(t, s));
            }

            inline void pushOrSet(const KeyT & key, const ScoreT & newScore) {
                if (!contains(key)) {
                    push(key, newScore);
                } else {
                    setScore(key, newScore);
                }
            }

            inline void clear(){ _data.clear(); _keyToId.clear(); }
            inline void swap(MaxHeap<KeyT, ScoreT> & h){ 
                _data.swap(h._data); 
                _keyToId.swap(h._keyToId);
            }

            inline bool contains(const KeyT & k) const { return _keyToId.find(k) != _keyToId.end(); }

        private:
            static inline int parentId(int id) { return (id - 1) / 2; }
            static inline int leftId(int id) { return id * 2 + 1; }
            static inline int rightId(int id) { return id * 2 + 2; }

            inline void swapKeys(int id1, int id2){
                std::swap(_keyToId[_data[id1].component], _keyToId[_data[id2].component]);
                std::swap(_data[id1], _data[id2]);
            }

            void maxHeapify(int id) {
                auto l = leftId(id);
                auto r = rightId(id);
                int largest = id;
                if (l < _data.size() && _scoreCompare(_data[id].score, _data[l].score)){
                    largest = l;
                }
                if (r < _data.size() && _scoreCompare(_data[largest].score, _data[r].score)) {
                    largest = r;
                }
                if (largest != id){
                    swapKeys(id, largest);
                    maxHeapify(largest);
                }
            }

            inline void makeMaxHeap() {
                for (int i = _data.size() / 2 - 1; i >= 0; --i){
                    maxHeapify(i);
                }
            }

        private:
            std::vector<Scored<KeyT, ScoreT>> _data;
            KeyToIdT _keyToId;
            ScoreCompareT _scoreCompare;
        };

        template <class KeyT, class ScoreT, class ScoreCompareT, class KeyToIdT>
        inline bool Contains(const MaxHeap<KeyT, ScoreT, ScoreCompareT, KeyToIdT> & h, const KeyT & k) {
            return h.contains(k);
        }



        // sparse dict
        template <class T>
        class Dictionary {
            struct Node {
                std::vector<std::unique_ptr<Node>> children;
                std::unique_ptr<T> val;
                bool isLeaf() const { return children.empty(); }
                Node() : val(nullptr){}
                explicit Node(size_t nc) : children(nc), val(nullptr) {}
                explicit Node(std::unique_ptr<T> v) : val(std::move(v)) {}
                Node(Node && n) : children(std::move(n.children)), val(std::move(n.val)) {}
                Node(const Node &) = delete;
                Node & operator = (Node && n) { children = std::move(n.children), val = std::move(n.val); return *this; }
                Node & operator = (const Node &) = delete;
                template <class Archive>
                void serialize(Archive & ar) {
                    ar(children, val);
                }
            };

        public:
            Dictionary() : _root(nullptr), _size(0){}
            explicit Dictionary(const std::vector<size_t> & ncs) : _nchildren(ncs), _root(nullptr), _size(0) {}
            explicit Dictionary(std::vector<size_t> && ncs) : _nchildren(std::move(ncs)), _root(nullptr), _size(0) {}
            Dictionary(Dictionary && t) : _nchildren(std::move(t._nchildren)), _root(std::move(t._root)), _size(t._size) { t._size = 0; }
            Dictionary(const Dictionary &) = delete;

            Dictionary & operator = (Dictionary && t) {
                _nchildren = std::move(t._nchildren);
                _root = std::move(t._root);
                std::swap(_size, t._size);
                return *this;
            }
            Dictionary & operator = (const Dictionary &) = delete;

        public:
            template <class IndexT, class TT, class = std::enable_if_t<std::is_integral<IndexT>::value>>
            void insert(const IndexT * inds, TT && val) {
                Node * leaf = create(inds);
                if (!leaf->val) {
                    leaf->val = std::make_unique<T>(std::forward<TT>(val));
                    _size++;
                } else {
                    leaf->val = std::make_unique<T>(std::forward<TT>(val));
                }
            }
            template <class IndexT, class TT>
            void insert(std::initializer_list<IndexT> inds, TT && val) {
                assert(inds.size() == _nchildren.size());
                return insert(inds.begin(), std::forward<TT>(val));
            }
            template <class IndexT, class = std::enable_if_t<std::is_integral<IndexT>::value>>
            bool contains(const IndexT * inds) const {
                return locate(inds) != nullptr;
            }
            template <class IndexT>
            bool contains(std::initializer_list<IndexT> inds) const {
                assert(inds.size() == _nchildren.size());
                return contains(inds.begin());
            }
            template <class IndexT, class = std::enable_if_t<std::is_integral<IndexT>::value>>
            const T & at(const IndexT * inds) const {
                auto n = locate(inds);
                assert(n->isLeaf());
                return *(n->val);
            }
            template <class IndexT, class = std::enable_if_t<std::is_integral<IndexT>::value>>
            T & at(const IndexT * inds) {
                auto n = locate(inds);
                assert(n->isLeaf());
                return *(n->val);
            }
            template <class IndexT>
            const T & at(std::initializer_list<IndexT> inds) const {
                assert(inds.size() == _nchildren.size()); 
                return at(inds.begin());
            }
            template <class IndexT>
            T & at(std::initializer_list<IndexT> inds) {
                assert(inds.size() == _nchildren.size()); 
                return at(inds.begin()); 
            }

            size_t size() const { return _size; }

            template <class Archive>
            void serialize(Archive & ar) {
                ar(_nchildren, _size, _root);
            }

        private:
            template <class IndexT = size_t>
            Node * create(const IndexT * inds) {
                if (_nchildren.empty()) {
                    return nullptr;
                }
                if (!_root) {
                    _root = std::make_unique<Node>(_nchildren[0]);
                }
                Node * parent = _root.get();
                int lastInd = inds[0];
                for (int i = 1; i < _nchildren.size(); i++) {
                    auto ind = inds[i];
                    size_t nc = _nchildren[i];
                    if (!parent->children[lastInd]) {
                        parent->children[lastInd] = std::make_unique<Node>(nc);
                    }
                    parent = parent->children[lastInd].get();
                    lastInd = ind;
                }
                if (!parent->children[lastInd]) {
                    parent->children[lastInd] = std::make_unique<Node>(nullptr);
                }
                return parent->children[lastInd].get();
            }

            template <class IndexT = size_t>
            Node * locate(const IndexT * inds) const {
                Node * curNode = _root.get();
                for (int i = 0; i < _nchildren.size(); i++) {
                    auto ind = inds[i];
                    size_t nc = _nchildren[i];
                    if (!curNode) {
                        return nullptr;
                    }
                    assert(ind < nc);
                    curNode = curNode->children[ind].get();
                }
                return curNode;
            }


        private:
            std::vector<size_t> _nchildren;
            size_t _size;
            std::unique_ptr<Node> _root;
        };


    }
}
 
