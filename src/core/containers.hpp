#ifndef PANORAMIX_CORE_CONTAINERS_HPP
#define PANORAMIX_CORE_CONTAINERS_HPP

#include "basic_types.hpp"
#include "utilities.hpp"
 
namespace panoramix {
    namespace core {

        
        // RTree Wrapper
        template <class T, class BoundingBoxFunctorT = DefaultBoundingBoxFunctor<T>>
        class RTreeWrapper {
        public:
            using BoxType = decltype(std::declval<BoundingBoxFunctorT>()(std::declval<T>()));
            using ValueType = typename BoxType::Type;
            static const int Dimension = BoxType::Dimension;

            inline explicit RTreeWrapper(BoundingBoxFunctorT getBB = BoundingBoxFunctorT())
                : _rtree(std::make_shared<third_party::RTree<T, ValueType, Dimension>>()),
                _getBoundingBox(getBB) {}

            template <class IteratorT>
            inline RTreeWrapper(IteratorT begin, IteratorT end, BoundingBoxFunctorT getBB = BoundingBoxFunctorT())
                : _rtree(std::make_shared<third_party::RTree<T, ValueType, Dimension>>()),
                _getBoundingBox(getBB){
                insert(begin, end);
            }

        public:
            inline size_t size() const { return _rtree->Count(); }
            inline bool empty() const { return size() == 0; }

            inline void clear() { return _rtree->RemoveAll(); }
            inline const BoundingBoxFunctorT & getBoundingBox() const { return _getBoundingBox; }

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
            inline int search(const BoxType & b, CallbackFunctorT callback) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, callback);
            }

            template <class CallbackFunctorT>
            inline int searchNear(const T & t, CallbackFunctorT callback) const {
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
            std::shared_ptr<third_party::RTree<T, ValueType, Dimension>> _rtree;
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
        template <class KeyT, class ScoreT = double, class KeyToIdT = std::unordered_map<KeyT, int>>
        class MaxHeap {            
        public:
            using iterator = typename std::vector<Scored<KeyT, ScoreT>>::iterator;
            using const_iterator = typename std::vector<Scored<KeyT, ScoreT>>::const_iterator;
            using value_type = Scored<KeyT, ScoreT>;

            inline MaxHeap() {}
            
            template <class IteratorT, 
            class = std::enable_if_t<
                std::is_same<std::iterator_traits<IteratorT>::value_type, core::Scored<KeyT, ScoreT>>::value>>
            inline MaxHeap(IteratorT begin, IteratorT end) {
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
            inline MaxHeap(IteratorT vbegin, IteratorT vend, const ScoreT & defaultScore = ScoreT()) {
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
                    while (id > 0 && _data[parentId(id)].score < _data[id].score){
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
                while (id > 0 && _data[parentId(id)].score < _data[id].score){
                    swapKeys(id, parentId(id));
                    id = parentId(id);
                }
            }
            
            inline void push(const KeyT & t, const ScoreT & s){
                push(core::ScoreAs(t, s));
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
                if (l < _data.size() && _data[id].score < _data[l].score){
                    largest = l;
                }
                if (r < _data.size() && _data[largest].score < _data[r].score){
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
        };

        template <class KeyT, class ScoreT, class KeyToIdT>
        inline bool Contains(const MaxHeap<KeyT, ScoreT, KeyToIdT> & h, const KeyT & k){
            return h.contains(k);
        }


    }
}
 
#endif