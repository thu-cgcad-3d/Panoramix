#ifndef PANORAMIX_CORE_UTILITIES_HPP
#define PANORAMIX_CORE_UTILITIES_HPP

#include <iterator>
#include <cassert>
 
namespace panoramix {
	namespace core {

		template <class ValueT, int dim>
		inline Eigen::Matrix<ValueT, dim, 1> EigenVec(const cv::Vec<ValueT, dim> & v) {
			Eigen::Matrix<ValueT, dim, 1> m;
			for(int i = 0; i < dim; i++)
				m(i) = v(i);
			return m;
		}

		template <class ValueT, int dim>
		inline cv::Vec<ValueT, dim> CVVec(const Eigen::Matrix<ValueT, dim, 1> & v) {
			cv::Vec<ValueT, dim> m;
			for(int i = 0; i < dim; i++)
				m(i) = v(i);
			return m;
		}

 
		// elements of container MUST support .exists
		template <class BoolExistIteratorT>
		class JumpIterator {

		public:
			using Iterator = BoolExistIteratorT;
			using value_type = typename std::iterator_traits<Iterator>::value_type;
			using reference = typename std::iterator_traits<Iterator>::reference;
			using pointer = typename std::iterator_traits<Iterator>::pointer;
			using difference_type = typename std::iterator_traits<Iterator>::difference_type;
			using iterator_category = std::forward_iterator_tag;

			/**
			* @brief JumpIterator
			* @param it_
			* @param end_
			*/
			inline JumpIterator(Iterator it_, Iterator end_)
				: _it(it_), _end(end_) {
				if (_it != _end && !_it->exists)
					++ (*this);
			}

			/**
			* @brief operator ++
			* @return
			*/
			inline JumpIterator & operator++() {
				assert(_it != _end);
				++_it;
				while (_it != _end && !_it->exists)
					++_it;
				return *this;
			}

			/**
			* @brief operator *
			* @return
			*/
			inline reference operator * () const {
				return *_it;
			}

			/**
			* @brief operator ->
			* @return
			*/
			inline pointer operator -> () const {
				return &(*_it);
			}

			/**
			* @brief operator ==
			* @param i
			* @return
			*/
			inline bool operator == (const JumpIterator & i) const {
				return _it == i._it;
			}

			/**
			* @brief operator !=
			* @param i
			* @return
			*/
			inline bool operator != (const JumpIterator & i) const {
				return !(*this == i);
			}

			/**
			* @brief internalIterator
			* @return
			*/
			inline Iterator internalIterator() const {
				return _it;
			}

		private:
			BoolExistIteratorT _it;
			BoolExistIteratorT _end;
		};

		/**
		* @brief class JumpContainerWrapper
		*/
		template <class BoolExistContainerT>
		class JumpContainerWrapper {
		public:
			using BoolExistIteratorT = typename BoolExistContainerT::iterator;
			using iterator = JumpIterator<BoolExistIteratorT>;
			using value_type = typename std::iterator_traits<iterator>::value_type;

			inline JumpContainerWrapper(BoolExistContainerT * cont_) : _cont(cont_){}
			inline iterator begin() { return iterator(std::begin(*_cont), std::end(*_cont)); }
			inline iterator end() { return iterator(std::end(*_cont), std::end(*_cont)); }
			inline iterator begin() const { return iterator(std::begin(*_cont), std::end(*_cont)); }
			inline iterator end() const { return iterator(std::end(*_cont), std::end(*_cont)); }

		private:
			BoolExistContainerT * _cont;
		};

		/**
		* @brief class ConstJumpContainerWrapper
		*/
		template <class BoolExistContainerT>
		class ConstJumpContainerWrapper {
		public:
			using BoolExistIteratorT = typename BoolExistContainerT::const_iterator;
			using iterator = JumpIterator<BoolExistIteratorT>;
			using value_type = typename std::iterator_traits<iterator>::value_type;

			inline ConstJumpContainerWrapper(const BoolExistContainerT * cont_) : _cont(cont_){}
			inline iterator begin() const { return iterator(std::begin(*_cont), std::end(*_cont)); }
			inline iterator end() const { return iterator(std::end(*_cont), std::end(*_cont)); }

		private:
			const BoolExistContainerT * _cont;
		};
 
	}
}
 
#endif