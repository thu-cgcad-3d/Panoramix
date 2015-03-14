#ifndef PANORAMIX_VIS_PROJECT_HPP
#define PANORAMIX_VIS_PROJECT_HPP

#include "../core/meta.hpp"
#include "../core/basic_types.hpp"
#include "../core/generic_topo.hpp"


namespace panoramix {
    namespace vis {

        struct Step;
        using StepPtr = std::shared_ptr<Step>;
        template <class ResultT> struct StepWithTypedResult;

        struct Step {
            std::string name;
            virtual void update(const std::vector<const Step*> & dependencies) = 0;
            template <class ResultT>
            inline bool has() const { 
                return dynamic_cast<StepWithTypedResult<ResultT>*>(this) != nullptr; 
            }
            template <class ResultT>
            inline const ResultT & storage() const {
                return dynamic_cast<StepWithTypedResult<ResultT>*>(this)->storage;
            }
            template <class ResultT>
            inline ResultT & storage() {
                return dynamic_cast<StepWithTypedResult<ResultT>*>(this)->storage;
            }
        };       


        template <class ResultT>
        struct StepWithTypedResult : Step {
            ResultT storage;
        };
        
        template <class UpdateFunctionT>
        struct StepWithTypedUpdater
            : StepWithTypedResult<typename core::FunctionTraits<UpdateFunctionT>::ResultType> {
            UpdateFunctionT updater;

            inline explicit StepWithTypedUpdater(UpdateFunctionT && fun) : updater(std::move(fun)) {}

            virtual void update(const std::vector<const Step*> & dependencies) override {
                assert(core::FunctionTraits<UpdateFunctionT>::ArgumentsNum == dependencies.size());
                updateUsingSequence(dependencies, 
                    typename core::SequenceGenerator<core::FunctionTraits<UpdateFunctionT>::ArgumentsNum>::type());
            }        
            template <int ...Idx>
            void updateUsingSequence(const std::vector<const Step*> & dependencies, core::Sequence<Idx...>){
                using ArgumentsTupleType =
                    typename core::FunctionTraits<UpdateFunctionT>::ArgumentsTupleType;
                //bool typeChecks[] = { dependencies.at(Idx)->has<typename std::tuple_element<Idx, ArgumentsTupleType>::type>() ... };
                //assert(std::all_of(std::begin(typeChecks), std::end(typeChecks), [](bool b){return b; }));
                storage = updater(dependencies.at(Idx)->storage<typename std::tuple_element<Idx, ArgumentsTupleType>::type>() ...);
            }
        };



        class ProjectCore {
        public:
            size_t size() const { return _steps.size(); }

            template <class UpdateFunctionT>
            inline int addStep(const std::string & name, UpdateFunctionT && updater, 
                const std::vector<int> & dependencies = std::vector<int>()) {
                auto step = std::make_shared<StepWithTypedUpdater<std::decay_t<UpdateFunctionT>>>(std::forward<UpdateFunctionT>(updater));
                step->name = name;
                for (int d : dependencies){
                    assert(d >= 0 && d < _steps.size());
                }
                _steps.push_back(step);
                _timeStamps.push_back(core::CurrentTime());
                _dependencies.push_back(dependencies);
                _stepIds[step.get()] = _steps.size() - 1;
                return _steps.size() - 1;
            }

            inline int id(const StepPtr & s) const { return _stepIds.at(s.get()); }
            inline int id(const Step * s) const { return _stepIds.at(s); }

            inline const std::string & stepNameAt(int id) const { return _steps[id]->name; }
            
            template <class ResultT>
            inline const ResultT & stepStorageAt(int id) const { 
                return _steps[id]->storage<ResultT>(); 
            }
            template <class ResultT>
            inline ResultT & stepStorageAt(int id) {
                return _steps[id]->storage<ResultT>();
            }

            inline const std::vector<int> & stepDependenciesAt(int id) const { return _dependencies[id]; }
            inline void setStepModified(int id) { _timeStamps[id] = core::CurrentTime(); }
            bool needsUpdate(int id) const {
                for (int d : _dependencies[id]){
                    if (_timeStamps[d] >= _timeStamps[id]){
                        return true;
                    }
                }
                return false;
            }
            
            /*using StepUpdateCallback = std::function<void(int stepId, int percent)>;*/
            template <class CallbackT>
            void updateAll(CallbackT && callback) {
                for (int i = 0; i < _steps.size(); i++){
                    if (needsUpdate(i)){
                        std::cout << "updating [" << _steps[i]->name << "]" << std::endl;
                        callback(i);
                        std::vector<const Step *> deps(_dependencies[i].size());
                        for (int k = 0; k < deps.size(); k++){
                            deps[k] = _steps[_dependencies[i][k]].get();
                        }
                        _steps[i]->update(deps);
                        _timeStamps[i] = core::CurrentTime();
                    }
                }
            }

        private:
            std::vector<StepPtr> _steps;
            std::vector<core::TimeStamp> _timeStamps;
            std::vector<std::vector<int>> _dependencies;
            std::unordered_map<Step const *, int> _stepIds;
        };
        

    }
}
 
#endif