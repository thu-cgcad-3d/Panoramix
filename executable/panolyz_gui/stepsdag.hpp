#pragma once
 
#include <QtWidgets>

#include "../../src/core/meta.hpp"
#include "../../src/core/basic_types.hpp"


struct StepWidgetInterface {
    virtual void refreshDataAsync() = 0; // on back-end
    virtual void refreshData() = 0; // on gui thread
    virtual void updatePainting() = 0; // on gui thread
    virtual void showWidget() = 0; // on gui thread
    virtual void hideWidget() = 0; // on gui thread
};


struct Lockable {
    QReadWriteLock lock;
    inline void lockForRead() { lock.lockForRead(); }
    inline void lockForWrite() { lock.lockForWrite(); }
    inline void unlock() { lock.unlock(); }
};

template <class T>
struct LockableType : Lockable {
    T component;
};


struct StepData : Lockable {
    panoramix::core::TimeStamp timeStamp;
    inline void setModified() { timeStamp = panoramix::core::CurrentTime(); }
    inline bool isOlderThan(const StepData & d) const { return timeStamp < d.timeStamp; }
    virtual StepWidgetInterface * createBindingWidgetAndActions(
        QList<QAction*> & actions, QWidget * parent = nullptr) {
        return nullptr;
    }
};

using DataPtr = std::shared_ptr<StepData>;

template <class T> struct DataOfType;

namespace {
    template <class T>
    struct HasBindingWidgetAndActionsImpl {
        template <class TT>
        static auto test(int) -> decltype(
            CreateBindingWidgetAndActions(std::declval<DataOfType<T>>(), std::declval<QList<QAction*>&>(), (QWidget*)nullptr),
            std::true_type()
            );
        template <class>
        static std::false_type test(...);
        static const bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value;
    };
}
template <class T>
struct HasBindingWidgetAndActions : std::integral_constant<bool, HasBindingWidgetAndActionsImpl<T>::value> {};


template <class T>
struct DataOfType : StepData {
    T content;
    DataOfType() {}
    DataOfType(const T & c) : content(c) {}
    DataOfType(T && c) : content(std::move(c)) {}
    DataOfType(const DataOfType &) = delete;

    virtual StepWidgetInterface * createBindingWidgetAndActions(
        QList<QAction*> & actions, QWidget * parent) {
        return createBindingWidgetAndActionsImpl(actions, parent, HasBindingWidgetAndActions<T>());
    }

    inline StepWidgetInterface * createBindingWidgetAndActionsImpl(
        QList<QAction*> & actions, QWidget * parent, std::true_type &){
        return CreateBindingWidgetAndActions(*this, actions, parent);
    }
    inline StepWidgetInterface * createBindingWidgetAndActionsImpl(
        QList<QAction*> & actions, QWidget * parent, std::false_type &){
        return nullptr;
    }
};



struct Step;
using StepPtr = std::shared_ptr<Step>;
template <class ResultT> struct StepWithTypedResult;


struct Step {
    DataPtr data;
    virtual void update(const std::vector<StepData *> & dependencies) = 0;
    inline bool null() const { return !data; }

    template <class T>
    inline T & contentAs() { return (dynamic_cast<DataOfType<T>*>(data.get()))->content; }
    template <class T>
    inline const T & contentAs() const { (dynamic_cast<const DataOfType<T>*>(data.get()))->content; }
};

template <class UpdateFunctionT>
struct StepWithTypedUpdater : Step {
    using ResultContentType = std::decay_t<typename panoramix::core::FunctionTraits<UpdateFunctionT>::ResultType>;
    using ResultDataType = DataOfType<ResultContentType>;
    using ArgumentsTupleType =
        typename panoramix::core::FunctionTraits<UpdateFunctionT>::ArgumentsTupleType;

    UpdateFunctionT updater;
    inline explicit StepWithTypedUpdater(UpdateFunctionT && fun)
        : updater(std::move(fun)) {
        data = std::make_shared<ResultDataType>();
    }
    inline explicit StepWithTypedUpdater(const UpdateFunctionT & fun)
        : updater(fun){
        data = std::make_shared<ResultDataType>();
    }

    virtual void update(const std::vector<StepData *> & dependencies) override {
        assert(panoramix::core::FunctionTraits<UpdateFunctionT>::ArgumentsNum == dependencies.size());
        updateUsingSequence(dependencies,
            typename panoramix::core::SequenceGenerator<panoramix::core::FunctionTraits<UpdateFunctionT>::ArgumentsNum>::type());
    }
    template <int ...Idx>
    void updateUsingSequence(const std::vector<StepData *> & dependencies, panoramix::core::Sequence<Idx...>){
        static_assert(panoramix::core::Sequence<std::is_base_of<StepData, std::decay_t<typename std::tuple_element<Idx, ArgumentsTupleType>::type>>::value ...>::All,
            "all input types must be derived from StepData");
        ResultContentType d =
            updater(*dynamic_cast<std::decay_t<typename std::tuple_element<Idx, ArgumentsTupleType>::type>*>(dependencies.at(Idx)) ...);
        data->lock.lockForWrite();
        contentAs<ResultContentType>() = std::move(d);
        data->lock.unlock();
    }
};






class StepsDAG : public QObject {    
    Q_OBJECT

public:
    explicit StepsDAG(QObject * parent = 0);

public:
    size_t size() const { return _steps.size(); }

    template <class UpdateFunctionT>
    inline int addStep(const QString & name, UpdateFunctionT && updater,
        const std::vector<int> & dependencies = std::vector<int>()) {
        auto step = std::make_shared<StepWithTypedUpdater<std::decay_t<UpdateFunctionT>>>(std::forward<UpdateFunctionT>(updater));
        for (int d : dependencies){
            assert(d >= 0 && d < _steps.size());
        }
        step->data->setModified();
        _names.push_back(name);
        _steps.push_back(step);
        _widgets.push_back(step->data->createBindingWidgetAndActions(_actions));
        _dependencies.push_back(dependencies);
        _stepIds[step.get()] = _steps.size() - 1;
        return _steps.size() - 1;
    }

    inline int id(const StepPtr & s) const { return _stepIds.value(s.get(), -1); }
    inline int id(const Step * s) const { return _stepIds.value(s, -1); }

    inline const QString & stepNameAt(int id) const { return _names[id]; }

    inline const QList<StepWidgetInterface *> widgets() const { return _widgets; }
    inline const QList<QAction*> actions() const { return _actions; }

    inline const std::vector<int> & stepDependenciesAt(int id) const { return _dependencies[id]; }

    bool needsUpdate(int id) const;    
    struct UpdateCallback {
        virtual void beforeUpdate(int stepId) const = 0;
        virtual void afterUpdate(int stepId) const = 0;
    };

    void updateAll(UpdateCallback const * callback = nullptr, bool forceSourceStepUpdate = false);

    Q_SIGNAL void dataUpdated(int index);
    Q_SIGNAL void messageUpdated(QString msg);

private:
    Q_SLOT void updateWidget(int index);

private:
    QList<QString> _names;
    QList<StepPtr> _steps;

    QList<StepWidgetInterface *> _widgets;
    QList<QAction*> _actions;

    QList<std::vector<int>> _dependencies;
    QHash<Step const *, int> _stepIds;
};
