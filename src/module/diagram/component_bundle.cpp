//
//  component_bundle.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component_bundle.h"
#include "utility/abort.h"
#include "utility/rng.h"
using namespace std;
using namespace diag;

template <typename T>
Bundle<T>::Bundle(string bundle_name)
{
    _bundle_name = bundle_name;
    for (int i = 0; i < MAX_BUNDLE; i++) {
        _component_bundle[i].Name = i;
        _component_name[i] = _component_bundle + i;
    }
    _available_space = 0;
}

template <typename T>
string Bundle<T>::BundleName()
{
    return _bundle_name;
}

template <typename T>
T *Bundle<T>::Add()
{
    if (_available_space > MAX_BUNDLE)
        ABORT("Too many objects >=" << MAX_BUNDLE);
    T *address = _component_name[_available_space];

    //need to be checked
    address->Name = _available_space;
    //////////

    _available_space++;
    return address;
}

template <typename T>
void Bundle<T>::Add(T *Target)
{
    if (_available_space >= MAX_BUNDLE)
        ABORT("Too many objects >=" << MAX_BUNDLE);
    _component_name[_available_space] = Target;
    _available_space++;
}

template <typename T>
void Bundle<T>::Remove(name name)
{
    if (DEBUGMODE && _available_space == 0)
        ABORT("Nothing to delete!");
    _available_space--;
    T *last = _component_name[_available_space];
    T *target = _component_name[name];

    //switch name
    target->Name = _available_space;
    last->Name = name;

    _component_name[name] = last;
    _component_name[_available_space] = target;
    return;
}

template <typename T>
void Bundle<T>::Remove(T *target)
{
    if (DEBUGMODE) {
        if (!Exist(target))
            ABORT("The component to remove do not belong to the bundle!");
        if (_available_space == 0)
            ABORT("Nothing to delete!");
    }
    _available_space--;
    int name = target->Name;
    T *last = _component_name[_available_space];

    //switch name
    target->Name = _available_space;
    last->Name = name;

    _component_name[name] = last;
    _component_name[_available_space] = target;
    return;
}

template <typename T>
T &Bundle<T>::operator[](name name)
{
    if (DEBUGMODE && name >= _available_space)
        ABORT(name << " exceeds the storage!");
    return *_component_name[name];
}

template <typename T>
T *Bundle<T>::operator()(name name)
{
    if (DEBUGMODE && name >= _available_space)
        ABORT(name << " exceeds the storage!");
    return _component_name[name];
}

template <typename T>
int Bundle<T>::HowMany()
{
    return _available_space;
}

template <typename T>
void Bundle<T>::Recover(int step)
{
    _available_space += step;
    //TODO: please take a careful thought on Recover, does the component has the correct name?
}

template <typename T>
bool Bundle<T>::Exist(T *target)
{
    if (target != nullptr && target->Name < _available_space && target >= _component_bundle && target < _component_bundle + MAX_BUNDLE)
        return true;
    else
        return false;
}

template <typename T>
T *Bundle<T>::RandomPick(RandomFactory &RNG)
{
    return _component_name[RNG.irn(0, _available_space - 1)];
}

namespace diag {
template class Bundle<GLine>;
template class Bundle<WLine>;
template class Bundle<Vertex>;
}
