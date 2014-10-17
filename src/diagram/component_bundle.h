//
//  component_bundle.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__component_bundle__
#define __Fermion_Simulator__component_bundle__

#include "component.h"
#include "utility.h"
#include "definition_global.h"
#include <string>

const int MAX_BUNDLE = 1024;
template <typename T>
class Bundle {
  private:
    T *_component_name[MAX_BUNDLE];
    T _component_bundle[MAX_BUNDLE];
    std::string _bundle_name;
    int _available_space;

  public:
    Bundle(std::string);
    T &Add();      //CAUTION: This function returns the address of the added object
    void Add(T &); //Add a copy of component to bundle
    std::string Name();
    void Remove(int name);
    void Remove(T &target);
    void Recover(int step);
    //step>0: Recover successive "step" steps of removing T*
    //step<0: Recover successive "step" steps of adding T*
    T &operator[](int name);
    int HowMany();
    bool Exist(int name);
    bool Exist(T &target);
};

int TestDiagram();
#endif /* defined(__Fermion_Simulator__component_bundle__) */
