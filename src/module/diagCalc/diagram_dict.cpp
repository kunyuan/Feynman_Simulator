//
// Created by kun on 10/5/17.
//
#include "diag_calculator.h"
#include "utility/dictionary.h"
#include <iostream>
using namespace diagCalc;

DiagramDict::DiagramDict(){

}

bool DiagramDict::FromDict(const Dictionary & dict){
    std::cout<<dict.HasKey("Sigma")<<std::endl;
    auto _Sigma=dict.Get<Dictionary>("Sigma");
    auto value=_Sigma.Get<int>("2");
    std::cout<<value<<std::endl;
//    _Sigma.PrettyString();
//    for(auto key=_Sigma.begin();key!=_Sigma.end();key++){
//        std::cout<<"hello"<<std::endl;
//        std::cout<<key<<std::endl;

//    }

}
