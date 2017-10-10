//
// Created by kun on 10/5/17.
//
#include "diag_calculator.h"
#include "utility/dictionary.h"
#include <iostream>
#include <array>
using namespace diagCalc;

DiagramDict::DiagramDict(){

}

bool DiagramDict::FromDict(const Dictionary & dict){
//    std::cout<<dict.HasKey("Sigma")<<std::endl;
//    auto _Sigma=dict.Get<Dictionary>("Sigma");
//    auto value=_Sigma.Get<int>("2");
//    std::cout<<value<<std::endl;
//    _Sigma.PrettyString();
//    for(auto key=_Sigma.begin();key!=_Sigma.end();key++){
//        std::cout<<"hello"<<std::endl;
//        std::cout<<key<<std::endl;

//    }
    auto diagram1=std::array<int, 2*MAX_ORDER>({2,1,3,4});
    auto diagram2=std::array<int, 2*MAX_ORDER>({1,3,2,4});
    auto spin1=std::array<int, 2*MAX_ORDER>({1,0,1,0});
    auto spin2=std::array<int, 2*MAX_ORDER>({0,1,0,0});
    auto FermiSign1=-1;
    auto FermiSign2=1;

    AllDiagramConfig.push_back(std::vector<std::array<int, 2*MAX_ORDER>>({diagram1, diagram2}));
    AllSpinConfig.push_back(std::vector<std::array<int, 2*MAX_ORDER>>({spin1, spin2}));
    AllFermiSignConfig.push_back(std::vector<int>(FermiSign1, FermiSign2));
//#define DiagramsList std::vector<std::array<int, 2*MAX_ORDER>>
//#define DiagramsList std::vector<std::array<int, 2*MAX_ORDER>>

    return true;
}
