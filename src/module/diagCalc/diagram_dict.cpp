//
// Created by kun on 10/5/17.
//
#include "diag_calculator.h"
#include "utility/dictionary.h"
#include <iostream>
#include <array>
#include <string>
using namespace diagCalc;
using namespace std;

DiagramDict::DiagramDict(){

}

bool DiagramDict::FromDict(const Dictionary & dict){
//    auto diagram1=std::array<int, 2*MAX_ORDER>({2,1,3,4});
//    auto diagram2=std::array<int, 2*MAX_ORDER>({1,3,2,4});
//    auto spin1=std::array<int, 2*MAX_ORDER>({1,0,1,0});
//    auto spin2=std::array<int, 2*MAX_ORDER>({0,1,0,0});
//    auto FermiSign1=-1;
//    auto FermiSign2=1;
//
//    AllDiagramConfig.push_back(std::vector<std::array<int, 2*MAX_ORDER>>({diagram1, diagram2}));
//    AllSpinConfig.push_back(std::vector<std::array<int, 2*MAX_ORDER>>({spin1, spin2}));
//    AllFermiSignConfig.push_back(std::vector<int>({FermiSign1, FermiSign2}));

    auto _Sigma=dict.Get<Dictionary>("Sigma");
    for(auto order=0;order<=MAX_ORDER;order++){
        if(_Sigma.HasKey(to_string(order))){
            auto _Dict=_Sigma.Get<Dictionary>(std::to_string(order));
            DiagramsList _diagramList;
            for (auto d: _Dict.Get<vector<vector<int>>>("Permutations")){
                array<int, 2*MAX_ORDER> _diagram;
                for(auto i=0;i<d.size();i++){
                    _diagram[i]=d[i];
                }
                _diagramList.push_back(_diagram);
            }
            AllDiagramConfig.push_back(_diagramList);

            SpinsList _spinsList;
            for (auto d: _Dict.Get<vector<vector<int>>>("Spins")){
                array<int, 2*MAX_ORDER> _spins;
                for(auto i=0;i<d.size();i++){
                    _spins[i]=d[i];
                }
                _spinsList.push_back(_spins);
            }
            AllSpinConfig.push_back(_spinsList);

            FermiSignList _signList;
            for (auto d: _Dict.Get<vector<int>>("FermiSigns")){
                _signList.push_back(d);
            }
            AllFermiSignConfig.push_back(_signList);
        }
    }
    cout<<AllDiagramConfig[0][0][0]<<","<<AllSpinConfig[0][0][1]<<endl;

    return true;
}

