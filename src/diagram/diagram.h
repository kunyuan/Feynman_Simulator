//
//  diagram_global.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__diagram_global__
#define __Fermion_Simulator__diagram_global__

#include <iostream>
#include "logger.h"
#include "component_bundle.h"
#include "rng.h"

const int MAX_K=10000;

int RandomPickK();
int RandomPickdSpin();

class Worm{
public:
    bool Exist;
    int Ira, Masha;    //extra line: Ira---"k,dSpin"--->Masha
    int K;
    int dSpin;
    
    Worm():Exist(false), Ira(0),Masha(0),K(0),dSpin(0){}
    Worm(int ira, int masha, int dk, int s)
    : Exist(true), Ira(ira), Masha(masha), K(dk), dSpin(s)
    {
    }
};

class Diagram {
  public:
    Diagram();
    int Order;
    Complex Phase, Weight;
    Bundle<GLine> G;
    Bundle<WLine> W;
    Bundle<Vertex> Ver;
    
    Worm Worm;
    bool IsWorm(const Vertex&);
    bool CanNotMoveWorm(int, const Vertex&);

    spin Spin(GLine &, const int&);
    spin Spin(WLine &, const int&, const int&);
    spin Spin(Vertex&, const int&);
    
    int Sublattice(GLine &, const int&);
    int Sublattice(WLine &, const int&);
    int Sublattice(Vertex &);
    
    string PrettyString(GLine &);
    string PrettyString(WLine &);
    string PrettyString(Vertex &);

    //Randomly Pick
    GLine& RandomPickG();
    WLine& RandomPickW();
    Vertex& RandomPickVer();
    
    //Diagram
    Vertex &NeighVer(GLine &, int dir);
    Vertex &NeighVer(WLine &, int dir);
    GLine &NeighG(Vertex &, int dir);
    WLine &NeighW(Vertex &);
    bool FixDiagram();

    //Diagram Check
    bool CheckG();
    bool CheckW();
    bool CheckVer();
    bool CheckDiagram();

    //Diagram IO
    bool ReadDiagram(std::string);
    bool ReadDiagram(std::istream &);
    void WriteDiagram(std::string);
    void WriteDiagram(std::ostream &);
    void WriteDiagram2gv(std::string);
    void WriteDiagram2gv(std::ostream &);
    
private:
    std::ostream& Component2gv(std::ostream &, GLine &);
    std::ostream& Component2gv(std::ostream &, WLine &);
    std::ostream& Component2gv(std::ostream &, Vertex &);
    
    template <typename T>
    std::ostream& Bundle2gv(std::ostream &, Bundle<T> &);
};
#endif /* defined(__Fermion_Simulator__diagram_global__) */
