//
//  parameter.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "job.h"
//using namespace std;

void Jobs::Read()
{
    Beta = 1.0;
    Order = 1;
    Seed = 1;
    Type = MC;
    return;
}
//6    #pid
//8,8    #L
//1.0    #Jcp
//1.0    #iniBeta
//0.0    #dBeta
//1.0    #finalBeta
//4    #Order
//false    #IsLoad
//2    #Type: MC
//10000    #Toss
//5000000    #Sample
//10    #Sweep
//-573165127    #Seed
//1.00_4_coll
//0.1    #Worm/Norm
//1.5,1.0,3.0,4.0    #Reweight
