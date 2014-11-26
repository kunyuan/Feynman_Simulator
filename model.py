PI=3.141592653589793238

def J1J2(beta):
    return {"Model": "J1J2",
            "Hopping": 0.0,
            "ChemicalPotential": [PI/2.0/beta*1j,PI/2.0/beta*1j],
            }
