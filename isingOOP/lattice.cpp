// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023
//test
#include "lattice.h"

namespace ising {
    int Lattice::getLength(){
        return length_;
    }
    void Lattice::setLength(int length){
        length_ = length;
    }
}