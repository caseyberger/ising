// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023
#include <cstdlib.h>
#include <iostream>

#include "lattice.h"

namespace ising {
    Lattice::Lattice(int length){
        Lattice::setLength(length); //set length
    }
    void Lattice::setLength(int length){
        length_ = length;
    }
    int Lattice::getLength(){
        return length_;
    }
    void Lattice::initialize(){
        int len = Lattice::getLength();
        int ** grid = new int*[len];
        for(int i = 0; i < len; i++){
            grid[i] = new int[len];
            for (int j = 0; j<len; j++){
                double r = ((double) std::rand())/((double) RAND_MAX);
                if (r<0.5)
                    {grid[i][j] = -1;}
                else
                    {grid[i][j] = 1;}
            }
        }
        grid_ = grid;
    }
    void Lattice::printLattice(){
        int len = Lattice::getLength();
        for(int i = 0; i < len; i++){
            for (int j = 0; j<len; j++){
                std::cout << grid_[i][j] << ",";
            }
            std::endl;
        }
    }
    
}