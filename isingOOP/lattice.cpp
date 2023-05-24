// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023

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
    void initialize(){
        double ** grid = new double*[length_];
        for(int i = 0; i < length_; i++){
            grid[i] = new double[length_];
            for (int j = 0; j<length_; j++){
                grid[i][j] = new int;
                double r = ((double)rand())/((double)RAND_MAX);
                if (r<0.5)
                    {grid[i][j] = -1;}
                else
                    {grid[i][j] = 1;}
            }
        }
        grid_ = grid;
    }
}