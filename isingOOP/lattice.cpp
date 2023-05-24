// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023
#include <iostream> //cout, endl
#include <iomanip> //setw

#include "lattice.h"

namespace ising {
    //public functions
    Lattice::Lattice(int length, double J){
        Lattice::setLength(length); //set length
        interactionJ_ = J; //set J (do I need to do this the same way I do length?)
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
                std::cout << std::setw(4) << grid_[i][j] << ",";
            }
            std::cout << std::endl;
        }
    }
    
    double Lattice::localEnergy(int i, int j){
        int *nnSpins = Lattice::getNeighborSpins_(i,j);
        int **nn = Lattice::getNeighbors_(i,j);
        std::cout << "i,j = " << i << "," << j << std::endl;
        std::cout << "nn = " << std::endl;
        std::cout<< nn[0][0] << "," << nn[0][1] << std::endl;
        std::cout<< nn[1][0] << "," << nn[1][1] << std::endl;
        std::cout<< nn[2][0] << "," << nn[2][1] << std::endl;
        std::cout<< nn[3][0] << "," << nn[3][1] << std::endl;
        
        double nnsum = nnSpins[0] + nnSpins[1] + nnSpins[2] + nnSpins[3];
        return - interactionJ * grid_[i][j] * nnsum;
    }
    
    
    //private functions
    int Lattice::plusOne_(int i){
        //int len = Lattice::getLength();
        if(i == length_){ return 0;}
        else{return i+1;}
    }
    
    int Lattice::minusOne_(int i){
        if(i == 0){ return length_-1;}
        else{return i-1;}
    }
    
    int** Lattice::getNeighbors_(int i, int j){
        int nn[4][2];
        nn[0][0] = Lattice::plusOne(i);
        nn[0][1] = j;
        nn[1][0] = i;
        nn[1][1] = Lattice::plusOne(j);
        nn[2][0] = Lattice::minusOne(i);
        nn[2][1] = j;
        nn[3][0] = i;
        nn[3][1] = Lattice::minusOne(j);
        return nn;
    }
    
    int* Lattice::getNeighborSpins_(int i, int j){
        int nnSpins[4];
        nnSpins[0] = grid_[Lattice::plusOne(i)][j];
        nnSpins[1] = grid_[i][Lattice::plusOne(j)];
        nnSpins[2] = grid_[Lattice::minusOne(i)][j];
        nnSpins[3] = grid[i][Lattice::minusOne(j)];
        return nnSpins;
    }
}