// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023
#include <iostream> //cout, endl
#include <iomanip> //setw
#include <vector> //vector
#include <cmath> //exp
#include <numeric> // iota
#include <list> //for lists needed in iota
#include <algorithm>  // shuffle
#include <random> //default_random_engine

#include "lattice.h"

namespace ising {
    //public functions
    Lattice::Lattice(int length, double J, double kBT){
        Lattice::setLength(length); //set length
        Lattice::setInteractionJ(J); //set J (do I need to do this the same way I do length?)
        Lattice::setTemperature(kBT);//see above note
    }
    void Lattice::setLength(int length){
        length_ = length;
    }
    
    void Lattice::setTemperature(double kBT){
        kBT_ = kBT;
    }
    
    void Lattice::setInteractionJ(double J){
        interactionJ_ = J;
    }
    
    int Lattice::getLength(){
        return length_;
    }
    
    double Lattice::getTemperature(){
        return kBT_;
    }
    double Lattice::getInteractionJ(){
        return interactionJ_;
    }
    void Lattice::initialize(){
        int ** grid = new int*[length_];
        for(int i = 0; i < length_; i++){
            grid[i] = new int[length_];
            for (int j = 0; j<length_; j++){
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
        double T = Lattice::getTemperature();
        double J = Lattice::getInteractionJ();
        std::cout << len << " x " << len << " lattice with T = " << T << " and J = " << J << std::endl;
        for(int i = 0; i < len; i++){
            for (int j = 0; j<len; j++){
                std::cout << std::setw(4) << grid_[i][j] << ",";
            }
            std::cout << std::endl;
        }
    }
    
    int* Lattice::getNeighbors(int i, int j){
        static int nn[8];
        nn[0] = Lattice::plusOne_(i);
        nn[1] = j;
        nn[2] = i;
        nn[3] = Lattice::plusOne_(j);
        nn[4] = Lattice::minusOne_(i);
        nn[5] = j;
        nn[6] = i;
        nn[7] = Lattice::minusOne_(j);
        return nn;
    }
    
    int* Lattice::getNeighborSpins(int i, int j){
        static int nnSpins[4];
        nnSpins[0] = grid_[Lattice::plusOne_(i)][j];
        nnSpins[1] = grid_[i][Lattice::plusOne_(j)];
        nnSpins[2] = grid_[Lattice::minusOne_(i)][j];
        nnSpins[3] = grid_[i][Lattice::minusOne_(j)];
        return nnSpins;
    }
    
    double Lattice::localEnergy(int i, int j){
        int *nnSpins = Lattice::getNeighborSpins(i,j);
        double nnsum = nnSpins[0] + nnSpins[1] + nnSpins[2] + nnSpins[3];
        double Eloc = (-1.*interactionJ_ * grid_[i][j] * nnsum);
        //std::cout << "Eloc = " << Eloc <<  std::endl;
        return Eloc;
    }
    
    double Lattice::Energy(){
        double Energy;
        for(int i = 0; i < length_; i++){
            for (int j = 0; j<length_; j++){
                Energy += Lattice::localEnergy(i,j);
            }
        }
        return Energy;
    }
    
    void Lattice::metropolisLoop(int nMC, std::vector<double> &mc_E){
        //adds energy for individual configurations to a vector so we can do stats for that temperature
        for (int n = 0; n < nMC; n++){
#ifdef TESTING_MODE
            std::cout << "Random sweep over lattice." << std::endl;
#endif 
            Lattice::sweepLattice_(); //iterate over the whole lattice, flipping spins if favorable
#ifdef TESTING_MODE
            std::cout << "Storing energy of current config in vector." << std::endl;
#endif
            mc_E.push_back(Lattice::Energy());
        }
    }
    
    //private functions
    int Lattice::plusOne_(int i){
        //int len = Lattice::getLength();
        if(i+1 == length_){ return 0;}
        else{return i+1;}
    }
    
    int Lattice::minusOne_(int i){
        if(i == 0){ return length_-1;}
        else{return i-1;}
    }
        
    void Lattice::flipSpin_(int i, int j){
        double E_init = 1.0*Lattice::localEnergy(i, j);
        double dE = -2.*E_init;
        double r = ((double) std::rand())/((double) RAND_MAX);
        if (dE < 0){
            grid_[i][j] = -1 * grid_[i][j];
        }
        else{
            if (r <= std::exp(-1.0*dE/kBT_)){
                grid_[i][j] = -1 * grid_[i][j];
            }
            else{
                dE = 0.0;
            }
        }
    }
    
    void Lattice::sweepLattice_(){
        //consider making two arrays for i and j, which are shuffled orders of the indices
        
        /*std::list<int> i_arr(length_);
        std::list<int> j_arr[length_];
        std::iota(i_arr.begin(), i_arr.end(), 1);
        std::iota(j_arr.begin(), j_arr.end(), 1);
        
        shuffle(i_arr.begin(), i_arr.end(), std::default_random_engine(1232));
        shuffle(j_arr.begin(), j_arr.end(), std::default_random_engine(7812));
        */
        
#ifdef TESTING_MODE
        std::cout << "Creating vectors to store i and j sites." << std::endl;
#endif
        std::vector<int> i_arr(length_);
        std::vector<int> j_arr(length_);
        std::iota(i_arr.begin(), i_arr.end(), 1);
        std::iota(j_arr.begin(), j_arr.end(), 1);
        
#ifdef TESTING_MODE
        std::cout << "Shuffling vectors with random default engine." << std::endl;
#endif      
        shuffle(i_arr.begin(), i_arr.end(), std::default_random_engine(1232));
        shuffle(j_arr.begin(), j_arr.end(), std::default_random_engine(8723));

#ifdef TESTING_MODE
        std::cout << "Iterating over i and j values and flipping spin." << std::endl;
#endif
        for(unsigned int i = 0; i < i_arr.size(); i++){
            for(unsigned int j = 0; i < j_arr.size(); j++){
                Lattice::flipSpin_(i_arr[i],j_arr[j]);
            }
        }
#ifdef TESTING_MODE
        std::cout << "End iteration." << std::endl;
#endif
#ifdef TESTING_MODE
        std::cout << "Clearing vectors to free up memory." << std::endl;
#endif
        i_arr.clear();
        j_arr.clear();
    }
}