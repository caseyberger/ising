// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023

#pragma once

namespace ising {
    class Lattice {
        public:
        //constructor
        Lattice(int length, double J, double kBT_);
        //public members (can be accessed in other functions)
        void setLength(int length);
        int getLength();
        void initialize();
        void printLattice();
        int* getNeighbors_(int i, int j);
        int* getNeighborSpins_(int i, int j);
        double localEnergy(int i, int j);
        double Energy();
        void flip_spin(int i, int j);
        void metropolisLoop(int nMC, vector<double> &mc_E);
        
        private: //used within the class itself only
        int length_;
        int **grid_;
        double interactionJ_;
        double kBT_;
        int plusOne_(int i);
        int minusOne_(int i);
        void sweepLattice_();
    };
}
