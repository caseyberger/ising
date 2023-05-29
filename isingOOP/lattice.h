// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023
//#include <vector>

#pragma once

namespace ising {
    class Lattice {
        public:
        //constructor
        Lattice(int length, double J, double kBT);
        //public members (can be accessed in other functions)
        void setLength(int length);
        void setTemperature(double kBT);
        void setInteractionJ(double J);
        int getLength();
        double getTemperature();
        double getInteractionJ();
        void initialize();
        void printLattice();
        int* getNeighbors(int i, int j);
        int* getNeighborSpins(int i, int j);
        double localEnergy(int i, int j);
        double Energy();
        double Magnetization();
        void metropolisLoop(int nMC, std::vector<double> &mc_E, std::vector<double> &mc_M);
        
        private: //used within the class itself only
        int length_;
        int **grid_;
        double interactionJ_;
        double kBT_;
        int plusOne_(int i);
        int minusOne_(int i);
        void flipSpin_(int i, int j);
        void sweepLattice_();
    };
}
