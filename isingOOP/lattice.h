// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023

#pragma once

namespace ising {
    class Lattice {
        public:
        //constructor
        Lattice(int length);
        //public members (can be accessed in other functions)
        void setLength(int length);
        int getLength();
        void initialize();
        void printLattice();
        
        private: //used within the class itself only
        int length_;
        int **grid_;
    };
}
