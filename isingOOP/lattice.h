// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023
// using ising::Lattice

#pragma once

namespace ising {
    class Lattice {
        //constructor to come
        //public members (can be accessed in other functions)
        public:
            int getLength();
            void setLength(int length);
        
        private: //used within the class itself only
            int length_;
    };
}
