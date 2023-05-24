// Casey Berger
// Created: May 24 2023
// Last edited: May 24, 2023


namespace ising {
    int Lattice::getLength(){
        return length_ * length_;
    }
    void Lattice::setLength(int length){
        length_ = length;
    }
}