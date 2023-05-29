// Casey Berger
// Created: May 7 2014
// Last edited: May 24, 2023 (OOP)

// include files
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdio.h>
#include <sstream>
#include <time.h>

//custom headers
#include "lattice.h"

using namespace std;
using ising::Lattice;

//function declaration
void read_in_inputs(int argc, char *argv[],int &len, double &J, double &Tmax, double &Tmin, double &dT, int &nMC);
void write_to_file(vector<double> &mc_E, vector<double> &mc_M, double T, double J, int len);
int main (int argc, char *argv[])
{
    int len,nMc;//length of lattice, number of monte carlo iterations
    double J, Tmax, Tmin, dT; //interaction strength, max temp, min temp, temperature iterator
    
    read_in_inputs(int argc, char *argv[],int &len, double &J, double &Tmax, double &Tmin, double &dT, int &nMC);
    
    Lattice L(len, J, Tmax);
    srand(7); //seed random number
    double T = Tmax; //current temp
    //temperature loop
    while (T>=Tmin){
        L.setTemperature(T);
        L.initialize();
#ifdef TESTING_MODE
        cout << "Lattice initialized"<<endl;
        L.printLattice();
#endif 
        vector<double> mc_E; //stores energies for monte carlo loop at one T
        vector<double> mc_M; //stores magnetizations for monte carlo loop at one T
        L.metropolisLoop(nMC, mc_E, mc_M);
        write_to_file(mc_E, mc_M, T, J, len);
        T -= dT;
    }
    return 0;
}

void read_in_inputs(int argc, char *argv[],int &len, double &J, double &Tmax, double &Tmin, double &dT, int &nMC)
{
    //read in parameters
#ifdef TESTING_MODE
    cout << "Function: read_in_inputs" << endl;
#endif
    string str, filename;
    const int n_params = 6;
    string inputs [n_params] = {"L","interactionJ", "Tmax","Tmin","dT","nMC"};//read in keywords for parameters
    if (argc != 2){ //exits if input file is not given
        cerr << "Usage: ./ising input.txt"<< endl << "Exiting program" << endl;
        exit(10);
    }
    else{
        ifstream input_file(argv[1]);
        if (!input_file.is_open()){
            cerr << "input file cannot be opened";
            exit(10);
        }
        else{
            int count = 0;
#ifdef TESTING_MODE
            cout << "Starting param search in file: ";
            for (int n=0; n<n_params; n++){
                cout << inputs[n] << ',';
            }
            cout << endl;
#endif  
            while (count < n_params) {
                while (getline(input_file, str)) {
                    //search for params in input
                    size_t found = str.find(inputs[count]);
                    size_t start;
                    if (found != string::npos) {
                        start = str.find_last_of(' ');
                        inputs[count] = str.substr(start + 1);
                        count++;
                    }
                    else{
                        //if your inputs file doesn't have that parameter listed 
                        cerr << "parameter "<< inputs[count] << " not in input file.";
                        exit(10);
                    }
                }
            }
            len = stod(inputs[0]);
            J = stod(inputs[1]);
            Tmax = stod(inputs[2]);
            Tmin = stod(inputs[3]);
            dT = stod(inputs[4]);
            nMC = stod(inputs[5]);
#ifdef TESTING_MODE
            cout << "parameters acquired: ";
            for (int n=0; n<n_params; n++){
                cout << inputs[n] << ',';
            }
            cout << endl;
#endif  
        }
    }
}

void write_to_file(vector<double> &mc_E, vector<double> &mc_M, double T, double J, int len)
{
    //cout << "write_to_file" << endl;
    //output both solutions to a .csv
    
    std::stringstream T_stream, J_stream, L_stream;
    T_stream << std::fixed << std::setprecision(3) << T; //truncate T for filename
    J_stream << std::fixed << std::setprecision(3) << J; //truncate J for filename
    L_stream << std::fixed << std::setprecision(1) << len; //truncate len for filename
    std::string str_T = T_stream.str();
    std::string str_J = J_stream.str();
    std::string str_L = L_stream.str();
    std::string fname = "data/ising_MC_T_" + str_T + "_J_" + str_J +"_L_" + str_L + ".csv";//output filename
    ofstream fout; //output stream
    fout.open(fname.c_str(),ios::out);
    
    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << "step,E,M,T,J,L" << endl;
    for (unsigned int n = 0; n<mc_E.size(); n++)
    {
        fout.setf(ios::fixed);
        fout << n << "," << mc_E[n] << ","<< mc_M[n] << "," << T << "," << J << "," << len << endl;
    }
    fout.close();
    mc_E.clear();
    mc_M.clear();
}