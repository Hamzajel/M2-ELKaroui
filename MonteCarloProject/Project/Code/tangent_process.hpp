#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

//#include "mc.hpp"
//#include "payoff.hpp"
//#include "pseudo_CEV.hpp"
//#include "BS.hpp"

using namespace std;

//Classe pour le calcul des greeks via malliavin, définition d'un fonction pour chaque greek 
//qui sort (sur base d'un bruit normal) le payoff simulé à prendre en compte dans le Monte Carlo 

// Greeks via méthode du flot tangent 

template <typename T_sim>
struct tangent_process{   
    tangent_process(double r,unsigned sample_size, double T, double K, T_sim sim, double size_h):
        r(r), sample_size(sample_size), T(T), K(K), sim(sim), size_h(size_h) {}

    double tan(double S_T, double Y_T, const string& type_payoff){
        return (der_payoff(S_T,K,type_payoff)*Y_T);
    }

    double delta(string type_payoff, vector<double> brow){
        if(type_payoff=="bin_put" || type_payoff=="bin_call"){
            cout<<"Erreur: méthode de flot tangent n'est pas disponible pour les options binaires"<<endl;
            exit(0);
        }
        vector<double> St_tan;
        St_tan=sim.St_tangent(brow);

        return (exp(-r*T)*(*this).tan(St_tan[0],St_tan[1],type_payoff));
    }

    double vega(string type_payoff, double* brow){
        if(type_payoff=="bin_put" || type_payoff=="bin_call"){
            cout<<"Erreur: méthode de flot tangent n'est pas disponible pour les options binaires"<<endl;
            exit(0);
        }
        double* St_tan;
        St_tan=sim.St_tangent(brow);

        return exp(-r*T)*(*this).tan(St_tan[0],St_tan[2],type_payoff);

    }    
    double lambda_delta(string type_payoff, double size_MC){
        if(type_payoff=="bin_put" || type_payoff=="bin_call"){
            cout<<"Erreur: méthode de flot tangent n'est pas disponible pour les options binaires"<<endl;
            exit(0);
        }
        vector<double> brow;
        vector<double> St_tan;

        double sum_v_delta=0;
        double sum_c_delta=0;
        double lambda_delta;

        int sign=1;
        if (type_payoff=="put"){
            sign=-1;
        } // si put alors c'est -1

        for(unsigned i=0;i<size_MC;i++){
            brow=normal_generation(size_h); 
            St_tan=sim.St_tangent(brow);

            sum_v_delta+=pow(sign*(exp(-r*T)*St_tan[1]-1),2);
            sum_c_delta+=exp(-r*T)*(*this).tan(St_tan[0],St_tan[1],type_payoff)*sign*(exp(-r*T)*St_tan[1]-1);


        }
        lambda_delta=sum_c_delta/sum_v_delta;

        return lambda_delta;
    }

    double delta_con_var(string type_payoff, vector<double> brow, double lambda_delta){
        if(type_payoff=="bin_put" || type_payoff=="bin_call"){
            cout<<"Erreur: méthode de flot tangent n'est pas disponible pour les options binaires"<<endl;
            exit(0);
        }
        int sign=1;
        if (type_payoff=="put"){
            sign=-1;
        } // si put alors c'est -1

        vector<double> St_tan;
        St_tan=sim.St_tangent(brow);

        return (exp(-r*T)*(*this).tan(St_tan[0],St_tan[1],type_payoff)-lambda_delta*sign*(exp(-r*T)*St_tan[1]-1));
    }

    private:
        unsigned sample_size;
        double r, T, K, size_h;
        std::normal_distribution<double> X;  
        T_sim sim;
        
};