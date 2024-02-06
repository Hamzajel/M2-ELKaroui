
#pragma once

#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>

//#include "mc.hpp"
//#include "payoff.hpp"
//#include "pseudo_CEV.hpp"
//#include "BS.hpp"

//Classe pour le calcul des greeks via malliavin, définition d'un fonction pour chaque greek 
//qui sort (sur base d'un bruit normal) le payoff simulé à prendre en compte dans le Monte Carlo 

using namespace std;

template <typename T_sim>
struct malliavin{
    malliavin() = default;
    malliavin(double r,unsigned sample_size, double T, double K, T_sim sim, double size_h):
        r(r), sample_size(sample_size), T(T), K(K), sim(sim), size_h(size_h) {}



    double delta(string type_payoff, vector<double> brow){
        vector<double> mal;
        mal=sim.malliavin_weight_delta(brow);

        return (exp(-r*T)*payoff(mal[1],K,std::move(type_payoff))*mal[0]);
    }
    double vega(string type_payoff, vector<double> brow){
        vector<double> mal;
        mal=sim.malliavin_weight_vega_new_TOUZI01(brow);
        return (exp(-r*T)*payoff(mal[1],K,std::move(type_payoff))*mal[0]);
    }

    double gamma(string type_payoff, vector<double> brow){
        vector<double> mal;
        mal=sim.malliavin_weight_gamma(brow);
        return (exp(-r*T)*payoff(mal[1],K,std::move(type_payoff))*mal[0]);
    }

    double lambda_delta(const string& type_payoff, double size_MC){
        vector<double> mal;
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
            mal=sim.malliavin_weight_delta(brow);

            sum_v_delta+=pow(sign*(exp(-r*T)*St_tan[1]-1),2);
            sum_c_delta+=exp(-r*T)*(payoff(mal[1],K,type_payoff)*mal[0])*sign*(exp(-r*T)*St_tan[1]-1);

        }
        lambda_delta=sum_c_delta/sum_v_delta;

        return lambda_delta;
    }

    double delta_con_var(string type_payoff, vector<double> brow, double lambda_delta){
        int sign=1;
        if (type_payoff=="put"){
            sign=-1;
        } // si put alors c'est -1

        vector<double> mal;
        vector<double> St_tan;

        mal=sim.malliavin_weight_delta(brow);
        St_tan=sim.St_tangent(brow);

        return (exp(-r*T)*(payoff(mal[1],K,type_payoff)*mal[0])-lambda_delta*sign*(exp(-r*T)*St_tan[1]-1));
    }
    private:
        unsigned sample_size{};
        double r{}, T{}, K{}, size_h{};
        std::normal_distribution<double> X;  
        T_sim sim;
};