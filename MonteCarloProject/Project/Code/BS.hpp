#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>


//#include "payoff.hpp"

//Classe pour le modéle Black and Scholes qui contient différentes fonctions utilisées dans les classes des méthodes 
//pour calculer les greeks, cette classe contient ainsi des fonctions qui vont permettre (avec un bruit Normal en entrée) 
//de simuler les payoffs nécéssaires pour le calcul Monte Carlo 

struct simu_payoff_BS{
    simu_payoff_BS(double S, double r, double sigma, double T, double K):// aujouter un arg pour le choix de la dynamique de S (B&S,...)
    S(S), r(r), sigma(sigma), T(T), K(K){}

    [[nodiscard]] double St(double X0,double s,double brow) const{
        return (X0*exp(s*sqrt(T)*brow+(r-(pow(s,2))/2)*T));
    }
    [[nodiscard]] vector<double> St_tangent(vector<double> brow) const{
        vector<double> ST_tan(3);
        ST_tan[0]=S*exp(sigma*sqrt(T)*brow[0]+(r-(pow(sigma,2))/2)*T);
        ST_tan[1]=ST_tan[0]/S;
        ST_tan[2]=ST_tan[0]*(brow[0]*sqrt(T)-sigma*T);
        
        return ST_tan;

    }

    double* St_bis(const vector<double>& brow,double eps){
            ST[0]=St(S,sigma,brow[0]);
            ST[1]=St(S+eps,sigma,brow[0]); 
            ST[2]=St(S-eps,sigma,brow[0]);
            ST[3]=St(S,sigma+eps,brow[0]);
            ST[4]=St(S,sigma-eps,brow[0]);
        
        return ST;
    }
    vector<double> payoff_diff_fin(const string& type_payoff, double eps, const vector<double>& brow){
        vector<double> pay(5);
        double* St=(*this).St_bis(brow,eps);
        pay[0]=payoff(St[0],K,type_payoff)*exp(-r*T);
        pay[1]=payoff(St[1],K,type_payoff)*exp(-r*T);
        pay[2]=payoff(St[2],K,type_payoff)*exp(-r*T);
        pay[3]=payoff(St[3],K,type_payoff)*exp(-r*T);
        pay[4]=payoff(St[4],K,type_payoff)*exp(-r*T);

        return pay;
    }
    [[nodiscard]] vector<double> malliavin_weight_delta(const vector<double>& brow) const{
        vector<double> out(2);
        out[0]=(sqrt(T)*brow[0])/(S*sigma*T);
        out[1]=(*this).St(S,sigma,brow[0]);
        return out;
    }
    [[nodiscard]] vector<double> malliavin_weight_gamma(const vector<double>& brow) const{
        vector<double> out(2);
        out[0]=((pow(sqrt(T)*brow[0],2)/(sigma*T))-(sqrt(T)*brow[0])-(1/sigma))/(pow(S,2)*sigma*T);
        out[1]=(*this).St(S,sigma,brow[0]);
        return out;
    }
    [[nodiscard]] vector<double> malliavin_weight_vega_new_TOUZI01(const vector<double>& brow) const{
        vector<double> out(2);
        out[0]=((pow(sqrt(T)*brow[0],2)/(sigma*T))-(sqrt(T)*brow[0])-(1/sigma));
        out[1]=(*this).St(S,sigma,brow[0]);
        return out;
    }
    protected:
    double S,r,sigma,T,K;
    double ST[5]{};
};