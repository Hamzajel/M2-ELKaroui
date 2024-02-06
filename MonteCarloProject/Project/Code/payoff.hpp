
#pragma once

using namespace std;

double payoff(double S, double K, const string& type_payoff){
    if(type_payoff=="call"){
        if(S>K){
            return S-K;
        }
        else return 0;
    }
    else if (type_payoff=="put"){
        if(S>K){
            return 0;
        }
        else return K-S;
        
    }
    else if (type_payoff=="bin_call"){
        if (S>K){
            return 1;
        }
        else return 0;
    }
    else if(type_payoff=="bin_put"){

        if (S>K){
            return 0;
        }
        else return 1;
    }
    return 0;
}
double der_payoff(double S, double K, const string& type_payoff){
    if(type_payoff=="call"){
        if(S>K){
            return 1;
        }
        else return 0;
    }
    else if (type_payoff=="put"){
        if(S>K){
            return 0;
        }
        else return -1;
        
    }
    else if (type_payoff=="bin_call"){
        if (S==K){
            return 1;
        }
        else return 0;
    }
    else if(type_payoff=="bin_put"){

        if (S==K){
            return -1;
        }
        else return 0;
    }
    return 0;
}
