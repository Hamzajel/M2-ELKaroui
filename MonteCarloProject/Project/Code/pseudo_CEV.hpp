
#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>


//#include "payoff.hpp"


//Classe pour le modéle Pseudo-CEV qui contient différentes fonctions utilisées dans les classes des méthodes 
//pour calculer les greeks, cette classe contient ainsi des fonctions qui vont permettre (avec un bruit Normal en entrée) 
//de simuler les payoffs nécéssaires pour le calcul Monte Carlo 


struct pseudo_CEV{

    pseudo_CEV(double X0,double K,double r,double sigma, double vega, double T, double h)
    : X0(X0),K(K), r(r), sigma(sigma),vega(vega),T(T),h(h){}

    [[nodiscard]] double value(double S, double sig,double v, double brow) const{
        return (S+r*S*h+(pow(S,(v+1))/(sqrt(1+pow(S,2))))*sig*sqrt(h)*brow);
        
    }
    [[nodiscard]] double value_vega(double S, double sig,double v, double brow,double eps) const{
        return (S+r*S*h+(pow(S,(v+1))/((sqrt(1+pow(S,2))))*(1)*sig+S*eps)*sqrt(h)*brow);
        
    }
    [[nodiscard]] double value_tan_delta(double S,double Y, double sig,double v, double brow) const{
        return(Y+r*Y*h+(((1+v)*pow(S,v)*sqrt(1+pow(S,2))-pow(S,2+v)/(sqrt(1+pow(S,2))))/(1+pow(S,2)))*Y*sig*sqrt(h)*brow);
        
    }
    [[nodiscard]] double value_tan_gamma(double S,double Y,double Z, double sig,double v, double brow) const{
        double first_der_sigma,second_der_sigma;
        first_der_sigma=(((1+v)*pow(S,v)*sqrt(1+pow(S,2))-pow(S,2+v)/(sqrt(1+pow(S,2))))/(1+pow(S,2)))*sig;
        second_der_sigma=sig*(1/(1+pow(S,2))*((1+v)*v*pow(S,v-1)*sqrt(1+pow(S,2))+(1+v)*pow(S,v+1)/(sqrt(1+pow(S,2)))-((2+v)*pow(S,v+1)*sqrt(1+pow(S,2))-pow(S,3+v)/(sqrt(1+pow(S,2))))/(1+pow(S,2)))+((1+v)*pow(S,v)*sqrt(1+pow(S,2))-pow(S,2+v)/(sqrt(1+pow(S,2))))*2*S/pow(1+pow(S,2),2));
        return (Z+r*Z*h+(first_der_sigma*Z+second_der_sigma*pow(Y,2))*sqrt(h)*brow);
    }
    /*[[nodiscard]] double inv_value_tan_delta(double S,double Y, double sig,double v, double brow) const{
        double sigma_prime=(((1+v)*pow(S,v)*sqrt(1+pow(S,2))-pow(S,2+v)/(sqrt(1+pow(S,2))))/(1+pow(S,2)))*sig;
        return Y+Y*(-r+pow(sigma_prime,2))*h-Y*sigma_prime*sqrt(h)*brow;
    }*/
    [[nodiscard]] double value_tan_vega(double S,double Y, double sig,double v, double brow) const{
        double gamma_prime=S;
        return Y+r*Y*h+(gamma_prime+(((1+v)*pow(S,v)*sqrt(1+pow(S,2))-pow(S,2+v)/(sqrt(1+pow(S,2))))/(1+pow(S,2)))*Y*sig)*sqrt(h)*brow;
        
    }

    [[nodiscard]] vector<double> St_tangent(vector<double> brow) const{
        vector<double> state(3);
        double size_d=T/h;
        int size= static_cast<int>(size_d);

        state[0]=X0,state[1]=1,state[2]=0;
            for(int i=1;i<size;++i){
                state[1]=(*this).value_tan_delta(state[0],state[1],sigma,vega,brow[i-1]);
                state[2]=(*this).value_tan_vega(state[0],state[2],sigma,vega,brow[i-1]);
                state[0]=(*this).value(state[0],sigma,vega,brow[i-1]);
            }
        return state;

    }

    [[nodiscard]] vector<double> St_bis(double eps, const vector<double>& brow) const{
        vector<double> state(5);
        double size_d=T/h;
        int size= static_cast<int>(size_d);

        state[0]=X0,state[1]=X0+eps,state[2]=X0-eps,state[3]=X0,state[4]=X0;
            for(int i=1;i<size;++i){
                state[0]=(*this).value(state[0],sigma,vega,brow[i]);
                state[1]=(*this).value(state[1],sigma,vega,brow[i]);
                state[2]=(*this).value(state[2],sigma,vega,brow[i]);
                state[3]=(*this).value_vega(state[3],sigma,vega,brow[i],eps);
                state[4]=(*this).value_vega(state[4],sigma,vega,brow[i],-eps);

            }
        return state;
    }
    
    
    [[nodiscard]] vector<double> payoff_diff_fin(const string& type_payoff, double eps, const vector<double>& brow) const{
        vector<double> pay(5);
        auto St=(*this).St_bis(eps,brow);
        pay[0]=payoff(St[0],K,type_payoff)*exp(-r*T);
        pay[1]=payoff(St[1],K,type_payoff)*exp(-r*T);
        pay[2]=payoff(St[2],K,type_payoff)*exp(-r*T);
        pay[3]=payoff(St[3],K,type_payoff)*exp(-r*T);
        pay[4]=payoff(St[4],K,type_payoff)*exp(-r*T);
        return pay;
    }
    [[nodiscard]] vector<double> malliavin_weight_delta(const vector<double>& brow) const{
        double size_d=T/h,Yt,St;
        int size= static_cast<int>(size_d);
        St=X0,Yt=1;
        double skorohod_int=0;
            for(int i=1;i<size;++i){
                Yt=(*this).value_tan_delta(St,Yt,sigma,vega,brow[i-1]);
                St=(*this).value(St,sigma,vega,brow[i-1]);
                skorohod_int+=(Yt/((pow(St,(vega+1))/(sqrt(1+pow(St,2))))*sigma))*sqrt(h)*brow[i-1];
                
            }
        skorohod_int*=(1/T);
        vector<double> out(2);
        out[0]=skorohod_int;
        out[1]=St;
        return out;
    }
    [[nodiscard]] vector<double> malliavin_weight_vega_new_TOUZI01(vector<double> brow) const{
        double size_d=T/h,skorohod_int=0,sigma_tilde;
        int size= static_cast<int>(size_d);
        double Yt[size],Zt[size],St[size],beta_bis[size],inv_Yt[size];

        St[0]=X0,Yt[0]=1,Zt[0]=0,inv_Yt[0]=1,brow[0]=0;
        beta_bis[0]=Zt[0]*inv_Yt[0];
        for(int i=1;i<size;i++){
            Zt[i]=(*this).value_tan_vega(St[i-1],Zt[i-1],sigma,vega,brow[i-1]); //derivate with respect to sigma
            St[i]=(*this).value(St[i-1],sigma,vega,brow[i-1]);
            sigma_tilde=(pow(St[i],(vega+1))*sigma)/(sqrt(1+pow(St[i],2)));
            skorohod_int+=(1/sigma_tilde)*Zt[i]*sqrt(h)*brow[i];
        }

        vector<double> out(2);
        out[0]=skorohod_int;
        out[1]=St[size-1];
        return out;
    }
    [[nodiscard]] vector<double> malliavin_weight_gamma(const vector<int>& brow) const{
        double size_d=T/h,Yt,St,Zt,sigma_tilde,sigma_prime;
        int size= static_cast<int>(size_d);
        double W_delta[size],der_W_delta[size];
        St=X0,Yt=1,Zt=0;
        W_delta[0]=0;
        double skorohod_int_delta=0,skorohod_int_gamma=0;
            for(int i=1;i<size;++i){
                //brow[i-1]=X(gen);
                Zt=(*this).value_tan_gamma(St,Yt,Zt,sigma,vega,brow[i-1]);
                Yt=(*this).value_tan_delta(St,Yt,sigma,vega,brow[i-1]);
                St=(*this).value(St,sigma,vega,brow[i-1]);
                sigma_tilde=(pow(St,(vega+1))*sigma)/(sqrt(1+pow(St,2))); 
                sigma_prime=(((1+vega)*pow(St,vega)*sqrt(1+pow(St,2))-pow(St,2+vega)/(sqrt(1+pow(St,2))))/(1+pow(St,2)))*sigma;
                //der_Yt=sigma_prime*Yt;
                W_delta[i]=(Yt/sigma_tilde*T);
                der_W_delta[i]=((Zt/sigma_tilde)-(Yt*sigma_prime)/pow(sigma_tilde,2))/T;
                skorohod_int_delta+=W_delta[i]*sqrt(h)*brow[i-1];
            }
        for(int i=1;i<size;i++){
            skorohod_int_gamma+=(der_W_delta[i])*sqrt(h)*brow[i-1];
        }

        skorohod_int_gamma+=pow(skorohod_int_delta,2);
        vector<double> out(2);
        out[0]=skorohod_int_gamma;
        out[1]=St;
        return out;
    }


    private:
    double X0,r,sigma,vega,T,h,K;
    //double ST_tan[3]{};


};