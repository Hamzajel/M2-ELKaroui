#pragma once

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

// Greeks via différence finie (pas constant et décroissant)
template <typename T_sim> 
struct diff_finie{
    diff_finie(T_sim sim,double eps, double size_h, double r, double T):
          sim(sim), eps(eps), size_h(size_h), r(r), T(T) {}

    double first_order(double F_up, double F_down,double epsilon){
        return((F_up-F_down)/(2*epsilon));
    }
    double second_order(double F_up, double F, double F_down, double epsilon){
        return((F_up+F_down-2*F)/(pow(epsilon,2)));

    }

    double delta(string type_payoff, vector<double> brow){
        vector<double> payoff_vec;
        payoff_vec= static_cast<const vector<double>>(sim.payoff_diff_fin(type_payoff, eps, brow));

        return (*this).first_order(payoff_vec[1],payoff_vec[2],eps);
    }
    double vega(string type_payoff, vector<double> brow){
        vector<double> payoff_vec;
        payoff_vec=sim.payoff_diff_fin(type_payoff,eps,brow);

        return (*this).first_order(payoff_vec[3],payoff_vec[4],eps);
    }
    double gamma(string type_payoff, vector<double> brow){
        vector<double> payoff_vec;
        payoff_vec=sim.payoff_diff_fin(type_payoff,eps,brow);

        return (*this).second_order(payoff_vec[1],payoff_vec[0],payoff_vec[2],eps);
    }

    double lambda_delta(string type_payoff, double size_MC){
        vector<double> brow;
        vector<double> St_tangent;
        vector<double> payoff_vec;

        double sum_v_delta=0;
        double sum_c_delta=0;
        double lambda_delta;

        int sign=1;
        if (type_payoff=="put"){
            sign=-1;
        } // si put alors c'est -1

        for(unsigned i=0;i<size_MC;i++){
            brow=normal_generation(size_h); 

            payoff_vec=sim.payoff_diff_fin(type_payoff,eps,brow);
            St_tangent=sim.St_tangent(brow);

            sum_v_delta+=pow(sign*(exp(-r*T)*St_tangent[1]-1),2);
            sum_c_delta+=(*this).first_order(payoff_vec[1],payoff_vec[2],eps)*sign*(exp(-r*T)*St_tangent[1]-1);

        }
        lambda_delta=sum_c_delta/sum_v_delta;

        return lambda_delta;
    }
    
    double delta_con_var(string type_payoff, vector<double> brow, double lambda_delta){
        int sign=1;
        if (type_payoff=="put"){
            sign=-1;
        } // si put alors c'est -1

        vector<double> St_tan;
        vector<double> payoff_vec;

        payoff_vec=sim.payoff_diff_fin(type_payoff,eps,brow);
        St_tan=sim.St_tangent(brow);

        return (*this).first_order(payoff_vec[1],payoff_vec[2],eps)-lambda_delta*sign*(exp(-r*T)*St_tan[1]-1);
    }
    
    vector<mean_var> MC_decrease(unsigned sample_size, string type_payoff, double init_eps, int h){
        double sum_x_delta=0, sum_xx_delta=0,sum_x_gamma=0,sum_xx_gamma=0,sum_x_vega=0,sum_xx_vega=0; 
        vector<double> payoff_vec;
        vector<mean_var> greek_dif(3);
        double epsilon;
        vector<double> brow;
        for(unsigned i=0;i<sample_size;++i){
            epsilon=pow((init_eps+i+1),(-0.25));
            brow=normal_generation(h);
            //delta
            payoff_vec=sim.payoff_diff_fin(type_payoff,epsilon,brow);
            sum_x_delta+=(*this).first_order(payoff_vec[1],payoff_vec[2],epsilon);
            sum_xx_delta+=pow((*this).first_order(payoff_vec[1],payoff_vec[2],epsilon),2);
            //gamma
            sum_x_gamma+=(*this).second_order(payoff_vec[1],payoff_vec[0],payoff_vec[2],epsilon);
            sum_xx_gamma+=pow((*this).second_order(payoff_vec[1],payoff_vec[0],payoff_vec[2],epsilon),2);
            //vega
            sum_x_vega+=(*this).first_order(payoff_vec[3],payoff_vec[4],epsilon);
            sum_xx_vega+=pow((*this).first_order(payoff_vec[3],payoff_vec[4],epsilon),2);
            }
        greek_dif[0]=mean_var(sample_size,sum_x_delta,sum_xx_delta);
        greek_dif[1]=mean_var(sample_size,sum_x_gamma,sum_xx_gamma);
        greek_dif[2]=mean_var(sample_size,sum_x_vega,sum_xx_vega);
        return greek_dif;
    }
    vector<mean_var> MC_decrease_anti(unsigned sample_size, string type_payoff, double init_eps, int h){
        double sum_x_delta=0, sum_xx_delta=0,sum_x_gamma=0,sum_xx_gamma=0,sum_x_vega=0,sum_xx_vega=0; 
        vector<double> payoff_vec;
        vector<double> payoff_vec_anti;

        vector<mean_var> greek_dif(3);

        double epsilon;
        vector<double> brow;
        for(unsigned i=0;i<sample_size;++i){
            brow=normal_generation(h);
            epsilon=pow((init_eps+i+1),(-0.25));
            //delta
            payoff_vec=sim.payoff_diff_fin(type_payoff,epsilon,brow);

            payoff_vec_anti=sim.payoff_diff_fin(type_payoff,epsilon,negative_tab(brow));

            
            sum_x_delta+=0.5*((*this).first_order(payoff_vec[1],payoff_vec[2],epsilon)+(*this).first_order(payoff_vec_anti[1],payoff_vec_anti[2],epsilon));
            sum_xx_delta+=pow(0.5*((*this).first_order(payoff_vec[1],payoff_vec[2],epsilon)+(*this).first_order(payoff_vec_anti[1],payoff_vec_anti[2],epsilon)),2);
            
            //gamma
            sum_x_gamma+=0.5*((*this).second_order(payoff_vec[1],payoff_vec[0],payoff_vec[2],epsilon)+(*this).second_order(payoff_vec_anti[1],payoff_vec_anti[0],payoff_vec_anti[2],epsilon));
            sum_xx_gamma+=pow(0.5*((*this).second_order(payoff_vec[1],payoff_vec[0],payoff_vec[2],epsilon)+(*this).second_order(payoff_vec_anti[1],payoff_vec_anti[0],payoff_vec_anti[2],epsilon)),2);
            //vega
            sum_x_vega+=0.5*((*this).first_order(payoff_vec[3],payoff_vec[4],epsilon)+(*this).first_order(payoff_vec_anti[3],payoff_vec_anti[4],epsilon));
            sum_xx_vega+=pow(0.5*((*this).first_order(payoff_vec[3],payoff_vec[4],epsilon)+(*this).first_order(payoff_vec_anti[3],payoff_vec_anti[4],epsilon)),2);
            }
        greek_dif[0]=mean_var(sample_size,sum_x_delta,sum_xx_delta);
        greek_dif[1]=mean_var(sample_size,sum_x_gamma,sum_xx_gamma);
        greek_dif[2]=mean_var(sample_size,sum_x_vega,sum_xx_vega);
        return greek_dif;
    }

    mean_var MC_decrease_cont_var(double rr, double Tt,unsigned sample_size,string type_payoff,double init_eps,int h){
        double sum_x_delta=0, sum_xx_delta=0;
        vector<double> payoff_vec;
        vector<double> payoff_vec_anti;

        mean_var greek_dif;
        vector<double> brow;
        vector<double> St_tangent;

        // calcul de lambda
        double sum_v_delta=0;
        double sum_c_delta=0;
        double lambda_delta=0;
        int sign=1;
        if (type_payoff=="put"){
            sign=-1;
        } // si put alors c'est -1


        double epsilon;

        for(unsigned i=0;i<sample_size*0.1;i++){
            epsilon=pow((init_eps+i+1),(-0.25));
            brow=normal_generation(h);

            payoff_vec=sim.payoff_diff_fin(type_payoff,epsilon,brow);
            St_tangent=sim.St_tangent(brow);

            sum_v_delta+=pow(sign*(exp(-r*T)*St_tangent[1]-1),2);
            sum_c_delta+=(*this).first_order(payoff_vec[1],payoff_vec[2],epsilon)*sign*(exp(-r*T)*St_tangent[1]-1);

        }
        lambda_delta=sum_c_delta/sum_v_delta;
        //calcul des greek avec variables de controle 
        for(unsigned i=0;i<0.9*sample_size;++i){
            epsilon=pow((init_eps+i+(sample_size*0.1)),(-0.25));
            brow=normal_generation(size_h); 
            //delta
            payoff_vec=sim.payoff_diff_fin(type_payoff,epsilon,brow);
            St_tangent=sim.St_tangent(brow);
            
            sum_x_delta+=(*this).first_order(payoff_vec[1],payoff_vec[2],epsilon)-lambda_delta*sign*(exp(-rr*Tt)*St_tangent[1]-1);
            sum_xx_delta+=pow((*this).first_order(payoff_vec[1],payoff_vec[2],epsilon)-lambda_delta*sign*(exp(-rr*Tt)*St_tangent[1]-1),2);
            }

        greek_dif=mean_var(sample_size*0.9,sum_x_delta,sum_xx_delta);
        return greek_dif;
    }   

    private:
        unsigned test{};
        double eps, size_h,r,T;
        T_sim sim;
        std::normal_distribution<double> X;  
};
