# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats as sps

def Put_BS_core(tau, K, DF, F, sigma):
    """
    Prix du put Black-Scholes en fonction de
    
    tau   : time to maturity
    K     : strike
    DF    : discount factor. DF = e^{-r*tau} si taux d'interet deterministe r.
    F     : prix forward du sous-jacent
    sigma : volatilite
    """
    sigma_sqrt_tau = sigma * np.sqrt(tau)
    
    d_1 = np.log(F/K)/sigma_sqrt_tau + sigma_sqrt_tau/2.
    
    d_2 = d_1 - sigma_sqrt_tau

    prix_put = DF * (K * sps.norm.cdf(-d_2) - F * sps.norm.cdf(-d_1))

    return prix_put


def Call_BS_core(tau, K, DF, F, sigma):
    """
    Prix du call Black-Scholes en fonction de
    
    tau   : time to maturity
    K     : strike
    DF    : discount factor. DF = e^{-r*tau} si taux d'interet deterministe r.
    F     : prix forward du sous-jacent
    sigma : volatilite
    """
    sigma_sqrt_tau = sigma * np.sqrt(tau)

    d_1 = np.log(F/K)/sigma_sqrt_tau + sigma_sqrt_tau/2.
    
    d_2 = d_1 - sigma_sqrt_tau

    prix_call = DF * (F * sps.norm.cdf(d_1) - K * sps.norm.cdf(d_2))

    return prix_call


def Put_BS(t, S, T, K, r, q, sigma):
    """
    Prix du put Black-Scholes en fonction des parametres BS usuels.
    """
    tau = T - t
    DF = np.exp(-r*tau)
    F = S * np.exp((r-q)*tau)
    
    return Put_BS_core(tau, K, DF, F, sigma)


def Vega_core(tau, K, DF, F, sigma):
    """
    Vega d'un put ou call Black-Scholes en fonction de
    
    tau   : time to maturity
    K     : strike
    DF    : discount factor. DF = e^{-r*tau} si taux d'interet deterministe r.
    F     : prix forward du sous-jacent
    sigma : volatilité
    """
    sigma_sqrt_tau = sigma * np.sqrt(tau)
    
    d_1 = np.log(F/K) / sigma_sqrt_tau + sigma_sqrt_tau / 2.

    vega = DF * F * np.sqrt(tau) * np.exp(-d_1**2 / 2) / np.sqrt(2*np.pi)
    
    return vega


def Vega(t, S, T, K, r, q, sigma):
    """
    Vega d'un put ou call Black-Scholes en fonction des parametres BS usuels.
    """
    tau = T - t
    DF = np.exp(-r*tau)
    F = S * np.exp((r-q)*tau)
    
    return Vega_core(tau, K, DF, F, sigma)


def volImplCore_Newton(tau, K, DF, F, price,
                          CallOrPutFlag = 1,
                          initial_point='automatic', prix_tol = 1.e-4, max_iter=50):
    """
    Volatilite implicite d'un call ou put de prix = price, lorsque les autres parametres sont:
    
    tau   : time to maturity
    K     : strike
    DF    : discount factor. DF = e^{-r*tau} si taux d'intéret détérministe r.
    F     : prix forward du sous-jacent
    
    CallOrPutFlag: = 1 if call, 0 if put
    
    Methode: Newton.
    """
    if initial_point == 'automatic':
        vol = np.sqrt( 2/tau * np.abs(np.log(F/K)) )
    else:
        vol = initial_point
    
    if CallOrPutFlag:
        current_price = Call_BS_core(tau, K, DF, F, vol)
    else:
        current_price = Put_BS_core(tau, K, DF, F, vol)
    
    critere_arret = np.abs(current_price - price)
    iterations = 0
    
    while ( (critere_arret > prix_tol) & (iterations < max_iter) ):
        iterations = iterations + 1
        
        vol = vol - (current_price - price) / Vega_core(tau, K, DF, F, vol)
        
        if CallOrPutFlag:
            current_price = Call_BS_core(tau, K, DF, F, vol)
        else:
            current_price = Put_BS_core(tau, K, DF, F, vol)
        
        critere_arret = np.abs(current_price - price)
    
    return vol, iterations


def volImplPutCore_dichotomie(tau, K, DF, F, price, prix_tol = 1.e-3, max_iter=50, a = 0.001, b = 2.):
    """
    Volatilite implicite d'un put de prix = price, lorsque les autres parametres sont:
    
    tau   : time to maturity
    K     : strike
    DF    : discount factor. DF = e^{-r*tau} si taux d'interet deterministe r
    F     : prix forward du sous-jacent
    
    Méthode de calcul: dichotomie sur l'intervalle [a, b].
    """
    prix_min = Put_BS_core(tau, K, DF, F, a)
    prix_max = Put_BS_core(tau, K, DF, F, b)
    
    ######################################################
    ## On verifie si le prix donné est entre les bornes
    ## choisies pour appliquer la méthode de bissection
    ######################################################
    check = (prix_min < price < prix_max)
    
    if check == False:
        print("""ATTENTION: Prix en dehors des bornes de bissection.
           Calcul de vol impossible pour les bornes choisies.
           Cette fonction renvoie vol = 0.""")
        return 0
    
    else:
        vol_min = a
        vol_max = b
        
        vol = (vol_min + vol_max) / 2
        prix_milieu = Put_BS_core(tau, K, DF, F, vol)
        
        critere_arret = np.abs(prix_milieu - price)
        iterations = 0
        
        while ( (critere_arret > prix_tol) & (iterations < max_iter) ):
            iterations = iterations + 1
            
            if prix_milieu - price > 0:
                vol_max = vol
                vol = (vol_min + vol_max) / 2
            else:
                vol_min = vol
                vol = (vol_min + vol_max) / 2
            
            prix_milieu = Put_BS_core(tau, K, DF, F, vol)
            
            critere_arret = np.abs(prix_milieu - price)
        
        return vol, iterations