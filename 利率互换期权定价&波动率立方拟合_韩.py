
'''
Author: Jianghao Han (jianghao@bu.edu), 04/05/2023
Description: Swaption Pricing and Risk Management under the SABR Model

'''


import numpy as np
import pandas as pd
from scipy.stats import norm
# import math
# import matplotlib.pyplot as plt
from scipy.optimize import minimize, root




# function 1:    
def instantaneous_forward_rate(F0, delta_time):
    ''' returns (constant) instantaneous forward rate (float) of a swap
        
        input:
        F0 (float): the current forward rate, which will start at some time point in the future, i.e. forward starting 
        delta_time(float): time(in years) between payments of fixed leg, eg: 0.5 for semi-annual payment
      
    '''
    ins_forward_rate = (1 / delta_time) * np.log(1 + delta_time * F0) 

    return ins_forward_rate



# function 2:    
def annuity(maturity, tenor, delta_time, instantaneous_forward_rate):
    ''' returns current annuity value(float) of a swap
        
        input:
        maturity (float): maturity of option within a swaption contract
        tenor (float): expiry(time to maturity) of IRS(swap)        
        delta_time(float): time(in years) between payments of fixed leg, eg: 0.5 for semi-annual payment
        instantaneous_forward_rate (float): (constant) instantaneous forward rate (float) of a swap
        
    '''
    discount_rate = instantaneous_forward_rate

    annuity = 0    
    
    for i in range(1, round(tenor / delta_time) + 1):
        # the total length of time to be discounted in years
        time_interval = maturity +  delta_time * i
        discount_factor = np.exp(- discount_rate * time_interval)
    
        annuity += delta_time * discount_factor  
    

    return annuity



# function 3:    
def swaption_pricing(A, F0, sigma, K, tenor, model = "BS"):
    ''' returns present rate(premium) (float) of a payer swaption when entering into an IRS
        
        input:
        A (float): the annuity value(present value) of a stream of payments with a notional value of 1 basis point
        F0 (float): the current forward rate, which will start at some time point in the future, i.e. forward starting 
        sigma (float): the annualized standard deviation of underlying forward rate 
        K (float): exercise/strike rate of swaption
        tenor (float): expiry(time to maturity) for a dynamic of underlying forward rate 
        model (str): the model used to price option, could be "BS" or "Bachelier"
        
    '''
    # Black-Scholes model payoff
    if model == "BS":
        d1 = (np.log(F0 / K) + 0.5 * sigma**2 * tenor) / (sigma * np.sqrt(tenor))
        d2 = (np.log(F0 / K) - 0.5 * sigma**2 * tenor) / (sigma * np.sqrt(tenor))
        
        # expected payoff 
        payoff = (F0 * norm.cdf(d1) - K * norm.cdf(d2))


    # Bachelier model payoff
    if model == "Bachelier":
        d1 = (F0 - K) / (sigma * np.sqrt(tenor))
        
        # expected payoff 
        payoff = sigma * np.sqrt(tenor) * (d1 * norm.cdf(d1) + norm.pdf(d1))
 
        
    # swaption pricing: discount total payoff to get present value    
    premium =  A * payoff
        
    return premium  # analytic solution  



# function 4:    
def implied_vol_Bachelier(T, K, F0, sigma0, alpha, beta, rho):
    ''' returns Bachelier(Normal) volatility of an option under the SABR model
        
        input:
        T (float): maturity of option within a swaption contract
        K (float): exercise/strike rate of swaption        
        F0 (float): the current forward rate, which will start at some time point in the future, i.e. forward starting 
        sigma0 (float): initial volatility of underlying forward rate
        alpha (float): SABR model's alpha parameter, which represents the initial value of the implied volatility
        beta (float): SABR model's beta parameter, which determines the shape of the volatility smile
        rho (float): SABR model's rho parameter, which represents the correlation between the underlying forward rate and its volatility
       
    '''
    epsilon = T * alpha ** 2
    zeta = (alpha / (sigma0 * (1 - beta))) * (F0**(1 - beta) - K**(1 - beta))
    delta = np.log((np.sqrt(1 - 2 * rho * zeta + zeta**2) + zeta - rho) / (1 - rho))
    
    # the midpoint of the current forward rate and the strike rate
    F_mid = (F0 + K) / 2
    gamma1 = beta / F_mid
    gamma2 = (beta * (beta - 1)) / F_mid**2    
    
    # asymptotic approximation formula
    implied_vol = alpha * (F0 - K) / delta * \
                  (1 + epsilon * ( (2 * gamma2 - gamma1 ** 2) / 24 * (sigma0 * F_mid**beta / alpha)**2 + \
                                    rho * gamma1 / 4 * sigma0 * F_mid**beta / alpha + \
                                    (2 - 3 * rho**2) / 24 ))

    return implied_vol
    

# function 5:    
def RSS_vol(parameter, T, K_list, F0, market_vol_list):
    ''' returns Residual Sum of Squares(a deviation metric) when we estimate market volatility using SABR model
        
        input:
        parameter (vector of float): (alpha, beta, rho, sigma0)     
        T (float): maturity of option within a swaption contract
        K_list (float): a list/np-array/dataframe of exercise/strike rate of swaption, length = 6 in this problem        
        F0 (float): the current forward rate, which will start at some time point in the future, i.e. forward starting 
        market_vol_list: a list of normal swaption volatilities observed in market, length = 6 in this problem 
        
        alpha (float): SABR model's alpha parameter, which represents the initial value of the implied volatility
        beta (float): SABR model's beta parameter, which determines the shape of the volatility smile
        rho (float): SABR model's rho parameter, which represents the correlation between the underlying forward rate and its volatility
        sigma0 (float): initial volatility of underlying forward rate        
        
    '''
    alpha, beta, rho, sigma0 = parameter  # Sequence Unpacking; prepared for parameter calibration
    RSS = 0
    
    for i in range(len(market_vol_list)):
        # estimate implied Bachelier(Normal) volatility of an option under the SABR model 
        vol_estimated = implied_vol_Bachelier(T, K_list[i], F0, sigma0, alpha, beta, rho)
        RSS += (vol_estimated - market_vol_list[i])**2
        
    return RSS


# function 6:    
def delta_BS(A, F0, sigma, K, tenor):
    ''' returns greek delta (float) of call option under BS model
        
        input:
        A (float): the annuity value(present value) of a stream of payments with a notional value of 1 basis point
        F0 (float): the current forward rate, which will start at some time point in the future, i.e. forward starting 
        sigma (float): the annualized standard deviation of underlying forward rate 
        K (float): exercise/strike rate of swaption
        tenor (float): expiry(time to maturity) for a dynamic of underlying forward rate
        
    '''
    d1 = (np.log(F0 / K) + 0.5 * sigma**2 * tenor) / (sigma * np.sqrt(tenor))
    # compute delta
    delta =  A * norm.cdf(d1)
        
    return delta  # analytic solution  






# Use the __main__ section for test call 
# This section will automatically be executed when the file is run in Python
if __name__ == '__main__':

    # # sample test call


    # parameter setting
    tenor = 5  # expiry for IRS in swaption, in years
    delta_time = 0.5  # fixed leg pay semi-annually
    F0 = [117.45, 120.60, 133.03, 152.05, 171.85]  # par swap rate listed, reported in bps
    bp = 0.0001  # 1bp
    F0 = [x * bp for x in F0]  # real par swap rate
    maturity = [1, 2, 3, 4, 5]  # maturity/expiry for options in swaption
    
    F0_table = np.tile(F0, (6, 1)).T  # stack list vertically and then transpose ; shape(5, 6)
    strike_adjust = np.array([-50, -25, -5, 5, 25, 50]) * bp  # real strike adjust
    K = F0_table + strike_adjust  # broadcasting in numpy    
    
    market_vol = np.array([[58.31, 51.51, 49.28, 48.74, 41.46, 37.33],
                           [51.72, 46.87, 43.09, 42.63, 38.23, 34.55],
                           [46.29, 44.48, 43.61, 39.36, 35.95, 32.55],
                           [45.72, 41.80, 38.92, 38.19, 34.41, 31.15],
                           [44.92, 40.61, 37.69, 36.94, 33.36, 30.21]])    # volatility listed, reported in bps
    
    market_vol = market_vol * bp  # reak market volatility
    
    
    

    # (a) Calculate the constant instantaneous forward rate for each swap
    print("(a) Calculate the constant instantaneous forward rate for each swap\n")
    forward_rate = [instantaneous_forward_rate(F, delta_time) for F in F0]
    print(forward_rate)
    print()
    

    # (b) Calculate the current annuity value for each swap
    print("(b) Calculate the current annuity value for each swap\n")
    A = [annuity(maturity[i], tenor, delta_time, forward_rate[i]) for i in range(5)]
    print(A)
    print()


    # (c) Calculate a table of premiums for each swaption using the Bachelier model
    print("(c) Calculate a table of premiums for each swaption using the Bachelier model\n")
    premium_table = np.empty((5, 6))
    premium_table[:] = np.nan   # filled with na initially

    # option(swaption) pricing to get premium using market volatility
    for i in range(market_vol.shape[0]):
        for j in range(market_vol.shape[1]):  # traverse market vol
        
            premium_table[i][j] = swaption_pricing(A[i], F0[i], market_vol[i][j] , K[i][j], tenor, model = "Bachelier")
    
    premium_df = pd.DataFrame(premium_table, 
                              index = ['1Y', '2Y', '3Y', '4Y', '5Y'],
                              columns = ['ATM-50', 'ATM-25', 'ATM-5', 'ATM+5', 'ATM+25', 'ATM+50'])  
   
    print(premium_df)
    print()    
 

    # (d) Calculate implied Bachelier volatility using Asymptotic Approximation under SABR model
    # Calibrate SABR parameters to minimize the distance between market and model/implied volatility
    print("(d) Calculate implied Bachelier volatility using Asymptotic Approximation under SABR model\n")
    print("Calibrate SABR parameters to minimize the distance between market and model/implied volatility\n")
    print()
    
    # beta = 0.5  # a general case just for simplification
    # initial guess for alpha, beta, rho, and sigma0
    initial_guess = [0.05, 0.5, 0, 0.05]
    bounds = ((0, 1), (0, 1), (-1, 1), (0, 1))
    
    SABR_parameter = pd.DataFrame(np.nan,                             
                                  index = [0, 1, 2, 3, 4],
                                  columns = ['alpha', 'beta', 'rho', 'sigma0'])

    # traverse all the market vol to calibrate parameters: alpha, beta, rho, sigma0
    for i in range(len(maturity)):  # length = 5
        optimization = minimize(RSS_vol,  # objective function
                                args = (maturity[i], K[i], F0[i], market_vol[i]),  # parameters other than unknown parameters
                                x0 = initial_guess, 
                                method = 'SLSQP', 
                                bounds = bounds)
        
        SABR_parameter.loc[i] = optimization.x
        
        # the minimal value of obejective function 
        print(f"After parameter calibration, the RSS of implied and market valotility for the {i+1}Y swaption is {optimization.fun}")
        print()
    
    # set new index for df for showing result
    new_index = ['1Y', '2Y', '3Y', '4Y', '5Y']
    SABR_parameter['maturity'] = new_index
    SABR_parameter.set_index('maturity', inplace=True)

    print("SABR parameters after calibration:\n")    
    print(SABR_parameter)
    print()  

   
    # (e) Comment on the relationship of the calibrated parameters as a function of expiry(maturity of option)
    print("(e) Comment on the relationship of the calibrated parameters as a function of expiry\n")
    print("When expiry increases, parameter alpha, rho and sigma0 decrease slowly while parameter beta increases slowly.")
    print()        
        

    # (f) Using calibrated SABR parameters to calculate the price and normal volatility 
    # of swaptions with strikes equal to ATM - 75 and ATM + 75
    print("(f) Calculate the price and normal volatility of swaptions when strikes equal to ATM - 75 and ATM + 75\n")

    # slice 2 columns from F0_table
    F0_f = F0_table[:, :2]  # just need 2 F0 array to adjust strike rate: ATM - 75 and ATM + 75
    strike_adjust_f = np.array([-75, 75]) * bp  # real strike adjust
    K_f = F0_f + strike_adjust_f  # broadcasting in numpy 

    ''' Idea:
    Now we do not have market vlatility when strike rates equal to ATM - 75 and ATM + 75.
    So we cannot price swaptions without reliable volatility data.
    We firstly need to calculate implied Bachelier(normal) volatility using Asymptotic Approximation 
    under SABR model based on the parameter sigma0 calibrated in (d), then we secondly use the 
    implied Bachelier(normal) volatility to price swaptions.
    In this situation, we set implied Bachelier(normal) volatility as market volatility when pricing swaption 
    '''
    
    # firstly compute implied Bachelier(normal) volatility  
    vol_f = np.empty((5, 2))
    vol_f[:] = np.nan   # filled with na initially    
    
    for i in range(5):  # row number of F0_f
        for j in range(2):  # column number of F0_f
        
            # use calibrated SABR parameters in (d)
            alpha = SABR_parameter["alpha"][i]
            beta = SABR_parameter["beta"][i]
            rho = SABR_parameter["rho"][i]
            sigma0 = SABR_parameter["sigma0"][i]
            
            # estimate implied vol
            vol_f[i][j]  = implied_vol_Bachelier(maturity[i], K_f[i][j], F0[i], sigma0, alpha, beta, rho)
            
    
    vol_f_df = pd.DataFrame(vol_f, 
                            index = ['1Y', '2Y', '3Y', '4Y', '5Y'], 
                            columns = ['ATM-75', 'ATM+75'])

    print("Implied Bachelier(normal) volatility:\n")     
    print(vol_f_df)    
    print()   

    
    # secondly compute the price(premium) of swaption  
    premium_f = np.empty((5, 2))
    premium_f[:] = np.nan   # filled with na initially 
    
    for i in range(5):  # row number of F0_f
        for j in range(2):  # column number of F0_f
        
            premium_f[i][j] = swaption_pricing(A[i], F0[i], vol_f[i][j] , K_f[i][j], tenor, model = "Bachelier")        
 
    
    premium_f_df = pd.DataFrame(premium_f, 
                            index = ['1Y', '2Y', '3Y', '4Y', '5Y'], 
                            columns = ['ATM-75', 'ATM+75'])

    print("Price(premium) of swaption:\n")     
    print(premium_f_df)    
    print()   


    # (g) Calculate the equivalent Black volatilities for each option
    print("(g) Calculate the equivalent Black volatilities for each option\n")
    vol_BS = np.empty((5, 6))
    vol_BS[:] = np.nan   # filled with na initially

    for i in range(5):
        for j in range(6):
            
            # set unknown vol_B = x
            def f(x):
                premium_BS = swaption_pricing(A[i], F0[i], x, K[i][j], tenor, model = "BS")
            
                equation = premium_BS - premium_table[i][j]  # premium_table from Bachelier
                
                return equation
            
            # solve unknown vol
            vol_solution = root(f, x0 = 0.1)
            # store
            vol_BS[i][j] = vol_solution.x


    vol_BS_df = pd.DataFrame(vol_BS, 
                             index = ['1Y', '2Y', '3Y', '4Y', '5Y'],
                             columns = ['ATM-50', 'ATM-25', 'ATM-5', 'ATM+5', 'ATM+25', 'ATM+50'])

    print("Implied Black(log normal) volatility:\n")  
    print(vol_BS_df)
    print()


    # (h) Calculate the delta of each options under Black’s model
    print("(h) Calculate the delta of each options under Black’s model\n")
    delta_Black = np.empty((5, 6))
    delta_Black[:] = np.nan   # filled with na initially

    for i in range(5):
        for j in range(6):
            
            delta_Black[i][j] = delta_BS(A[i], F0[i], vol_BS[i][j], K[i][j], tenor)


    delta_Black_df = pd.DataFrame(delta_Black, 
                             index = ['1Y', '2Y', '3Y', '4Y', '5Y'],
                             columns = ['ATM-50', 'ATM-25', 'ATM-5', 'ATM+5', 'ATM+25', 'ATM+50'])

    print("Delta of each options under BS model:\n")  
    print(delta_Black_df)
    print()


    # (i) Estimate a SABR smile adjusted delta for each option
    print("(i) Estimate a SABR smile adjusted delta for each option\n")
    delta_adjust = np.empty((5, 6))
    delta_adjust[:] = np.nan   # filled with na initially

    for i in range(5):
        for j in range(6):
            
            # use calibrated SABR parameters in (d)
            alpha = SABR_parameter["alpha"][i]
            beta = SABR_parameter["beta"][i]
            rho = SABR_parameter["rho"][i]
            sigma0 = SABR_parameter["sigma0"][i]
            
            # shift forward rate F0
            F0_up = F0[i] + 1 * bp 
            F0_down = F0[i] - 1 * bp 

            # lead to the shift of vol
            vol_up  = implied_vol_Bachelier(maturity[i], K[i][j], F0_up, sigma0, alpha, beta, rho)
            vol_down  = implied_vol_Bachelier(maturity[i], K[i][j], F0_down, sigma0, alpha, beta, rho)

            # then lead to the shift of premium
            premium_up = swaption_pricing(A[i], F0_up, vol_up, K[i][j], tenor, model = "Bachelier")
            premium_down = swaption_pricing(A[i], F0_down, vol_down, K[i][j], tenor, model = "Bachelier")

            # then lead to the adjusted delta
            delta_adjust[i][j] = (premium_up - premium_down) / (F0_up - F0_down)


    delta_adjust_df = pd.DataFrame(delta_adjust, 
                             index = ['1Y', '2Y', '3Y', '4Y', '5Y'],
                             columns = ['ATM-50', 'ATM-25', 'ATM-5', 'ATM+5', 'ATM+25', 'ATM+50'])

    print("SABR smile adjusted delta for each option:\n")  
    print(delta_adjust_df)
    print()

    print("Comparing SABR smile adjusted delta with BS delta, we can find that", end = " ")
    print("adjusted delta is smaller than BS delta in most cases.")
    print()
    print("The reason might be that SABR model captures the volatility smile while BS does not.", end = " ")
    print("Then the implied volatility smile absorbs(reduces) part of the sensitivity of the premium with respect to forward rate.")












