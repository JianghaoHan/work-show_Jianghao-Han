
'''
Author: Jianghao Han (jianghao@bu.edu), 02/06/2023
Description: build CDS survival curve

'''

import numpy as np
import pandas as pd

from sympy import symbols, solve  # solve unknown variable in equation



# function 1: 
def survival_prob(maturity):
    """ returns the survival probability of a CDS at maturity(in years) 
    
        maturity: maturity(in years) of CDS which cannot be 4 (int)

    """
    prob = symbols("prob")
    i = maturity # just for convenience

    # compute RPV01 and contingent leg 
    if i == 1:        
        RPV01 = 0.5 * (time_list[i] - time_list[i-1]) * np.exp(-rf * i) * (prob_list[i-1] + prob)
        contingent_leg_pv = (1-R) * np.exp(-rf * i) * (prob_list[i-1] - prob)
        
    if i == 2:        
        RPV01 = (0.5 * (time_list[1] - time_list[0]) * np.exp(-rf * 1) * (prob_list[0] + prob_list[1]) +
                 0.5 * (time_list[i] - time_list[i-1]) * np.exp(-rf * i) * (prob_list[i-1] + prob))

        contingent_leg_pv = ((1-R) * np.exp(-rf * 1) * (prob_list[0] - prob_list[1]) + 
                             (1-R) * np.exp(-rf * i) * (prob_list[i-1] - prob))

    if i == 3:        
        RPV01 = (0.5 * (time_list[1] - time_list[0]) * np.exp(-rf * 1) * (prob_list[0] + prob_list[1]) +
                 0.5 * (time_list[2] - time_list[1]) * np.exp(-rf * 2) * (prob_list[1] + prob_list[2]) +
                 0.5 * (time_list[i] - time_list[i-1]) * np.exp(-rf * i) * (prob_list[i-1] + prob))
    
        contingent_leg_pv = ((1-R) * np.exp(-rf * 1) * (prob_list[0] - prob_list[1]) +
                             (1-R) * np.exp(-rf * 2) * (prob_list[1] - prob_list[2]) +
                             (1-R) * np.exp(-rf * i) * (prob_list[i-1] - prob))

    if i == 5:        
        RPV01 = (0.5 * (time_list[1] - time_list[0]) * np.exp(-rf * 1) * (prob_list[0] + prob_list[1]) +
                 0.5 * (time_list[2] - time_list[1]) * np.exp(-rf * 2) * (prob_list[1] + prob_list[2]) +
                 0.5 * (time_list[3] - time_list[2]) * np.exp(-rf * 3) * (prob_list[2] + prob_list[3]) +
                 0.5 * (time_list[i] - time_list[i-1]) * np.exp(-rf * i) * (prob_list[i-2] + prob))
    
        contingent_leg_pv = ((1-R) * np.exp(-rf * 1) * (prob_list[0] - prob_list[1]) +
                             (1-R) * np.exp(-rf * 2) * (prob_list[1] - prob_list[2]) +
                             (1-R) * np.exp(-rf * 3) * (prob_list[2] - prob_list[3]) +
                             (1-R) * np.exp(-rf * i) * (prob_list[i-2] - prob))        
        
    # compute premium leg    
    premium_leg_pv = premium_list[i] * RPV01        
        
    # slove prob
    prob_slove = float(solve(premium_leg_pv - contingent_leg_pv, prob)[0])
    
    return prob_slove
    

# function 2: 
def CDS_present_value(maturity, rf, spread, recovery_rate):
    """ returns the present value a CDS whose maturity is 4  
    
        maturity: maturity(in years) of CDS which can only be 4 (int)
        rf: risk free interest rate 
        spread: contractual spread of CDS at inception
        recovery_rate: expected recovery rate after credit event

    """
    
    RPV01 = 0.5 * ((1 - 0) * np.exp(-rf * 1) * (prob_list[0] + prob_list[1]) +
                     (2 - 1) * np.exp(-rf * 2) * (prob_list[1] + prob_list[2]) +
                     (3 - 2) * np.exp(-rf * 3) * (prob_list[2] + prob_list[3]) +
                     (4 - 3) * np.exp(-rf * 4) * (prob_list[3] + prob_4))

    premium_leg_pv = spread * RPV01
    
    contingent_leg_pv = (1 - recovery_rate) * (np.exp(-rf * 1) * (prob_list[0] - prob_list[1]) +
                               np.exp(-rf * 2) * (prob_list[1] - prob_list[2]) +
                               np.exp(-rf * 3) * (prob_list[2] - prob_list[3]) +
                               np.exp(-rf * 4) * (prob_list[3] - prob_4))     

    CDS_pv = premium_leg_pv - contingent_leg_pv  
    
    return CDS_pv




# Use the __main__ section for test call 
# This section will automatically be executed when the file is run in Python
if __name__ == '__main__':

    # # sample test call
    
    
    # (a) bootstrap the CDS survival curve 
    ''' Idea:
    At time 0, the probability of CDS survivaing to time 0 is prob_0 = 1.
    We know the present value of premium leg equals to the present value of 
    contingent leg. Based on this formula and the information about 1-year CDS,
    we can compute the probability of CDS survivaing to time 1, say prob_1.
    
    Continuing to use this Bootstrap method, we can get prob_2, prob_3 and prob_5.
    With survival probability, we can compute piecewise-constant hazard rate
    h_0_1, h_1_2, h_2_3 and h_3_5.
    Note:h_0_1 represents constant hazard rate between time 0 and time 1.
    
    With all the piecewise-constant hazard rate, we can build CDS survival curve.
    '''
    
    print("(a) bootstrap the CDS survival curve\n")
    rf = 0.03  # risk free interest rate
    R = 0.4 # revovery rate after credit event
    bp = 0.0001 # 1 base point
    prob_0 = 1
    
    time_list = [0, 1, 2, 3, 3, 5]  # different maturity of sample CDS
    premium_list = [0, 100*bp, 110*bp, 120*bp, 0, 140*bp] # sample premium setted at inception
    prob_list = [prob_0]  # prob_0 = 1, waiting to add new probabilities

    # store information about CDS survival curve
    CDS_info = pd.DataFrame()

    # compute survival probablity of CDS at different maturity
    for maturity in [1, 2, 3, 5]:  # # survive to 1,2,3 and 5 year
        prob_survival = survival_prob(maturity)
        prob_list.append(prob_survival)

    # compute piecewise-constant hazard rate  
    h_0_1 = (np.log(prob_list[0]) - np.log(prob_list[1])) / (1 - 0)
    h_1_2 = (np.log(prob_list[1]) - np.log(prob_list[2])) / (2 - 1)
    h_2_3 = (np.log(prob_list[2]) - np.log(prob_list[3])) / (3 - 2)
    h_3_5 = (np.log(prob_list[3]) - np.log(prob_list[4])) / (5 - 3)
        
    print("survival probability of CDS at year 1, 2, 3 and 5: " + str(prob_list[1:])+"\n")
    print(f"hazard rate between year 0 to year 1: {h_0_1}")
    print(f"hazard rate between year 1 to year 2: {h_1_2}")
    print(f"hazard rate between year 2 to year 3: {h_2_3}")
    print(f"hazard rate between year 3 to year 5: {h_3_5}\n")
    
    
    # (b) get fair spread for a 4y CDS that starts at time 0(today)
    print("(b) get fair spread for a 4y CDS that starts at time 0(today)\n")
    prob_4 = np.exp(np.log(prob_list[3]) - h_3_5 * (4 - 3)) # survival probablity at year 4
    
    premium4 = symbols("premium4")
    equation = premium4 * (0.5 * (time_list[1] - time_list[0]) * np.exp(-rf * 1) * (prob_list[0] + prob_list[1]) +
                           0.5 * (time_list[2] - time_list[1]) * np.exp(-rf * 2) * (prob_list[1] + prob_list[2]) +
                           0.5 * (time_list[3] - time_list[2]) * np.exp(-rf * 3) * (prob_list[2] + prob_list[3]) +
                           0.5 * (4 - 3) * np.exp(-rf * 4) * (prob_list[3] + prob_4)) -\
                            ((1-R) * np.exp(-rf * 1) * (prob_list[0] - prob_list[1]) +
                             (1-R) * np.exp(-rf * 2) * (prob_list[1] - prob_list[2]) +
                             (1-R) * np.exp(-rf * 3) * (prob_list[2] - prob_list[3]) +
                             (1-R) * np.exp(-rf * 4) * (prob_list[3] - prob_4)) 
    
    # slove premium_4
    premium_4 = float(solve(equation, premium4)[0])        
    print(f"fair spread for a 4y CDS (in bps): {premium_4 / bp}\n")


    # (c) compute the mark-to-market value of 5y CDS bought 1 year ago with spread 80bps
    print("(c) compute the mark-to-market value of 5y CDS bought 1 year ago with spread 80bps\n")          
    spread_0 = 80 * bp
    
    # for a 5y CDS, spread hypothetically issued at year 1 should equal to premium_4 for 4y CDS
    spread_1 = premium_4  # we can get this constractual from (b)
    
    # at year 1
    RPV01_1 = 0.5 * ((1 - 0) * np.exp(-rf * 1) * (prob_list[0] + prob_list[1]) +
                     (2 - 1) * np.exp(-rf * 2) * (prob_list[1] + prob_list[2]) +
                     (3 - 2) * np.exp(-rf * 3) * (prob_list[2] + prob_list[3]) +
                     (4 - 3) * np.exp(-rf * 4) * (prob_list[3] + prob_4))
    
    # mark-to-market value
    MTM_value = (spread_1 - spread_0) * RPV01_1
    print(f"pseudo contractual spread issued at year 1 for a 5y CDS (in bps): {spread_1 / bp}\n")
    print(f"Mark-to-Market value at year 1 for a 5y CDS (in bps): {MTM_value / bp}\n")


    # (d) compute the dv01 of 4y CDS in (b) with respect to spread from on SDC Curve in (a)
    print("(d) compute the dv01 of 4y CDS with respect to spread\n")

    # get spread of 4y CDS and set delta         
    spread = premium_4
    delta = 1 * bp  # small change of spread
    
    # compute the present value of 4y CDS
    premium_leg_pv = spread * RPV01_1  # RPV01_1 of 5y CDS = RPV01_0 of 4y CDS
  
    contingent_leg_pv = (1 - R) * (np.exp(-rf * 1) * (prob_list[0] - prob_list[1]) +
                                   np.exp(-rf * 2) * (prob_list[1] - prob_list[2]) +
                                   np.exp(-rf * 3) * (prob_list[2] - prob_list[3]) +
                                   np.exp(-rf * 4) * (prob_list[3] - prob_4))   
    
    
    CDS_pv = premium_leg_pv - contingent_leg_pv  # should be 0
    
    CDS_pv_plus_1bp = (spread + delta) * RPV01_1 - contingent_leg_pv
    CDS_pv_minus_1bp = (spread - delta) * RPV01_1 - contingent_leg_pv
    
    DV01_spread = 0.5 * (abs(CDS_pv_plus_1bp - CDS_pv) + abs(CDS_pv_minus_1bp - CDS_pv))
    print(f"DV01 of 4y CDS with respect to spread: {DV01_spread}\n")
    

    # (e) compute the dv01 of 4y CDS in (b) wrt the interest rate curve(risk free rate).
    print("(e) compute the dv01 of 4y CDS wrt risk free rate\n")    
    rf = 0.03  # risk free interest rate            
    delta = 1 * bp  # small change of risk free rate  

    CDS_pv_plus_1bp = CDS_present_value(4, rf + delta, spread, R) 
    CDS_pv_minus_1bp = CDS_present_value(4, rf - delta, spread, R) 

    DV01_rf = 0.5 * (abs(CDS_pv_plus_1bp - CDS_pv) + abs(CDS_pv_minus_1bp - CDS_pv))
    print(f"DV01 of 4y CDS with respect to risk free rate: {DV01_rf}\n")    


    # (f) compute sensitivity of 4y CDS wrt R.
    print("(e) compute sensitivity of 4y CDS wrt revovery rate R\n")
    R = 0.4 # revovery rate after credit event
    delta = 0.01  # small change of revovery rate

    CDS_pv_plus = CDS_present_value(4, rf, spread, R + delta) 
    CDS_pv_minus = CDS_present_value(4, rf, spread, R - delta) 

    sensitivity_R = 0.5 * (abs(CDS_pv_plus - CDS_pv) + abs(CDS_pv_minus - CDS_pv))
    print(f"Sensitivity of 4y CDS wrt R (when R change 1%): {sensitivity_R}\n") 

    # Note: par value of bond associated with CDS
    print("【Note】: we assume par value of bond is 1")
    
    

