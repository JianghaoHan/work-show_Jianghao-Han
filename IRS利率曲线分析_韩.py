
'''
Author: Jianghao Han (jianghao@bu.edu), 02/12/2023
Description: IRS curve

'''

import numpy as np
import pandas as pd
from copy import deepcopy

from sympy import symbols, solve  # solve unknown variable in equation
import matplotlib.pyplot as plt



# function 1: 
def DiscountFactor_to_ForwardRate(discount_factor, delta_time):
    ''' returns forward rate(float) given discount factor during the same period
    

        input:
        discount_factor: discount factor from time s to time t(s and t are future time points)
        delta_time: time(in years) between payments of fixed leg of IRS, eg: 0.5 for semi-annual payment

    '''
    forward_rate = (1 / delta_time) * (discount_factor**-1 - 1) 
    
    return forward_rate
    

# function 1:     
def ForwardRate_to_DiscountFactor(forward_rate, delta_time):
    ''' returns discount factor(float) given forward rate during the same period
    

        input:
        forward_rate: forward rate from time s to time t(s and t are future time points)
        delta_time: time(in years) between payments of fixed leg of IRS, eg: 0.5 for semi-annual payment

    '''
    discount_factor = (delta_time * forward_rate + 1)**-1 
    
    return discount_factor    


# function 1:    
def DiscountFactor_to_DiscountRate(discount_factor, delta_time):
    ''' returns discount rate(float) given discount factor during the same period
    

        input:
        discount_factor: discount factor from time s to time t(s and t are future time points)
        delta_time: time(in years) between payments of fixed leg of IRS, eg: 0.5 for semi-annual payment

    '''
    discount_rate = - (1 / delta_time) * np.log(discount_factor)
    
    return discount_rate    
    
    
    


# function 1: 
def build_IRS_curve(maturity, swap_rate, delta_time):
    """ returns some information of IRS as following (numpy array): 
            forward_rate: piecewise-constant forward rate just setted at maturity                        
            discount_factor: discount factor matched with forward_rate 
            discount_rate: discount rate for zero coupon bond whose survival time equals to maturity

    
        input:
        maturity(list): maturity of IRS contract
        swap_rate(list): annualized and fixed return rate of fixed leg of IRS
        delta_time(float): time(in years) between payments of fixed leg, eg: 0.5 for semi-annual payment
               
    """
    # compute time interval between maturity; get time_interval list
    time_interval = [maturity[0]]  # first time interval
    
    for i in range(1, len(maturity)):
        maturity_diff = maturity[i] - maturity[i - 1]
        time_interval.append(maturity_diff)        

    # compute forward rate list
    forward_rate = [swap_rate[0]]   # one period IRS
    discount_factor = [ForwardRate_to_DiscountFactor(forward_rate[0], time_interval[0])] # ZCB
    discount_rate = [DiscountFactor_to_DiscountRate(discount_factor[0], time_interval[0])] # zero rate

    # just check    
    # fixed_leg_pv = time_interval[0] * swap_rate[0] * discount_factor[0]
    # floating_leg_pv = time_interval[0] * forward_rate[0] * discount_factor[0]

    # bootstrap to solve new forward rate    
    for i in range(1, len(maturity)):
        # denominator of swap rate
        fixed_leg_pv_deno = 0
        
        floating_leg_pv = 0
        
        
        for j in range(i):
                        
            fixed_leg_pv_deno += time_interval[j] * discount_factor[j]          
            floating_leg_pv += time_interval[j] * forward_rate[j] * discount_factor[j]

        # solve new discount_factor, then transfer it into new forward_rate and new discount_rate
        new_discount_factor = symbols("new_discount_factor")
        fixed_leg_pv_deno += time_interval[i] * new_discount_factor        
        fixed_leg_pv = swap_rate[i] * fixed_leg_pv_deno
        
        # need simplification firstly
        floating_leg_pv += 1 - new_discount_factor
        
        IRS_pv = fixed_leg_pv - floating_leg_pv  # should be zero
        
        # slove new_discount_factor 
        new_disc_factor = float(solve(IRS_pv, new_discount_factor)[0])
        new_forward_rate = DiscountFactor_to_ForwardRate(new_disc_factor, time_interval[i])
        new_discount_rate = DiscountFactor_to_DiscountRate(new_disc_factor, time_interval[i])
        
        # update new information for IRS curve
        forward_rate.append(new_forward_rate)
        discount_factor.append(new_disc_factor)
        discount_rate.append(new_discount_rate)
        
        
    return forward_rate, discount_factor, discount_rate
        
        
        



# Use the __main__ section for test call 
# This section will automatically be executed when the file is run in Python
if __name__ == '__main__':

    # # sample test call
    
    
    # bootstrap the IRS curve to get piecewise-constant forward rate
    maturity = [1,2,3,4,5,7,10,30]
    swap_rate = [0.028438,0.0306,0.03126,0.03144,0.03150,0.03169,0.03210,0.03237]  
    delta_time = 0.5  # fixed leg pay semi-annually


    forward_rate, discount_factor, discount_rate = build_IRS_curve(maturity, swap_rate, delta_time)

    data = np.array([maturity, swap_rate, forward_rate, discount_factor, discount_rate]).transpose()
    
    IRS_info = pd.DataFrame(data, columns = ["maturity", "swap_rate", "forward_rate", "discount_factor", "discount_rate"])


    # (a) extract the constant forward rate for the first year 
    print("(a) constant forward rate for the first year\n")
    print(IRS_info["forward_rate"][0])
    print()
    

    # (b) find the forward rate from one year to two years
    print("(b) find the forward rate from one year to two years\n")
    print(IRS_info["forward_rate"][1])
    print()


    # (c) extract piecewise constant forward rates for the entire curve
    # Comment on the forward rates vs. the swap rates
    print("(c) extract piecewise constant forward rates for the entire curve\n")
    print(IRS_info[["forward_rate"]])
    print()
    
    # visualize
    plt.plot(maturity, forward_rate, label = "forward rate")
    plt.plot(maturity, swap_rate, label = "swap rate")
    plt.title("forward rates v.s. the swap rates")
    plt.legend()
    plt.show()
    print("Forward rates are higher than swap rates all the time, especially for short term rate.\n")
    print("The reason for the positive spread of forward rate is that forward rates ", end = "") 
    print("are expectations of future floating rates, which are likely to rise in the future.\n")
          

    # (d) compute the breakeven swap rate of a 15Y swap
    print("(d) compute the breakeven swap rate of a 15Y swap\n")  
    time_interval = [1,1,1,1,1,2,3,5] # len(time_interval) = 8
    
    # element wise computation
    fixed_leg_pv_deno = np.sum(np.multiply(time_interval, IRS_info["discount_factor"]))
    floating_leg_pv = np.sum(np.array(time_interval) * np.array(IRS_info["forward_rate"]) * np.array(IRS_info["discount_factor"]))

    swap_rate_15Y = floating_leg_pv / fixed_leg_pv_deno

    print(swap_rate_15Y)
    print()


    # (e) compute discount factors and zero rates 
    # Comment on the differences in the zero rates and swap rates
    print("(e) compute discount factors and zero rates\n") 
    print(IRS_info[["maturity", "swap_rate", "forward_rate", "discount_factor", "discount_rate"]])
    print()


    # (f) all forward rates up 100 basis points and re-calculate swap rates
    # Are these rates equivalent to having shifted the swap rates directly?
    print("(f) all forward rates up 100 basis points and re-calculate swap rates\n") 
    IRS_info["forward_rate_new"] = IRS_info["forward_rate"] + 0.01
    time_interval = [1,1,1,1,1,2,3,20] # original interval, len(time_interval) = 8

    # new forward rate leads to new discount factor
    IRS_info["discount_rate_new"] = (1/np.array(time_interval))  * np.log(np.multiply(time_interval, IRS_info["forward_rate_new"]) + 1)
    IRS_info["discount_factor_new"] = np.exp(- IRS_info["discount_rate_new"] * np.array(time_interval))

    # compute new swap rate
    IRS_info["swap_rate_new"] = np.nan
    IRS_info["swap_rate_new"][0] = IRS_info["forward_rate_new"][0] # one period


    for i in range(1, len(time_interval)):
        # element wise computation
        fixed_leg_pv_deno = np.sum(np.multiply(time_interval[:i], IRS_info["discount_factor"][:i]))
        floating_leg_pv = np.sum(np.array(time_interval[:i]) * np.array(IRS_info["forward_rate_new"][:i]) * np.array(IRS_info["discount_factor_new"][:i]))
    
        swap_rate = floating_leg_pv / fixed_leg_pv_deno  
        
        IRS_info["swap_rate_new"][i] = swap_rate

    IRS_info["swap_change"] = IRS_info["swap_rate_new"] - IRS_info["swap_rate"]
        
    print(IRS_info[["maturity", "swap_rate", "swap_rate_new", "swap_change"]])
    print()

    # visualize
    plt.plot(maturity, IRS_info["swap_rate"], label = "swap_rate")
    plt.plot(maturity, IRS_info["swap_rate_new"], label = "swap_rate_new")
    plt.plot(maturity, IRS_info["swap_change"], label = "swap_change")
    plt.title("original swap rates v.s. new swap rates")
    plt.legend()
    plt.show()
  
    print("Shifting forward rates are not equivalent to shifting swap rates directly!")
    print("The change of swap rates is smaller than 100bp.\n")


    # (g) the new bearish swap rates
    print("(g) the new bearish swap rates\n") 
    shift = np.array([0,0,0,5,10,15,25,50]) / 10000
    IRS_info["bear_swap_rate"] = IRS_info["swap_rate"] + shift
    print(IRS_info[["bear_swap_rate"]])
    print()
    
    # (h) re-run bootstrapping with this new bearish swap rates to get forward rates
    # Comment on the changes to the forward rates
    print("(h) re-run bootstrapping with this new bearish swap rates to get forward rates\n")     
        
    # bootstrap the IRS curve to get piecewise-constant forward rate
    forward_rate_bear = build_IRS_curve(maturity, IRS_info["bear_swap_rate"], delta_time)[0]

    data_bear = np.array([maturity, IRS_info["forward_rate"], forward_rate_bear]).transpose()    
    IRS_info_bear = pd.DataFrame(data_bear, columns = ["maturity", "forward_rate", "forward_rate_bear"])
    IRS_info_bear["forward_rate_change"] = IRS_info_bear["forward_rate_bear"] - IRS_info_bear["forward_rate"]
    IRS_info_bear["swap_rate_shift"] = shift  
    
    print(IRS_info_bear)  
    print()
    
    # visualize
    plt.plot(maturity, IRS_info_bear["forward_rate"], label = "forward_rate")
    plt.plot(maturity, IRS_info_bear["forward_rate_bear"], label = "forward_rate_bear")
    plt.plot(maturity, IRS_info_bear["forward_rate_change"], label = "forward_rate_change")
    plt.plot(maturity, IRS_info_bear["swap_rate_shift"], label = "swap_rate_shift")
    plt.title("forward rates change when it comes to bear market")
    plt.legend()
    plt.show()
    
    print("When it comes to bear market, the rise of forward rates is greater than the rise of swap rates.")
    print("And long-term forward rates increase faster than short-term because investors panic about the future.")      
    print("When it comes to a short-term bear market, investors will be more pessimistic about the long-term market.")      
    print()

    
    # (i) the new bullish swap rates
    print("(i) the new bullish swap rates\n") 
    shift = np.array([-50,-25,-15,-10,-5,0,0,0]) / 10000
    IRS_info["bull_swap_rate"] = IRS_info["swap_rate"] + shift
    print(IRS_info[["bull_swap_rate"]])
    print()
    
    # (j) re-run bootstrapping with this new bullish swap rates to get forward rates
    # Comment on the changes to the forward rates
    print("(j) re-run bootstrapping with this new bullish swap rates to get forward rates\n")     
        
    # bootstrap the IRS curve to get piecewise-constant forward rate
    forward_rate_bull = build_IRS_curve(maturity, IRS_info["bull_swap_rate"], delta_time)[0]

    data_bull = np.array([maturity, IRS_info["forward_rate"], forward_rate_bull]).transpose()    
    IRS_info_bull = pd.DataFrame(data_bull, columns = ["maturity", "forward_rate", "forward_rate_bull"])
    IRS_info_bull["forward_rate_change"] = IRS_info_bull["forward_rate_bull"] - IRS_info_bull["forward_rate"]
    IRS_info_bull["swap_rate_shift"] = shift  
    
    print(IRS_info_bull)  
    print()
    
    # visualize
    plt.plot(maturity, IRS_info_bull["forward_rate"], label = "forward_rate")
    plt.plot(maturity, IRS_info_bull["forward_rate_bull"], label = "forward_rate_bull")
    plt.plot(maturity, IRS_info_bull["forward_rate_change"], label = "forward_rate_change")
    plt.plot(maturity, IRS_info_bull["swap_rate_shift"], label = "swap_rate_shift")
    plt.title("forward rates change when it comes to bull market")
    plt.legend()
    plt.show()
    
    print("When it comes to bull market, the decline of forward rates is smaller than the decline of swap rates.")
    print("And short-term forward rates decrease faster than long-term because investors optimistic about the near future.")      
    print("However investors are not optimistic about the long-term market just because of the immediate bull market.")      

    
    
    
    
    
    
    


