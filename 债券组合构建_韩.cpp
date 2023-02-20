
/*

Author: Jianghao Han (jianghao@bu.edu), 10/10/2022
Description: Bond Pricing in C++
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <iomanip>  // control output
#include <cstdlib>


// function 1: bond pricing without coupon using DCF model
double bond_price_without_coupon (double fv, int n, double ytm) {
	/*
        fv is the future (maturity) value of the bond;
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond

        input fv, ytm: any positive number (double)
        input n: any positive number (int)
	*/

	double pv;  // local variable, pv for present value

	pv = fv * pow(1 + ytm, -n);

	return pv;

}


// function 2: bond duration(sensitivity of the bond price to a change in bond yield)
double bond_duration_without_coupon (double fv, int n, double ytm, double delta) {
	/*
        fv is the future (maturity) value of the bond;
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond
		delta is the change of YTM

        input fv, ytm: any positive number (double)
		input delta: any number (double)
        input n: any positive number (int)
	*/

	double pv;  // local variable, pv for present value
	double pv_plus;
	double pv_minus; 
	double duration;

	pv = bond_price_without_coupon (fv, n, ytm);
	pv_minus = bond_price_without_coupon (fv, n, ytm + delta);
	pv_plus = bond_price_without_coupon (fv, n, ytm - delta);

	duration = (((pv_plus - pv_minus) / 2) / pv) / delta;

	return duration;

}


// function 3: bond pricing with coupon using DCF model
double bond_price_with_coupon (double fv, int n, double ytm, int m, double c) {
	/*
        fv is the future (maturity) value of the bond;
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond;
		m is the number of coupon payments per year;
		c is the annual coupon rate expressed as percentage per year

        input fv, ytm, c: any positive number (double)
        input n, m: any positive number (int)		
	*/

	double pv_fv;  // local variable, pv for present value; pv of bond face value
	double pv_coupon; // pv of coupon flow
	double pv_sum;
	double coupon;

	pv_fv = fv * pow(1 + ytm / m, - (n * m));

	coupon = fv * c * (1 / m);  // coupon per payment
	pv_coupon = coupon * (1 - pow(1 + ytm / m, - (n * m))) / (ytm / m);

	pv_sum = pv_fv + pv_coupon;

	return pv_sum;

}


// function 4: coupon bond duration(sensitivity of the bond price to a change in bond yield)
double bond_duration_with_coupon (double fv, int n, double ytm, int m, double c, double delta) {
	/*
        fv is the future (maturity) value of the bond;
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond;
		m is the number of coupon payments per year;
		c is the annual coupon rate expressed as percentage per year
		delta is the change of YTM

        input fv, ytm, c: any positive number (double)
		input delta: any number (double)
        input n, m: any positive number (int)
	*/

	double pv;  // local variable, pv for present value
	double pv_plus;
	double pv_minus; 
	double duration;

	pv = bond_price_with_coupon (fv, n, ytm, m, c);
	pv_minus = bond_price_with_coupon (fv, n, ytm + delta, m, c);
	pv_plus = bond_price_with_coupon (fv, n, ytm - delta, m, c);

	duration = (((pv_plus - pv_minus) / 2) / pv) / delta;

	return duration;

}


// function 5: bond convexity witout coupon
double bond_convexity_without_coupon (double fv, int n, double ytm, double delta) {
	/*
        fv is the future (maturity) value of the bond;
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond
		delta is the change of YTM

        input fv, ytm: any positive number (double)
		input delta: any number (double)
        input n: any positive number (int)
	*/

	double pv;  // local variable, pv for present value
	double pv_plus;
	double pv_minus; 
	double convexity;

	pv = bond_price_without_coupon (fv, n, ytm);
	pv_minus = bond_price_without_coupon (fv, n, ytm + delta);
	pv_plus = bond_price_without_coupon (fv, n, ytm - delta);

	convexity = (pv_plus + pv_minus - 2 * pv) / (pv * delta * delta);

	return convexity;

}


// function 6: coupon bond convexity
double bond_convexity_with_coupon (double fv, int n, double ytm, int m, double c, double delta) {
	/*
        fv is the future (maturity) value of the bond;
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond;
		m is the number of coupon payments per year;
		c is the annual coupon rate expressed as percentage per year
		delta is the change of YTM

        input fv, ytm, c: any positive number (double)
		input delta: any number (double)
        input n, m: any positive number (int)
	*/

	double pv;  // local variable, pv for present value
	double pv_plus;
	double pv_minus; 
	double convexity;

	pv = bond_price_with_coupon (fv, n, ytm, m, c);
	pv_minus = bond_price_with_coupon (fv, n, ytm + delta, m, c);
	pv_plus = bond_price_with_coupon (fv, n, ytm - delta, m, c);

	convexity = (pv_plus + pv_minus - 2 * pv) / (pv * delta * delta);

	return convexity;

}


// function 7: portfolio pricing without coupon; portfolio consisted of 3 bonds
double portfolio_price_without_coupon (double fv1, int n1, double ytm1, double position1,
								  	   double fv2, int n2, double ytm2, double position2,
									   double fv3, int n3, double ytm3, double position3) {
	/*
        fvi is the future (maturity) value of the bond, i = 1,2,3;
        ni is the number of years until maturity, i = 1,2,3;
        ytmi is the yield to maturity of a bond, i = 1,2,3;
		positioni is the share hold of a asset, i = 1,2,3

        input fv, ytm: any positive number (double)
        input n: any positive number (int)
		input position: any number (double) positive for long and negative for short
	*/

	double pv1, pv2, pv3;  // local variable, pv for present value
	double portfolio_pv;
	
	// create iterator
	double fv_array[] = {fv1, fv2, fv3};
	int n_array[] = {n1, n2, n3};
	double ytm_array[] = {ytm1, ytm2, ytm3};
	double position_array[] = {position1, position2, position3};
	
	// store the sum of pricing 
	portfolio_pv = 0;

	for (int i = 0; i < (sizeof(n_array) / sizeof(n_array[0])); i++) {
		portfolio_pv += position_array[i] * (fv_array[i] * pow(1 + ytm_array[i], -n_array[i]));
	}

	return portfolio_pv;

}


// function 8: bond cashflow
double * bond_cashflow (double fv, int n, double ytm, int m, double c, double principal_rate) {
	/*
        fv is the future (maturity) value of the bond;
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond;
		m is the number of coupon payments per year;
		c is the annual coupon rate expressed as percentage per year;
		principal_rate is repayment rate of its principal annually

        input fv, ytm, c, principal_rate: any positive number (double)
		input delta: any number (double)
        input n, m: any positive number (int)	
	*/
	double principal_pay;
	double coupon_pay;
	double annual_pay;
	double * cashflow;
	cashflow = new double [n * m];  // empty array for store

	principal_pay = fv * principal_rate;

	for (int i = 0; i < n * m; i++) {
		coupon_pay = fv * c;
		annual_pay = principal_pay + coupon_pay;
		fv -= principal_pay;

		cashflow[i] = annual_pay;
	}

	return cashflow;

}


// function 9: bond pricing for amortizing bond using DCF
double amortizing_bond_pricing (int n, double ytm, int m, double cashflow[5]) {
	/*
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond;
		m is the number of coupon payments per year;	
		cashflow is the array containing bond cashflows of payment

        input ytm: any positive number (double)
        input n, m: any positive number (int)	
		input cashflow: array[n * m]
	
	*/
	double bond_pv;
	double cashflow_discounted;

	bond_pv = 0;

	for (int i = 0; i < n * m; i++) {
		cashflow_discounted = cashflow[i] * pow(1 + ytm / m , -(i + 1));
		bond_pv += cashflow_discounted;

	}

	return bond_pv;

}


// function 10: amortizing bond duration
double amortizing_bond_duration (int n, double ytm, int m, double cashflow[5], double delta) {
	/*
        n is the number of years until maturity;
        ytm is the yield to maturity of a bond;
		m is the number of coupon payments per year;
		cashflow is the array containing bond cashflows of payment
		delta is the change of YTM

        input fv, ytm, c: any positive number (double)
		input delta: any number (double)
        input n, m: any positive number (int)
	*/

	double pv;  // local variable, pv for present value
	double pv_plus;
	double pv_minus; 
	double duration;

	pv = amortizing_bond_pricing (n, ytm, m, cashflow);
	pv_minus = amortizing_bond_pricing (n, ytm + delta, m, cashflow);
	pv_plus = amortizing_bond_pricing (n, ytm - delta, m, cashflow);

	duration = (((pv_plus - pv_minus) / 2) / pv) / delta;

	return duration;

}














/*
Use the main() section for test call 
This section will automatically be executed when the file is run
*/

// cout = c out; endl = end line

using namespace std;

int main()
{ 

	// (a) Calculate prices of a zero coupon bond that pays $100 at maturity
	cout << setiosflags(ios::fixed) << setprecision(2) << endl;  // set the decimal format of the output
	printf("(a)\n");
	double bond_price_test1 = bond_price_without_coupon (100, 1, 0.025);
	cout << "bond_price_without_coupon (100, 1, 0.025) returns " << bond_price_test1 << endl;

	double bond_price_test2 = bond_price_without_coupon (100, 2, 0.026);
	cout << "bond_price_without_coupon (100, 2, 0.026) returns " << bond_price_test2 << endl;

	double bond_price_test3 = bond_price_without_coupon (100, 3, 0.027);
	cout << "bond_price_without_coupon (100, 3, 0.027) returns " << bond_price_test3 << endl;

	double bond_price_test4 = bond_price_without_coupon (100, 5, 0.03);
	cout << "bond_price_without_coupon (100, 5, 0.03) returns " << bond_price_test4 << endl;

	double bond_price_test5 = bond_price_without_coupon (100, 10, 0.035);
	cout << "bond_price_without_coupon (100, 10, 0.035) returns " << bond_price_test5 << endl;

	double bond_price_test6 = bond_price_without_coupon (100, 30, 0.04);
	cout << "bond_price_without_coupon (100, 30, 0.04) returns " << bond_price_test6 << endl;

	cout << endl;
	printf("Comment: The pricing of 1-year bond is the highest. It is reasonable because it discounts the least number of times.\n");
	cout << endl;


	// (b) Calculate the duration of each zero coupon bond (sensitivity of the bond price to a change in bond yield)
	printf("(b)\n");
	double delta = 0.0001;  // 1 bp
	printf("delta = %.4f \n", delta);  // string format

	double duration_test1 = bond_duration_without_coupon (100, 1, 0.025, delta);
	cout << "bond_duration_without_coupon (100, 1, 0.025, delta) returns " << duration_test1 << endl; 

	double duration_test2 = bond_duration_without_coupon (100, 2, 0.026, delta);
	cout << "bond_duration_without_coupon (100, 2, 0.026, delta) returns " << duration_test2 << endl;

	double duration_test3 = bond_duration_without_coupon (100, 3, 0.027, delta);
	cout << "bond_duration_without_coupon (100, 3, 0.027, delta) returns " << duration_test3 << endl;

	double duration_test4 = bond_duration_without_coupon (100, 5, 0.03, delta);
	cout << "bond_duration_without_coupon (100, 5, 0.03, delta) returns " << duration_test4 << endl;

	double duration_test5 = bond_duration_without_coupon (100, 10, 0.035, delta);
	cout << "bond_duration_without_coupon (100, 10, 0.035, delta) returns " << duration_test5 << endl;

	double duration_test6 = bond_duration_without_coupon (100, 30, 0.04, delta);
	cout << "bond_duration_without_coupon (100, 30, 0.04, delta) returns " << duration_test6 << endl;

	cout << endl;
	printf("Comment: Bond prices and bond yields are negatively related.\n");
	cout << endl;

	// (c) Calculate prices of coupon bonds that pay $100 at maturity at 3% annually coupon rate until maturity.
	printf("(c)\n");
	double bond_price_coupon1 = bond_price_with_coupon (100, 1, 0.025, 1, 0.03);
	cout << "bond_price_with_coupon (100, 1, 0.025, 1, 0.03) returns " << bond_price_coupon1 << endl;

	double bond_price_coupon2 = bond_price_with_coupon (100, 2, 0.026, 1, 0.03);
	cout << "bond_price_with_coupon (100, 2, 0.026, 1, 0.03) returns " << bond_price_coupon2 << endl;

	double bond_price_coupon3 = bond_price_with_coupon (100, 3, 0.027, 1, 0.03);
	cout << "bond_price_with_coupon (100, 3, 0.027, 1, 0.03) returns " << bond_price_coupon3 << endl;

	double bond_price_coupon4 = bond_price_with_coupon (100, 5, 0.03, 1, 0.03);
	cout << "bond_price_with_coupon (100, 5, 0.03, 1, 0.03) returns " << bond_price_coupon4 << endl;

	double bond_price_coupon5 = bond_price_with_coupon (100, 10, 0.035, 1, 0.03);
	cout << "bond_price_with_coupon (100, 10, 0.035, 1, 0.03) returns " << bond_price_coupon5 << endl;

	double bond_price_coupon6 = bond_price_with_coupon (100, 30, 0.04, 1, 0.03);
	cout << "bond_price_with_coupon (100, 30, 0.04, 1, 0.03) returns " << bond_price_coupon6 << endl;

	cout << endl;
	printf("Comment: Prices of 1-year, 2-year and 3-year bonds are above $100; Prices of 10-year and 30-year bonds are below. $100. The reason is that the YTM of short term bonds is lower, so the discount is less when using DCF. \n");
	cout << endl;


	// (d) Calculate the duration of each coupon bond using finite differences.
	// double delta = 0.0001;  // 1 bp
	printf("(d)\n");
	printf("delta = %.4f \n", delta);  // string format

	double duration_coupon1 = bond_duration_with_coupon (100, 1, 0.025, 1, 0.03,delta);
	cout << "bond_duration_with_coupon (100, 1, 0.025, 1, 0.03,delta) returns " << duration_coupon1 << endl; 

	double duration_coupon2 = bond_duration_with_coupon (100, 2, 0.026, 1, 0.03, delta);
	cout << "bond_duration_with_coupon (100, 2, 0.026, 1, 0.03, delta) returns " << duration_coupon2 << endl;

	double duration_coupon3 = bond_duration_with_coupon (100, 3, 0.027, 1, 0.03, delta);
	cout << "bond_duration_with_coupon (100, 3, 0.027, 1, 0.03, delta) returns " << duration_coupon3 << endl;

	double duration_coupon4 = bond_duration_with_coupon (100, 5, 0.03, 1, 0.03, delta);
	cout << "bond_duration_with_coupon (100, 5, 0.03, 1, 0.03, delta) returns " << duration_coupon4 << endl;

	double duration_coupon5 = bond_duration_with_coupon (100, 10, 0.035, 1, 0.03, delta);
	cout << "bond_duration_with_coupon (100, 10, 0.035, 1, 0.03, delta) returns " << duration_coupon5 << endl;

	double duration_coupon6 = bond_duration_with_coupon (100, 30, 0.04, 1, 0.03, delta);
	cout << "bond_duration_with_coupon (100, 30, 0.04, 1, 0.03, delta) returns " << duration_coupon6 << endl;

	cout << endl;
	printf("Comment: Zero-coupon bonds have higher duration because coupon bonds have more cashflows, which will result in earlier bond payoffs. \n");
	cout << endl;


	// (e) Calculate the convexity of each bond
	// double delta = 0.0001;  // 1 bp
	printf("(e)\n");
	printf("delta = %.4f \n", delta);  // string format

	double convexity_test1 = bond_convexity_without_coupon (100, 1, 0.025, delta);
	cout << "bond_convexity_without_coupon (100, 1, 0.025, delta) returns " << convexity_test1 << endl; 

	double convexity_test2 = bond_convexity_without_coupon (100, 2, 0.026, delta);
	cout << "bond_convexity_without_coupon (100, 2, 0.026, delta) returns " << convexity_test2 << endl;

	double convexity_test3 = bond_convexity_without_coupon (100, 3, 0.027, delta);
	cout << "bond_convexity_without_coupon (100, 3, 0.027, delta) returns " << convexity_test3 << endl;

	double convexity_test4 = bond_convexity_without_coupon (100, 5, 0.03, delta);
	cout << "bond_convexity_without_coupon (100, 5, 0.03, delta) returns " << convexity_test4 << endl;

	double convexity_test5 = bond_convexity_without_coupon (100, 10, 0.035, delta);
	cout << "bond_convexity_without_coupon (100, 10, 0.035, delta) returns " << convexity_test5 << endl;

	double convexity_test6 = bond_convexity_without_coupon (100, 30, 0.04, delta);
	cout << "bond_convexity_without_coupon (100, 30, 0.04, delta) returns " << convexity_test6 << endl;

	cout << endl;

	double convexity_coupon1 = bond_convexity_with_coupon (100, 1, 0.025, 1, 0.03, delta);
	cout << "bond_convexity_with_coupon (100, 1, 0.025, 1, 0.03, delta) returns " << convexity_coupon1 << endl; 

	double convexity_coupon2 = bond_convexity_with_coupon (100, 2, 0.026, 1, 0.03, delta);
	cout << "bond_convexity_with_coupon (100, 2, 0.026, 1, 0.03, delta) returns " << convexity_coupon2 << endl;

	double convexity_coupon3 = bond_convexity_with_coupon (100, 3, 0.027, 1, 0.03, delta);
	cout << "bond_convexity_with_coupon (100, 3, 0.027, 1, 0.03, delta) returns " << convexity_coupon3 << endl;

	double convexity_coupon4 = bond_convexity_with_coupon (100, 5, 0.03, 1, 0.03, delta);
	cout << "bond_convexity_with_coupon (100, 5, 0.03, 1, 0.03, delta) returns " << convexity_coupon4 << endl;

	double convexity_coupon5 = bond_convexity_with_coupon (100, 10, 0.035, 1, 0.03, delta);
	cout << "bond_convexity_with_coupon (100, 10, 0.035, 1, 0.03, delta) returns " << convexity_coupon5 << endl;

	double convexity_coupon6 = bond_convexity_with_coupon (100, 30, 0.04, 1, 0.03, delta);
	cout << "bond_convexity_with_coupon (100, 30, 0.04, 1, 0.03, delta) returns " << convexity_coupon6 << endl;

	cout << endl;
	printf("Comment: The second derivatives are positive. \n");
	cout << endl;


	// (f) Calculate the initial value of the portfolio (portfolio pricing)
	printf("(f)\n");
	double portfolio_pv = portfolio_price_without_coupon (100, 1, 0.025, 1,
														  100, 2, 0.026, -2,
														  100, 3, 0.027, 1);
	cout << "The initial value of the portfolio is " << portfolio_pv << endl;

	cout << endl;


	// (g) Calculate the duration and convexity of this portfolio.
	printf("(g)\n");
	double bond_price_abs_sum = 1 * bond_price_test1 + 2 * bond_price_test2 + 1 * bond_price_test3;
	double portfolio_duration_without_coupon = duration_test1 * (1 * bond_price_test1 / bond_price_abs_sum) + \
											   duration_test2 * (-2 * bond_price_test2 / bond_price_abs_sum) + \	
											   duration_test3 * (1 * bond_price_test3 / bond_price_abs_sum);
	printf("delta = %.4f \n", delta);  // string format
	cout << "The duration of this portfolio is " << abs(portfolio_duration_without_coupon) << endl;
	// according to fomula above, we can derive the expression of position of bond 2 when duration neutral (portfolio_duration = 0)

	double portfolio_convexity_without_coupon = convexity_test1 * (1 * bond_price_test1 / bond_price_abs_sum) + \
											   convexity_test2 * (-2 * bond_price_test2 / bond_price_abs_sum) + \	
											   convexity_test3 * (1 * bond_price_test3 / bond_price_abs_sum);
	cout << "The convexity of this portfolio is " << portfolio_convexity_without_coupon << endl;

	cout << endl;
	printf("Comment: The convexity is bigger. Actually the convexity is bigger than the square of the duration. \n");
	cout << endl;


	// (h) Adjust the number of units of the short position in the two year zero-coupon bond so that the portfolio is duration neutral
	// adjusted_unit is positive
	printf("(h)\n");
	double adjusted_unit = (duration_test1 * 1 * bond_price_test1 + duration_test3 * 1 * bond_price_test3) / (duration_test2 * bond_price_test2);
	cout << "In order to reach duration neutral, the adjusted number of units of the short position in the two year zero-coupon bond is " << adjusted_unit << endl;
	// check duration
	double portfolio_duration_without_coupon_adjusted = duration_test1 * (1 * bond_price_test1 / bond_price_abs_sum) + \
											   duration_test2 * (-adjusted_unit * bond_price_test2 / bond_price_abs_sum) + \	
											   duration_test3 * (1 * bond_price_test3 / bond_price_abs_sum);
	cout << "The adjusted duration of this portfolio is " << abs(portfolio_duration_without_coupon_adjusted) << endl;	

	cout << endl;


	// (i) adjusted portfolio present value when delta = 100 bps 
	printf("(i)\n");
	double portfolio_pv_adjusted = portfolio_price_without_coupon (100, 1, 0.025, 1,
																   100, 2, 0.026, -adjusted_unit,
																   100, 3, 0.027, 1);
	cout << "The adjusted portfolio (duration netral) present value is " << portfolio_pv_adjusted << endl;	
	delta = 0.01;  // 100 bps
	printf("delta = %.2f \n", delta);  // string format	
	double portfolio_pv_adjusted_100bp_up = portfolio_price_without_coupon (100, 1, 0.025 + delta, 1,
														  					100, 2, 0.026 + delta, -adjusted_unit,
														  					100, 3, 0.027 + delta, 1);
	cout << "The adjusted portfolio present value when delta = 100 bps is " << portfolio_pv_adjusted_100bp_up << endl;
	cout << endl;

	// (j) adjusted portfolio present value when delta = -100 bps 
	printf("(j)\n");
	delta = -0.01;  // -100 bps
	printf("delta = %.2f \n", delta);  // string format	
	double portfolio_pv_adjusted_100bp_down = portfolio_price_without_coupon (100, 1, 0.025 - delta, 1,
														  					100, 2, 0.026 - delta, -adjusted_unit,
														  					100, 3, 0.027 - delta, 1);
	cout << "The adjusted portfolio present value when delta = -100 bps is " << portfolio_pv_adjusted_100bp_down << endl;
	cout << endl;
	printf("Portfolio present value barely changes, proving duration-neutral !!! This is a portfolio we would want to own because there is no Interest Rate Risk. \n");

	cout << endl;


	// (k) Print the cashflows of bond annul payment
	printf("(k)\n");
	double * cashflow;
	cashflow = bond_cashflow (100, 5, 0.03, 1, 0.03, 0.2);
	cout << "cashflows of a 5-year amortizing bond that repays 20% of its principal annually and pays a 3% coupon (annually): " << endl;	

	for (int i = 0; i < 5; i++){
		cout << cashflow[i] << " ";
		cout << endl;
	}

 	cout << endl;


	// (i) Calculate the price and duration of the amortizing bond using finite differences.
	printf("(i)\n");
	double amortizing_bond_pv =	amortizing_bond_pricing (5, 0.03, 1, cashflow);
	cout << "bond pricing for amortizing bond is " << amortizing_bond_pv << endl;

	delta = 0.0001;  // 1 bp
	double amortizing_duration = amortizing_bond_duration (5, 0.03, 1, cashflow, delta);
	cout << "duration for amortizing bond is " << amortizing_duration << endl;
	cout << endl;
	cout << "Bond prcing: amortizing bond > coupon bond > zero-coupon bond." << endl;
	cout << "Bond duration: amortizing bond < coupon bond < zero-coupon bond." << endl;	


}


















