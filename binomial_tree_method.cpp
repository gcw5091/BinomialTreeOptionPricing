// Black-Scholes Method

#include <iostream>
#include <cmath>

using namespace std;

// Black-Scholes formula for European call option
double blackScholesCall(double S, double K, double r, double sigma, double T) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);

    double N_d1 = 0.5 * (1 + erf(d1 / sqrt(2)));
    double N_d2 = 0.5 * (1 + erf(d2 / sqrt(2)));

    return S * N_d1 - K * exp(-r * T) * N_d2;
}

int main() {
    double S, K, r, sigma, T;

    // Input parameters
    cout << "Enter current stock price (S): ";
    cin >> S;
    cout << "Enter strike price (K): ";
    cin >> K;
    cout << "Enter risk-free interest rate (r): ";
    cin >> r;
    cout << "Enter volatility (sigma): ";
    cin >> sigma;
    cout << "Enter time to expiration in years (T): ";
    cin >> T;

    // Calculate and output option price
    double optionPrice = blackScholesCall(S, K, r, sigma, T);
    cout << "Option price: " << optionPrice << endl;

    return 0;
}