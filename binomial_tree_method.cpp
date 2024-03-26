#include <iostream>
#include <cmath>
#include <vector>

// CDF of the standard normal distribution
double cdf(double x) {
    return 0.5 * (1 + erf(x / sqrt(2)));
}

// Black Scholes Model for Call Options
double blackScholesCall(double S0, double K, double T, double r, double sigma) {
    double d1 = (log(S0 / K) + (r + sigma * sigma * 0.5) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);

    double callPrice = S0 * cdf(d1) - K * exp(-r * T) * cdf(d2);
    return callPrice;
}

double BinomialTree(double S0, double K, double T, double r, double sigma, int N) {
    double delta_t = T / N;
    double u = exp(sigma * sqrt(delta_t));
    double d = 1 / u;
    double p = (exp(r * delta_t) - d)/(u - d);
    // Initialize the price array to store option prices at each node
    std::vector<double> optionPrices(N + 1); // Adding 1 to include index 0 to N
    // Calculate option prices at expiration
    for (int i = 0; i <= N; ++i) {
        double ST = S0 * pow(u, N - i) * pow(d, i);
        optionPrices[i] = std::max(ST - K, 0.0);  // Payoff of a call option at expiration
    }

    // Calculate option prices at previous nodes. Overwrite the vector to save memory
    for (int j = N - 1; j >= 0; --j) {
        for (int i = 0; i <= j; ++i) {
            optionPrices[i] = std::max(exp(-r * delta_t) * (optionPrices[i] * p + optionPrices[i+1] * (1 - p)), 0.0);
        }
    }

    // The option price at time 0 is at index 0
    return optionPrices[0];
}

int main() {
    double S0 = 10;            // Current stock price
    double K = 10;            // Strike price
    double T = 0.25;         // Time to expiration (in years)
    double r = 0.02;        // Risk-free interest rate
    double sigma = 0.25;   // Volatility
    // Number of Time Steps for Binomial Tree
    int N_1 = 10;      
    int N_2 = 100;    
    int N_3 = 1000;    
    int N_4 = 10000;    
    int N_5 = 100000;    

    double callPrice_BS = blackScholesCall(S0, K, T, r, sigma);
    std::cout << "Black Scholes Model Call Option Price:              " << callPrice_BS << std::endl;
    double callPrice_1 = BinomialTree(S0, K, T, r, sigma, N_1);
    double callPrice_2 = BinomialTree(S0, K, T, r, sigma, N_2);
    double callPrice_3 = BinomialTree(S0, K, T, r, sigma, N_3);
    double callPrice_4 = BinomialTree(S0, K, T, r, sigma, N_4);
    double callError_1 = log10(callPrice_BS - callPrice_1);
    double callError_2 = log10(callPrice_BS - callPrice_2);
    double callError_3 = log10(callPrice_BS - callPrice_3);
    double callError_4 = log10(callPrice_BS - callPrice_4);

    std::cout << "Binomial Tree Method Call Option Price N = 10:      " << callPrice_1 << std::endl;
    std::cout << "Binomial Tree Method Call Option Price N = 100:     " << callPrice_2 << std::endl;
    std::cout << "Binomial Tree Method Call Option Price N = 1,000:   " << callPrice_3 << std::endl;
    std::cout << "Binomial Tree Method Call Option Price N = 10,000:  " << callPrice_4 << std::endl;

    // Provide E for each value of N
    std::cout << std::endl;
    std::cout << "Binomial Tree Method Call Option Log(Error) N = 10:      " << callError_1 << std::endl;
    std::cout << "Binomial Tree Method Call Option Log(Error) N = 100:     " << callError_2 << std::endl;
    std::cout << "Binomial Tree Method Call Option Log(Error) N = 1,000:   " << callError_3 << std::endl;
    std::cout << "Binomial Tree Method Call Option Log(Error) N = 10,000:  " << callError_4 << std::endl;

    return 0;
}
