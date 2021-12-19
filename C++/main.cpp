#include <iostream>
#include "bessel.h"
#include "utils.h"
#include "exact.h"
#include "parallel.h"
#include <boost/random.hpp>
#include <complex>
#include <random>
#include <fstream>


// Runs on a single core (standard).
int main()
{

    // Model parameters
    double kappa = 1;
    double theta = 0.2*0.2;
    double sigma = 0.3;
    double v0 = 0.2*0.2;
    double r = 0.1;
    double S0 = 100.;
    double K = 100.;
    double rho = -0.7;
    double T = 30.0/365.0;

    // Simulation parameters
    const long numberPaths = 10000;
    const int steps = 30;

    boost::mt19937 rng(200000);
    boost::uniform_real<> uni_dist(0, 1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);

    double dt = T / ((double)(steps));

    vector<double> times;
    for (int i = 0; i <= steps; i++) {
        times.push_back(i * dt);
    }
        
    double S[steps+1]{};
    double V[steps+1]{};
    double intV[steps+1]{};

    double intsqrvdW;
    double m;
    double s;

    std::ofstream write_output1;
    std::ofstream write_output2;
    std::ofstream write_output3;

    write_output1.open("S1_paths.txt",std::ios_base::app);
    write_output2.open("V_paths.txt", std::ios_base::app);
    write_output3.open("intV_paths.txt", std::ios_base::app);


    assert(write_output1.is_open());
    assert(write_output2.is_open());
    assert(write_output3.is_open());

    for (int i = 0; i < numberPaths; i++)
    {
        S[0] = S0;
        V[0] = v0;
        intV[0] = 0;
        for (long j = 1; j <= steps; j++)
        {
            V[j] = quantileSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V[j - 1], times[j]);
            intV[j] = quantileIntSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V[j - 1], times[j], V[j]);
            intsqrvdW = (V[j] - V[j - 1] - kappa * theta * dt + kappa * intV[j]) / sigma;
            m = (r * dt - 0.5 * intV[j] + rho * intsqrvdW);
            s = sqrt((1 - rho * rho) * intV[j]);
            S[j] = S[j - 1] * exp(m + s * qnorm(uni()));
            write_output1 << S[j - 1] << " ";
            write_output2 << V[j - 1] << " ";
            write_output3 << intV[j - 1] << " ";
        }
        write_output1 << S[steps] << "\n";
        write_output2 << V[steps] << "\n";
        write_output3 << intV[steps] << "\n";

        if (i % 100 == 0) { std::cout << "Finished " << i << " of " << numberPaths << " steps. \n"; }
    }

    write_output1.close();
    write_output2.close();
    write_output3.close();


	return 0;
}



// Runs on 4 cores!
//int main()
//{
//
//     Model parameters
//    double kappa = 1;
//    double theta = 0.2*0.2;
//    double sigma = 1;
//    double v0 = 0.2*0.2;
//    double r = 0.1;
//    double S0 = 100.;
//    double K = 100.;
//    double rho = -0.7;
//    double T = 30.0/365.0;
//
//     Simulation parameters
//    const long numberPaths = 10000;
//    const int steps = 30;
//
//     Generating random variables
//    boost::mt19937 rng(3213217);
//    boost::uniform_real<> uni_dist(0, 1);
//    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
//
//    double dt = T / ((double)(steps));
//
//    vector<double> times;
//    for (int i = 0; i <= steps; i++) {
//        times.push_back(i * dt);
//    }
//    
//     Creating one container for each of the 4 cores, since we want to run loops entirely parallel.
//    double S1[steps+1]{};
//    double V1[steps+1]{};
//    double intV1[steps+1]{};
//
//    double intsqrvdW1;
//    double m1;
//    double s1;
//
//    double S2[steps + 1]{};
//    double V2[steps + 1]{};
//    double intV2[steps + 1]{};
//
//    double intsqrvdW2;
//    double m2;
//    double s2;
//
//    double S3[steps + 1]{};
//    double V3[steps + 1]{};
//    double intV3[steps + 1]{};
//
//    double intsqrvdW3;
//    double m3;
//    double s3;
//
//    double S4[steps + 1]{};
//    double V4[steps + 1]{};
//    double intV4[steps + 1]{};
//
//    double intsqrvdW4;
//    double m4;
//    double s4;
//
//    std::ofstream write_output11;
//    std::ofstream write_output21;
//    std::ofstream write_output31;
//    std::ofstream write_output12;
//    std::ofstream write_output22;
//    std::ofstream write_output32;
//    std::ofstream write_output13;
//    std::ofstream write_output23;
//    std::ofstream write_output33;
//    std::ofstream write_output14;
//    std::ofstream write_output24;
//    std::ofstream write_output34;
//
//    write_output11.open("S1_high1.txt",std::ios_base::app);
//    write_output21.open("V_high1.txt", std::ios_base::app);
//    write_output31.open("intV_high1.txt", std::ios_base::app);
//    write_output12.open("S1_high2.txt", std::ios_base::app);
//    write_output22.open("V_high2.txt", std::ios_base::app);
//    write_output32.open("intV_high2.txt", std::ios_base::app);
//    write_output13.open("S1_high3.txt", std::ios_base::app);
//    write_output23.open("V_high3.txt", std::ios_base::app);
//    write_output33.open("intV_high3.txt", std::ios_base::app);
//    write_output14.open("S1_high4.txt", std::ios_base::app);
//    write_output24.open("V_high4.txt", std::ios_base::app);
//    write_output34.open("intV_high4.txt", std::ios_base::app);
//
//    assert(write_output11.is_open());
//    assert(write_output21.is_open());
//    assert(write_output31.is_open());
//
//    assert(write_output12.is_open());
//    assert(write_output22.is_open());
//    assert(write_output32.is_open());
//
//    assert(write_output13.is_open());
//    assert(write_output23.is_open());
//    assert(write_output33.is_open());
//
//    assert(write_output14.is_open());
//    assert(write_output24.is_open());
//    assert(write_output34.is_open());
//  
//
//    parallel_for(numberPaths, [&](int start, int end) {
//        for (int i = start; i < end; ++i)
//        {
//             The parallelized loop splits the loop into 4 (one for each core)
//             To avoid problems with saving data, we need each core to save into their own file
//            if (i >= numberPaths * 3 / 4) {
//                S1[0] = S0;
//                V1[0] = v0;
//                intV1[0] = 0;
//                for (long j = 1; j <= steps; j++)
//                {
//                    V1[j] = quantileSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V1[j - 1], times[j]);
//                    intV1[j] = quantileIntSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V1[j - 1], times[j], V1[j]);
//                    intsqrvdW1 = (V1[j] - V1[j - 1] - kappa * theta * dt + kappa * intV1[j]) / sigma;
//                    m1 = (r * dt - 0.5 * intV1[j] + rho * intsqrvdW1);
//                    s1 = sqrt((1 - rho * rho) * intV1[j]);
//                    S1[j] = S1[j - 1] * exp(m1 + s1 * qnorm(uni()));
//                    write_output11 << S1[j-1] << " ";
//                    write_output21 << V1[j-1] << " ";
//                    write_output31 << intV1[j-1] << " ";
//                }
//                write_output11 << S1[steps] << "\n";
//                write_output21 << V1[steps] << "\n";
//                write_output31 << intV1[steps] << "\n";
//
//            }
//            if (i >= numberPaths * 2 / 4 && i < numberPaths * 3 / 4) {
//                S2[0] = S0;
//                V2[0] = v0;
//                intV2[0] = 0;
//                for (long j = 1; j <= steps; j++)
//                {
//                    V2[j] = quantileSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V2[j - 1], times[j]);
//                    intV2[j] = quantileIntSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V2[j - 1], times[j], V2[j]);
//                    intsqrvdW2 = (V2[j] - V2[j - 1] - kappa * theta * dt + kappa * intV2[j]) / sigma;
//                    m2 = (r * dt - 0.5 * intV2[j] + rho * intsqrvdW2);
//                    s2 = sqrt((1 - rho * rho) * intV2[j]);
//                    S2[j] = S2[j - 1] * exp(m2 + s2 * qnorm(uni()));
//                    write_output12 << S2[j - 1] << " ";
//                    write_output22 << V2[j - 1] << " ";
//                    write_output32 << intV2[j - 1] << " ";
//                }
//                write_output12 << S2[steps] << "\n";
//                write_output22 << V2[steps] << "\n";
//                write_output32 << intV2[steps] << "\n";
//            }
//            if (i >= numberPaths * 1 / 4 && i < numberPaths * 2 / 4) {
//                S3[0] = S0;
//                V3[0] = v0;
//                intV3[0] = 0;
//                for (long j = 1; j <= steps; j++)
//                {
//                    V3[j] = quantileSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V3[j - 1], times[j]);
//                    intV3[j] = quantileIntSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V3[j - 1], times[j], V3[j]);
//                    intsqrvdW3 = (V3[j] - V3[j - 1] - kappa * theta * dt + kappa * intV3[j]) / sigma;
//                    m3 = (r * dt - 0.5 * intV3[j] + rho * intsqrvdW3);
//                    s3 = sqrt((1 - rho * rho) * intV3[j]);
//                    S3[j] = S3[j - 1] * exp(m3 + s3 * qnorm(uni()));
//                    write_output13 << S3[j - 1] << " ";
//                    write_output23 << V3[j - 1] << " ";
//                    write_output33 << intV3[j - 1] << " ";
//                }
//                write_output13 << S3[steps] << "\n";
//                write_output23 << V3[steps] << "\n";
//                write_output33 << intV3[steps] << "\n";
//            }
//            if (i < numberPaths * 1 / 4) {
//                S4[0] = S0;
//                V4[0] = v0;
//                intV4[0] = 0;
//                for (long j = 1; j <= steps; j++)
//                {
//                    V4[j] = quantileSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V4[j - 1], times[j]);
//                    intV4[j] = quantileIntSquareRootProcess(uni(), kappa, theta, sigma, times[j - 1], V4[j - 1], times[j], V4[j]);
//                    intsqrvdW4 = (V4[j] - V4[j - 1] - kappa * theta * dt + kappa * intV4[j]) / sigma;
//                    m4 = (r * dt - 0.5 * intV4[j] + rho * intsqrvdW4);
//                    s4 = sqrt((1 - rho * rho) * intV4[j]);
//                    S4[j] = S4[j - 1] * exp(m4 + s4 * qnorm(uni()));
//                    write_output14 << S4[j - 1] << " ";
//                    write_output24 << V4[j - 1] << " ";
//                    write_output34 << intV4[j - 1] << " ";
//                }
//                write_output14 << S4[steps] << "\n";
//                write_output24 << V4[steps] << "\n";
//                write_output34 << intV4[steps] << "\n";
//            }
//            if (i % 100 == 0) { std::cout << "Finished " << i << " of " << numberPaths << " steps. \n";}
//        }
//    });
//
//    write_output11.close();
//    write_output21.close();
//    write_output31.close();
//    write_output12.close();
//    write_output22.close();
//    write_output32.close();
//    write_output13.close();
//    write_output23.close();
//    write_output33.close();
//    write_output14.close();
//    write_output24.close();
//    write_output34.close();
//    
//
//	return 0;
//}

