#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <tools/filter.hpp>
#include <tools/Logger.hpp>
#include <chrono>
#include <fstream>

loglevel_e loglevel = logWARNING;

using namespace std;
using namespace wbc;

BOOST_AUTO_TEST_CASE(test_moving_median_filter){
    MovingMedianFilter filter(5);
    BOOST_CHECK_EQUAL(filter.apply(5.0), 5.0);
    BOOST_CHECK_EQUAL(filter.apply(9.0), 7.0);
    BOOST_CHECK_EQUAL(filter.apply(3.0), 5.0);
    BOOST_CHECK_EQUAL(filter.apply(4.0), 4.5);
    BOOST_CHECK_EQUAL(filter.apply(1.0), 4.0);
}

BOOST_AUTO_TEST_CASE(test_low_pass_filter){
    const double cutoff_freq = 10.0; // Hz
    const double sample_time = 1e-3; // s
    LowPassFilter filter(cutoff_freq, sample_time);

    Eigen::VectorXd x0 = Eigen::VectorXd::Constant(3, 2.0);
    Eigen::VectorXd x1 = Eigen::VectorXd::Constant(3, 4.0);

    // First sample initializes the filter state
    BOOST_CHECK_SMALL((filter.apply(x0) - x0).norm(), 1e-12);

    // Single step: y_1 = alpha*x_1 + (1-alpha)*y_0
    double alpha = sample_time / (sample_time + 1.0 / (2.0 * M_PI * cutoff_freq));
    Eigen::VectorXd y1_expected = alpha * x1 + (1.0 - alpha) * x0;
    BOOST_CHECK_SMALL((filter.apply(x1) - y1_expected).norm(), 1e-12);

    // Step response converges to the input
    for(int i = 0; i < 10000; i++)
        filter.apply(x1);
    BOOST_CHECK_SMALL((filter.apply(x1) - x1).norm(), 1e-6);

    // reset() re-initializes the state with the next sample
    filter.reset();
    BOOST_CHECK_SMALL((filter.apply(x0) - x0).norm(), 1e-12);
}

BOOST_AUTO_TEST_CASE(test_multi_dim){

    const int N_JOINT = 20;
    const int WINDOWSIZE = 50;

    std::vector<MovingAverageFilter> filters;
    for(int i = 0; i < N_JOINT; i++)
        filters.push_back(MovingAverageFilter(WINDOWSIZE));
    
    ofstream myfile;
    myfile.open ("filter_benchmark.txt");

    double avg_time = 0;
    int count = 0;
    for(double t = 0; t < 2*M_PI; t+=0.001){
        Eigen::VectorXd input = Eigen::VectorXd::Random(N_JOINT) / 10.0;
        for(int j = 0; j < N_JOINT; j++)
            input[j]+=sin(t);
        if(count % 1000 == 0)
            input += Eigen::VectorXd::Random(N_JOINT);
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for(int j = 0; j < N_JOINT; j++){
            double res = filters[j].apply(input[j]);
            myfile <<input[j]<<" "<<res<<" ";
        }
        myfile << sin(t);
        myfile << std::endl;        
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        avg_time += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        count++;
    }
    avg_time /= count;
    std::cout<<"Avg time: "<<avg_time<<" [us]"<<std::endl;
}
