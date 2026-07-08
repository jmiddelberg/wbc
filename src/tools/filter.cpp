#include "filter.hpp"
#include <iostream>
#include <cmath>
#include <cassert>

namespace wbc {

double MovingMedianFilter::apply(double new_value) 
{
    if (values_.size() < window_size_)
        values_.conservativeResize(values_.size()+1);
    else
        values_.segment(0,window_size_-1) = values_.segment(1,window_size_-1);
    values_[values_.size()-1] = new_value;
    sorted_values_ = values_;

    std::sort(std::begin(sorted_values_), std::end(sorted_values_));
    size_t mid = sorted_values_.size() / 2;

    if (sorted_values_.size() % 2 == 0) {
        return (sorted_values_[mid - 1] + sorted_values_[mid]) / 2.0;
    } else {
        return sorted_values_[mid];
    }
}

double MovingAverageFilter::apply(double new_value) {
    if (values_.size() < window_size_)
        values_.conservativeResize(values_.size()+1);
    else
        values_.segment(0,window_size_-1) = values_.segment(1,window_size_-1);
    values_[values_.size()-1] = new_value;
    return values_.mean();
}

LowPassFilter::LowPassFilter(double cutoff_freq, double sample_time) : initialized_(false) {
    assert(cutoff_freq > 0);
    assert(sample_time > 0);
    double rc = 1.0 / (2.0 * M_PI * cutoff_freq);
    alpha_ = sample_time / (sample_time + rc);
}

const Eigen::VectorXd& LowPassFilter::apply(const Eigen::VectorXd& new_value) {
    if (!initialized_ || y_.size() != new_value.size()) {
        y_ = new_value;
        initialized_ = true;
    }
    else
        y_ = alpha_ * new_value + (1.0 - alpha_) * y_;
    return y_;
}


}  // namespace wbc