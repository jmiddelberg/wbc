#ifndef WBC_TOOLS_FILTER_HPP
#define WBC_TOOLS_FILTER_HPP  

#include <vector>
#include <Eigen/Core>

namespace wbc {

class Filter{
protected:
    size_t window_size_;
public:
    Filter(size_t window_size) : window_size_(window_size) {}
    virtual ~Filter() {}
    virtual double apply(double new_value) = 0;
};

class MovingMedianFilter : public Filter {
protected:
    Eigen::VectorXd values_;
    Eigen::VectorXd sorted_values_;
public:
    MovingMedianFilter(size_t window_size) : Filter(window_size) {}

    virtual double apply(double new_value);

};

class MovingAverageFilter : public Filter {
protected:
    Eigen::VectorXd values_;
public:
    MovingAverageFilter(size_t window_size) : Filter(window_size) {}

    virtual double apply(double new_value);

};

/**
 * First order (exponential) low-pass filter for vector-valued signals:
 * \f[
 *    y_k = \alpha x_k + (1 - \alpha) y_{k-1}, \qquad \alpha = \frac{\Delta t}{\Delta t + \frac{1}{2\pi f_c}}
 * \f]
 * The filter state is initialized with the first input sample, so there is no transient from zero.
 */
class LowPassFilter {
protected:
    double alpha_;
    bool initialized_;
    Eigen::VectorXd y_;
public:
    /** @param cutoff_freq Cutoff frequency \f$f_c\f$ in Hz. Has to be > 0.
      * @param sample_time Time \f$\Delta t\f$ in seconds between two consecutive calls of apply(). Has to be > 0.*/
    LowPassFilter(double cutoff_freq, double sample_time);

    /** Add a new sample and return the filtered signal*/
    const Eigen::VectorXd& apply(const Eigen::VectorXd& new_value);

    /** Re-initialize the filter state with the next input sample*/
    void reset(){initialized_ = false;}
};

}

#endif  // WBC_TOOLS_FILTER_HPP