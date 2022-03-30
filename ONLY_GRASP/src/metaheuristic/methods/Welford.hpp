/**
 * Welford's online algorithm to compute mean and standard deviation.
 * https://gist.github.com/alexalemi/2151722
 */

class Welford {

 public:

    Welford() = default;

    void update(float x) {
        ++k;
        const auto new_mean = mean + (x - mean)*1.0f/static_cast<float>(k);
        //const auto new_standard_deviation = standard_deviation + (x - mean)*(x - new_mean);
        mean = new_mean;
        //standard_deviation = new_standard_deviation;
    }

    float get_mean() { return static_cast<float>(mean); }

    //float get_standard_deviation() {
    //    if (k==1) { return 0; }
    //    return std::sqrt(standard_deviation/(static_cast<float>(k) - 1.0f));
    //}

    void reset() {
        k = 0;
        mean = 0.0;
        //standard_deviation = 0.0;
    }

 private:

    unsigned long k = 0;
    double mean = 0.0;
    //double standard_deviation = 0.0;


};

