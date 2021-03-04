#include <vector>
#include <algorithm>

#include "WeightClass.hh"


WeightClass::WeightClass() {}

double WeightClass::normalize_dose() {
    double max_ap_dose = *std::max_element(this->dose.begin(), this->dose.end());
    std::for_each(this->dose.begin(), this->dose.end(), [max_ap_dose](double &d){ d /= max_ap_dose; });
    this->norm = max_ap_dose;
    return max_ap_dose;
}

double WeightClass::average_dose() {
    double avg_dose = 0.0;
    int num_nonzero = 0;
    for (size_t i = 0; i < this->dose.size(); i++) {
        if (this->dose[i] > 0) {
            avg_dose += this->dose[i];
            num_nonzero += 1;
        }
    }

    avg_dose /= num_nonzero;

	std::for_each(this->dose.begin(), this->dose.end(), [avg_dose](double &d){ d /= avg_dose; });
    return avg_dose;
}

