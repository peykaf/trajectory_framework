#include "RecalcAperture.hh"
#include "json.hh"
#include "DenseDose.hh"

using json = nlohmann::json;

RecalcAperture::RecalcAperture(json ap_json) {
    int row_index = 0;
    for (auto row : ap_json["rows"]) {
        RecalcApertureRow ap_row;
        ap_row.l_bound = row[0];
        ap_row.r_bound = row[1];
        ap_row.row_index = row_index;
        this->rows.push_back(ap_row);
        row_index += 1;
    }

    DenseDose dose;
    dose.from_file(ap_json["dose_filename"]);
    this->dose = dose.grid;
}

json RecalcAperture::to_json() {
    std::vector<std::vector<double> > ap_vector;
    for (auto row : this->rows) {
        std::vector<double> ap_row = {row.l_bound, row.r_bound};
        ap_vector.push_back(ap_row);
    }

    json ap_json;
    ap_json["rows"] = ap_vector;
    ap_json["weight"] = this->weight;

    return ap_json;
}