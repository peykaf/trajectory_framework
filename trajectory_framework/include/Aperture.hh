#ifndef APERTURE_H
#define APERTURE_H 1
#include <vector>
#include "json.hh"
#include "WeightClass.hh"

class ApertureRow {
public:
	int row_index;
    size_t l_bound;
    size_t r_bound;
};

class Aperture : public WeightClass {
public:
	Aperture() {}
	Aperture(size_t num_rows, size_t num_columns) {
		this->init_aperture(num_rows);
		this->num_beamlet_columns = num_columns;
	}

	void init_aperture(size_t num_rows, size_t num_columns) {
		this->init_aperture(num_rows);
		this->num_beamlet_columns = num_columns;
	}

	void init_aperture(size_t num_rows) {
		this->rows.resize(num_rows);
		for (size_t row_i = 0; row_i < this->rows.size(); row_i++) {
			this->rows[row_i].row_index = row_i;
			this->rows[row_i].l_bound = 0;
			this->rows[row_i].r_bound = 0;
		}
	}

	void print_aperture() {
		for (auto &row: this->rows) {
			std::cout << "|";
			for (size_t i = 0; i < row.l_bound; i++) {
				std::cout << "-";
			}
			for (size_t i = row.l_bound; i < row.r_bound; i++) {
				std::cout << "O";
			}
			for (size_t i = row.r_bound; i < this->num_beamlet_columns; i++) {
				std::cout << "-";
			}
			std::cout << "|";
			std::cout << std::endl;
		}
	}

	std::vector<ApertureRow> rows;
	size_t num_beamlet_columns;
	double size;
	nlohmann::json to_json() {return nullptr;};
};
#endif
