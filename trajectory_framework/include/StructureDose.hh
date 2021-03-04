#ifndef StructureDose_H
#define StructureDose_H 1
#include <vector>
#include <string>
#include <fstream>
#include "json.hh"
#include "HardConstraint.hh"
#include "DVConstraint.hh"
#include "MeanConstraint.hh"
#include "SparseDose.hh"

class StructureDose {
public:
	StructureDose();
	std::vector<std::vector<double> > doses;
	int structure_index;
}


#endif