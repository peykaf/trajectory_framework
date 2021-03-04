#ifndef KVATBeamlet_H
#define KVATBeamlet_H 1

#include "WeightClass.hh"
#include "json.hh"
#include <vector>
#include <fstream>
#include <algorithm>


class KVATBeamlet : public WeightClass {
public:
	KVATBeamlet() {};

	nlohmann::json to_json() {return nullptr;};
	int index;
};

#endif
