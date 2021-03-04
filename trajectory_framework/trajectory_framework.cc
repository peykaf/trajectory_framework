// Copyright Marc-Andre Renaud 2017

#include <string>
#include "json.hh"

#include "TrajectoryGenerator.hh"
#include "IpoptOptimiser.hh"
#include "PricingIPOPTWrapper.hh"
#include "PricingGenerator.hh"
#include "RobustIPOPTWrapper.hh"
#include "WeightedRobustIPOPTWrapper.hh"
#include "RecalcIPOPTWrapper.hh"

#include "SAOptimisation.hh"
#include "IPSAOptimisation.hh"
#include "FluenceMapOptimisation.hh"
#include "FMO.hh"
#include "BrachyOptimisation.hh"
#include "RecalcOptimisation.hh"
#include "MixedOptimisation.hh"
#include "MixedArcOptimisation.hh"
#include "KVATOptimisation.hh"
#include "RobustOptimisation.hh"
#include "RobustArcOptimisation.hh"
#include "WeightedRobustOptimisation.hh"
#include "WeightedRobustArcOptimisation.hh"
#include "WeightedRobustRecalc.hh"

using json = nlohmann::json;

void launch_optimisation(TrajectoryGenerator &traj_obj)
{
    IpoptOptimiser *opt_obj = new IpoptOptimiser(&traj_obj);
    opt_obj->run();
}

void launch_optimisation(PricingGenerator &traj_obj)
{
    PricingIPOPTWrapper *opt_obj = new PricingIPOPTWrapper(&traj_obj);
    opt_obj->run();
}

void launch_optimisation(RobustOptimisation &traj_obj)
{
    RobustIPOPTWrapper *opt_obj = new RobustIPOPTWrapper(&traj_obj);
    opt_obj->run();
}

void launch_optimisation(WeightedRobustOptimisation &traj_obj)
{
    WeightedRobustIPOPTWrapper *opt_obj = new WeightedRobustIPOPTWrapper(&traj_obj);
    opt_obj->run();
}

void launch_optimisation(WeightedRobustRecalc &traj_obj)
{
    RecalcIPOPTWrapper *opt_obj = new RecalcIPOPTWrapper(&traj_obj);
    opt_obj->run();
}

void read_input_file(std::string filename)
{
    std::string buffer;
    std::vector<std::string> temp_vector;

    std::cout << "Opening: " << filename << "\n";
    std::ifstream input_file(filename);

    try
    {
        if (!input_file.is_open())
            throw "Could not open file";

        std::string extension = find_extension(filename);
        // Read whole file into string
        std::string input_string((std::istreambuf_iterator<char>(input_file)),
                                 (std::istreambuf_iterator<char>()));
        input_file.close();

        auto input_json = json::parse(input_string);
        std::string name = remove_extension(filename);
        input_json["name"] = name;

        if (input_json["opt_type"] == "brachy")
        {
            BrachyOptimisation traj_obj = BrachyOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "FMO")
        {
            FMO traj_obj = FMO(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "sparse_fluence_map")
        {
            FluenceMapOptimisation traj_obj = FluenceMapOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "recalc")
        {
            RecalcOptimisation traj_obj = RecalcOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "mixed")
        {
            MixedOptimisation traj_obj = MixedOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "mixed_arc")
        {
            MixedArcOptimisation traj_obj = MixedArcOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "kvat")
        {
            KVATOptimisation traj_obj = KVATOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "IPSA")
        {
            IPSAOptimisation traj_obj = IPSAOptimisation(input_json);
            traj_obj.launch();
        }
        else if (input_json["opt_type"] == "robust")
        {
            RobustOptimisation traj_obj = RobustOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "robust_arc")
        {
            RobustArcOptimisation traj_obj = RobustArcOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "weighted_robust")
        {
            WeightedRobustOptimisation traj_obj = WeightedRobustOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "weighted_robust_arc")
        {
            WeightedRobustArcOptimisation traj_obj = WeightedRobustArcOptimisation(input_json);
            launch_optimisation(traj_obj);
        }
        else if (input_json["opt_type"] == "weighted_robust_recalc")
        {
            WeightedRobustRecalc traj_obj = WeightedRobustRecalc(input_json);
            launch_optimisation(traj_obj);
        }
    }
    catch (const std::string err)
    {
        std::cout << "Error!\n";
        std::cout << err << "\n";
        exit(1);
    }
    catch (char const *err)
    {
        std::cout << "Error!\n";
        std::cout << err << "\n";
        exit(1);
    }
}

int main(int argc, char *argv[])
{
    if (argc <= 1)
    {
        std::cout << "syntax: trajectory_framework <filename>" << std::endl;
        return -1;
    }

    std::string filename(argv[1]);
    read_input_file(filename);
    return 0;
}
