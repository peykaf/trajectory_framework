import json
import sys
import glob
import matplotlib
import numpy
import argparse
from matplotlib import pyplot
import matplotlib.patches as mpatches

class OptOutput():
    def __init__(self, filename):
        self.filename = filename
        self.name = ".".join(self.filename.split(".")[:-1]).replace("_", " ")

        self.modalities = {
            "photons": {
                "beams": {},
            },
            "electrons": {
                "beams": {},
            }
        }

        self.overall_data = {}

        self.read_file(filename)


    def read_file(self, filename):
        with open(filename) as myfile:
            data = json.load(myfile)
            data = data[1:]
            num_apertures = []
            cost_functions = []
            aperture_sizes = []
            total_weights = []
            running_times = []
            pruned_apertures = []

            electron_beams = {}
            photon_beams = {}

            for iteration in data:
                if iteration["active_apertures"] <= 0:
                    continue

                num_apertures.append(iteration["active_apertures"])
                cost_functions.append(iteration["cost_function_value"])
                aperture_sizes.append(iteration["avg_aperture_size"])
                running_times.append(iteration.get("running_time", 0))
                pruned_apertures.append(iteration.get("pruned_apertures", 0))

                total_weight = 0

                if "electron_cpt_stats" in iteration or "cpt_stats" in iteration:
                    silly = "electron_cpt_stats" if "electron_cpt_stats" in iteration else "cpt_stats"
                    if iteration[silly]:
                        for cpt_id, cpt in enumerate(iteration[silly]):
                            if cpt["energy"] <= 0:
                                continue

                            if cpt_id not in electron_beams:
                                electron_beams[cpt_id] = {}

                            if cpt["energy"] not in electron_beams[cpt_id]:
                                electron_beams[cpt_id][cpt["energy"]] = {"apertures": [], "weights": []}

                            electron_beams[cpt_id][cpt["energy"]]["apertures"].append(cpt["num_apertures"])
                            electron_beams[cpt_id][cpt["energy"]]["weights"].append(cpt["weight_sum"])
                            total_weight += cpt["weight_sum"]

                if "photon_cpt_stats" in iteration and type(iteration["photon_cpt_stats"]) is list:
                    for cpt_id, cpt in enumerate(iteration["photon_cpt_stats"]):
                        if cpt["energy"] <= 0:
                            continue

                        if cpt_id not in photon_beams:
                            photon_beams[cpt_id] = {}

                        if cpt["energy"] not in photon_beams[cpt_id]:
                            photon_beams[cpt_id][cpt["energy"]] = {"apertures": [], "weights": []}

                        photon_beams[cpt_id][cpt["energy"]]["apertures"].append(cpt["num_apertures"])
                        photon_beams[cpt_id][cpt["energy"]]["weights"].append(cpt["weight_sum"])
                        total_weight += cpt["weight_sum"]

                total_weights.append(total_weight)

            self.overall_data["num_apertures"] = num_apertures
            self.overall_data["cost_functions"] = cost_functions
            self.overall_data["aperture_sizes"] = aperture_sizes
            self.overall_data["total_weights"] = total_weights
            self.overall_data["running_times"] = running_times
            self.overall_data["pruned_apertures"] = pruned_apertures

            self.modalities["electrons"]["beams"] = electron_beams
            self.modalities["photons"]["beams"] = photon_beams


    def get_apertures_per_energy(self, apertures=None):
        modalities = {"photons": {}, "electrons": {}}

        for modality in modalities.keys():
            energy_dict = {}
            for angle, energies in self.modalities[modality]["beams"].iteritems():
                for energy, iterations in energies.iteritems():
                    if energy not in energy_dict:
                        energy_dict[energy] = []

                    energy_dict[energy].append(numpy.array(iterations["apertures"]))

            for energy, ap in energy_dict.iteritems():
                modalities[modality][energy] = sum(ap)

        if apertures:
            for modality in modalities.keys():
                for energy in modalities[modality]:
                    aps = modalities[modality][energy]
                    e_dict = {}
                    for index, ap in enumerate(apertures):
                        if ap not in e_dict:
                            e_dict[ap] = []
                        e_dict[ap].append(aps[index])

                    new_aps = numpy.array([max(e_dict[ap]) for ap in sorted(e_dict.keys())])
                    modalities[modality][energy] = new_aps


        return modalities


    def get_weights_per_energy(self, apertures=None):
        modalities = {"photons": {}, "electrons": {}}

        for modality in modalities.keys():
            energy_dict = {}
            for angle, energies in self.modalities[modality]["beams"].iteritems():
                for energy, iterations in energies.iteritems():
                    if energy not in energy_dict:
                        energy_dict[energy] = []

                    energy_dict[energy].append(numpy.array(iterations["weights"]))

            for energy, weights in energy_dict.iteritems():
                print energy
                modalities[modality][energy] = sum(weights)

        if apertures:
            for modality in modalities.keys():
                for energy in modalities[modality]:
                    weights = modalities[modality][energy]
                    e_dict = {}
                    for index, ap in enumerate(apertures):
                        if ap not in e_dict:
                            e_dict[ap] = []
                        e_dict[ap].append(weights[index])

                    weights = numpy.array([max(e_dict[ap]) for ap in sorted(e_dict.keys())])
                    modalities[modality][energy] = weights

        return modalities


    def get_apertures_per_beam(self, angle):
        if angle not in self.modalities["electrons"]["beams"] and angle not in self.modalities["photons"]["beams"]:
            raise Exception("Angle not in data.")

        beam_apertures = {"photons": {}, "electrons": {}}

        # Electrons
        if angle in self.modalities["electrons"]["beams"]:
            for energy, iterations in self.modalities["electrons"]["beams"][angle].iteritems():
                beam_apertures["electrons"][energy] = numpy.array(iterations["apertures"])

        # Repeat same thing for photons
        if angle in self.modalities["photons"]["beams"]:
            for energy, iterations in self.modalities["photons"]["beams"][angle].iteritems():
                beam_apertures["photons"][energy] = numpy.array(iterations["apertures"])

        return beam_apertures

    def get_weights_per_beam(self, angle):
        if angle not in self.modalities["electrons"]["beams"] and angle not in self.modalities["photons"]["beams"]:
            raise Exception("Angle not in data.")

        beam_weights = {"photons": {}, "electrons": {}}

        # Electrons
        if angle in self.modalities["electrons"]["beams"]:
            for energy, iterations in self.modalities["electrons"]["beams"][angle].iteritems():
                beam_weights["electrons"][energy] = numpy.array(iterations["weights"])

        # Repeat same thing for photons
        if angle in self.modalities["photons"]["beams"]:
            for energy, iterations in self.modalities["photons"]["beams"][angle].iteritems():
                beam_weights["photons"][energy] = numpy.array(iterations["weights"])

        return beam_weights


class OptVisualisation():
    def __init__(self, filenames):
        if isinstance(filenames, str):
            filenames = [filenames]

        self.filenames = filenames
        self.optData = [OptOutput(filename) for filename in filenames]

        self.linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
        self.color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        self.hatches = ['/', '\\' , '|', '-', '+', 'x', 'o', 'O', '.', '*']
        self.current_max_weights = 0.0

    def print_beam_angles(self):
        for opt in self.optData:
            print "For %s:" % opt.name
            if len(opt.modalities["electrons"]["beams"].keys()) > 0:
                print "Electron angles: " + ", ".join([str(x) for x in sorted(opt.modalities["electrons"]["beams"].keys())])

            if len(opt.modalities["photons"]["beams"].keys()) > 0:
                print "Photon angles: " + ", ".join([str(x) for x in sorted(opt.modalities["photons"]["beams"].keys())])

    def plot_cost(self, x_axis="iterations"):
        for opt in self.optData:
            apertures = opt.overall_data["num_apertures"]
            cost_functions = opt.overall_data["cost_functions"]
            if x_axis == "iterations":
                x_data = numpy.arange(len(apertures)) + 1
                y_data = cost_functions
            else:
                cost_dict = {}
                for index, ap in enumerate(apertures):
                    if ap not in cost_dict:
                        cost_dict[ap] = []
                    cost_dict[ap].append(cost_functions[index])

                x_data = sorted(cost_dict.keys())
                y_data = [min(cost_dict[ap]) for ap in x_data]

            pyplot.semilogy(x_data, y_data, label=opt.name)
            #pyplot.plot(x_data, y_data, label=opt.name)

        if x_axis == "iterations":
            pyplot.xlabel("num. iterations")
        else:
            pyplot.xlabel("num. apertures")
        pyplot.ylabel("cost function value")

        pyplot.legend()
        pyplot.show()

    def plot_time(self, x_axis="iterations"):
        x_datas = []
        y_datas = []
        for opt in self.optData:
            apertures = opt.overall_data["num_apertures"]
            running_times = opt.overall_data["running_times"]
            if x_axis == "iterations":
                x_data = numpy.arange(len(apertures)) + 1
                y_data = running_times
            else:
                time_dict = {}
                for index, ap in enumerate(apertures):
                    if ap not in time_dict:
                        time_dict[ap] = []
                    time_dict[ap].append(running_times[index])

                x_data = sorted(time_dict.keys())
                y_data = [min(time_dict[ap]) for ap in x_data]

            x_datas.append(x_data)
            y_datas.append(y_data)

            #pyplot.plot(x_data, y_data, label=opt.name)

        if x_axis == "iterations":
            pyplot.xlabel("num. iterations")
        else:
            pyplot.xlabel("num. apertures")

        max_time = max([y_data[-1] for y_data in y_datas])
        print max_time

        for x_data, y_data, opt in zip(x_datas, y_datas, self.optData):
            pyplot.plot(x_data, numpy.array(y_data), label=opt.name)

        #pyplot.ylabel("relative running time")
        pyplot.ylabel("running time (s)")

        pyplot.legend(loc=2)
        pyplot.show()

    def plot_pruned(self, x_axis="iterations"):
        for opt in self.optData:
            apertures = opt.overall_data["num_apertures"]
            pruned_aps = opt.overall_data["pruned_apertures"]
            if x_axis == "iterations":
                x_data = numpy.arange(len(apertures)) + 1
                y_data = pruned_aps
            else:
                pruned_dict = {}
                for index, ap in enumerate(apertures):
                    if ap not in pruned_dict:
                        pruned_dict[ap] = []
                    pruned_dict[ap].append(pruned_aps[index])

                x_data = sorted(pruned_dict.keys())
                y_data = [min(pruned_dict[ap]) for ap in x_data]

            pyplot.plot(x_data, y_data, label=opt.name)

        if x_axis == "iterations":
            pyplot.xlabel("num. iterations")
        else:
            pyplot.xlabel("num. apertures")
        pyplot.ylabel("pruned apertures")

        pyplot.legend()
        pyplot.show()

    def plot_size(self, x_axis="iterations"):
        for opt in self.optData:
            apertures = opt.overall_data["num_apertures"]
            cost_functions = opt.overall_data["aperture_sizes"]
            if x_axis == "iterations":
                x_data = numpy.arange(len(apertures))
            else:
                x_data = apertures

            pyplot.plot(x_data, cost_functions, label=opt.name)

        if x_axis == "iterations":
            pyplot.xlabel("num. iterations")
        else:
            pyplot.xlabel("num. apertures")
        pyplot.ylabel("aperture size / mm$^2$")

        pyplot.legend()
        pyplot.show()

    def plot_total_weight(self, x_axis="iterations"):
        for opt in self.optData:
            apertures = opt.overall_data["num_apertures"]
            cost_functions = opt.overall_data["total_weights"]
            if x_axis == "iterations":
                x_data = numpy.arange(len(apertures))
            else:
                x_data = apertures

            pyplot.plot(x_data, cost_functions, label=opt.name)

        if x_axis == "iterations":
            pyplot.xlabel("num. iterations")
        else:
            pyplot.xlabel("num. apertures")
        pyplot.ylabel("total weight")

        pyplot.legend()
        pyplot.show()


    def plot_apertures(self, beam=None, relative=False, stacked=False, x_axis="iterations"):
        if stacked:
            if beam:
                self._plot_beam_energies_stacked(beam, relative, x_axis)
            else:
                self._plot_energies_stacked(relative, x_axis)
        else:
            if beam:
                self._plot_beam_energies(beam, relative, x_axis)
            else:
                self._plot_energies(relative, x_axis)

    def plot_weights(self, beam=None, relative=False, stacked=False, x_axis="iterations"):
        if stacked:
            if beam:
                self._plot_beam_weights_stacked(beam, relative, x_axis)
            else:
                self._plot_weights_stacked(relative, x_axis)
        else:
            if beam:
                self._plot_beam_weights(beam, relative, x_axis)
            else:
                self._plot_weights(relative, x_axis)


    def _plot_energies(self, relative, x_axis):
        for opt_index, opt in enumerate(self.optData):
            apertures = opt.overall_data["num_apertures"]

            max_apertures = max(apertures)

            if x_axis == "iterations":
                modalities = opt.get_apertures_per_energy()
            else:
                modalities = opt.get_apertures_per_energy(apertures)


            data = [numpy.array(modalities["electrons"][energy], dtype=float) for energy in sorted(modalities["electrons"].keys())]
            data += [numpy.array(modalities["photons"][energy], dtype=float) for energy in sorted(modalities["photons"].keys())]

            #if len(self.optData) > 1:
            #    energy_list = [opt.name + " " + str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
            #    energy_list += [opt.name + " " + str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]
            #else:
            energy_list = [str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
            energy_list += [str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]

            if relative:
                pyplot.ylim([0, 1])
                for energy_index in range(len(data)):
                    for iteration in range(len(data[energy_index])):
                        data[energy_index][iteration] /= float(sorted(list(set(apertures)))[iteration])

            if x_axis == "iterations":
                pyplot.xlabel("iteration number")
                pyplot.ylabel("apertures")

                pyplot.xlim([0, len(apertures)])
                x_data = numpy.arange(len(apertures)) + 1
            else:
                pyplot.xlabel("total apertures")
                pyplot.ylabel("apertures")
                pyplot.xlim([0, max(apertures)])
                x_data = sorted(list(set(apertures)))

            for index, energy_data in enumerate(data):
                if opt_index == 0:
                    label = energy_list[index]
                else:
                    # Disable the label
                    label = "_" + energy_list[index]

                pyplot.plot(x_data,
                            energy_data,
                            label=label,
                            color=self.color_cycle[index % len(self.color_cycle)],
                            ls=self.linestyles[opt_index % len(self.linestyles)])

        ax = pyplot.gca()
        handles, labels = ax.get_legend_handles_labels()

        # reverse the order
        ax.legend(handles[::-1], labels[::-1], loc=2)

        pyplot.show()


    def _plot_energies_stacked(self, relative, x_axis):
        opt = self.optData[0]

        apertures = opt.overall_data["num_apertures"]

        max_apertures = max(apertures)
        if x_axis == "iterations":
            modalities = opt.get_apertures_per_energy()
        else:
            modalities = opt.get_apertures_per_energy(apertures)

        data = [numpy.array(modalities["electrons"][energy], dtype=float) for energy in sorted(modalities["electrons"].keys())]
        data += [numpy.array(modalities["photons"][energy], dtype=float) for energy in sorted(modalities["photons"].keys())]

        energy_list = [str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
        energy_list += [str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]

        if relative:
            pyplot.ylim([0, 1])
            for energy_index in range(len(data)):
                for iteration in range(len(data[energy_index])):
                    data[energy_index][iteration] /= float(apertures[iteration])

        colors = [self.color_cycle[i] for i in range(len(energy_list))]
        hatches = [self.hatches[i % len(self.hatches)] for i in range(len(energy_list))]

        patches = [mpatches.Patch(color=colors[index], hatch=hatches[index]) for index, val in enumerate(energy_list)]
        if x_axis == "iterations":
            pyplot.xlabel("iteration number")
            pyplot.ylabel("apertures")

            pyplot.xlim([0, len(apertures)])
            x_data = numpy.arange(len(apertures)) + 1
        else:
            pyplot.xlabel("total apertures")
            pyplot.ylabel("apertures")
            pyplot.xlim([0, max(apertures)])
            x_data = sorted(list(set(apertures)))


        polys = pyplot.stackplot(x_data, data, colors=colors)

        #for index, poly in enumerate(polys):
        #    poly.set_hatch(hatches[index])

        pyplot.legend(patches[::-1], energy_list[::-1], loc=2)

        pyplot.show()

    def _plot_beam_energies(self, angle, relative, x_axis):
        for opt_index, opt in enumerate(self.optData):
            apertures = opt.overall_data["num_apertures"]

            max_apertures = max(apertures)

            beam_energies = opt.get_apertures_per_beam(angle)

            data = [numpy.array(beam_energies["electrons"][energy], dtype=float) for energy in sorted(beam_energies["electrons"].keys())]
            data += [numpy.array(beam_energies["photons"][energy], dtype=float) for energy in sorted(beam_energies["photons"].keys())]

            energy_list = [opt.name + ", %i$^\circ$," % angle + str(energy) + " MeV" for energy in sorted(beam_energies["electrons"].keys())]
            energy_list += [opt.name + ", %i$^\circ$," % angle + str(energy) + " MV" for energy in sorted(beam_energies["photons"].keys())]

            if relative:
                pyplot.ylim([0, 1])
                for energy_index in range(len(data)):
                    for iteration in range(len(data[energy_index])):
                        data[energy_index][iteration] /= float(apertures[iteration])

            if x_axis == "iterations":
                pyplot.xlabel("iteration number")
                pyplot.ylabel("apertures")

                pyplot.xlim([0, len(apertures)])
                x_data = numpy.arange(len(apertures)) + 1
            else:
                pyplot.xlabel("total apertures")
                pyplot.ylabel("apertures")
                pyplot.xlim([0, max(apertures)])
                x_data = apertures

            for index, energy_data in enumerate(data):
                pyplot.plot(x_data,
                            energy_data,
                            label=energy_list[index],
                            color=self.color_cycle[index % len(self.color_cycle)],
                            ls=self.linestyles[opt_index % len(self.linestyles)])

        pyplot.legend(loc=2)
        pyplot.show()


    def _plot_beam_energies_stacked(self, angle, relative, x_axis):
        opt = self.optData[0]

        apertures = opt.overall_data["num_apertures"]

        max_apertures = max(apertures)

        modalities = opt.get_apertures_per_beam(angle)

        data = [numpy.array(modalities["electrons"][energy], dtype=float) for energy in sorted(modalities["electrons"].keys())]
        data += [numpy.array(modalities["photons"][energy], dtype=float) for energy in sorted(modalities["photons"].keys())]

        energy_list = ["%i$^\circ$," % angle + str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
        energy_list += ["%i$^\circ$," % angle + str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]

        if relative:
            pyplot.ylim([0, 1])
            for energy_index in range(len(data)):
                for iteration in range(len(data[energy_index])):
                    data[energy_index][iteration] /= float(apertures[iteration])

        colors = [self.color_cycle[i] for i in range(len(energy_list))]

        patches = [mpatches.Patch(color=colors[index]) for index, val in enumerate(energy_list)]

        if x_axis == "iterations":
            pyplot.xlabel("iteration number")
            pyplot.ylabel("apertures")

            pyplot.xlim([0, len(apertures)])
            x_data = numpy.arange(len(apertures)) + 1
        else:
            pyplot.xlabel("total apertures")
            pyplot.ylabel("apertures")
            pyplot.xlim([0, max(apertures)])
            x_data = apertures

        pyplot.stackplot(x_data, data, colors=colors)

        pyplot.legend(patches[::-1], energy_list[::-1], loc=2)

        pyplot.show()


    def _plot_weights(self, relative, x_axis):
        for opt_index, opt in enumerate(self.optData):
            apertures = opt.overall_data["num_apertures"]

            max_apertures = max(apertures)

            if x_axis == "iterations":
                modalities = opt.get_weights_per_energy()
            else:
                modalities = opt.get_weights_per_energy(apertures)

            data = [modalities["electrons"][energy] for energy in sorted(modalities["electrons"].keys())]
            data += [modalities["photons"][energy] for energy in sorted(modalities["photons"].keys())]

            max_weights = max([max(d) for d in data])
            if max_weights > self.current_max_weights:
                self.current_max_weights = max_weights
            total_weights = sum(data)

            #if len(self.optData) > 1:
            #    energy_list = [opt.name + " " + str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
            #    energy_list += [opt.name + " " + str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]
            #else:
            energy_list = [str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
            energy_list += [str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]

            if relative:
                pyplot.ylim([0, 1.05])
                for energy_index in range(len(data)):
                    for iteration in range(len(data[energy_index])):
                        data[energy_index][iteration] /= float(total_weights[iteration])
            else:
                pyplot.ylim([0, self.current_max_weights + 0.1 * self.current_max_weights])

            if x_axis == "iterations":
                pyplot.xlabel("iteration number")
                pyplot.ylabel("apertures")

                pyplot.xlim([0, len(apertures)])
                x_data = numpy.arange(len(apertures)) + 1
            else:
                pyplot.xlabel("total apertures")
                pyplot.xlim([0, max(apertures)])
                x_data = sorted(list(set(apertures)))

            pyplot.ylabel("weight")

            for index, energy_data in enumerate(data):
                if opt_index == 0:
                    label = energy_list[index]
                else:
                    label = "_" + energy_list[index]
                pyplot.plot(x_data,
                            energy_data,
                            label=label,
                            color=self.color_cycle[index % len(self.color_cycle)],
                            ls=self.linestyles[opt_index % len(self.linestyles)])

        ax = pyplot.gca()
        handles, labels = ax.get_legend_handles_labels()

        # reverse the order
        ax.legend(handles[::-1], labels[::-1], loc=2)

        pyplot.legend(loc=2)
        pyplot.show()


    def _plot_weights_stacked(self, relative, x_axis):
        opt = self.optData[0]

        apertures = opt.overall_data["num_apertures"]

        max_apertures = max(apertures)

        if x_axis == "iterations":
            modalities = opt.get_weights_per_energy()
        else:
            modalities = opt.get_weights_per_energy(apertures)

        data = [modalities["electrons"][energy] for energy in sorted(modalities["electrons"].keys())]
        data += [modalities["photons"][energy] for energy in sorted(modalities["photons"].keys())]

        data_length = [len(x) for x in data]
        print data_length

        energy_list = [str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
        energy_list += [str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]

        print energy_list
        max_weights = [max(d) for d in data]

        total_weights = sum(data)
        max_total_weight = max(total_weights)

        if relative:
            pyplot.ylim([0, 1])
            for energy_index in range(len(data)):
                for iteration in range(len(data[energy_index])):
                    data[energy_index][iteration] /= float(total_weights[iteration])
        else:
            pyplot.ylim([0, max_total_weight])

        colors = [self.color_cycle[i] for i in range(len(energy_list))]

        patches = [mpatches.Patch(color=colors[index]) for index, val in enumerate(energy_list)]

        if x_axis == "iterations":
            pyplot.xlabel("iteration number")
            pyplot.xlim([0, len(apertures)])
            x_data = numpy.arange(len(apertures)) + 1
        else:
            pyplot.xlabel("total apertures")
            pyplot.xlim([0, max(apertures)])
            x_data = sorted(list(set(apertures)))
        pyplot.ylabel("weight")

        pyplot.stackplot(x_data, data, colors=colors)

        pyplot.legend(patches[::-1], energy_list[::-1], loc=2)
        pyplot.show()

    def _plot_beam_weights(self, angle, relative, x_axis):
        for opt_index, opt in enumerate(self.optData):
            apertures = opt.overall_data["num_apertures"]

            max_apertures = max(apertures)

            modalities = opt.get_weights_per_beam(angle)

            data = [modalities["electrons"][energy] for energy in sorted(modalities["electrons"].keys())]
            data += [modalities["photons"][energy] for energy in sorted(modalities["photons"].keys())]

            max_weights = max([max(d) for d in data])

            total_weights = sum(data)

            energy_list = [opt.name + ", %i$^\circ$ " % angle + str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
            energy_list += [opt.name + ", %i$^\circ$ " % angle + str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]

            if relative:
                pyplot.ylim([0, 1])
                for energy_index in range(len(data)):
                    for iteration in range(len(data[energy_index])):
                        data[energy_index][iteration] /= float(total_weights[iteration])
            else:
                pyplot.ylim([0, max_weights + 0.1 * max_weights])

            if x_axis == "iterations":
                pyplot.xlabel("iteration number")
                pyplot.xlim([0, len(apertures)])
                x_data = numpy.arange(len(apertures)) + 1
            else:
                pyplot.xlabel("total apertures")
                pyplot.xlim([0, max(apertures)])
                x_data = apertures
            pyplot.ylabel("weight")

            for index, energy_data in enumerate(data):
                pyplot.plot(x_data,
                            energy_data,
                            label=energy_list[index],
                            color=self.color_cycle[index % len(self.color_cycle)],
                            ls=self.linestyles[opt_index % len(self.linestyles)])

        pyplot.legend(loc=2)
        pyplot.show()

    def _plot_beam_weights_stacked(self, angle, relative, x_axis):
        opt = self.optData[0]

        apertures = opt.overall_data["num_apertures"]

        max_apertures = max(apertures)

        modalities = opt.get_weights_per_beam(angle)

        data = [modalities["electrons"][energy] for energy in sorted(modalities["electrons"].keys())]
        data += [modalities["photons"][energy] for energy in sorted(modalities["photons"].keys())]

        energy_list = ["%i$^\circ$," % angle + str(energy) + " MeV" for energy in sorted(modalities["electrons"].keys())]
        energy_list += ["%i$^\circ$," % angle + str(energy) + " MV" for energy in sorted(modalities["photons"].keys())]

        max_weights = [max(d) for d in data]

        total_weights = sum(data)
        max_total_weight = max(total_weights)

        if relative:
            pyplot.ylim([0, 1])
            for energy_index in range(len(data)):
                for iteration in range(len(data[energy_index])):
                    data[energy_index][iteration] /= float(total_weights[iteration])
        else:
            pyplot.ylim([0, max_total_weight])

        colors = [self.color_cycle[i] for i in range(len(energy_list))]

        patches = [mpatches.Patch(color=colors[index]) for index, val in enumerate(energy_list)]

        if x_axis == "iterations":
            pyplot.xlabel("iteration number")

            pyplot.xlim([0, len(apertures)])
            x_data = numpy.arange(len(apertures)) + 1
        else:
            pyplot.xlabel("total apertures")
            pyplot.xlim([0, max(apertures)])
            x_data = apertures

        pyplot.ylabel("weight")
        pyplot.stackplot(x_data, data, colors=colors)

        pyplot.legend(patches[::-1], energy_list[::-1], loc=2)
        pyplot.show()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot data about optimisation run')
    parser.add_argument('filenames', help='The data filename for optimisation', type=str, nargs='+')
    parser.add_argument('-c', '--cost', action="store_true", help='plot cost function as a function of iteration number')
    parser.add_argument('-p', '--pruned', action="store_true", help='plot pruned apertures as a function of iteration number')
    parser.add_argument('-i', '--size', action="store_true", help='plot aperture size as a function of iteration number')
    parser.add_argument('-tt', '--time', action="store_true", help='plot running time as a function of iteration number')
    parser.add_argument('-ww', '--total_weight', action="store_true", help='plot total weights as a function of iteration number')
    parser.add_argument('-a', '--apertures', action="store_true", help='plot apertures as a function of iteration number')
    parser.add_argument('-w', '--weights', action="store_true", help='plot weights as a function of iteration number')
    parser.add_argument('-r', '--relative', action="store_true", help='relative plot')
    parser.add_argument('-s', '--stacked', action="store_true", help='stacked area plot, only one dataset will be plotted')
    parser.add_argument('-t', '--totalap', action="store_true", help='plot num. apertures on x axis instead of iteration number')

    parser.add_argument('-b', '--beam', type=int, help='specify a beam angle to plot')

    args = parser.parse_args()

    opt_obj = OptVisualisation(args.filenames)

    if args.totalap:
        x_axis = "apertures"
    else:
        x_axis = "iterations"

    if not args.apertures and not args.weights and not args.cost and not args.size and not args.total_weight and not args.pruned and not args.time:
        opt_obj.print_beam_angles()
        exit()

    if args.stacked and len(args.filenames) > 1:
        exit("Error: only one dataset can be plotted in a stacked area plot.")

    if args.cost:
        opt_obj.plot_cost(x_axis=x_axis)

    if args.size:
        opt_obj.plot_size(x_axis=x_axis)

    if args.pruned:
        opt_obj.plot_pruned(x_axis=x_axis)

    if args.time:
        opt_obj.plot_time(x_axis=x_axis)

    if args.total_weight:
        opt_obj.plot_total_weight(x_axis=x_axis)

    if args.apertures:
        opt_obj.plot_apertures(beam=args.beam, relative=args.relative, stacked=args.stacked, x_axis=x_axis)

    if args.weights:
        opt_obj.plot_weights(beam=args.beam, relative=args.relative, stacked=args.stacked, x_axis=x_axis)


