import argparse
import os


def get_options(args=None):

    parser = argparse.ArgumentParser(description='creating nodes')

    # data
    parser.add_argument('--num_beams', default=300, help='number of beams')
    parser.add_argument('--type', default='all', help='type of the nodes, "all" points or current "ck" path')

    # plot
    parser.add_argument('--plot_graph', default=False, help='whether to plot the graph or not')

    opts = parser.parse_args(args)

    return opts
