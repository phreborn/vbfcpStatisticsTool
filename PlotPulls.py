#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import yaml
import sys
import os
from collections import OrderedDict


def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True, help="Steering file to configure everything")
    p.add_argument("--top", type=str, help="Plot top X parameter")
    p.add_argument("--relative", action='store_true', help="Plot relative impact")
    p.add_argument("--poi", type=str, default="mu_XS_ttH", help="define the POI to plot")
    args = p.parse_args()
    with open(args.config, 'r') as inf:
        config = ordered_load(inf, yaml.SafeLoader)
    return config, args

if __name__ == '__main__':
    config, args = parse_args(sys.argv[1:])
    for ws, attribute in config["workspaces"].iteritems():
        for type, attribute_list in attribute.iteritems():
            if not isinstance(attribute_list, list):
                print "attribute list needs to be a list to plot different POIs for the same output"
                raise SystemExit
            for this_attribute in attribute_list:
                cmd = "./bin/plot_pulls.exe"
                this_input = os.path.join("root-files", ws + "_" + type, "pulls/")
                cmd += " --input " + this_input
                if "poi" in this_attribute.keys():
                    cmd += " --poi " + str(this_attribute["poi"])
                    args.poi = this_attribute["poi"]
                else:
                    cmd += " --poi " + args.poi
                cmd += " --scale_poi " + str(this_attribute["scale_poi"])
                cmd += " --scale_theta " + str(this_attribute["scale_theta"])
                if "map" in this_attribute.keys():
                    cmd += " --map " + str(this_attribute["map"])

                if args.relative:
                    cmd += " --relative true"

                if "legend" in this_attribute.keys():
                    cmd += " --label " + str("\"" + this_attribute["legend"] + "\"")
                if args.top is not None:
                    cmd += " --top " + str(args.top)
                print cmd
                os.system(cmd)
                if args.top is not None:
                    if int(args.top) < 100:
                        outputfile = os.path.join("pdf-files", "ranking_" + args.poi + "_rank_0001_to_00" + args.top + ".pdf")
                    else:
                        outputfile = os.path.join("pdf-files", "ranking_" + args.poi + "_rank_0001_to_0" + args.top + ".pdf")
                    new_outputfile = os.path.join("pdf-files", "ranking_" + args.poi + "_" + ws + "_" + type + "_top_" + args.top + ".pdf")
                else:
                    outputfile = os.path.join("pdf-files", "ranking_" + args.poi + "_rank_0001_to_00-1.pdf")
                    new_outputfile = os.path.join("pdf-files", "ranking_" + args.poi + "_" + ws + "_" + type + ".pdf")

                if os.path.exists(outputfile):
                    #  print "Move " + outputfile + " --> " + new_outputfile
                    #  out_dir = os.path.join("pdf-files", ws)
                    #  if not os.path.exists(out_dir):
                    #      os.makedirs(out_dir)
                    cmd = "mv " + outputfile + " " + new_outputfile
                    print cmd
                    os.system(cmd)
                    new_cmd = cmd.replace("pdf-files", "png-files")
                    new_cmd = new_cmd.replace(".pdf", ".png")
                    print new_cmd
                    os.system(new_cmd)
