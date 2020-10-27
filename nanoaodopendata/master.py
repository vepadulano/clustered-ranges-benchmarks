import csv
from collections import namedtuple
from time import time


import ROOT

import numpy


Range = namedtuple(
    "Range", ["start", "end", "filelist", "friend_info"], defaults=[None, None])


def get_clusters_master(treename, filelist):

    clusters = []
    offset = 0

    for filename in filelist:
        f = ROOT.TFile.Open(str(filename))
        t = f.Get(treename)

        entries = t.GetEntriesFast()
        it = t.GetClusterIterator(0)
        start = it()
        end = 0

        while start < entries:
            end = it()
            cluster = (start + offset, end + offset, offset, filename)
            clusters.append(cluster)
            start = end

        offset += entries

    return clusters


def get_clustered_ranges_master(treename, filelist, npartitions,
                                friend_info=None):
    clusters = get_clusters_master(treename, filelist)
    numclusters = len(clusters)

    partSize = numclusters // npartitions
    remainder = numclusters % npartitions

    i = 0  # Iterator
    ranges = []
    entries_to_process = 0

    while i < numclusters:
        index_start = i
        start = clusters[i][0]
        i = i + partSize
        if remainder > 0:
            i += 1
            remainder -= 1
        index_end = i
        if i == numclusters:
            end = clusters[-1][1]
        else:
            end = clusters[i - 1][1]

        range_files = []
        for idx in range(index_start, index_end):
            current_file = clusters[idx][3]
            if range_files and range_files[-1] == current_file:
                continue
            range_files.append(clusters[idx][3])

        offset_first_cluster = clusters[index_start][2]
        ranges.append(Range(start - offset_first_cluster,
                            end - offset_first_cluster,
                            range_files,
                            friend_info))
        entries_to_process += (end - start)

    return ranges


if __name__ == "__main__":
    treename = "Events"
    filelist = [
        "root://eospublic.cern.ch//eos/opendata/cms/derived-data/"
        "AOD2NanoAODOutreachTool/{}.root".format(filename)
        for filename in ["DYJetsToLL",
                         "GluGluToHToTauTau",
                         "Run2012BC_DoubleMuParked_Muons",
                         "Run2012B_TauPlusX",
                         "Run2012C_TauPlusX",
                         "TTbar",
                         "VBF_HToTauTau",
                         "W1JetsToLNu",
                         "W2JetsToLNu",
                         "W3JetsToLNu",
                         ]
    ]

    npartitions = 5

    times_master = []
    for _ in range(1000):
        start = time()
        get_clustered_ranges_master(treename, filelist, npartitions)
        end = time() - start
        times_master.append(end)

    with open("master.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows([[element] for element in times_master])

    with open("clustered_ranges_benchmark.txt", "+a") as f:
        f.write("\nMaster clustered ranges:\n{}\n"
                .format("\n\n".join(
                    map(str,
                        get_clustered_ranges_master(treename,
                                                    filelist,
                                                    npartitions)
                        )))
                )
        f.write("\nMaster benchmark over 1000 executions: {} +- {} seconds\n"
                .format(
                    round(numpy.mean(times_master), 4),
                    round(numpy.std(times_master), 4)
                )
                )
