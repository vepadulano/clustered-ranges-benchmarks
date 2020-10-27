import csv
from collections import namedtuple
from time import time

import ROOT

import numpy


Range = namedtuple(
    "Range", ["start", "end", "filelist", "friend_info"], defaults=[None, None])


def get_clusters_PR(treename, filelist):

    clusters = []
    cluster = namedtuple(
        "cluster", ["start", "end", "offset", "filetuple"])
    fileandindex = namedtuple("fileandindex", ["filename", "fileindex"])
    offset = 0
    fileindex = 0

    for filename in filelist:
        f = ROOT.TFile.Open(str(filename))
        t = f.Get(treename)

        entries = t.GetEntriesFast()
        it = t.GetClusterIterator(0)
        start = it()
        end = 0

        while start < entries:
            end = it()
            clusters.append(cluster(start + offset, end + offset, offset,
                                    fileandindex(filename, fileindex)))
            start = end

        fileindex += 1
        offset += entries

    return clusters


def get_clustered_ranges_Enric(treename, filelist, npartitions,
                               friend_info=None):
    clusters = get_clusters_PR(treename, filelist)
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
            filename, fileindex = clusters[idx][3]
            if (range_files and clusters[idx - 1][3].fileindex == fileindex):
                continue
            range_files.append(filename)  # here we add the file to range_files

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

    times_enric_pr = []
    for _ in range(1000):
        start = time()
        get_clustered_ranges_Enric(treename, filelist, npartitions)
        end = time() - start
        times_enric_pr.append(end)

    with open("enric_pr.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows([[element] for element in times_enric_pr])

    with open("clustered_ranges_benchmark.txt", "+a") as f:
        f.write("\nEnric PR clustered ranges:\n{}\n"
                .format("\n\n".join(
                    map(str,
                        get_clustered_ranges_Enric(treename,
                                                   filelist,
                                                   npartitions)
                        )))
                )
        f.write("\nEnric PR benchmark over 1000 executions: {} +- {} seconds\n"
                .format(
                    round(numpy.mean(times_enric_pr), 4),
                    round(numpy.std(times_enric_pr), 4)
                )
                )
