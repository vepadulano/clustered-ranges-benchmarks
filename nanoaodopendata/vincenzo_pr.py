import csv
from collections import namedtuple
from time import time

import ROOT

import numpy


Range = namedtuple(
    "Range", ["start", "end", "filelist", "friend_info"], defaults=[None, None])


def _n_even_chunks(iterable, n_chunks):
    """Yield `n_chunks` as even chunks as possible from `iterable`."""
    last = 0
    for i in range(1, n_chunks + 1):
        cur = int(round(i * (len(iterable) / n_chunks)))
        yield iterable[last:cur]
        last = cur


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


def get_clustered_ranges_PR(treename, filelist, npartitions,
                            friend_info=None):
    filesandclusters = get_clusters_PR(treename, filelist)

    clusters_split_by_partition = list(
        _n_even_chunks(filesandclusters, npartitions))

    partitions_startentries = [min(clusters)[0]
                               for clusters
                               in clusters_split_by_partition]
    partitions_endentries = [max(clusters)[1]
                             for clusters in clusters_split_by_partition]
    partitions_offset = [clusters[0].offset
                         for clusters in clusters_split_by_partition]

    partitions_filelist = [
        [
            filetuple.filename
            for filetuple in sorted(set([  # set to take unique file indexes
                cluster.filetuple for cluster in clusters
            ]))
        ]
        for clusters in clusters_split_by_partition
    ]

    clustered_ranges = [
        Range(start - offset, end - offset, range_files, friend_info)
        for start, end, offset, range_files
        in zip(partitions_startentries, partitions_endentries,
               partitions_offset, partitions_filelist)
    ]

    return clustered_ranges


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

    times_vincenzo_pr = []
    for _ in range(1000):
        start = time()
        get_clustered_ranges_PR(treename, filelist, npartitions)
        end = time() - start
        times_vincenzo_pr.append(end)

    with open("vincenzo_pr.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows([[element] for element in times_vincenzo_pr])

    with open("clustered_ranges_benchmark.txt", "+a") as f:
        f.write("\nVincenzo PR clustered ranges:\n{}\n"
                .format("\n\n".join(
                    map(str,
                        get_clustered_ranges_PR(treename,
                                                filelist,
                                                npartitions)
                        )))
                )
        f.write("\nVincenzo PR benchmark over 1000 executions: {} +- {} seconds\n"
                .format(
                    round(numpy.mean(times_vincenzo_pr), 4),
                    round(numpy.std(times_vincenzo_pr), 4)
                )
                )
