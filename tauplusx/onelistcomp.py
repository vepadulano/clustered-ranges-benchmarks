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


def get_clustered_ranges_onelistcomp(treename, filelist, npartitions,
                                     friend_info=None):
    filesandclusters = get_clusters_PR(treename, filelist)

    clusters_split_by_partition = list(
        _n_even_chunks(filesandclusters, npartitions))

    clustered_ranges = [
        Range(
            min(clusters)[0] - clusters[0].offset,
            max(clusters)[1] - clusters[0].offset,
            [
                filetuple.filename
                for filetuple in sorted(set([  # set to take unique file indexes
                    cluster.filetuple for cluster in clusters
                ]))
            ],
            friend_info
        )
        for clusters in clusters_split_by_partition
    ]

    return clustered_ranges


if __name__ == "__main__":
    treename = "Events"
    filelist = [
        "root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012B_TauPlusX.root",
        "root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012C_TauPlusX.root"
    ]
    npartitions = 2

    times_onelistcomp = []
    for _ in range(1000):
        start = time()
        get_clustered_ranges_onelistcomp(treename, filelist, npartitions)
        end = time() - start
        times_onelistcomp.append(end)

    with open("onelistcomp.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows([[element] for element in times_onelistcomp])

    with open("clustered_ranges_benchmark.txt", "+a") as f:
        f.write("\nOneListComp clustered ranges:\n{}\n"
                .format("\n\n".join(
                    map(str,
                        get_clustered_ranges_onelistcomp(treename,
                                                         filelist,
                                                         npartitions)
                        )))
                )
        f.write("\nOneListComp benchmark over 1000 executions: {} +- {} seconds\n"
                .format(
                    round(numpy.mean(times_onelistcomp), 4),
                    round(numpy.std(times_onelistcomp), 4)
                )
                )
