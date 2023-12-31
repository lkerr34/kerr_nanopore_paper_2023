#! /usr/bin/env python3

import sys
import csv
import argparse
import gzip

class SiteStats:
    def __init__(self, g_size, g_seq):
        self.num_reads = 0
        self.called_sites = 0
        self.called_sites_methylated = 0
        self.uncalled_sites = 0
        self.group_size = g_size
        self.sequence = g_seq
        
def update_call_stats(key, num_cpg_sites, is_methylated, is_uncalled, sequence):
    if key not in sites:
        sites[key] = SiteStats(num_cpg_sites, sequence)

    sites[key].num_reads += 1
    if is_uncalled > 0:
        sites[key].uncalled_sites += num_cpg_sites
    else:
        sites[key].called_sites += num_cpg_sites
        if is_methylated > 0:
            sites[key].called_sites_methylated += num_cpg_sites  

parser = argparse.ArgumentParser( description='Calculate methylation frequency at genomic CpG sites')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2.0)
parser.add_argument('-s', '--split-groups', action='store_true')
args, input_files = parser.parse_known_args()
assert(args.call_threshold is not None)

sites = dict()
# iterate over input files and collect per-site stats
for f in input_files:
    if f[-3:] == ".gz":
        in_fh = gzip.open(f, 'rt')
    else:
        in_fh = open(f)
    csv_reader = csv.DictReader(in_fh, delimiter='\t')
    for record in csv_reader:

        num_sites = int(record['num_motifs'])
        llr = float(record['log_lik_ratio'])

        # Save ambiguous CpGs as uncalled
        if abs(llr) < args.call_threshold:
            is_uncalled = 1
            is_methylated = 0
            # if this is a multi-cpg group and split-groups is set, break up these sites
            if args.split_groups and num_sites > 1:
                c = str(record['chromosome'])
                s = int(record['start'])
                e = int(record['end'])
            
                # find the position of the first CG dinucleotide
                sequence = record['sequence']
                cg_pos = sequence.find("CG")
                first_cg_pos = cg_pos
                while cg_pos != -1:
                    key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                    update_call_stats(key, 1, is_methylated, is_uncalled, "split-group")
                    cg_pos = sequence.find("CG", cg_pos + 1)
            else:
                sequence = record['sequence']
                key = (str(record['chromosome']), int(record['start']), int(record['end']))
                update_call_stats(key, num_sites, is_methylated, is_uncalled, sequence)
        
        else:
        
            sequence = record['sequence']

            is_methylated = llr > 0
            is_uncalled = 0

            # if this is a multi-cpg group and split_groups is set, break up these sites
            if args.split_groups and num_sites > 1:
                c = str(record['chromosome'])
                s = int(record['start'])
                e = int(record['end'])

                # find the position of the first CG dinucleotide
                sequence = record['sequence']
                cg_pos = sequence.find("CG")
                first_cg_pos = cg_pos
                while cg_pos != -1:
                    key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                    update_call_stats(key, 1, is_methylated, is_uncalled, "split-group")
                    cg_pos = sequence.find("CG", cg_pos + 1)
            else:
                key = (str(record['chromosome']), int(record['start']), int(record['end']))
                update_call_stats(key, num_sites, is_methylated, is_uncalled, sequence)

# header
print("\t".join(["chromosome", "start", "end", "num_motifs_in_group", "num_sites", "called_sites_methylated", "methylated_frequency", "num_uncalled", "total_coverage", "group_sequence"]))

sorted_keys = sorted(list(sites.keys()), key = lambda x: x)

for key in sorted_keys:
    (c, s, e) = key
    T = float(sites[key].called_sites) + float(sites[key].uncalled_sites)
    if sites[key].called_sites > 0:
        f = float(sites[key].called_sites_methylated) / float(sites[key].called_sites)
        print("%s\t%s\t%s\t%d\t%d\t%d\t%.3f\t%d\t%d\t%s" % (c, s, e, sites[key].group_size, sites[key].called_sites, sites[key].called_sites_methylated, f, sites[key].uncalled_sites, T, sites[key].sequence))
    else:
        f = "NA"
        print("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%s" % (c, s, e, sites[key].group_size, sites[key].called_sites, sites[key].called_sites_methylated, f, sites[key].uncalled_sites, T, sites[key].sequence))
       
