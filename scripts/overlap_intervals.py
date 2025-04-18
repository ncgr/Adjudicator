import os
import argparse


class ProcessWAO:
    def __init__(self, wao):
        self.wao = os.path.abspath(wao)
        self.trees = {}
        self.query_trees = {}
        self.found = {}
        self.found_query = {}
        self.lost = {}

    def read_bed_file(self):
        with open(self.wao, "r") as file:
            for line in file:
                fields = line.strip().split()
                ref = fields[0]
                query = fields[9]
                if ref not in self.trees:
                    self.trees[ref] = intervaltree.IntervalTree()
                start = int(fields[1])
                end = int(fields[2])
                if len(fields) > 3:
                    query = fields[-1]
                    if query not in self.query_trees[ref]:
                        self.query_trees[ref][query] = intervaltree.IntervalTree()
                    self.query_trees[ref][query][start:end] = (start, end)
                self.trees[ref][start:end] = (start, end)

    def check_intervals(self):
        with open(self.filepath2, "r") as file:
            for line in file:
                fields = line.strip().split()
                ref = fields[0]
                if ref not in self.found:
                    self.found[ref] = []
                if ref not in self.trees:
                    self.found[ref] = None
                    continue
                start = int(fields[1])
                end = int(fields[2])
                gene_number = int(fields[-1])
                overlap_intervals = self.trees[ref].overlap(start, end)
                if overlap_intervals:
                    self.found[ref].append(gene_number)
                else:
                    if ref not in self.lost:
                        self.lost[ref] = []
                    self.lost[ref].append(gene_number)
                if self.query_trees[ref]:
                    if ref not in self.found_query:
                        self.found_query[ref] = {}
                    for query in self.query_trees[ref]:
                        overlap_intervals_query = self.query_trees[ref][query].overlap(start, end)
                        if query not in self.found_query[ref]:
                            self.found_query[ref][query] = []
                        if overlap_intervals_query:
                            self.found_query[ref][query].append(gene_number)

    def make_intervals(self, sorted_list):
        intervals = []
        current_interval = [sorted_list[0]]
        for i in range(1, len(sorted_list)):
            if sorted_list[i] == current_interval[-1] + 1:
                current_interval.append(sorted_list[i])
            else:
                intervals.append((current_interval[0], current_interval[-1]))
                current_interval = [sorted_list[i]]
        if current_interval:
            intervals.append((current_interval[0], current_interval[-1]))
        return intervals

    def print_found_genes(self, output_file):
        with open(output_file, 'w') as found:
            for ref in self.found:
                if self.found[ref]:
                    intervals = '\t'.join(map(str, self.make_intervals(self.found[ref])))
                    found.write(f"{ref}\t{intervals}\n")
    
        with open(f"{output_file}.byquery.tsv", 'w') as found:
            for ref in self.found_query:
                if self.found_query[ref]:
                    for query in self.found_query[ref]:
                        if self.found_query[ref][query]:
                            intervals = '\t'.join(map(str, self.make_intervals(self.found_query[ref][query])))
                            found.write(f"{ref}\t{query}\t{intervals}\n")

        
    def print_lost_genes(self, output_file):
        with open(output_file, 'w') as lost:
            for ref in self.lost:
                if self.lost[ref]:
                    intervals = '\t'.join(map(str, self.make_intervals(self.lost[ref])))
                    lost.write(f"{ref}\t{intervals}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process the output of wao for overlap.")
    parser.add_argument("wao", type=str, help="wao output from bedtools")
    parser.add_argument("out_prefix", type=str, help="the output file prefix")

    args = parser.parse_args()

    wao = args.wao
    prefix = args.out_prefix
    ProcessWAO(wao)
