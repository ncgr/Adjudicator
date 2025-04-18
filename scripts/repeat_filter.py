import sys
import re


#awk -F"\t" '{if(l[$9]){c[$9] += $NF;r[$9]=$0}else{l[$9] = $5 - $4;c[$9] += $NF;r[$9]=$0;d[$9] = $1"\t"$4"\t"$5}}END{for (a in c){if(c[a]/l[a] <= 10){print a"\t"c[a]/l[a]"\t"d[a]}}}' all_repeats_exons.gff3

def filter_models(wao, coverage, black_list):
    """Filters wao overlap between exons and repeats"""
    get_id = r"ID=(.+?);"
    get_parent = r"Parent=(.+?);"
    overlaps = {}
    with open(wao) as w:
        for line in w:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            exon_ids = re.findall(get_id, f"{fields[8]};", re.IGNORECASE)
            parent_ids = re.findall(get_parent, f"{fields[8]};", re.IGNORECASE)
            exon_id = exon_ids[0]
            parent_id = parent_ids[0]
            if parent_id not in overlaps:
                overlaps[parent_id] = {"exons": {}, "total_length": 0, "total_coverage": 0}
#            else:
#                print(f"PROOOF: {overlaps[parent_id]} {exon_id}")
            if exon_id not in overlaps[parent_id]["exons"]:
                overlaps[parent_id]["exons"][exon_id] = {"length": int(fields[4]) - int(fields[3]) + 1, 
                                                          "overlap": int(fields[-1])
                                                        }
                overlaps[parent_id]["total_length"] += int(fields[4]) - int(fields[3]) + 1
            else:
                overlaps[parent_id]["exons"][exon_id]["overlap"] += int(fields[-1])
#            print(exon_id)
#            print(parent_id)
    for parent in overlaps:
        parent_coverage = 0
        exon_sum = overlaps[parent]["total_length"]
        for exon in overlaps[parent]["exons"]:
            parent_coverage += overlaps[parent]["exons"][exon]["overlap"]

#            print(f"{parent}\t{exon}\t{length}\t{coverage}\t{float(coverage/length)}\t{exon_sum}")
        if black_list:
            if float(parent_coverage/exon_sum) > coverage:
                print(f"{parent}\t{parent_coverage}\t{exon_sum}")
        else:
            if float(parent_coverage/exon_sum) <= coverage:
                print(f"{parent}\t{parent_coverage}\t{exon_sum}")


if __name__ == '__main__':
    intersect_wao = sys.argv[1]  # get bed file from first argument put click in here later
    maximum_repeat_coverage = float(sys.argv[2])  # minimum feature overlap
    black_list = False
    if len(sys.argv) == 4:
        black_list = True
    filter_models(intersect_wao, maximum_repeat_coverage, black_list) 
