import sys
import re


#awk -F"\t" '{if(l[$9]){c[$9] += $NF;r[$9]=$0}else{l[$9] = $5 - $4;c[$9] += $NF;r[$9]=$0;d[$9] = $1"\t"$4"\t"$5}}END{for (a in c){if(c[a]/l[a] <= 10){print a"\t"c[a]/l[a]"\t"d[a]}}}' all_repeats_exons.gff3

def sum_features_per_transcript(gff3):
    """Sums CDS features for each transcript in gff3 file."""
    sums = {}
    get_parent = r"Parent=(.+?);"
    with open(gff3) as g:
        for line in g:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] == "CDS":
                transcript = re.findall(get_parent, f"{fields[8]};", re.IGNORECASE)[0]
                if transcript not in sums:
                    sums[transcript] = 0
                sums[transcript] += int(fields[4]) - int(fields[3])
    return sums



def filter_models(wao, coverage, black_list, sums):
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
            if fields[8] == fields[17]:  # the feature compares itself
                continue
            exon_ids = re.findall(get_id, f"{fields[8]};", re.IGNORECASE)
            parent_ids = re.findall(get_parent, f"{fields[8]};", re.IGNORECASE)
            exon_id_left = exon_ids[0]
            parent_id_left = parent_ids[0]
            exon_ids = re.findall(get_id, f"{fields[17]};", re.IGNORECASE)
            parent_ids = re.findall(get_parent, f"{fields[17]};", re.IGNORECASE)
            exon_id_right = exon_ids[0]
            parent_id_right = parent_ids[0]
            parent_combo = f"{parent_id_left} {parent_id_right}"
            parent_combo_reverse = f"{parent_id_right} {parent_id_left}"
            exon_combo = f"{exon_id_left} {exon_id_right}"
            exon_combo_reverse = f"{exon_id_right} {exon_id_left}"
            if parent_combo_reverse in overlaps:  # don't care because I've seen the combo
                continue
            if parent_combo not in overlaps:
                overlaps[parent_combo] = {"exons": {}, "total_length_left": 0, "total_length_right": 0, "total_coverage": 0, "seen": {}, "all_cds_left": sums[parent_id_left], "all_cds_right": sums[parent_id_right]}
#            if parent_id not in overlaps:
#                overlaps[parent_id] = {"exons": {}, "total_length": 0, "total_coverage": 0}
#            else:
#                print(f"PROOOF: {overlaps[parent_id]} {exon_id}")
            if exon_combo_reverse in overlaps[parent_combo]["exons"]:  # don't care because I've seen the combo
                continue
            if exon_combo not in overlaps[parent_combo]["exons"]:
                overlaps[parent_combo]["exons"][exon_combo] = {"length_left": int(fields[4]) - int(fields[3]),
                                                               "length_right": int(fields[13]) - int(fields[12]), 
                                                               "overlap": int(fields[-1]) - 1
                                                              }
                if exon_id_left not in overlaps[parent_combo]["seen"]:
                    overlaps[parent_combo]["total_length_left"] += int(fields[4]) - int(fields[3])
                    overlaps[parent_combo]["seen"][exon_id_left] = 1
                if exon_id_right not in overlaps[parent_combo]["seen"]:
                    overlaps[parent_combo]["total_length_right"] += int(fields[13]) - int(fields[12])
                    overlaps[parent_combo]["seen"][exon_id_right] = 1
#            else:
#                overlaps[parent_combo]["exons"][exon_combo]["overlap"] += int(fields[-1])
#            print(exon_id)
#            print(parent_id)
    print(f"#parent_combo\tcombo_coverage\tsum_left\tsum_right\tleft_transcript_cds_sum\tright_transcript_cds_sum")
    for parent in overlaps:
        parent_coverage = 0
        exon_sum_left = overlaps[parent]["total_length_left"]
        exon_sum_right = overlaps[parent]["total_length_right"]
        all_cds_left = overlaps[parent]["all_cds_left"]
        all_cds_right = overlaps[parent]["all_cds_right"]
        for exon in overlaps[parent]["exons"]:
            parent_coverage += overlaps[parent]["exons"][exon]["overlap"]
            #print(parent, exon, exon_sum_left, exon_sum_right)
#            print(f"{parent}\t{exon}\t{length}\t{coverage}\t{float(coverage/length)}\t{exon_sum}")
        print(f"{parent}\t{parent_coverage}\t{exon_sum_left}\t{exon_sum_right}\t{all_cds_left}\t{all_cds_right}")
        #if black_list:
        #    if float(parent_coverage/exon_sum) > coverage:
        #        print(f"{parent}\t{parent_coverage}\t{exon_sum}")
        #else:
        #    if float(parent_coverage/exon_sum) <= coverage:
        #        print(f"{parent}\t{parent_coverage}\t{exon_sum}")


if __name__ == '__main__':
    intersect_wao = sys.argv[1]  # get gff file from first argument put click in here later
    maximum_overlap_coverage = float(sys.argv[2])  # minimum feature overlap
    cds_transcripts_gff3 = sys.argv[3]  # gff3 file to sum CDS for transcripts
    black_list = False
    if len(sys.argv) == 5:
        black_list = True
    sums = sum_features_per_transcript(cds_transcripts_gff3)
    filter_models(intersect_wao, maximum_overlap_coverage, black_list, sums)
