#!/usr/bin/env python3

import sys


def get_families(genefamilies):
    '''Read genefamilies file return families dict'''
    families = {}
    with open(genefamilies) as gopen:  # make class method for this
        for line in gopen:  # iterate through lines in genefamilies file
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            feature = fields[0]
            genefamily = fields[1]
            model = fields[2]
#            print(fields[-1])
            bscore = float(fields[5])
            score = float(fields[4])
            if score/bscore >= 1.5:  # disclude inflated scoring due to multiple copies or fusions
                pass
            families[feature] = {'gene_family': genefamily, 'score': bscore,
                                 'model': model}
    return families


def choose_model(overlaps, genefamiliesa, genefamiliesb, minoverlap, no_orphans, name):
    '''Choose feature with best gene family assignment between overlapping features'''
    families = {'a': get_families(genefamiliesa), 'b': get_families(genefamiliesb)}
    my_ids = {}
    alt_ids = {}
    eliminated_ids = {}
    with open(overlaps) as vopen:
        for line in vopen:
            line = line.rstrip()
            if not line:
                continue
            fields = line.split('\t')
#            print(line)
            overlap = int(fields[-1])
            featureid1 = fields[8].split(';')[0].split('=')[1]
#            print(featureid1)
            record1 = "\t".join(fields[:9])  # reconstruct original gff feature1
#            print(record1)
            if not int(overlap):  # feature has no overlap add to retain
                my_ids[featureid1] = {'family': None, 'score': 0, 'model': None,
                                      'parent': featureid1, 'record': record1}
                continue
            record2 = "\t".join(fields[9:-1])  # reconstruct original gff feature2
#            print(record2)
            featureid2 = fields[-2].split(';')[0].split('=')[1]
            start1 = int(fields[3])
            stop1 = int(fields[4])
            if overlap/(stop1 - start1) <= minoverlap:  # skip feature
                continue
            score1 = 0
            score2 = 0
            family1 = None
            model1 = None
            family2 = None
            model2 = None
#            print(featureid1, featureid2)
            if not featureid1 in my_ids:
                my_ids[featureid1] = {'family': None, 'score': 0, 'model': None,
                                      'parent': featureid1, 'record': record1}
            if featureid1 in families['a']:  # feature a has a family
                family1 = families['a'][featureid1]['gene_family']
                score1 = families['a'][featureid1]['score']
                model1 = families['a'][featureid1]['model']
            if featureid2 in families['b']:  # feature b has a family
                family2 = families['b'][featureid2]['gene_family']
                score2 = families['b'][featureid2]['score']
                model2 = families['b'][featureid2]['model']
            if not family1 and not family2:  # niether features have families
#                print(f'No Gene Family Assignments, {featureid1}, {featureid2}')
#                eliminated_ids[featureid1] = "no familiy1 or family2"
                eliminated_ids[featureid2] = f"{featureid1}\t{featureid2}\t{score1}\t{score2}\tno gene families"
                continue
            if family1 == family2:  # features have same family
#                print(f'Same families, {featureid1}, {featureid2}')
                if score1 == score2:  # features are equal
                    if my_ids[featureid1]['score'] <= score1:
                        my_ids[featureid1]['score'] = score1
                        my_ids[featureid1]['family'] = family1
                        my_ids[featureid1]['model'] = model1
                        my_ids[featureid1]['parent'] = featureid1
                        my_ids[featureid1]['record'] = record1
                        eliminated_ids[featureid2] = f"{featureid1}\t{featureid2}\t{score1}\t{score2}\tscore1 = score2"
#                        if featureid2 in my_ids:
#                            del my_ids[featureid2]
#                    print(f'Model1 and Model2 Same Score')
                elif score1 > score2:  # take score 1
                    if my_ids[featureid1]['score'] <= score1:
                        my_ids[featureid1]['score'] = score1
                        my_ids[featureid1]['family'] = family1
                        my_ids[featureid1]['model'] = model1
                        my_ids[featureid1]['parent'] = featureid1
                        my_ids[featureid1]['record'] = record1
                        eliminated_ids[featureid2] = f"{featureid1}\t{featureid2}\t{score1}\t{score2}\tscore1 > score2"
#                        if featureid2 in my_ids:
#                            del my_ids[featureid2]
#                    print(f'Model1 Higher Score, {featureid1}, {score1}, {score2}')
                else:  # take score 2
                    if my_ids[featureid1]['score'] < score2:
                        my_ids[featureid1]['score'] = score2
                        my_ids[featureid1]['family'] = family2
                        my_ids[featureid1]['model'] = model2
                        my_ids[featureid1]['parent'] = featureid2
                        my_ids[featureid1]['record'] = record2
                        eliminated_ids[featureid1] = f"{featureid1}\t{featureid2}\t{score1}\t{score2}\tscore1 < score2"
#                    print(f'Model2 Higher Score, {featureid2}, {featureid1}, {score1}, {score2}')
            else:
                if not family1:  # family 1 has no model but family 2 does
#                    print(f'Feature1 has no gene model, Feature2 does, {featureid1}, {featureid2}')
                    if my_ids[featureid1]['score'] < score2:
                        my_ids[featureid1]['score'] = score2
                        my_ids[featureid1]['family'] = family2
                        my_ids[featureid1]['model'] = model2
                        my_ids[featureid1]['parent'] = featureid2
                        my_ids[featureid1]['record'] = record2
                        eliminated_ids[featureid1] = f"{featureid1}\t{featureid2}\t{score1}\t{score2}\tscore1 NA"
                elif not family2:  # family 2 has no model but family 1 does
                    if my_ids[featureid1]['score'] <= score1:
                        my_ids[featureid1]['score'] = score1
                        my_ids[featureid1]['family'] = family1
                        my_ids[featureid1]['model'] = model1
                        my_ids[featureid1]['parent'] = featureid1
                        my_ids[featureid1]['record'] = record1
                        eliminated_ids[featureid2] = f"{featureid1}\t{featureid2}\t{score1}\t{score2}\tscore2 NA"
#                        if featureid2 in my_ids:
#                            del my_ids[featureid2]
#                    print(f'Feature2 has no gene model, Feature1 does, {featureid1}, {featureid2}')
                else:
                    if featureid2 not in eliminated_ids:
                        if not family2 and no_orphans:
                            continue
                        my_ids[featureid2] = {'family': family2, 'score': score2, 
                                              'model': model2, 'parent': featureid2,
                                              'record': record2}
#                    print(f'Different Families, {featureid1}, {featureid2}, {family1}, {family2}')

#            if featureid in intervals:
#                print(f'{line};genefamily={intervals[featureid]["gene_family"]};genefamily_score={intervals[featureid]["score"]}')
    for f in my_ids:
#        print(f)
#        print(my_ids[f])
#        if my_ids[f]['parent'] not in eliminated_ids and f not in eliminated_ids:
        if my_ids[f]['parent'] not in eliminated_ids:
            print(my_ids[f]['record'])
#        if my_ids[f]['model']:
#            print(my_ids[f]['record'])
#        else:
#            print(f)
    eliminated = open(f'./{name}.eliminated.tsv', 'w')
    for f in eliminated_ids:
        eliminated.write(f'{eliminated_ids[f]}' + '\n')
    eliminated.close()


if __name__ == '__main__':
    intersect = sys.argv[1]  # get bed file from first argument put click in here later
    genefamiliesa = sys.argv[2]  # get gene family assignments and scores for a
    genefamiliesb = sys.argv[3]  # get gene family assignments and scores for b
    minoverlap = float(sys.argv[4])  # minimum feature overlap
    no_orphans = sys.argv[5]
    name = sys.argv[6]
#    genefamilies = (genefamilies1, genefamilies2)
    choose_model(intersect, genefamiliesa, genefamiliesb, minoverlap, no_orphans, name)
