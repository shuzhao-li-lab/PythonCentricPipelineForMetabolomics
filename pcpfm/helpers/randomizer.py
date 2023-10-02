import pandas as pd
from copy import deepcopy
import random
import numpy as np
import sys

def bin_to_quantiles(values, quantiles):
    bins = []
    for value in values:
        for quantile, quant_cutoff in enumerate(quantiles):
            if quant_cutoff < value:
                pass
            elif quant_cutoff >= value:
                bins.append(quantile)
                break
    return bins

def draw_date_to_bin(values, cutoff):
    bins = []
    for value in values:
        if value < cutoff:
            bins.append(0)
        else:
            bins.append(1)
    return bins

def extract(row, columns):
    return dict(zip(columns, row[columns].values))

def create_id(sample):
    sample['id'] = sample['Sample ID']
    return sample

def cleanup_sample(sample):
    for x in ['JAX ID', 'DOB', 'Date Drawn']:
        del sample[x]
    return sample

def group_samples(samples):
    groups = {}
    for sample in samples:
        if sample['id'] not in groups:
            groups[sample['id']] = []
        groups[sample['id']].append(sample)
    return groups

def gen_expected(samples, num_batches=1):
    counts = {}
    for sample in samples:
        for field in ['Clone ID', 'KO strategy', 'Media']:
            if field not in counts:
                counts[field] = {}
            val = sample[field]
            if val not in counts[field]:
                counts[field][val] = 0
            if num_batches:
                counts[field][val] += 1/num_batches
            else:
                counts[field][val] += 1
    return counts

def gen_weights(samples):
    counts = {}
    for sample in samples:
        for field in ['Clone ID', 'KO strategy', 'Media']:
            val = sample[field]
            if (field, val) not in counts:
                counts[(field, val)] = 0
            counts[(field, val)] += 1
    max_count = max(counts.values())
    print({k: len(samples)/v for k,v in counts.items()}
)
    return {k: max_count/v for k,v in counts.items()}

def gen_weights2(samples):
    counts = {}
    for sample in samples:
        for field in ['Clone ID', 'KO strategy', 'Media']:
            val = sample[field]
            if (field, val) not in counts:
                counts[(field, val)] = 0
            counts[(field, val)] += 1
    min_per_field = {}
    for key, value in counts.items():
        if key[0] not in min_per_field:
            min_per_field[key[0]] = 0
        min_per_field[key[0]] = max(value, min_per_field[key[0]])
    weights = {k: min_per_field[k[0]]/v for k,v in counts.items()}
    print(weights)
    return weights

def calc_score(expected, observed, weights=None):
    delta = {}
    for field in expected.keys():
        for value, count in expected[field].items():
            if field in observed and value in observed[field]:
                if weights:
                    delta[field + value] = (expected[field][value] - observed[field][value]) * weights[(field, value)]
                else:
                    delta[field + value] = (expected[field][value] - observed[field][value])
            else:
                if weights:
                    delta[field + value] = count * weights[(field, value)]
                else:
                    delta[field + value] = count 
    score = 0
    for field, delta_value in delta.items():
        score += delta_value ** 2

    return score

def organize_locations(positions, by="JAX ID"):
    return {p[by]: p for p in positions}

def map_to_location(sample, positions):
    if sample["JAX ID"] in positions:
        sample["S.No"] = "Box: " + positions[sample["JAX ID"]]["Box"] + " ,Loc: " + positions[sample["JAX ID"]]["Location"] + " ,S.No: " + str(positions[sample["JAX ID"]]["S.No"])
    return sample
    
def get_d_values(d):
    values = []
    for _, value in d.items():
        if type(value) is not dict:
            values.append(value)
        else:
            values.extend(get_d_values(value))
    return values

def split_samples(samples, num_batches, weights, round_robin=False):
    expected = gen_expected(samples, num_batches)
    groups = [[] for _ in range(num_batches)]
    scores = [calc_score(expected, group, weights) for group in groups]
    working_samples = random.sample(samples, len(samples))
    used = [False for _ in groups]
    while working_samples:
        best_improved = (None, None)
        for j, already_used in enumerate(used):
            for i, sample in enumerate(working_samples):
                copies = [deepcopy(group) for group in groups]
                for copy in copies:
                    copy.append(sample)
                test_scores = [calc_score(expected, gen_expected(copy, num_batches), weights) for copy in copies]
                if not already_used:
                    improvement = scores[j] - test_scores[j]
                    if best_improved[1] is None or improvement > best_improved[1]:
                        best_improved = [[i,j], improvement]


        sample_to_add = working_samples.pop(best_improved[0][0])
        groups[best_improved[0][1]].append(sample_to_add)
        scores[best_improved[0][1]] = calc_score(expected, gen_expected(groups[best_improved[0][1]], num_batches), weights)
        used[best_improved[0][1]] = round_robin
        if np.all(used):
            used = [False for _ in groups]
        working_samples = random.sample(working_samples, len(working_samples))

    return [random.sample(group, len(group)) for group in random.sample(groups, len(groups))]

samples = pd.read_excel(sys.argv[1])
samples = samples.apply(extract, axis=1, args=(samples.columns,))
samples = [create_id(sample) for sample in samples]

#weights = gen_weights(samples)
weights = None
groups = group_samples(samples)
expected = gen_expected(samples)

all_samples = []
total_count = 0
for b_no, batch in enumerate(split_samples(samples, 1, weights)):
    smallest_subgroup = min(get_d_values(gen_expected(batch)))
    for sub_b_no, sub_batch in enumerate(split_samples(batch, int(round(smallest_subgroup)), weights)):
        print("\t", sub_b_no, len(sub_batch))
        for k,v in gen_expected(sub_batch).items():
            print("\t\t", k,v)
        for sample in random.sample(sub_batch, len(sub_batch)):
            total_count += 1
            sample["batch_no"] = b_no
            sample["sub_batch_no"] = sub_b_no
            all_samples.append(sample)
randomized = pd.DataFrame(all_samples)
randomized.drop(columns='Box label (Media)', inplace=True)
randomized['run_order'] = np.arange(len(randomized)) + 1
randomized.to_csv(sys.argv[2], sep=",")