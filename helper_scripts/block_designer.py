import pandas as pd
from copy import deepcopy
import random
import numpy as np
import sys


identity_field = sys.argv[3]
to_stratify = sys.argv[4:-1]
NUM_BATCHES = int(sys.argv[-1])

def extract(row, columns):
    return dict(zip(columns, row[columns].values))

def gen_expected(samples, num_batches):
    counts = {}
    for sample in samples:
        for field in to_stratify:
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

def calc_score(expected, observed, weights=None):
    delta = {}
    for field in expected.keys():
        for value, count in expected[field].items():
            value = str(value)
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

def get_value_counts(d):
    values = []
    for _, value in d.items():
        if type(value) is not dict:
            values.append(value)
        else:
            values.extend(get_value_counts(value))
    return values

def split_samples(samples, num_batches, weights, round_robin=False):
    expected = gen_expected(samples, num_batches)
    groups = [[] for _ in range(num_batches)]
    scores = [calc_score(expected, group, weights) for group in groups]
    working_samples = random.sample(samples, len(samples))
    used = [False for _ in groups]
    last_score = np.inf
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
        print("Total Penalty: ", np.sum(scores), " Improvement: ", last_score - np.sum(scores))
        last_score = np.sum(scores)
        if np.all(used):
            used = [False for _ in groups]
        working_samples = random.sample(working_samples, len(working_samples))
    return [random.sample(group, len(group)) for group in random.sample(groups, len(groups))]

if sys.argv[1].endswith("xlsx"):
    samples = pd.read_excel(sys.argv[1])
elif sys.argv[1].endswith("csv"):
    samples = pd.read_csv(sys.argv[1])

samples = samples.apply(extract, axis=1, args=(samples.columns,))
samples = list(samples)

new_samples = []
for sample in samples:
    new_sample = {}
    for k, v in sample.items():
        try:
            new_sample[k.rstrip()] = v.rstrip()
        except:
            new_sample[k] = v
    new_samples.append(new_sample)
samples = new_samples
fields = list(samples[0].keys())

for f in fields:
    if f not in to_stratify and f != identity_field:
        for sample in samples:
            del sample[f]

#all_samples = []
#_, sample_sets = split_samples(samples, NUM_BATCHES, None)
#while sample_sets:
#    sample_set = sample_sets.pop()
#    num_subgroups = len(sample_set)
#    smallest_subgroup = max(get_value_counts(gen_expected(sample_set, )))


# This function implements the Euclidean 
# algorithm to find H.C.F. of two number
import math

def find_gcd(iterable):
    gcd = int(0)
    for val in iterable:
        gcd = int(math.gcd(gcd, int(val)))
    return gcd


all_samples = []
for b_no, batch in enumerate(split_samples(samples, NUM_BATCHES, None)):
    num_subgroups = find_gcd(get_value_counts(gen_expected(batch, NUM_BATCHES)))        
    for sub_b_no, sub_batch in enumerate(split_samples(batch, int(num_subgroups), None)):
        for sample in random.sample(sub_batch, len(sub_batch)):
            sample["batch_no"] = b_no
            sample["sub_batch_no"] = sub_b_no
            all_samples.append(sample)
randomized = pd.DataFrame(all_samples)
randomized['run_order'] = np.arange(len(randomized)) + 1
randomized.to_csv(sys.argv[2], sep=",")