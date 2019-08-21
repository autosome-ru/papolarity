from pasio import *
import numpy as np
import logging

logger = logging.getLogger('pasio')
logger.setLevel(logging.WARNING)

alpha = 1
beta = 1
scorer_factory = ScorerFactory(alpha, beta)
length_regularization_multiplier = 0
length_regularization_function = lambda x: x
split_number_regularization_multiplier = 0
algorithm = 'rounds'
no_split_constant = True
window_size = 2500
window_shift = 1250
num_rounds = None
split_at_gaps = False

def pasio_splitter():
    square_splitter = SquareSplitter(scorer_factory,
        length_regularization_multiplier = 0,
        length_regularization_function = lambda x: x,
        split_number_regularization_multiplier = 0,
        split_number_regularization_function = lambda x: x)

    if algorithm == 'exact':
        splitter = square_splitter
    else:
        if no_split_constant:
            square_splitter = ReducerCombiner(NotConstantReducer(), square_splitter)
        else:
            square_splitter = ReducerCombiner(NotZeroReducer(), square_splitter)
        ###
        if algorithm == 'slidingwindow':
            sliding_window = SlidingWindow(window_size = window_size, window_shift = window_shift)
            reducer = SlidingWindowReducer(sliding_window = sliding_window, base_reducer = square_splitter)
            splitter = ReducerCombiner(reducer, square_splitter)
        elif algorithm == 'rounds':
            sliding_window = SlidingWindow(window_size = window_size, window_shift = window_shift)
            base_reducer = SlidingWindowReducer(sliding_window = sliding_window, base_reducer = square_splitter)
            reducer = RoundReducer(base_reducer = base_reducer, num_rounds = num_rounds)
            splitter = ReducerCombiner(reducer, NopSplitter(scorer_factory))
    return splitter

def stabile_segments(profile, splitter):
    counts = np.array(profile)
    split_candidates = np.arange(len(counts) + 1)
    score, splits = splitter.split(counts, split_candidates)
    scorer = splitter.scorer(counts, splits)

    for (start, stop, mean_count) in zip(splits, splits[1:], scorer.mean_counts()):
        yield (start, stop, mean_count)
