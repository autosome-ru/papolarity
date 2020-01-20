import numpy as np

def flatten(xs):
    result = []
    for x in xs:
        if isinstance(x, list):
            result.extend(flatten(x))
        else:
            result.append(x)
    return result

def drop_none(xs):
    return [x for x in xs if x is not None]

def tsv_string_empty_none(row):
    row_strings = [(str(value) if value is not None else '') for value in row]
    return '\t'.join(row_strings)

def take_the_only(arr):
    if len(arr) > 1:
        raise Exception('Several elements when the only one is expected')
    if len(arr) == 0:
        raise Exception('No elements when one is expected')
    return arr[0]

def common_subsequence(iterators, key=lambda x: x, check_sorted=False):
    sentinel = object()
    for (key, aligned_objects) in align_iterators(iterators, key=key, object_missing=sentinel, check_sorted=check_sorted):
        if sentinel not in aligned_objects:
            yield (key, aligned_objects)

_sentinel_key = object()
def align_iterators(iterators, key=lambda x: x, object_missing=None, check_sorted=False):
    '''
    Take a list of iterables yielding objects whose `key` values increase
    strictly monotonically. At each iteration yield a pair consisting of
    a key and a tuple of objects from these interators which have the same key.
    If one of iterators doesn't contain corresponding value or yields,
    then `object_missing` (which is `None` by default) is substituted 
    at that place.
    `key` can be either a callable object (same key for each iterator)
    or a list of callable objects (one key per iterator)
    '''
    iterators = [iter(iterator) for iterator in iterators]
    num_iters = len(iterators)
    exhausted = [False] * num_iters
    objects = [None] * num_iters
    if callable(key):
        key = [key for _ in iterators]
    else:
        assert len(iterators) == len(key)
    key_values = [_sentinel_key] * num_iters
    for idx in range(num_iters):
        objects[idx], key_values[idx], exhausted[idx] = _next_unless_exhausted(iterators[idx], key[idx], exhausted[idx], object_missing=object_missing)
    while not all(exhausted):
        min_key = min(k for k in key_values if k is not _sentinel_key)
        matching_idxs = {idx  for (idx, k) in enumerate(key_values)  if k == min_key}
        aligned_objects = [(objects[idx] if idx in matching_idxs else object_missing)  for idx in range(num_iters)]
        yield (min_key, aligned_objects)
        for idx in matching_idxs:
            old_key_value = key_values[idx]
            objects[idx], key_values[idx], exhausted[idx] = _next_unless_exhausted(iterators[idx], key[idx], exhausted[idx])
            if check_sorted and (key_values[idx] is not _sentinel_key) and (key_values[idx] <= old_key_value):
                raise ValueError(f"Iterator at index `{idx}` is not sorted: `{key_values[idx]}` follows `{old_key_value}`")

def _next_unless_exhausted(iterator, key, exhausted, object_missing=None):
    '''
    Utility function to take values from an iterator and transform them to keys
    until iterator is exhausted and report iterator's status
    '''
    if not exhausted:
        try:
            obj = next(iterator)
            k = key(obj)
            return (obj, k, False)
        except StopIteration:
            return (object_missing, _sentinel_key, True)
    else:
        return (object_missing, _sentinel_key, True)

def get_constant_intervals(profile):
    if len(profile) == 0:
        return []
    last_value_indices = np.nonzero(np.diff(profile))[0]
    start_indices_inclusive = np.concatenate(([0], last_value_indices + 1))
    end_indices_exclusive = np.concatenate( (last_value_indices + 1, [len(profile)]) )
    for (start, end) in zip(start_indices_inclusive, end_indices_exclusive):
        yield (start, end, profile[start])
