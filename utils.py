def take_the_only(arr):
    if len(arr) > 1:
        raise Exception('Several elements when the only one is expected')
    if len(arr) == 0:
        raise Exception('No elements when one is expected')
    return arr[0]

def common_subsequence(iterators, key=lambda x: x, check_sorted=False):
    for (key, aligned_objects) in align_iterators(iterators, key=key, check_sorted=check_sorted):
        if all(aligned_objects):
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
    '''
    iterators = [iter(iterator) for iterator in iterators]
    num_iters = len(iterators)
    exhausted = [False] * num_iters
    objects = [None] * num_iters
    keys = [_sentinel_key] * num_iters
    for idx in range(num_iters):
        objects[idx], keys[idx], exhausted[idx] = _next_unless_exhausted(iterators[idx], key, exhausted[idx], object_missing=object_missing)
    while not all(exhausted):
        min_key = min(k for k in keys if k is not _sentinel_key)
        matching_idxs = {idx  for (idx, k) in enumerate(keys)  if k == min_key}
        aligned_objects = [(objects[idx] if idx in matching_idxs else None)  for idx in range(num_iters)]
        yield (min_key, aligned_objects)
        for idx in matching_idxs:
            old_key = keys[idx]
            objects[idx], keys[idx], exhausted[idx] = _next_unless_exhausted(iterators[idx], key, exhausted[idx])
            if check_sorted and (keys[idx] is not _sentinel_key) and (keys[idx] <= old_key):
                raise ValueError(f"Iterator at index `{idx}` is not sorted: `{keys[idx]}` follows `{old_key}`")

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

def pool_profiles(profiles):
    return [sum(pos_vals) for pos_vals in zip(*profiles)]
