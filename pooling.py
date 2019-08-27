import sys

def pooling(profiles):
    return [sum(pos_vals) for pos_vals in zip(*profiles)]

def join_sorted(iter_1, iter_2, key=lambda x: x):
    iter_1 = iter(iter_1)
    iter_2 = iter(iter_2)
    try:
        obj_1 = next(iter_1)
        obj_2 = next(iter_2)
        k_1 = key(obj_1)
        k_2 = key(obj_2)
        while True:
            if k_1 == k_2:
                yield (k_1, obj_1, obj_2)
                obj_1 = next(iter_1)
                obj_2 = next(iter_2)
                k_1 = key(obj_1)
                k_2 = key(obj_2)
            elif k_1 < k_2:
                # print(f'join_sorted - skip {k_1} from first iterator', file=sys.stderr)
                obj_1 = next(iter_1)
                k_1 = key(obj_1)
            elif k_1 > k_2:
                # print(f'join_sorted - skip {k_2} from second iterator', file=sys.stderr)
                obj_2 = next(iter_2)
                k_2 = key(obj_2)
    except StopIteration:
        pass
