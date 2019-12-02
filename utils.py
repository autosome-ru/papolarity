def take_the_only(arr):
    if len(arr) > 1:
        raise Exception('Several elements when the only one is expected')
    if len(arr) == 0:
        raise Exception('No elements when one is expected')
    return arr[0]
