from utils import pool_profiles
import logging
import pasio

logger = logging.getLogger('pasio')
logger.setLevel(logging.WARNING)

def make_joint_segmentation(coverages, splitter):
    pooled_coverage = pool_profiles(coverages)
    if len(pooled_coverage) == 0:
        return None
    return list(pasio.segments_with_scores(pooled_coverage, splitter))

def stabilize_profile(profile, segments):
    stable_profile = np.zeros_like(profile)
    for segment in segments:
        start = segment.start
        stop = segment.stop
        stable_profile[start:stop] = np.mean(profile[start:stop])
    return stable_profile
