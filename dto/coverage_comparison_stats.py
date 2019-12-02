_coverage_stats_fields = [
    'transcript_info', 'slope',
    'control_mean_coverage', 'experiment_mean_coverage',
    'control_total_coverage', 'experiment_total_coverage',
    'control_polarity_score', 'experiment_polarity_score',
    'num_segments',
    'multipoint_slope', 'profile_difference',
]

class CoverageComparisonStats(namedtuple('CoverageComparisonStats', _coverage_stats_fields)):
    @property
    def gene_id(self):
        return self.transcript_info.gene_id

    @property
    def transcript_id(self):
        return self.transcript_info.transcript_id

    @property
    def transcript_length(self):
        return self.transcript_info.transcript_length

    @property
    def cds_start(self):
        return self.transcript_info.cds_start

    @property
    def cds_stop(self):
        return self.transcript_info.cds_stop

    def __repr__(self):
        fields = [
            self.gene_id, self.transcript_id, self.transcript_length, self.cds_start, self.cds_stop,
            self.slope,
            self.control_mean_coverage, self.experiment_mean_coverage,
            self.control_total_coverage, self.experiment_total_coverage,
            self.control_polarity_score, self.experiment_polarity_score,
            self.num_segments,
            self.multipoint_slope, self.profile_difference,
        ]
        return '\t'.join(map(str, fields))

    @classmethod
    def header(cls):
        fields = [
            'gene_id', 'transcript_id', 'transcript_length', 'cds_start', 'cds_stop',
            'slope',
            'control_mean_coverage', 'experiment_mean_coverage',
            'control_total_coverage', 'experiment_total_coverage',
            'control_polarity_score', 'experiment_polarity_score',
            'num_segments',
            'multipoint_slope', 'profile_difference',
        ]
        return '\t'.join(fields)

    @classmethod
    def from_string(cls, line):
        row = line.rstrip('\n').split('\t')
        gene_id, transcript_id, transcript_length, cds_start, cds_stop, \
            slope, \
            control_mean_coverage, experiment_mean_coverage, \
            control_total_coverage, experiment_total_coverage, \
            control_polarity_score, experiment_polarity_score, \
            num_segments, \
            multipoint_slope, profile_difference, \
            *rest = row
        info = {
            'transcript_info': CodingTranscriptInfo(gene_id, transcript_id, int(transcript_length), int(cds_start), int(cds_stop)),
            'slope': float(slope),
            'control_mean_coverage': float(control_mean_coverage), 'experiment_mean_coverage': float(experiment_mean_coverage),
            'control_total_coverage': int(control_total_coverage), 'experiment_total_coverage': int(experiment_total_coverage),
            'control_polarity_score': float(control_polarity_score), 'experiment_polarity_score': float(experiment_polarity_score),
            'num_segments': int(num_segments),
            'multipoint_slope': float(multipoint_slope),
            'profile_difference': float(profile_difference),
        }
        return cls(**info)

    @classmethod
    def make_from_profiles(cls, transcript_info, cds_profile_control, cds_profile_experiment, segments):
        info = {
            'transcript_info': transcript_info,
            'slope': slope_by_profiles(cds_profile_control, cds_profile_experiment, segments, mode='center'),
            'control_mean_coverage': np.mean(cds_profile_control), 'experiment_mean_coverage': np.mean(cds_profile_experiment),
            'control_total_coverage': np.sum(cds_profile_control), 'experiment_total_coverage': np.sum(cds_profile_experiment),
            'control_polarity_score': polarity_score(cds_profile_control), 'experiment_polarity_score': polarity_score(cds_profile_experiment),
            'num_segments': len(segments),
            'multipoint_slope': slope_by_profiles(cds_profile_control, cds_profile_experiment, segments, mode='every_point'),
            'profile_difference': profile_difference(cds_profile_control, cds_profile_experiment, segments),
        }
        return cls(**info)

    @classmethod
    def each_in_file(cls, filename):
        with open(filename) as f:
            header = f.readline()
            for line in f:
                yield cls.from_string(line)

    @classmethod
    def print(cls, infos, file=sys.stdout, extended=False):
        if extended:
            infos = list(infos)
            best_transcript_ids = set(info.transcript_id for info in choose_best_transcript(infos))

            print(cls.header() + '\t' + '\t'.join(['geom_mean_coverage', 'polarity_difference', 'is_best_transcript']), file=file)
            for info in infos:
                if info.transcript_id in best_transcript_ids:
                    is_best_transcript = '+'
                else:
                    is_best_transcript = '-'
                print(str(info) + '\t' + '\t'.join(map(str, [info.geom_mean_coverage(), info.polarity_difference(), is_best_transcript])), file=file)
        else:
            print(cls.header(), file=file)
            for info in infos:
                print(info, file=file)

    def geom_mean_coverage(self):
        return (self.control_mean_coverage * self.experiment_mean_coverage) ** 0.5

    def polarity_difference(self):
        return self.experiment_polarity_score - self.control_polarity_score
