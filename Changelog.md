# Release notes 
### Version 1.1.1
* Fixed bug which caused transcripts with zero CDS start coordinates (relative to Transcript start) as non-coding. It's typical for transcripts with non-annotated TSS. This bug resulted in broken `cds_features.tsv` and can cause problems at later stages of pipeline.

### Version 1.1.0
* Fixed problem for contigs written in both upper and lower cases. Tools now can sort transcripts in case-insensitive manner and check sortedness accordingly (and do it by default).
Old syntax `papolarity get_coverage sample.bam --sort` is analogous to new syntax `papolarity get_coverage sample.bam --sort case-sensitive`. But we recommend to use `--sort case-insensitive` option (see corresponding changes in `protocol_paper.sh`).

### Version 1.0.0
First stable version.
