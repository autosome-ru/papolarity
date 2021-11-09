# Papolarity
Papolarity is a Python package for analysis of transcript-level short read coverage profiles.

For a single sample, for each transcript papolarity allows for computing the classic polarity metric which, in the case of Ribo-Seq, reflects ribosome positional preferences.

For comparison versus a control sample, papolarity estimates an improved metric, the relative linear regression slope of coverage along transcript length. This involves de-noising by profile segmentation with a Poisson model (using [pasio](https://github.com/autosome-ru/pasio/)), and aggregation of Ribo-Seq coverage within segments, thus achieving reliable estimates of the regression slope.

*Publication:* [Assessing Ribosome Distribution Along Transcripts with Polarity Scores and Regression Slope Estimates](https://pubmed.ncbi.nlm.nih.gov/33765281/). 
Methods Mol Biol. 2021;2252:269-294. doi:10.1007/978-1-0716-1150-0_13

![papolarity logo](supplementary_files/papolarity_logo.png)

## Toolkit

Papolarity provide a toolkit to perform different tasks necessary for processing transcriptomic data such as Ribo-Seq alignments.

Installation:
```
python -m pip install papolarity
```

There're good chances that you'd also install pasio: `python -m pip install pasio`

The package is organized as a single entry point for a set of subcommands.

You can run it with one of these commands:

* `papolarity [arguments]`

* `python -m papolarity [arguments]` - if you need to specify a certain version of python to run a package.

Python 3.7+ is supported; probably this restriction will be relaxed later.

There are no conventions about a structure of folders and file names. All files that are used by tool are always specified in command line arguments.

Papolarity have a few conventions about file extensions: all files with `.gz` extension are treated as gzip archives. Input files with names ending with `.gz` will be automatically unpacked, output files will be automatically packed. Character `-` instead of filename will be treated as stdin or stdout. It can be useful to use papolarity in pipelined commands.

You can follow the protocol to get the idea how these tools are supposed to be used. If you need to customize pipelines, please reference to help for corresponding tools:
`papolarity --help` lists all available tools. `papolarity <cmd> --help` shows description of all arguments and options for a specified tool.

## Protocol

In our paper "Assessing Ribosome Distribution Along Transcripts with Polarity Scores and Regression Slope Estimates" ([doi:10.1007/978-1-0716-1150-0_13](https://doi.org/10.1007/978-1-0716-1150-0_13)) we describe a protocol for Ribo-Seq analysis. In a file [protocol-paper.sh](https://github.com/autosome-ru/papolarity/blob/master/protocol-paper.sh) you can find a script we used in a paper to process our datasets. It's slightly modified for better readability compared to a paper, and is more easily customizable. Also it has a few additional commands to generate plots which are absent in paper. Steps are named after paper sections.

You can use this protocol as is or change any parts you wish. As long as you comply with data formats and use consistent data (e.g. all files should be clipped in the same manner, or non-clipped at all), papolarity will work, order of commands, folder names, filenames and so on doesn't matter.

To run this pipeline, you should have several auxiliary tools installed: [csvtk](https://bioinf.shenwei.me/csvtk/), [GNU parallel](https://www.gnu.org/software/parallel/), and python package [pasio](https://github.com/autosome-ru/pasio/).
