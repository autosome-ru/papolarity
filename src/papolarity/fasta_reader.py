def fasta_from_file(filename):
    '''Yields pairs (header, sequence). Properly handles multiline FASTA.'''
    with open(filename) as fasta_file:
        header = None
        for line in fasta_file:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header != None:
                    yield (header, ''.join(sequences))
                header = line[1:].lstrip()
                sequences = []
            else:
                sequences.append(line.strip())
        if header != None:
            yield (header, ''.join(sequences))
