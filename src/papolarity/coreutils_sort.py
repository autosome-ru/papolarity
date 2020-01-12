import pprint
import six
import subprocess
import pybedtools
from pybedtools.logger import logger
from pybedtools.helpers import BUFSIZE, BEDToolsError
from pybedtools import BedTool

def coreutils_sort(self, **kwargs):
    if 'instream' not in kwargs:
        kwargs['instream'] = self.fn

    cmds, tmp, stdin = self.handle_coreutils_sort_kwargs(prog='sort', **kwargs)
    stream = call_coreutils_sort(cmds, tmp, stdin=stdin)
    return BedTool(stream)

def handle_coreutils_sort_kwargs(self, prog='sort', instream=None, **kwargs):
    """
    Handle coreutils sort program calls.

    *kwargs* are passed directly from the calling method (self.coreutils_sort).

    This method figures out, given how this BedTool was constructed, what
    to send to BEDTools programs -- for example, an open file to stdin with
    the `-` argument, or a filename with the `-a` argument.

    *instream* can be e.g., self.fn or 'a.bed' or an iterator.
    """
    pybedtools.logger.debug(
        'BedTool.handle_coreutils_sort_kwargs() got these kwargs:\n%s',
        pprint.pformat(kwargs))

    stdin = None

    # Decide how to send instream to sort.
    # If it's a BedTool, then get underlying stream
    if isinstance(instream, BedTool):
        instream = instream.fn

    # Filename? No pipe, just provide the file
    if isinstance(instream, six.string_types):
        stdin = None
        input_fn = instream
    # A generator or iterator: pipe it as a generator of lines
    else:
        stdin = (str(i) for i in instream)
        input_fn = '-'

    # If stream not specified, then a tempfile will be created
    if kwargs.pop('stream', None):
        tmp = None
    else:
        output = kwargs.pop('output', None)
        if output:
            tmp = output
        else:
            tmp = BedTool._tmp()

    additional_args = kwargs.pop('additional_args', None)

    # Parse the kwargs into BEDTools-ready args
    cmds = [prog]

    for key, value in sorted(list(kwargs.items()), reverse=True):
        if isinstance(value, bool):
            if value:
                cmds.append('--' + key)
            else:
                continue
        elif isinstance(value, list) or isinstance(value, tuple):
            value = list(map(str, value))

            # sort --key 1,1 --key 2,2r -k 5,5
            for val in value:
                if len(key) == 1:
                    cmds.append('-' + key)
                else:
                    cmds.append('--' + key)
                cmds.append(str(val))
        else:
            cmds.append('--' + key)
            cmds.append(str(value))

    if additional_args:
        cmds.append(additional_args)

    cmds.append(input_fn)
    return cmds, tmp, stdin



def call_coreutils_sort(cmds, tmpfn=None, stdin=None):
    """
    Use subprocess.Popen to call BEDTools and catch any errors.

    Output goes to *tmpfn*, or, if None, output stays in subprocess.PIPE and
    can be iterated over.

    *stdin* is an optional file-like object that will be sent to
    subprocess.Popen.

    Prints some useful help upon getting common errors.

    """
    input_is_stream = stdin is not None
    output_is_stream = tmpfn is None

    try:
        # coming from an iterator, sending as iterator
        if input_is_stream and output_is_stream:
            logger.debug('helpers.call_coreutils_sort(): input is stream, output is stream')
            logger.debug('helpers.call_coreutils_sort(): cmds=%s', ' '.join(cmds))
            p = subprocess.Popen(cmds,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 stdin=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            for line in stdin:
                p.stdin.write(line.encode())
            p.stdin.close()  # This is important to prevent deadlocks

            output = (i.decode('UTF-8') for i in p.stdout)
            stderr = None

        # coming from an iterator, writing to file
        if input_is_stream and not output_is_stream:
            logger.debug('helpers.call_coreutils_sort(): input is stream, output is file')
            logger.debug('helpers.call_coreutils_sort(): cmds=%s', ' '.join(cmds))
            outfile = open(tmpfn, 'wb')
            p = subprocess.Popen(cmds,
                                 stdout=outfile,
                                 stderr=subprocess.PIPE,
                                 stdin=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            if hasattr(stdin, 'read'):
                stdout, stderr = p.communicate(stdin.read())
                p.stdin.close()
            else:
                for item in stdin:
                    p.stdin.write(item.encode())
                stdout, stderr = p.communicate()
                p.stdin.close()
            output = tmpfn
            outfile.close()

        # coming from a file, sending as iterator
        if not input_is_stream and output_is_stream:
            logger.debug('helpers.call_coreutils_sort(): input is filename, output is stream')
            logger.debug('helpers.call_coreutils_sort(): cmds=%s', ' '.join(cmds))
            p = subprocess.Popen(cmds,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            output = (i.decode('UTF-8') for i in p.stdout)
            stderr = None

        # file-to-file
        if not input_is_stream and not output_is_stream:
            logger.debug('helpers.call_coreutils_sort(): input is filename, output is filename (%s)', tmpfn)
            logger.debug('helpers.call_coreutils_sort(): cmds=%s', ' '.join(cmds))
            outfile = open(tmpfn, 'wb')
            p = subprocess.Popen(cmds,
                                 stdout=outfile,
                                 stderr=subprocess.PIPE,
                                 bufsize=BUFSIZE)
            stdout, stderr = p.communicate()
            output = tmpfn
            outfile.close()

        if stderr:
            if isinstance(stderr, bytes):
                stderr = stderr.decode('UTF_8')
            raise BEDToolsError(subprocess.list2cmdline(cmds), stderr)

    except (OSError, IOError) as err:
        print('%s: %s' % (type(err), os.strerror(err.errno)))
        print('The command was:\n\n\t%s\n' % subprocess.list2cmdline(cmds))
        raise err

    return output

pybedtools.BedTool.coreutils_sort = coreutils_sort
pybedtools.BedTool.handle_coreutils_sort_kwargs = handle_coreutils_sort_kwargs
