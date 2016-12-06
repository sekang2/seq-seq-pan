class FormatError(Exception):
    pass


class InputError(Exception):
    pass


class ParameterError(Exception):
    def __init__(self, parameter, value, range_text):
        self.parameter = parameter
        self.value = value
        self.range_text = range_text


class ConsensusFastaFormatError(FormatError):
    def __init__(self):
        self.message = "ERROR: Wrong format of consensus fasta header. " \
                       "Please rebuild consensus fasta with current script version."


class ConsensusFastaError(InputError):
    def __init__(self):
        self.message = "ERROR: Consensus fasta contains none or more than one entry. " \
                       "Please rebuild consensus fasta with current script version."


class ConsensusFastaInputError(InputError):
    def __init__(self, character):
        self.message = "ERROR: Your consensus sequence contains a non-IUPAC character: " + character + "."


class XMFAHeaderFormatError(FormatError):
    def __init__(self, header):
        self.message = "ERROR: XMFA sequence header \"" + header + "\" does not apply to XMFA header format rule: " \
                                                                   "\"> seq:start-end strand\"."


class LcbInputError(InputError):
    def __init__(self, lcbnr):
        self.message = 'ERROR: Problem with LCB Nr. {0}: Entries must not be of different length.'.format(lcbnr)


class ConsensusXMFAInputError(InputError):
    def __init__(self):
        self.message = "ERROR: XMFA with more than 2 genomes provided for splitting or merging LCBs. " \
                       "Please align genomes to consensus sequence one by one, " \
                       "creating a new consensus sequence for every added genome."


class ConsensusFastaIdxFormatError(InputError):
    def __init__(self, message):
        self.message = "ERROR: Format of consensus fasta index file not correct:" + message


class CoordinateOutOfBoundsError(InputError):
    def __init__(self, coord, source):
        source = "consensus" if source == "c" else str(source)
        self.message = "ERROR: Position (" + str(coord) + ") not part of sequence (" + source + ".)"


class CoordinatesInputError(InputError):
    def __init__(self):
        self.message = 'ERROR: Please provide indices for mapping in correct format: ' \
                       'First line: source_seq\tdest_seq[,dest_seq2,...] using "c" or sequence number. ' \
                       'Then one coordinate per line.'


class ConsensusGenomeNumberError(InputError):
    def __init__(self):
        self.message = "ERROR: Could not assign XMFA sequence number to consensus file. " \
                       "Was consensus sequence used for aligning this XMFA file?"


class ConsensusCorruptError(InputError):
    def __init__(self, pos):
        self.message = (
            "ERROR: Could not find delimiter sequence at position " + str(pos)
            + ". Your consensus sequence is incorrect.")
