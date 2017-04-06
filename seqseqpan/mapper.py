import bisect
from collections import defaultdict

from seqseqpan.exception import CoordinateOutOfBoundsError


class Mapper:
    def map_coordinates(self, alignment, consensus, source, destinations, coordinates):

        if type(destinations) is not list:
            destinations = [destinations]

        destinations = sorted(destinations)
        coordinates = sorted(coordinates)
        coord_dict = defaultdict(dict)

        if min(coordinates) < 1:
            raise CoordinateOutOfBoundsError(min(coordinates),source)

        add_consensus = ("c" in destinations)

        # store and do not reorder!
        lcbs = alignment.lcbs

        # remove consensus-destination from list (calculation is different)
        if add_consensus:
            destinations = sorted(destinations)[:-1]

        if source == "c":

            cons_len = sum([l.length for l in lcbs])
            if max(coordinates) > cons_len:  # coordinates higher than alignment length can't be mapped
                raise CoordinateOutOfBoundsError(max(coordinates), source)

            i = 0
            for coord in coordinates:

                idx = bisect.bisect_left(consensus.block_start_indices, coord)
                idx -= 1

                lcb = lcbs[idx]

                # calculate position within current consensus block - there are no gaps in consensus sequence!
                pos_within_block = coord - consensus.block_start_indices[idx] - 1

                if add_consensus:
                    coord_dict[coord]["c"] = coord

                coord_dict[coord].update(self._get_coords_for_entries(lcb.entries, destinations, pos_within_block))
                i += 1
        else:

            source_blocks = {}
            # store ends of blocks for finding lcb for coords
            for i in range(len(lcbs)):
                e = lcbs[i].get_entry(int(source))
                if e is not None:
                    source_blocks[e.end] = {"lcb": i, "entry": e}

            if max(coordinates) > max(source_blocks.keys()):  # coordinates higher than genome length can't be mapped
                raise CoordinateOutOfBoundsError(max(coordinates), source)

            for coord in coordinates:
                source_ends = sorted(source_blocks.keys())
                idx_in_ends = bisect.bisect_left(source_ends, coord)
                lcb_idx = source_blocks[source_ends[idx_in_ends]]["lcb"]

                if lcb_idx > len(lcbs):
                    raise CoordinateOutOfBoundsError(coord, source)

                source_entry = source_blocks[source_ends[idx_in_ends]]["entry"]
                lcb = lcbs[lcb_idx]

                if source_entry.strand == "+":
                    pos_within_block_without_gaps = coord - source_entry.start
                else:
                    pos_within_block_without_gaps = source_entry.end - coord

                # add consensus coordinates to dict
                if add_consensus:
                    consensus_length = sum([lcb.length for lcb in lcbs[0:lcb_idx]])
                    coord_dict[coord]["c"] = consensus_length + pos_within_block_without_gaps + 1

                # check if dests other than consensus are needed
                if len(destinations) > 0:
                    # add 1 for get_position_within_entry_with_gaps as it works with one-based indices
                    # subtract 1 afterwards as we are working with zero-based indices here
                    pos_within_block = source_entry.get_position_within_entry_with_gaps(pos_within_block_without_gaps + 1 ) - 1

                    coord_dict[coord].update(self._get_coords_for_entries(lcb.entries, destinations, pos_within_block))

        return coord_dict

    def _get_coords_for_entries(self, entries, destinations, pos_within_block):
        coord_dict = {}

        for e in entries:  # faster than looping through entries everytime to get entry of genome x
            if str(e.genome_nr) in destinations:

                is_gap = False

                # gaps are always counted from the left, independent of strand!
                e_gap_sublist = e.get_gap_sublist(0, pos_within_block)
                e_sum_gaps = 0

                if len(e_gap_sublist) > 0:
                    e_sum_gaps = sum(end - start for start, end in e_gap_sublist.items())
                    last_gap = sorted(e_gap_sublist.items())[-1]
                    is_gap = (last_gap[0] <= pos_within_block < last_gap[1])

                if not is_gap:
                    if e.strand == "+": # only position calculation with known gaps is dependent on strand!
                        coord_dict[str(e.genome_nr)] = e.start + (pos_within_block - e_sum_gaps)
                    else:
                        coord_dict[str(e.genome_nr)] = (e.end - (pos_within_block - e_sum_gaps)) * -1

        return coord_dict
