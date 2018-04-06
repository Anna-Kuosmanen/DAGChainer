#!/usr/bin/env python
"""
Creates a graph from bam file

Author: aekuosma

"""

import pysam, sys, numpy, collections

MAX_ALTERNATIVE_SOURCES = 10
MAX_ALTERNATIVE_SINKS = 10
window_size = 20

def splice_overlapping_exons(exon_list):
    new_list = list()

    first_start = -1
    last_end = -1
    chrom = exon_list[0][0]

    exon_ends = list()

    for exon in exon_list:
        exon_ends.append(exon[2])
        if first_start == -1 or exon[1] < first_start:
            first_start = exon[1]
        if last_end == -1 or exon[2] > last_end:
            last_end = exon[2]

    positions = [0]*(last_end - first_start + 2)
	
    for exon in exon_list:
        for i in range(exon[1], exon[2]+1):
            positions[i-first_start] = positions[i-first_start] + 1
	
    start = -1
    prev_pos = -1
	
    for index, pos in enumerate(positions):
        if (index + first_start) in exon_ends:
            if start == -1:
                new_list.append((chrom, index+first_start, index+first_start))
            else:
                new_list.append((chrom, start+first_start, index+first_start))
            if pos == 0:
                start = -1
            else:
                start = index+1
            prev_pos = pos
            continue
        if start == -1 and pos == 0:
            continue
        elif start == -1 and pos != 0:
            start = index
            prev_pos = pos
            continue
        elif start != -1 and pos != prev_pos:
            if index > start:
                new_list.append((chrom, start+first_start, first_start+index-1))
            prev_pos = pos
            if pos == 0:
                start = -1
            else:
                start = index
	
    if start != -1:
        new_list.append((chrom, first_start + start, last_end))
        
    return new_list


def calculate_threshold(bamfile, exon_list):

    handle = pysam.Samfile(bamfile, "rb" )	
    threshold = 0
	
    for exon in exon_list:	
        coverage = 0
        it = handle.fetch(exon[0], exon[1], exon[2])
        for read in it:
            if read.is_secondary:
                continue
        
            ranges = list()
            pos = read.pos
            cigar = read.cigar
				
            for part in cigar:
                if part[0] == 0:
                    ranges.append((pos, pos + part[1]))
                if part[0] != 1:
                    pos = pos + part[1]
					
            for rang in ranges:
                if rang[0] > exon[2] or rang[1] < exon[1]:
                    continue
                else:
                    if exon[1] <= rang[0] and exon[2] >= rang[1]:
                        coverage = coverage + (rang[1] - rang[0] + 1)
                        continue
						
        coverage = 1.0*coverage/(exon[2] - exon[1]+1)
		
        threshold = threshold + coverage
	
    threshold = threshold/len(exon_list)
	
    handle.close()
    return threshold


def find_novel_exons(exon_list, bamfile, threshold):

    handle = pysam.Samfile(bamfile, "rb" )

    new_list = list()
    new_sources = list()
    new_sinks = list()

    for ind, exon in enumerate(exon_list):

        if ind == (len(exon_list) - 1):
            break

        second_exon = exon_list[ind + 1]

        # Next to each other
        if exon[2] + 1 == second_exon[1]:
            continue

        gap_start = (long)(exon[2])+1
        gap_end = (long)(second_exon[1]) - 1

        it = handle.fetch(exon[0], gap_start-1, gap_end)

        positions = [0]*(gap_end - gap_start + 1)

        for read in it:

            if read.is_secondary:
                continue
            if read.mapq == 0:
                continue

            for pos in read.positions:
                if pos + 1 >= gap_start and pos + 1 <= gap_end:
                    positions[pos-gap_start] = positions[pos-gap_start] + 1


        start = -1
        prev_pos = -1
        average_coverage = -1

        chrom = exon[0]

        for index, pos in enumerate(positions):
            if start == -1 and pos == 0:
                continue
            elif start == -1 and pos > 0:
                start = index
                prev_pos = pos
                continue
            elif start != -1 and pos > 0:
                average_coverage = average_coverage + pos
            elif start != -1 and pos == 0:
                average_coverage = average_coverage/(index-start + 1)
                if average_coverage >= threshold:
                    new_list.append((chrom, start+gap_start, gap_start+index))
                    if start == 0 and index != gap_end:
                        new_sinks.append((chrom, start+gap_start, gap_start+index))
                prev_pos = pos
                start = -1

        if start != -1 and gap_end - (gap_start + start) > 2:
            average_coverage = average_coverage/(gap_end - (gap_start + start)+1)
            if average_coverage >= threshold:
                new_list.append((chrom, gap_start + start, gap_end))
                if start != 0:
                    new_sources.append((chrom, gap_start + start, gap_end))

    handle.close()
    return sorted(new_list + exon_list, key=lambda exon: (exon[0], exon[1], exon[2])), new_sources, new_sinks



def find_exons(bamfile, chrom, area_start, area_end, threshold):
    
    handle = pysam.Samfile(bamfile, "rb" )
	
    # For ease of use, these are 0-based, end-exclusive (convert before appending to exon_list)
    exon_starts = list()
    exon_ends = list()
			
    it = handle.fetch(chrom, area_start-1, area_end)
		
		
    for read in it:

        if read.is_secondary:
            continue
        if read.mapq == 0:
            continue
 
        # Not spliced, not interesting
        if len(read.cigar) == 1:
            continue
        
        offset = 0
        
        # M = 0 I = 1 D = 2 N = 3 S = 4 H = 5 P = 6
        for part in read.cigar:
            # The splice gap
            if part[0] == 3:
                if read.pos + offset not in exon_ends:
                    if (read.pos + offset) not in exon_ends:
                        exon_ends.append(read.pos + offset)
                if read.pos + offset + part[1] not in exon_starts:
                    if (read.pos + offset + part[1]) not in exon_starts:
                        exon_starts.append(read.pos + offset + part[1])
            if part[0] in [0, 2, 3]:
                offset = offset + part[1]

        
    exon_starts = sorted(exon_starts)
    exon_ends = sorted(exon_ends)

#    print exon_starts
#    print exon_ends
    
    # Just one exon, do coverage search over the whole area
    if len(exon_starts) == 0 and len(exon_ends) == 0:
        it = handle.fetch(chrom, area_start - 1, area_end)
        start = -1
        end = -1
        
        for read in it:
            if read.is_secondary:
                continue
            if read.mapq == 0:
                continue
            if start == -1:
                start = read.pos
            elif start > read.pos:
                start = read.pos
            if end == -1:
                end = read.aend
            elif end < read.aend:
                end = read.aend

        if start != -1 and end != -1:
            return(((chrom, start+1, end+1), ["0"], ["0"]))
        else:
            return ((list(), list(), list()))


    exon_list = list()
    sources = list()
    sinks = list()

    # Pair up exons
    pairless_starts = list()
    pairless_ends = list()
	
    while len(exon_starts) > 0 and len(exon_ends) > 0:

        if exon_starts[0] <= exon_ends[0]:
            # there is no start[1] or start[1] is after end[0], therefore 0 and 0 are pair
            if len(exon_starts) == 1 or exon_starts[1] > exon_ends[0]:
                if exon_starts[0] != exon_ends[0]:
                    exon_list.append((chrom, exon_starts[0] + 1, exon_ends[0]))
                del exon_starts[0]
                del exon_ends[0]

            # start[1] is before end[0], start[0] is pairless
            else:
                pairless_starts.append(exon_starts[0])
                del exon_starts[0]

        # If end[0] is smaller than start[0], end[0] is pairless
        else:
            pairless_ends.append(exon_ends[0])
            del exon_ends[0]
		
	
    for start in exon_starts:
        pairless_starts.append(start)
    for end in exon_ends:
        pairless_ends.append(end)
        

    exon_list = sorted(exon_list, key=lambda exon: (exon[0], exon[1], exon[2]))
#    print "After pairing"
#    print exon_list
#    print pairless_starts
#    print pairless_ends

	# Find the starts and ends of exons without pair and mark source/sink if applicable
    for start_index, start in enumerate(pairless_starts):
        sink = True
        next_start = -1

        # Find the start of the exon following this one (that know when to stop searching)
        for exon in exon_list:

            if next_start == -1 and exon[1] > start:
                next_start = exon[1]
            if exon[1] > start and exon[1] < next_start:
                next_start = exon[1]

        # Check whether the next exon is also one with pairless start
        if start_index < (len(pairless_starts)-1) and next_start > pairless_starts[start_index+1]:
            next_start = pairless_starts[start_index + 1]+1
	if start_index < (len(pairless_starts) -1) and next_start == -1:
            next_start = pairless_starts[start_index + 1]+1


        pair_end = start + 1
        temp_end = -1
        it = handle.fetch(chrom, start, start+1)
        while temp_end != pair_end:
            # Ran into the next exon, stop
            if next_start != -1 and temp_end > next_start:
                pair_end = next_start - 1
                sink = False
                break
        
            # Move forward
            if temp_end != -1:
                pair_end = temp_end
            for read in it:
            
                if read.is_secondary:
                    continue
                if read.mapq == 0:
                    continue
                 
                if pair_end not in read.positions:
                    continue
                prev_read_pos = pair_end
                for pos in read.positions:
                    if pos == prev_read_pos + 1:
                        if pos > temp_end:
                            temp_end = pos
                        prev_read_pos = pos

            if temp_end == -1:
                temp_end = pair_end
            it = handle.fetch(chrom, temp_end, temp_end+1)
    
        if start+1 != pair_end:    
            exon_list.append((chrom, start + 1, pair_end))
            exon_list = sorted(exon_list, key=lambda exon: (exon[0], exon[1], exon[2]))

            if sink:
                sinks.append((chrom, start + 1, pair_end))
	
    for end_index, end in enumerate(pairless_ends):
        source = True
        prev_end = -1
    	
        # Find the exon before this one, that can check when to stop
        for exon in exon_list:
            if prev_end == -1 and exon[2] < end:
                prev_end = exon[2]
            elif prev_end < exon[2] and exon[2] < end:
                prev_end = exon[2]

        # No exon before this one
        if prev_end == -1:
            # Note: End-exclusive, -1
            pair_start = end - 1
            temp_start = -1
            it = handle.fetch(chrom, end-1, end)
            while temp_start != pair_start:
                if temp_start != -1:
                    pair_start = temp_start

                for read in it:
                    if read.is_secondary:
                        continue
                    if read.mapq == 0:
                        continue
                        
                    if pair_start not in read.positions:
                        continue

                    prev_read_pos = pair_start
                    for pos in reversed(read.positions):
                        if pos == prev_read_pos - 1:
                            if temp_start == -1 or pos < temp_start:
                                temp_start = pos
                            prev_read_pos = pos

                if temp_start == -1:
                    temp_start = pair_start
                it = handle.fetch(chrom, temp_start, temp_start+1)

            if pair_start + 1 != end:
                exon_list.append((chrom, pair_start + 1, end))
                exon_list = sorted(exon_list, key=lambda exon: (exon[0], exon[1], exon[2]))

                if source:
                    sources.append((chrom, pair_start + 1, end))
       
        # There is exon before this one        
        # Might be a source, might be the final fragment
        else:
            source = True
            pair_start = end -1
            temp_start = -1
            it = handle.fetch(chrom, end-1, end)
            while temp_start != pair_start:
                if temp_start != -1 and temp_start < prev_end:
                    pair_start = prev_end + 1
                    source = False
                    break
            
                if temp_start != -1:
                    pair_start = temp_start
                for read in it:
                
                    if read.is_secondary:
                        continue
                    if read.mapq == 0:
                        continue
                        
                    if pair_start not in read.positions:
                        continue
                
                    prev_read_pos = pair_start+1
                    for pos in reversed(read.positions):
                        if pos == prev_read_pos - 1:
                            if temp_start == -1 or pos < temp_start:
                                temp_start = pos
                            prev_read_pos = pos

                if temp_start == -1:
                    temp_start = pair_start
                it = handle.fetch(chrom, temp_start, temp_start+1)

            if pair_start != end:
                exon_list.append((chrom, pair_start, end))
                exon_list = sorted(exon_list, key=lambda exon: (exon[0], exon[1], exon[2]))
                if source:
                    sources.append((chrom, pair_start, end))
	
	# Check if there are any missed exons (for example bigger exon in one transcript, with "middle" dropping off in another)
	# Threshold requirement should filter out "junk"
    exon_list = sorted(exon_list, key=lambda exon: (exon[0], exon[1], exon[2]))
#    print "After filling up pairs"
#    print exon_list
    exon_list = splice_overlapping_exons(exon_list)

    exon_list = sorted(exon_list, key=lambda exon: (exon[0], exon[1], exon[2]))
#    print "After splicing"
#    print exon_list

    cov_threshold = calculate_threshold(bamfile, exon_list)/10
    if cov_threshold < 10.0:
        cov_threshold = 10.0

    # The new exons can also be sources or sinks
    exon_list, new_sources, new_sinks = find_novel_exons(exon_list, bamfile, cov_threshold)
    sources = sources + new_sources
    sinks = sinks + new_sinks
    #print "After novel"
    #print exon_list
    alternative_exons = list()
    alternative_sinks = list()
    alternative_sources = list()
    deletable_indices = list()

	# Check for alternative 5' and 3'
    for index, exon in enumerate(exon_list):
        chrom = exon[0]
        start = exon[1]
        end = exon[2]
	    
	    # Check if there's at least threshold difference between coverage of first and coverage of last base, investigate further if so (assume there's no both source and sink in one exon)
        first_cov = 0
        it = handle.fetch(chrom, start-1, start)
        for read in it:
            if start-1 in read.positions:
                first_cov = first_cov + 1
	            
        second_cov = 0
        it = handle.fetch(chrom, end-1, end)
        for read in it:
            if end-1 in read.positions:
                second_cov = second_cov + 1
                
        if first_cov == 0 or second_cov == 0:
            continue
	            
	    # Possible sink
        if 1.0*(first_cov - second_cov)/first_cov > threshold:            
            coverage_table = [0]*(end-start + 1)
            prev_start = start
            for alignment in handle.fetch(chrom, start-1, end):
                for pos in alignment.positions:
                    if (pos > end or pos < start):
                        continue
                    if alignment.is_paired and alignment.pos < alignment.mpos:
                        continue
                    coverage_table[pos - start] = coverage_table[pos - start] +1
            slope_found = False
            descending_count = 0
            level_count = 0

            candidate_end = -1

            for ind, item in enumerate(coverage_table):
                if ind >= len(coverage_table) - window_size-1 or coverage_table[ind-window_size] == 0:
                    break
            
                if item > coverage_table[ind + window_size]:
                    descending_count = descending_count + 1
                # Likely that this isn't a sink, but U-shaped instead
                elif item * 1.1 < coverage_table[ind + window_size]:
                    break
                else:
                    level_count = level_count + 1
                    
                if descending_count > window_size and 1.0*(coverage_table[ind-window_size] - coverage_table[ind])/(coverage_table[ind-window_size]+1) > threshold and 1.0*(coverage_table[ind-window_size] - coverage_table[ind])/(coverage_table[ind-window_size]+1) <= 0.8:
                    level_count = 0
                    slope_found = True
                
                if level_count > window_size and descending_count > window_size and 1.0*(coverage_table[ind-window_size] - coverage_table[ind])/(coverage_table[ind-window_size]+1) > threshold and 1.0*(coverage_table[ind-window_size] - coverage_table[ind])/(coverage_table[ind-window_size]+1) <= 0.8:
                    candidate_end = start + ind - window_size/6
                    alternative_exons.append((chrom, prev_start, candidate_end))
                    alternative_sinks.append((chrom, prev_start, candidate_end))
                    descending_count = 0
                    level_count = 0
                    slope_found = False
                    prev_start = candidate_end + 1
                    candidate_end = -1

                    
            # There's probably alternative end
            if candidate_end != -1:
                deletable_indices.append(index)
                alternative_exons.append((chrom, candidate_end + 1, end))
                alternative_sinks.append((chrom, candidate_end + 1, end))
                
            else:
                if slope_found and exon not in sinks:
                    alternative_sinks.append(exon)
                # Check if there's big difference between first and last position (could be very small slope)
                elif not slope_found and exon not in sinks:
                    if 1.0*(coverage_table[0] - coverage_table[-1])/(coverage_table[0]+1) > threshold and 1.0*(coverage_table[0] - coverage_table[-1])/(coverage_table[0]+1) <= 0.8:
                        alternative_sinks.append(exon)
                        
        
        # Possible source
        if 1.0*(second_cov - first_cov)/second_cov > threshold:
            coverage_table = [0]*(end-start + 1)
            prev_start = start
            for alignment in handle.fetch(chrom, start-1, end):
                for pos in alignment.positions:
                    if (pos > end-1 or pos < start):
                        continue
                    if alignment.is_paired and alignment.pos > alignment.mpos:
                        continue
                    coverage_table[pos - start] = coverage_table[pos - start] +1
            slope_found = False
            ascending_count = 0
            level_count = 0
            
            candidate_start = -1
            for ind, item in enumerate(coverage_table):
                if ind >= len(coverage_table) - window_size-1:
                    break
            
                if item < coverage_table[ind + window_size]:
                    ascending_count = ascending_count + 1
                # Likely that this isn't a source after all, but instead U-shaped
                elif item > 1.1*coverage_table[ind + window_size]:
                    break
                else:
                    level_count = level_count + 1
                    
                    
                # Slope found at start of exon
                if level_count < window_size/2 and ascending_count > window_size and 1.0*(coverage_table[ind] - coverage_table[ind-window_size])/(coverage_table[ind]+1) > threshold and 1.0*(coverage_table[ind] - coverage_table[ind-window_size])/(coverage_table[ind]+1) <= 0.8:
                    slope_found = True
                    level_count = 0
                    ascending_count = 0
                
                # Slope found after there's been some part of exon, mark alternative start
                elif ascending_count > window_size and 1.0*(coverage_table[ind] - coverage_table[ind-window_size])/(coverage_table[ind]+1) > threshold and 1.0*(coverage_table[ind] - coverage_table[ind-window_size])/(coverage_table[ind]+1) <= 0.8:
                    candidate_start = start + ind-(window_size/6)
                    alternative_exons.append((chrom, prev_start, candidate_start - 1))      
                    if slope_found:
                        alternative_sources.append((chrom, prev_start, candidate_start - 1))
                    slope_found = False
                    level_count = 0
                    ascending_count = 0
                    prev_start = candidate_start
                    continue
                    
            if candidate_start != -1:
                deletable_indices.append(index)
                alternative_exons.append((chrom, candidate_start, end))
                alternative_sources.append((chrom, candidate_start, end))
                
            else:
                if slope_found and exon not in sources:
                    alternative_sources.append(exon)
                # Check if there's big difference between first and last position (could be very small slope)
                elif not slope_found and exon not in sources:
                    if 1.0*(coverage_table[-1] - coverage_table[0])/(coverage_table[-1]+1) > threshold and 1.0*(coverage_table[-1] - coverage_table[0])/(coverage_table[-1]+1) <= 0.8:
                        alternative_sources.append(exon)

#    print "Alternative exons"
#    print alternative_exons

    if len(alternative_sources) < MAX_ALTERNATIVE_SOURCES and len(alternative_sinks) < MAX_ALTERNATIVE_SINKS:
        offset = 0              
        for ind in deletable_indices:
            if exon_list[ind-offset] in sources:
                del sources[sources.index(exon_list[ind-offset])]
            if exon_list[ind-offset] in sinks:
                del sinks[sinks.index(exon_list[ind-offset])]
            del(exon_list[ind-offset])
            offset = offset +1

        exon_list = exon_list + alternative_exons
        exon_list = sorted(exon_list, key=lambda exon: (exon[0], exon[1], exon[2]))

        sources = sources + alternative_sources
        sinks = sinks + alternative_sinks

    if exon_list[0] not in sources:
        sources.append(exon_list[0])
    if exon_list[-1] not in sinks:
        sinks.append(exon_list[-1])
	       
    sources = sorted(sources, key=lambda exon: (exon[0], exon[1], exon[2]))
    sinks = sorted(sinks, key=lambda exon: (exon[0], exon[1], exon[2]))
	
    source_indices = list()
    sink_indices = list()

    for source in sources:
        if source in exon_list:
            source_indices.append(str(exon_list.index(source)))
    for sink in sinks:
        if sink in exon_list:
            sink_indices.append(str(exon_list.index(sink)))
	
    handle.close()

    return exon_list, source_indices, sink_indices

# Find arcs and calculate node and arc coverage from bam file
def parse_bam(bamfile, exon_list):

    handle = pysam.Samfile(bamfile, "rb" )
	
    exon_weights = list()
    adjacency_list = list()
    edge_weights = list()
    
    # Only one exon, python drops the list structure
    if not isinstance(exon_list, list):
        coverage = 0
        it = handle.fetch(exon_list[0], exon_list[1]-1, exon_list[2])
        for read in it:
		
            if read.is_secondary:
                continue
            if read.mapq == 0:
                continue		
            for pos in read.positions:
                if pos+1 >= exon_list[1] and pos+1 <= exon_list[2]:
                    coverage = coverage + 1
					
        coverage = 1.0*coverage/(exon_list[2] - exon_list[1]+1)
        exon_weights.append(str(coverage))
        adjacency_list.append(list())
        edge_weights.append(list())
        return exon_weights, adjacency_list, edge_weights
    

    for index, exon in enumerate(exon_list):
        # Count exon weight
        coverage = 0
        it = handle.fetch(exon[0], exon[1]-1, exon[2])
        for read in it:
		
            if read.is_secondary:
                continue
            if read.mapq == 0:
                continue		
            for pos in read.positions:
                if pos+1 >= exon[1] and pos+1 <= exon[2]:
                    coverage = coverage + 1
					
        coverage = 1.0*coverage/(exon[2] - exon[1]+1)
        exon_weights.append(str(coverage))
		
        adjacencies = list()
        adjacency_weights = list()
		
        for second_index in range(index+1, len(exon_list)):
            second_exon = exon_list[second_index]			
            it = handle.fetch(exon[0], exon[2]-1, exon[2])

            # Created from splitting overlapping exons			
            if exon[2] + 1 == second_exon[1]:
                shared_reads = 0
			
                for read in it:
				
                    if read.is_secondary:
                        continue
	            if read.mapq == 0:
                        continue			

                    if exon[2]-1 in read.positions and second_exon[1]-1 in read.positions:
                        shared_reads = shared_reads + 1
				
                if shared_reads > 0:
                    adjacencies.append(str(second_index))
                    adjacency_weights.append(str(shared_reads))
				
			
            # Otherwise
            else:
			
                common_reads = 0
			
                # Cigar operations: M = 0 I = 1 D = 2 N = 3 S = 4 H = 5 P = 6
                for read in it:
                    if read.is_secondary:
                        continue
                    if read.mapq == 0:
                        continue
                    cigar = read.cigar

                    # areas that the actual parts of read cover
                    ranges = list()
                    pos = read.pos
                    for part in cigar:
                        if part[0] == 0:
                            ranges.append((pos+1, pos + part[1]))
                        if part[0] != 1 and part[0] != 4 and part[0] != 5:
                            pos = pos + part[1]
		
                    for ind, rang in enumerate(ranges):
                        if rang[0] > exon[2] or rang[1] < exon[1]:
                            continue
                        # If last, end
                        if rang == ranges[-1]:
                            continue
                        else:
                            if exon[2] == rang[1]: 
                                # Check if the second exon is in the next range (blocks false arcs between first and other than second exon for reads spanning more than 2 exons)
                                candidate_second_range = ranges[ind+1]
                                if second_exon[1] == candidate_second_range[0]:
                                    common_reads = common_reads + 1
					
								
				
                if common_reads > 0:			
                    adjacencies.append(str(second_index))
                    adjacency_weights.append(str(common_reads))
		

        adjacency_list.append(adjacencies)
        edge_weights.append(adjacency_weights)

    handle.close()
    
    return exon_weights, adjacency_list, edge_weights


def print_graph(exon_list, exon_weights, adjacency_list, edge_weights, sources, sinks, output):
	handle = open(output, "w")
	
	# Number of vertices
	if isinstance(exon_list, list):
	    handle.write(str(len(exon_list)))
	else:
	    handle.write("1")
	handle.write("\n")
	
	# List of adjacencies
	for adjacencies in adjacency_list:
		handle.write(" ".join(adjacencies))
		handle.write("\n")
	
	# Vertex weights
	handle.write(" ".join(exon_weights))
	handle.write("\n")
	
	# Edge weights
	for weights in edge_weights:
		handle.write(" ".join(weights))
		handle.write("\n")
		
    # Sources
	handle.write(" ".join(sources))
	handle.write("\n")
    
    # Sinks
	handle.write(" ".join(sinks))
	handle.write("\n")
		
	handle.close()

def print_nodes(exon_list, output):
    handle = open(output, "w")

    if isinstance(exon_list, list):
        for exon in exon_list:
            handle.write(str(exon[0]) + "\t" + str(exon[1]) + "\t" + str(exon[2]) + "\n")
    else:
        handle.write(str(exon_list[0]) + "\t" + str(exon_list[1]) + "\t" + str(exon_list[2]) + "\n")

    handle.close()

def validate_graph(adjacencies):
   new_sources = list()
   new_sinks = list()

   no_of_exons = len(adjacencies)

   lonely = [True]*no_of_exons;

   for index, adjacency in enumerate(adjacencies):
      # If it has outneighbors, it's not lonely
      if len(adjacency) > 0:
          lonely[index] = False
      # If someone points to it (inneighbor), it's not lonely
      for item in adjacency:
          lonely[int(item)] = False

   for index, item in enumerate(lonely):
       if(item):
           new_sources.append(str(index))
           new_sinks.append(str(index))

   return new_sources, new_sinks

# Find subpath constraints created by reads spanning more than two exons
def find_subpath_constraints(exon_list, bamfile):

   # Subpath as key, time it appeared as value
   subpath_coverages = collections.OrderedDict()
   handle = pysam.Samfile(bamfile, "rb" )

   # Select all the reads between start of first exon and end of last exon
   it = handle.fetch(exon_list[0][0], exon_list[0][1], exon_list[-1][2])

   # Check all the reads in this area
   for read in it:
      covered_exons = list()
      for pos in read.positions:
         for index, exon in enumerate(exon_list):
            if (pos+1) >= exon[1] and (pos+1) <= exon[2]:
               if str(index) not in covered_exons:
                  covered_exons.append(str(index))
               break

      if len(covered_exons) > 2 and " ".join(covered_exons) not in subpath_coverages.keys():
         subpath_coverages[" ".join(covered_exons)] = 1
      elif len(covered_exons) > 2:
         subpath_coverages[" ".join(covered_exons)] = subpath_coverages[" ".join(covered_exons)] + 1

   handle.close()

   # Normalize by the length
   for key in subpath_coverages:
      path = key.split(" ")
      length = 0
      for node in path:
         length += (exon_list[int(node)][2]-exon_list[int(node)][1]+1)
      subpath_coverages[key] = 1.0*subpath_coverages[key]/length

   # Lists to return stuff sorted
   subpath_list = list()
   subpath_cov = list()

   for key in subpath_coverages:
      subpath_list.append(key)
      subpath_cov.append(subpath_coverages[key])

   return subpath_list, subpath_cov

def add_subpath_information(subpath_list, subpath_coverages, output):
   handle = open(output, "a")
   # If there are no constraints, write 0
   if len(subpath_list) == 0:
	handle.write("0\n")
   else:
       if isinstance(subpath_list, list):
          handle.write("%s\n" % (len(subpath_list)))
          for subpath in subpath_list:
             handle.write("%s\n" % (subpath))
             # The next lines are the coverages
          for coverage in subpath_coverages:
		handle.write("%s\n" % (coverage))
       else:
          handle.write("1\n")
          handle.write("%s\n" % (subpath_list))
          handle.write("%s\n" % (subpath_coverages))    

   handle.close()

# Caller gives chromosome, start and end within which to search
def create_graph_without_annotation(bamfile, chrom, start, end, graph, nodes, threshold, debug):

    exon_list, sources, sinks = find_exons(bamfile, chrom, start, end, threshold)
    if isinstance(exon_list, list):
    	exon_list = sorted(exon_list, key=lambda exon: (exon[0], exon[1], exon[2]))
    if len(exon_list) == 0:
        return 0
    exon_weights, adjacency_list, edge_weights = parse_bam(bamfile, exon_list)
    #	exon_weights = adjust_for_slope(bamfile, exon_list, exon_weights, sources, sinks)
    if debug:
        print exon_list
        print exon_weights
        print adjacency_list
        print edge_weights

    # Checks for exons disconnected from graph and marks them as both source and sink
    new_sources, new_sinks = validate_graph(adjacency_list)
    for source in new_sources:
        if source not in sources:
            sources.append(source)
    for sink in new_sinks:
        if sink not in sinks:
            sinks.append(sink)

    sources = sorted(sources)
    sinks = sorted(sinks)

    if debug:
        print sources
        print sinks
	

    print_graph(exon_list, exon_weights, adjacency_list, edge_weights, sources, sinks, graph)
    print_nodes(exon_list, nodes)
    if isinstance(exon_list, list):
        subpath_list, subpath_coverages = find_subpath_constraints(exon_list, bamfile)
        if debug:
            print subpath_list
        add_subpath_information(subpath_list, subpath_coverages, graph)

    return 1

