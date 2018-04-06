/* createGraphWithSubpathConstraints.h
 *
 *  Created on: November 17, 2014
 *      Author: aekuosma
 */
#ifndef CREATEGRAPH_H_
#define CREATEGRAPH_H_

#include <iostream>
#include <algorithm>
#include "../../SamReader.h"
#include "../../SamAlignment.h"
#include "MPCUtil.h"
#include "MPCHeaders.h"

// Used for alternative primes, not in use currently
//static int MAX_ALTERNATIVE_SOURCES = 10;
//static int MAX_ALTERNATIVE_SINKS = 10;
//static int window_size = 20;
// Set this as something suitably large if you don't want to cut
// Colinear aligner requires exons to be at most 1024 bases long
static int MAX_EXON_LENGTH = 1024;

class createGraphWithSubpathConstraints {


// Represents an exon
struct Exon {

	std::string chrom;
	long start;
	long end;

	bool operator< (const Exon& str) const{
		if(chrom != str.chrom)
			return(chrom < str.chrom);
		else
			return (start < str.start);
	}

	bool operator==(const Exon& str) const{
		return (chrom == str.chrom && start == str.start && end == str.end);
	}

	bool operator!=(const Exon& str) const{
		return !(chrom == str.chrom && start == str.start && end == str.end);
	}

	std::string toString() {
		std::stringstream sstm;
		sstm << "(" << chrom << ", " << start << ", " << end << ")";
		return sstm.str();
	}
};

// Represents a range of continuous alignment (i.e. spliced alignment creates 2 or more ranges)
struct AlignmentRange {

	// Name of the read, required to parse the different parts of same read
	std::string name;
	std::string chrom;
	long start;
	long end;
	// Used for checking exon boundaries
	bool first_fragment;
	bool last_fragment;
	// How many'th fragment is this, 0-based
	int index;

	std::string toString() {
		std::stringstream sstm;
  		sstm << "(" << name << ", " << chrom << ", " << start << ", " << end << "(first: " << first_fragment << ", last: " << last_fragment << ", fragment index: " << index << "))";
		return sstm.str();
	}


};

// Sort AlignmentRanges by name
struct less_than_name {
	inline bool operator() (const AlignmentRange& range1, const AlignmentRange& range2) {
		return (range1.name < range2.name);
	}
};

// Sort AlignmentRanges by start position
struct less_than_start{
	inline bool operator() (const AlignmentRange& range1, const AlignmentRange& range2) {
	if(range1.start != range2.start)
        	return (range1.start < range2.start);
	else
		return (range1.end < range2.end);
	}
};

// Sort AlignmentRanges by end position
struct less_than_end {
	inline bool operator() (const AlignmentRange& range1, const AlignmentRange& range2) {
		if(range1.end != range2.end)
        		return (range1.end < range2.end);
		else
			return(range1.start < range2.start);
	}
};

private:

	static Exon make_exon(std::string chrom_, long start_, long end_) {
		Exon myexon;
		myexon.chrom = chrom_;
		myexon.start = start_;
		myexon.end = end_;

		return myexon;
	}

	static AlignmentRange make_range(std::string name_, std::string chrom_, long start_, long end_, bool first_, bool last_, int index_) {
		AlignmentRange range;
		range.name = name_;
	        range.chrom = chrom_;
		range.start = start_;
		range.end = end_;
		range.first_fragment = first_;
		range.last_fragment = last_;
	        range.index = index_;

	        return range;
	}

	// Returns the index of the exon that contains the given position
	static int find_exon_index(std::vector<Exon> vec, long pos) {

		size_t mid, left = 0 ;
		size_t right = vec.size();
		while (left < right) {
			mid = left + (right - left)/2;
			if (pos > vec[mid].end){
				left = mid+1;
			}
			else if (pos < vec[mid].start){                                 			right = mid;
			}
			else {                                                          			return mid;
			}
		}
		return -1;      
	}

	// Returns the ranges that contains the given position
	static vector<AlignmentRange> find_ranges(std::vector<AlignmentRange> ranges, long pos) {

		//size_t mid, left = 0 ;
		//size_t right = ranges.size();

		std::vector<AlignmentRange> results;
		// TODO binary search (Not as easy as you'd think)
		for(unsigned i=0;i<ranges.size();i++)
	 		if(ranges[i].start <= pos && ranges[i].end >= pos)
	     			results.push_back(ranges[i]);

		return results;
	}

	// Returns the ranges that end at given position
	static std::vector<AlignmentRange> find_ranges_end(std::vector<AlignmentRange> ranges, long pos) {
 		size_t mid, left = 0 ;
		size_t right = ranges.size();

		std::vector<AlignmentRange> results;

		while (left < right) {
			mid = left + (right - left)/2;
			if (pos > ranges[mid].end){
	 			left = mid+1;
			}
			else if (pos < ranges[mid].end){                        				right = mid;
			}
			// Found one range, look around it for other range
			else {
				// scroll backward till find the first range
				while(mid > 0 && ranges[mid-1].end >= pos) {
					if(mid >= 1)
						mid--;
					else
						break;
				}
				while(ranges[mid].end == pos) {                 					results.push_back(ranges[mid]);
					mid++;
					if(mid == ranges.size())
						break;
				}
				return results;
			}                                                       		}
		return results;
	}

	// Returns the ranges that start at given position
	static std::vector<AlignmentRange> find_ranges_start(std::vector<AlignmentRange> ranges, long pos) {
		size_t mid, left = 0 ;
		size_t right = ranges.size();

		std::vector<AlignmentRange> results;

		while (left < right) {
			mid = left + (right - left)/2;
			if (pos > ranges[mid].start){
				left = mid+1;
			}
			else if (pos < ranges[mid].start){                      				right = mid;
			}
			// Found one range, look around it for other range
			else {
				// scroll backward till find the first range
				while(mid > 0 && ranges[mid-1].start >= pos) {
					if(mid >= 1)
						 mid--;
					else
						 break;
				}        
				while(ranges[mid].start == pos) {                                                         
					results.push_back(ranges[mid]);
					mid++;
					if(mid == ranges.size())
						break;

				}
				return results;
			}                                                                                                               
		}

		return results;
	}


	// Returns the coverage of the given area
	// TODO binary search (not as easy as you'd think, glitches because the coordinates that are not used for sorting can be out of order)
	static double count_coverage(std::vector<AlignmentRange> ranges, long start, long end) {
		int coverage = 0;
		for(unsigned i=0; i<ranges.size();i++) {
			// Whole range is within the area                                   
			if(ranges[i].start >= start && ranges[i].end <= end)
				coverage += ranges[i].end-ranges[i].start + 1;
				// Whole area is within the range
			else if(start >= ranges[i].start && end <= ranges[i].end)
				 coverage += (end - start + 1);
			// Overlap from start
			else if(ranges[i].end >= start && ranges[i].end <= end)
				 coverage += ranges[i].end - start + 1;
			// Overlap from end
			else if(ranges[i].start >=start && ranges[i].start <= end)
				coverage += end - ranges[i].start + 1;
		}
		return 1.0*coverage/(end-start+1);
	}


	// Parses ranges from SAM file
	// 0-based
	static std::vector<AlignmentRange> parse_ranges(SamReader* reader, std::string chrom, long start, long end) {
		std::vector<AlignmentRange> ranges;

		std::vector<SamAlignment> alignments = reader->readAlignmentsRegion(chrom, start, end);

		for(unsigned al=0;al<alignments.size();al++) {
		
			SamAlignment alignment = alignments.at(al);
			if(alignment.isSecondary())
				continue;
			if(alignment.getMappingQuality() == 0)
				continue;

			long pos = alignment.getReferenceStart();

			long prev_start = alignment.getReferenceStart();
			bool first_fragment = true;
			int counter = 0;

			std::vector<std::pair<int,char> > cigar = alignment.getCigarTuples();
			for(unsigned j=0;j<cigar.size();j++) {
				if(cigar.at(j).second == 'N') {
					AlignmentRange tmp = make_range(alignment.getQueryName(), chrom, prev_start, pos-1, first_fragment, false, counter);
					ranges.push_back(tmp);
					first_fragment = false;
					counter++;
					prev_start = pos+cigar.at(j).first;
				}
				// If not insert or soft/hard clip, update pos
				if(cigar.at(j).second != 'I' && cigar.at(j).second != 'S' && cigar.at(j).second != 'H')
					pos += cigar.at(j).first;
			}
			AlignmentRange tmp = make_range(alignment.getQueryName(), chrom, prev_start, pos-1, first_fragment, true, counter);
			ranges.push_back(tmp);

		}
		alignments.clear();
		return ranges;
	}

	// When parsing exons, there can be overlaps. Splice into pseudoexons.
	static std::vector<Exon> splice_overlapping_exons(std::vector<Exon> &exon_list) {
		std::vector<Exon> new_list;

		int first_start = -1;
		int last_end = -1;
		std::string chrom = exon_list.at(0).chrom;

		std::vector<int> exon_ends;

		for(unsigned i=0;i<exon_list.size();i++) {
			exon_ends.push_back(exon_list.at(i).end);
			if(first_start == -1 || exon_list.at(i).start < first_start)
				first_start = exon_list.at(i).start;
			if(last_end == -1 or exon_list.at(i).end > last_end)
				last_end = exon_list.at(i).end;
		}

		// Holds the "exon coverage" of each position, that is, number of exons that cover each position
		int *positions = new int[last_end-first_start+2];

		for(int i=0;i<last_end-first_start+2;i++)
			positions[i] = 0;

		for(unsigned i=0;i<exon_list.size();i++)
			for(int j=exon_list.at(i).start;j<=exon_list.at(i).end;j++)
				positions[j-first_start]++;


		int start = -1;
		int prev_pos = -1;
		for(int i=0;i<last_end-first_start+2;i++) {
			bool in_exon_ends = false;
			for(unsigned j=0;j<exon_ends.size();j++){
				if(i+first_start == exon_ends.at(j)) {
					in_exon_ends = true;
					break;
				}
			}
			if(in_exon_ends) {
				if(start == -1)
					new_list.push_back(make_exon(chrom, i+first_start, i+first_start));
				else
					new_list.push_back(make_exon(chrom, start+first_start, i+first_start));

				if(positions[i] == 0)
					start = -1;
				else
					start = i+1;
				prev_pos = positions[i];
				continue;
			}
			if(start == -1 && positions[i] == 0)
				continue;
			else if(start == -1 && positions[i] != 0) {
				start = i;
				prev_pos = positions[i];
				continue;
			}
			else if(start != -1 && positions[i] != prev_pos) {
				if(i > start)
					new_list.push_back(make_exon(chrom, start+first_start, first_start+i-1));
				prev_pos = positions[i];
				if(positions[i] == 0)
					start = -1;
				else
					start = i;            
			}

			// Cut if too long
			if(start != -1 && i-start+1 > MAX_EXON_LENGTH) {

				new_list.push_back(make_exon(chrom, start+first_start, i+first_start-1));
				start = i;
			}
		}

		if(start != -1)
			new_list.push_back(make_exon(chrom, first_start + start, last_end));

		delete[] positions;

		return new_list;
	}


	static std::vector<Exon> find_exons(std::vector<AlignmentRange> &ranges, std::vector<AlignmentRange> &ranges_by_end, std::string chrom, long area_start, long area_end, double alt_prime_threshold) {

		std::vector<Exon> exon_list;
		long prev_start = ranges.at(0).start;
		long prev_end = ranges.at(0).end;
		// Prevents exons from likely mismappings at ends
		long prev_exon_end = -1;
		// Prevents exons from likely mismappings at starts
		bool last_was_first_fragment = false;

		for(unsigned i=0;i<ranges.size();i++) {
			AlignmentRange range = ranges.at(i);
			// Add the previous exon if there's no overlap and it hasn't been added yet
			if(range.start > prev_end) {
				Exon exon = make_exon(range.chrom, prev_start, prev_end);
				if(find(exon_list.begin(), exon_list.end(), exon) == exon_list.end()) {
					// Check that it's not a likely mismapping
					if(prev_exon_end < (prev_end-5)) {
						exon_list.push_back(exon);
					}
				}
				prev_start = range.start;
				if(range.end > prev_end)
					prev_end = range.end;
			}
			// Check that the previous exon isn't too close (if it was created from the first fragment)
			// If that's the case, might need to modify start
			if(exon_list.size() > 0 && last_was_first_fragment && !range.first_fragment && exon_list.at(exon_list.size()-1).start >= (range.start-5)) {
				Exon exon = exon_list.at(exon_list.size()-1);
				exon_list.pop_back();
				exon.start = range.start;

				exon_list.push_back(exon);
			}


			// If it's neither first nor last, then it's a whole exon by itself
			if(!range.first_fragment && !range.last_fragment) {
				Exon exon = make_exon(range.chrom, range.start, range.end);
				if(find(exon_list.begin(), exon_list.end(), exon) == exon_list.end()) {
					exon_list.push_back(exon);

					prev_exon_end = range.end;
				}

				if(range.end > prev_end)
					prev_end = range.end;

				// If there is a reason to believe that the previous start was mismapping, correct it
				if(last_was_first_fragment && (range.start-prev_start) <= 5)
					prev_start = range.start;

	
			}
			// Range ends at end of exon, add this exon if it hasn't been added already
			else if(!range.last_fragment) {
				Exon exon = make_exon(range.chrom, prev_start, range.end);

				if(find(exon_list.begin(), exon_list.end(), exon) == exon_list.end()) {
					exon_list.push_back(exon);

					prev_exon_end = range.end;
				}

				if(range.end > prev_end)
					prev_end = range.end;

			}
			// An exon starts here and there's overlap (if there wasn't overlap, previous if set prev_start = range.start), mark previous exon and move prev_start
			else if(!range.first_fragment && prev_start != range.start) {
				Exon exon = make_exon(range.chrom, prev_start, range.start-1);

				if(!last_was_first_fragment || (range.start-prev_start) > 5) {
					if(find(exon_list.begin(), exon_list.end(), exon) == exon_list.end()) {
						exon_list.push_back(exon);

						prev_exon_end = range.end;
					}
				}
	    
				prev_start = range.start;
				if(range.end > prev_end)
					prev_end = range.end;
			}

			// Overlap by last fragment (so we know start, but not end of this exon, update candidate end if it's bigger than the current end)
			else if(range.last_fragment && range.start <= prev_end && prev_end < range.end) {
				prev_end = range.end;
			}
			if(range.first_fragment)
				last_was_first_fragment = true;
			else
				last_was_first_fragment = false;

		}
		// Add last exon
		exon_list.push_back(make_exon(ranges.at(0).chrom, prev_start, prev_end));


		// Splice overlaps
		exon_list = splice_overlapping_exons(exon_list);

		return exon_list;
	}

	/* Find arcs and calculate node and arc coverage */
	static std::tuple<std::vector<double>, std::vector<std::vector<int> >, std::vector<std::vector<long> > > parse_sam(std::vector<AlignmentRange> &ranges_by_start, std::vector<AlignmentRange> &ranges_by_end, std::vector<Exon> &exon_list){

		// One pass over the ranges to link starts and ends to names

		std::vector<double> exon_weights;
		std::vector<vector<int> > adjacency_list;
		std::vector<vector<long> > edge_weights;
		// The ranges that start at this position
		std::vector<vector<AlignmentRange> > all_start_alignments;
		// The ranges that end at this position
		std::vector<vector<AlignmentRange> > all_end_alignments;

		// Count exon weights and fill the ranges for arcs
		for(unsigned i=0;i<exon_list.size();i++) {
			exon_weights.push_back(count_coverage(ranges_by_start, exon_list.at(i).start, exon_list.at(i).end));
			std::vector<AlignmentRange> start_alignments = find_ranges_start(ranges_by_start, exon_list.at(i).start);
			sort(start_alignments.begin(), start_alignments.end(), less_than_name());
			all_start_alignments.push_back(start_alignments);

			std::vector<AlignmentRange> end_alignments = find_ranges_end(ranges_by_end, exon_list.at(i).end);
			sort(end_alignments.begin(), end_alignments.end(), less_than_name());
			all_end_alignments.push_back(end_alignments);
		}

		// Find adjacencies and their weights
		for(unsigned i=0;i<exon_list.size();i++) {
			std::vector<int> adjacencies;
			std::vector<long> adjacency_weights;


				for(unsigned j=i+1;j<exon_list.size();j++ ) {
					if(exon_list.at(i).end < exon_list.at(j).start-1) {

						std::vector<AlignmentRange> end_alignments = all_end_alignments.at(i);
						std::vector<AlignmentRange> start_alignments = all_start_alignments.at(j);

						int shared = 0;
						unsigned start_counter = 0;
						unsigned end_counter = 0;

						while(true) {

							if(start_counter == start_alignments.size() || end_counter == end_alignments.size())
								break;
							if(start_alignments.at(start_counter).name == end_alignments.at(end_counter).name) {
								// Checks that this is actually next piece
								if(end_alignments.at(end_counter).index == start_alignments.at(start_counter).index-1)
									shared++;

									start_counter++;
									end_counter++;
							}
							else if(start_alignments.at(start_counter).name > end_alignments.at(end_counter).name)
								end_counter++;

							else
								start_counter++;
						}


						if(shared > 0) {
							adjacencies.push_back(j);
							adjacency_weights.push_back(shared);
						}
				}
				// The exons are actually pseudoexons, they're possibly from same range, check with find_range function
				// This should not increase the running time too much
				else {
					std::vector<AlignmentRange> start_alignments = find_ranges(ranges_by_start, exon_list.at(j).start);
					std::vector<AlignmentRange> end_alignments = find_ranges(ranges_by_start, exon_list.at(i).end);

					sort(start_alignments.begin(), start_alignments.end(), less_than_name());
					sort(end_alignments.begin(), end_alignments.end(), less_than_name());

					int shared = 0;
					unsigned start_counter = 0;
					unsigned end_counter = 0;

					while(true) {

						if(start_counter == start_alignments.size() || end_counter == end_alignments.size())
							break;
						if(start_alignments.at(start_counter).name == end_alignments.at(end_counter).name) {
							shared++;

							start_counter++;
							end_counter++;
						}
						else if(start_alignments.at(start_counter).name > end_alignments.at(end_counter).name)
							end_counter++;

						else
							start_counter++;
					}


					if(shared > 0) {
						adjacencies.push_back(j);
						adjacency_weights.push_back(shared);
					}

				}
			}
			adjacency_list.push_back(adjacencies);
			edge_weights.push_back(adjacency_weights);
		}

		return make_tuple(exon_weights, adjacency_list, edge_weights);

	}


	// Print the graph file
	static void print_graph(std::vector<Exon> &exon_list, std::vector<double> &exon_weights, std::vector<std::vector<int> > &adjacencies, std::vector<std::vector<long> > &edge_weights, std::vector<int> &sources, std::vector<int> &sinks, std::string output) {
		std::ofstream handle(output.c_str());

		// Number of nodes
		handle << exon_list.size() << std::endl;
		// Adjacencies
		for(unsigned i=0;i<adjacencies.size();i++) {
			if(adjacencies.at(i).size() == 0) {
				handle << std::endl;
				continue;
			}

			for(unsigned j=0;j<adjacencies.at(i).size()-1;j++) {
				handle << adjacencies.at(i).at(j) << " ";
			}
			handle << adjacencies.at(i).at(adjacencies.at(i).size()-1) << std::endl;
		}
		// Vertex weights
		for(unsigned i=0;i<exon_weights.size()-1;i++)
			handle << exon_weights.at(i) << " ";
		handle << exon_weights.at(exon_weights.size()-1) << std::endl;

		// Edge weights
		for(unsigned i=0;i<edge_weights.size();i++) {
			if(edge_weights.at(i).size() == 0) {
				handle << std::endl;
				continue;
			}
			for(unsigned j=0;j<edge_weights.at(i).size()-1;j++) {
				handle << edge_weights.at(i).at(j) << " ";
			}
			handle << edge_weights.at(i).at(edge_weights.at(i).size()-1) << std::endl;
		}

		// Sources
		for(unsigned i=0;i<sources.size()-1;i++)
			handle << sources.at(i) << " ";
		handle << sources.at(sources.size()-1) << std::endl;

		// Sinks
		for(unsigned i=0;i<sinks.size()-1;i++)
			handle << sinks.at(i) << " ";
		handle << sinks.at(sinks.size()-1) << std::endl;

		handle.close();
	}

	// Print the nodes file
	// Exons are carried as 0 based, +1 them to be consistent with gtf format
	static void print_nodes(vector<Exon> &exon_list, string output) {
		std::ofstream handle(output.c_str());

		for(unsigned i=0;i<exon_list.size();i++) {
			handle << exon_list.at(i).chrom << "\t" << (exon_list.at(i).start +1) << "\t" << (exon_list.at(i).end+1) << std::endl;
		}

		handle.close();
	}

	// Find sources and sinks
	// If it has no in-neighbors, it's a source, and if it has no out-neighbors, it's a sink
	static std::tuple<std::vector<int>, std::vector<int> > find_sources_and_sinks(std::vector<vector<int> > &adjacencies){
		std::vector<int> new_sources;
		std::vector<int> new_sinks;

		int no_of_exons = adjacencies.size();
		bool* source = new bool[no_of_exons];
		bool* sink = new bool[no_of_exons];

		for(int i=0;i<no_of_exons;i++){
			source[i] = true;
			sink[i] = true;
		}

		for(unsigned i=0;i<adjacencies.size();i++) {
		// If it has outneighbors, it's not sink
		if(adjacencies.at(i).size() > 0)
			sink[i] = false;

		}
		// If someone points to it (inneighbor), it's not a source
		for(unsigned i=0;i<adjacencies.size();i++) {
			for(unsigned j=0;j<adjacencies.at(i).size();j++)
				source[adjacencies.at(i).at(j)] = false;
		}

		for(int i=0;i<no_of_exons;i++) {
			if(source[i]) 
				new_sources.push_back(i);
			if(sink[i])
				new_sinks.push_back(i);
		}

		delete[] source;
		delete[] sink;

		return std::make_tuple(new_sources, new_sinks);
	}


	/* Find subpath constraints created by reads spanning more than two exons */
	// TODO use ranges for this? Would it be faster?
	static std::tuple<std::vector<std::string>,std::vector<double> > find_subpath_constraints(std::vector<Exon> &exon_list, std::vector<vector<int> > adjacency_list, SamReader* reader){

		std::vector<string> subpath_list;
		std::vector<double> subpath_coverages;

		std::string chrom = exon_list.at(0).chrom;

		std::stringstream sstm;

		// Select all the reads between start of first exon and end of last exon
		std::vector<SamAlignment> alignments = reader->readAlignmentsRegion(chrom, exon_list.at(0).start-1, exon_list.at(exon_list.size()-1).end);

		// Check all the reads in this area
		for(unsigned al=0;al<alignments.size();al++) {
			SamAlignment alignment = alignments.at(al);
			std::vector<int> covered_exons;

			long pos = alignment.getReferenceStart();
			std::vector<long> positions;

			std::vector<std::pair<int,char> > cigar = alignment.getCigarTuples();
			for(unsigned i=0;i<cigar.size();i++) {

				if(cigar.at(i).second == 'M' || cigar.at(i).second == 'X' || cigar.at(i).second == '=' || cigar.at(i).second == 'D') {
					for(int j=0;j<cigar.at(i).first;j++)
						positions.push_back(pos+j);
				}
				// If not insert or soft/hard clip, update pos
				if(cigar.at(i).second != 'I' && cigar.at(i).second != 'S' && cigar.at(i).second != 'H')
				pos += cigar.at(i).first;
			}

			for(unsigned i=0;i<exon_list.size();i++) {
				// To get viable subpath constraints covering at least 2 exons, either start or end of the exon or both have to be covered
					if(find(positions.begin(), positions.end(), exon_list.at(i).start) != positions.end() || find(positions.begin(), positions.end(), exon_list.at(i).end) != positions.end()) {
						covered_exons.push_back(i);
					}
			}
			if(covered_exons.size() >= 2) {
				for(unsigned i=0;i<covered_exons.size()-1;i++)
					sstm << covered_exons.at(i) << " ";
				sstm << covered_exons.at(covered_exons.size()-1);
			}
			else
				continue;
			if(subpath_list.size() == 0) {
				subpath_list.push_back(sstm.str());
				subpath_coverages.push_back(1.0);
			}
			else {
				for(unsigned i=0;i<subpath_list.size();i++) {
					if(subpath_list.at(i) == sstm.str()) {
						subpath_coverages.at(i)++;
						break;
					}
					// Go here if we're at last element, and the string hasn't been found
					else if(i == subpath_list.size()-1) {
						subpath_list.push_back(sstm.str());
						subpath_coverages.push_back(1.0);
						// Need to break here, otherwise it'll read this again
						break;
					}
				}
			}
			sstm.str("");
		}

		// Normalize by the length
		for(unsigned i=0;i<subpath_list.size();i++) {
			long total_length = 0;
			std::vector<string> subpath;
			split(subpath, subpath_list.at(i), is_any_of(MPCUtil::SPACE_SEPARATOR));
			for(unsigned j=0;j<subpath.size();j++) {
				total_length += (exon_list.at(MPCUtil::getIntValue(subpath.at(j))).end-exon_list.at(MPCUtil::getIntValue(subpath.at(j))).start + 1);
			}
			subpath_coverages.at(i) = subpath_coverages.at(i)/total_length;
		}

		alignments.clear();

		return make_tuple(subpath_list, subpath_coverages);
	}


	// Add subpath information to the graph file
	static void add_subpath_information(std::vector<std::string> &subpath_list, std::vector<double> &subpath_coverages, std::string output) {
		std::ofstream handle(output.c_str(), ios::app);
		// If there are no constraints, write 0
		if(subpath_list.size() == 0)
			handle << "0" << std::endl;
		else {
			handle << subpath_list.size() << std::endl;
			for(unsigned i=0;i<subpath_list.size();i++)
				handle << subpath_list.at(i) << std::endl;
			for(unsigned i=0;i<subpath_coverages.size();i++)
				handle << subpath_coverages.at(i) << std::endl;
		}
		handle.close();
	}



public:

	/* Caller gives chromosome, start and end within which to search */
	static int create_graph_with_subpath_constraints_without_annotation(SamReader* reader, std::string chrom, long start, long end, std::string graph, std::string nodes, double threshold, bool debug) {

		std::vector<Exon> exon_list;
		std::vector<int> sources;
		std::vector<int> sinks;

		std::vector<AlignmentRange> alignment_ranges_by_start;

		alignment_ranges_by_start = parse_ranges(reader, chrom, start, end);

		if(alignment_ranges_by_start.size() == 0) {
			std::cerr << "Fatal error, failed to parse any contiguous blocks from the alignments. Exiting." << std::endl;
			exit(1);
		}


		sort(alignment_ranges_by_start.begin(), alignment_ranges_by_start.end(), less_than_start());

		// Copy
		vector<AlignmentRange> alignment_ranges_by_end = alignment_ranges_by_start;

		sort(alignment_ranges_by_end.begin(), alignment_ranges_by_end.end(), less_than_end());

		exon_list = find_exons(alignment_ranges_by_start, alignment_ranges_by_end, chrom, start, end, threshold);

		sort(exon_list.begin(), exon_list.end());
		if(exon_list.size() == 0) {
			// Open the files just to create them
			ofstream handle(graph.c_str());
			handle.close();
			handle.open(nodes.c_str());
			handle.close();
			return 0;
		}
		std::vector<double> exon_weights;
		std::vector<std::vector<int> > adjacency_list;
		std::vector<std::vector<long> > edge_weights;
		tie(exon_weights, adjacency_list, edge_weights) = parse_sam(alignment_ranges_by_start, alignment_ranges_by_end, exon_list);
		if(debug) {
			std::cerr << "Exon list:" << std::endl;
			for(unsigned i=0;i<exon_list.size();i++)
				std::cerr << exon_list.at(i).toString() << ", ";
			std::cerr << std::endl;
			std::cerr << "Exon weights:" << std::endl;
			for(unsigned i=0;i<exon_weights.size();i++)
				std::cerr << exon_weights.at(i) << ", ";
			std::cerr << std::endl;
			std::cerr << "Adjacency list:" << std::endl;
			for(unsigned i=0;i<adjacency_list.size();i++) {
				std::cerr << "[";
				for(unsigned j=0;j<adjacency_list.at(i).size();j++) {
					std::cerr << adjacency_list.at(i).at(j) << ", ";
				}
				std::cerr << "], ";
			}
			std::cerr << std::endl;
			for(unsigned i=0;i<edge_weights.size();i++) {
				std::cerr << "[";
				for(unsigned j=0;j<edge_weights.at(i).size();j++) {
					std::cerr << edge_weights.at(i).at(j) << ", ";
				}
				std::cerr << "], ";
			}
			std::cerr << std::endl;

		}

		tie(sources, sinks) = find_sources_and_sinks(adjacency_list);

		if(debug){
			std::cerr << "Sources:" << std::endl;
			for(unsigned i=0;i<sources.size();i++)
				std::cerr << sources.at(i) << ", ";
			std::cerr << std::endl;
			std::cerr << "Sinks:" << std::endl;
			for(unsigned i=0;i<sinks.size();i++)
				std::cerr << sinks.at(i) << ", ";
			std::cerr << std::endl;
		}

		print_graph(exon_list, exon_weights, adjacency_list, edge_weights, sources, sinks, graph);
		print_nodes(exon_list, nodes);
		if(exon_list.size() > 1) {
			std::vector<std::string> subpath_list;
			std::vector<double> subpath_coverages;
			// TODO: Would this be faster to use ranges? It is just one pass over the sam file
			tie(subpath_list, subpath_coverages) = find_subpath_constraints(exon_list, adjacency_list, reader);
			if(debug) {
				std::cerr << "Subpaths:" << std::endl;
				for(unsigned i=0;i<subpath_list.size();i++)
					std::cerr << subpath_list.at(i) << std::endl;
				std::cerr << "Subpath coverages:" << std::endl;
				for(unsigned i=0;i<subpath_coverages.size();i++)
					std::cerr << subpath_coverages.at(i) << std::endl;
			}
			add_subpath_information(subpath_list, subpath_coverages, graph);
		}

		alignment_ranges_by_start.clear();
		alignment_ranges_by_end.clear();

		return 1;
	}

	/* Caller gives chromosome, start and end within which to search */
	/* Difference to above is that this doesn't add subpath constraints */
	static int create_graph_without_annotation(SamReader* reader, std::string chrom, long start, long end, std::string graph, std::string nodes, double threshold, bool debug) {

		std::vector<Exon> exon_list;
		std::vector<int> sources;
		std::vector<int> sinks;

		std::vector<AlignmentRange> alignment_ranges_by_start;

		alignment_ranges_by_start = parse_ranges(reader, chrom, start, end);

		if(alignment_ranges_by_start.size() == 0) {
			std::cerr << "Fatal error, failed to parse any contiguous blocks from the alignments. Exiting." << std::endl;
			exit(1);
		}


		sort(alignment_ranges_by_start.begin(), alignment_ranges_by_start.end(), less_than_start());

		// Copy
		vector<AlignmentRange> alignment_ranges_by_end = alignment_ranges_by_start;

		sort(alignment_ranges_by_end.begin(), alignment_ranges_by_end.end(), less_than_end());

		exon_list = find_exons(alignment_ranges_by_start, alignment_ranges_by_end, chrom, start, end, threshold);

		sort(exon_list.begin(), exon_list.end());
		if(exon_list.size() == 0) {
			// Open the files just to create them
			ofstream handle(graph.c_str());
			handle.close();
			handle.open(nodes.c_str());
			handle.close();
			return 0;
		}
		std::vector<double> exon_weights;
		std::vector<std::vector<int> > adjacency_list;
		std::vector<std::vector<long> > edge_weights;
		tie(exon_weights, adjacency_list, edge_weights) = parse_sam(alignment_ranges_by_start, alignment_ranges_by_end, exon_list);
		if(debug) {
			std::cerr << "Exon list:" << std::endl;
			for(unsigned i=0;i<exon_list.size();i++)
				std::cerr << exon_list.at(i).toString() << ", ";
			std::cerr << std::endl;
			std::cerr << "Exon weights:" << std::endl;
			for(unsigned i=0;i<exon_weights.size();i++)
				std::cerr << exon_weights.at(i) << ", ";
			std::cerr << std::endl;
			std::cerr << "Adjacency list:" << std::endl;
			for(unsigned i=0;i<adjacency_list.size();i++) {
				std::cerr << "[";
				for(unsigned j=0;j<adjacency_list.at(i).size();j++) {
					std::cerr << adjacency_list.at(i).at(j) << ", ";
				}
				std::cerr << "], ";
			}
			std::cerr << std::endl;
			for(unsigned i=0;i<edge_weights.size();i++) {
				std::cerr << "[";
				for(unsigned j=0;j<edge_weights.at(i).size();j++) {
					std::cerr << edge_weights.at(i).at(j) << ", ";
				}
				std::cerr << "], ";
			}
			std::cerr << std::endl;

		}

		tie(sources, sinks) = find_sources_and_sinks(adjacency_list);

		if(debug){
			std::cerr << "Sources:" << std::endl;
			for(unsigned i=0;i<sources.size();i++)
				std::cerr << sources.at(i) << ", ";
			std::cerr << std::endl;
			std::cerr << "Sinks:" << std::endl;
			for(unsigned i=0;i<sinks.size();i++)
				std::cerr << sinks.at(i) << ", ";
			std::cerr << std::endl;
		}

		print_graph(exon_list, exon_weights, adjacency_list, edge_weights, sources, sinks, graph);
		print_nodes(exon_list, nodes);

		alignment_ranges_by_start.clear();
		alignment_ranges_by_end.clear();

		return 1;
	}



};
#endif /* CREATEGRAPH_H */
