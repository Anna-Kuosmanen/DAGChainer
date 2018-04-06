/**
 * Call for one path file at a time, results are returned as a vector
 *
 * Path files of format:
 * Path 1
 * Path 2
 * ...
 *
 * Author: aekuosma
 *
**/
#include <iostream>
#include <tuple>

#include "MPCUtil.h"
#include "../../SamReader.h"

class convertPaths {

public:

	// samfile required to extract the strand
	static std::vector<std::string> convert_path_no_annotation(std::string nodesfile, std::string pathsfile, int gene_id, std::string samfile, bool debug) {
		if(debug) {
			std::cerr << "Nodes file: " << nodesfile << std::endl;
			std::cerr << "Paths file: " << pathsfile << std::endl;
		}

		std::vector<std::string> transcripts;

		// Test for the existence of the files
		if (FILE *file = fopen(pathsfile.c_str(), "r")) {
			fclose(file);
		} else {
			std::cerr << "convertPaths internal error: Paths file " << pathsfile << " doesn't exist." << endl << "Some transcripts might be missing from the results file." << std::endl;
			return transcripts;
		}
		if (FILE *file = fopen(nodesfile.c_str(), "r")) {
			fclose(file);
		} else {
			std::cerr << "convertPaths internal error: Nodes file " << nodesfile << " doesn't exist." << endl << "Some transcripts might be missing from the results file." << std::endl;
			return transcripts;
		}
		if (FILE *file = fopen(samfile.c_str(), "r")) {
			fclose(file);
		} else {
			std::cerr << "convertPaths: SAM file doesn't exist." << endl << "Exiting." << std::endl;
			exit(2);
		}

		std::ifstream in_nodes(nodesfile.c_str());

		std::ifstream in_paths(pathsfile.c_str());

		std::vector<std::tuple<std::string, int, int> > nodes;

		std::string line;

		// Read the nodes
		while(true) {
			getline(in_nodes, line);

			if(line == "")
				break;

			std::vector<std::string> parts;

			split(parts, line, is_any_of(MPCUtil::TAB_SEPARATOR));

			nodes.push_back(make_tuple(parts.at(0), MPCUtil::getIntValue(parts.at(1)), MPCUtil::getIntValue(parts.at(2))));

		}

		in_nodes.close();

		std::vector<std::vector<int> > paths;

		// Read the paths
		while(true) {
			getline(in_paths, line);
			if(line == "")
				break;
			std::vector<std::string> parts;

			split(parts, line, is_any_of(MPCUtil::SPACE_SEPARATOR));

			// Get rid of the empty entry caused by the space at the end of path
			if(parts.at(parts.size()-1) == "")
				parts.pop_back();

			std::vector<int> nodes;

			for(unsigned i=0;i<parts.size();i++)
				nodes.push_back(MPCUtil::getIntValue(parts.at(i)));

			paths.push_back(nodes);

		}	

		in_paths.close();
	    
		if(paths.size() == 0) {
			return transcripts;
		}

		// Sort the vector of paths that the transcripts are in order
		sort(paths.begin(), paths.end());
	    
		int transcript_count = 1;
	    
		SamReader* sam_handle = new SamReader();

		sam_handle->Open(samfile);

		std::stringstream sstm;

		// For every path, convert nodes into exons
		for(unsigned i=0;i<paths.size();i++) {

			std::string source = "Traphlor";
			std::string frame = ".";
			int score = 1;
			char strand = '.';

			bool forward = false;
			bool reverse = false;

			// Iterate over the reads crossing on the path to check for strand tags
			for(unsigned j=0;j<paths.at(i).size();j++) {
				std::tuple<std::string, int, int> selected_exon = nodes[paths.at(i).at(j)];
				std::vector<SamAlignment> alignments = sam_handle->readAlignmentsRegion(get<0>(selected_exon), get<1>(selected_exon), get<2>(selected_exon));

				for(unsigned al=0;al< alignments.size();al++){
					SamAlignment alignment = alignments.at(al);
					if(alignment.getCigarTuples().size() == 1)
						continue;

					if(alignment.hasTag("XS:A")) {
						std::string strand_tag = alignment.getValueForTag("XS:A");
						if(strand_tag == "+")
							forward = true;
						else if(strand_tag == "-")
							reverse = true;
					}
				}

				alignments.clear();
			}


			if(forward && !reverse)
				strand = '+';
			else if(!forward && reverse)
				strand = '-';

			std::vector<std::string> exons;
			int last_end = -1;

			std::string seqname;
			int transcript_start = -1;
			int transcript_end = -1;

			for(unsigned j=0;j<paths.at(i).size();j++) {
				std::tuple<std::string, int, int> current_node = nodes[paths.at(i).at(j)];
				// Get the reference and transcript start from first node
				if(j==0) {
					seqname = get<0>(current_node);
					transcript_start = get<1>(current_node);
				}
				if(j == paths.at(i).size()-1)
					transcript_end = get<2>(current_node);

				int start = get<1>(current_node);
				int end = get<2>(current_node);

				// split exon, merge this and the previous one
				if(last_end != -1 && last_end == start - 1) {
					std::vector<std::string> last_exon;
					split(last_exon, exons.at(exons.size()-1), is_any_of(MPCUtil::TAB_SEPARATOR));
					exons.pop_back();
					sstm << last_exon.at(0) << "\t" << last_exon.at(1) << "\t" << last_exon.at(2) << "\t" << last_exon.at(3) << "\t" << end << "\t"
						<< last_exon.at(5) << "\t" << last_exon.at(6) << "\t" << last_exon.at(7);
					exons.push_back(sstm.str());
					sstm.str("");
				}
				else {
					sstm << seqname << "\t" << source << "\texon\t" << start << "\t" << end << "\t" << score << "\t" << strand << "\t" << frame;
					exons.push_back(sstm.str());
					sstm.str("");
				}
				last_end = end;
			}

			sstm << seqname << "\t" << source << "\ttranscript\t" << transcript_start << "\t" << transcript_end << "\t" << score << "\t" << strand << "\t" 
				<< frame << "\tgene_id \"Traph." << gene_id << "\"; transcript_id \"Traph." << gene_id << "." << transcript_count << "\"" << endl;

			transcripts.push_back(sstm.str());
			sstm.str("");

			for(unsigned j=0;j<exons.size();j++) {
				sstm << exons.at(j) << "\tgene_id \"Traph." << gene_id << "\"; transcript_id \"Traph." << gene_id << "." << transcript_count << "\"; exon_number \""<< j+1 << "\"" << endl;
				transcripts.push_back(sstm.str());
				sstm.str("");
			}

			transcript_count++;
		}

		in_paths.close();
		sam_handle->Close();
		delete sam_handle;
		return transcripts;
	}

};
