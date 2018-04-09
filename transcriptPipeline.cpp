/*
 * Transcript prediction pipeline that takes as input:
 *   - SAM file of alignments for short reads
 *   - FASTA file of long reads
 *   (- FASTA file of the genome to extract the splicing graph sequence)
 *
 * The pipeline
 *   - creates a splicing graph from the short read alignments
 *   - aligns long reads to the splicing graph using colinear chaining
 *   - predicts transcripts using Traphlor
 *
 * Created: Feb 27, 2018
 * Author: aekuosma
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <getopt.h>

#include "SamReader.h"
#include "GenomeReader.h"
#include "FastaReader.h"
#include "SequenceGraph.h"
#include "ColinearSolver.h"
#include "utils.h"

#include "traphlor/RNA_MPC_SC/createGraphWithSubpathConstraints.h"
#include "traphlor/RNA_MPC_SC/convertPaths.h"
#include "traphlor/RNA_MPC_SC/MPCGraph.h"
#include "traphlor/RNA_MPC_SC/MPCPrintUtil.h"



double ALT_PRIME_THRESHOLD = 0.3;
std::string GRAPH_FILE_PREFIX = "gene.graph";
std::string NODES_FILE_PREFIX = "gene.nodes";
std::string GFA_FILE_SUFFIX = ".gfa";
std::string PATHS_FILE_SUFFIX = ".sol";
std::string TRANSCRIPT_FILE = "transcripts.gtf";

std::string SAMTOOLS = "/cs/work/home/aekuosma/software/samtools/samtools";
std::string VG = "/cs/work/home/aekuosma/software/vg/bin/vg";

void print_help(char const *name) {

	std::cerr << "usage: " << name << " [options]" << std::endl 
		<< std::endl;

	std::cerr << "Required options:" << std::endl
		<< "    -s SHORT --short SHORT        Short reads alignment (SAM) file." << std::endl
		<< "    -l LONG --long LONG           Long reads (FASTA) file." << std::endl
		<< "    -g GENOME --genome GENOME     Genome (FASTA) file." << std::endl
		<< "    -o OUTPUT --output OUTPUT     Output directory." << std::endl
		<< "Other options:" << std::endl
		<< "    -m SEEDLEN --minimum SEEDLEN  Minimum seed length (default 5)." << std::endl
		<< "    -t INT --stringency INT       How strict the subpath reporting is. Higher values allow for more distant" << std::endl
		<< "                                  anchors to form a chain (Range: 0-5. Default: 0). If using very high minimum" << std::endl
		<< "                                  seed length, it is adviced to adjust this parameter, as no anchors might map to smaller exons." << std::endl
		<< "    -d  --debug                   Debug mode on." << std::endl;


	std::cerr << std::endl;
}


void process_region(SamReader* samreader, GenomeReader* genomereader, std::string reference, long start, long end, int file_tally, int mem_length_threshold, bool debug, std::string long_reads_file, std::string odir, int stringency) {

	std::stringstream sstm;

	// Create the splicing graph
	sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << (file_tally + 1);
	std::string graphFile = sstm.str();
	sstm.str("");
	sstm << odir << "/tmp/" << NODES_FILE_PREFIX << "_" << (file_tally + 1);
	std::string nodesFile = sstm.str();
	sstm.str("");

	createGraphWithSubpathConstraints::create_graph_without_annotation(samreader, reference, start, end, graphFile, nodesFile, ALT_PRIME_THRESHOLD, debug);

	// Write the nodes file with sequence instead of coordinates for SequenceGraph
	sstm << odir << "/tmp/" << NODES_FILE_PREFIX << "_" << (file_tally + 1) << ".fa";

	std::string nodesFaFile = sstm.str();
	sstm.str("");

	std::ifstream nodes_in(nodesFile.c_str());
	std::ofstream nodes_out(nodesFaFile.c_str());

	int nodenumbers = 0;

	std::string line;

	while(getline(nodes_in, line)) {
		if(line == "")
			break;
		// Remember that graph creation writes gtf type
		// -1 coordinates to make 0-based
		std::vector<std::string> parts = split(line, '\t');
		std::string seq = genomereader->readSequenceUpperCase(parts.at(0), atol(parts.at(1).c_str())-1, atol(parts.at(2).c_str())-1);

		nodes_out << nodenumbers << std::endl;
		nodes_out << seq << std::endl;

		nodenumbers++;

	}

	nodes_in.close();

	// Create the SequenceGraph object from the splicing graph
	SequenceGraph* SGraph = new SequenceGraph();

	// Note: reading the splicing graph supports the naive anchor finding
	// which requires as last argument pattern length.
	// Since we're not using that with GCSA2, to save space it can be set as 1
	bool readsuccess = SGraph->readSplicingGraphFile(graphFile, nodesFaFile, 1);

	if(!readsuccess) {
		std::cerr << "Failed to create the sequence graph, check the existance of the graph file and nodes file." << std::endl;
		return;

	}


	ColinearSolver* solver = new ColinearSolver(SGraph);

	// Write the SequenceGraph object as GFA for MEM finding

	sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << (file_tally + 1) << GFA_FILE_SUFFIX;

	std::string gfaout = sstm.str();
	sstm.str("");

	SGraph->outputGFA(gfaout, true);

	// Create the GCSA2 index by calling VG

	// Convert GFA to VG
	sstm << VG << " view --vg --gfa-in " << gfaout << " > " << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << (file_tally +1) << ".vg";

	system(sstm.str().c_str());
	sstm.str("");

	// Index

	sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << (file_tally +1) <<  ".gcsa";

	std::string gcsa_name = sstm.str();

	sstm << ".lcp";
	std::string lcp_name = sstm.str();
	sstm.str("");

	// TODO Should be able to do this without outside calls
	sstm << VG << " index -g " << gcsa_name << " -k 16 " << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << (file_tally +1) << ".vg";

	system(sstm.str().c_str());
	sstm.str("");


	ColinearChain bestchain;

	FastaReader* fastareader = new FastaReader();

	bool opensuccess = fastareader->Open(long_reads_file);

	if(!opensuccess) {
		std::cerr << "Failure to open the long reads file." << std::endl;
		return;
	}


	FastaEntry entry;

	std::vector<std::string> subpath_list;
	std::vector<std::vector<int> > subpath_vector_list;
	std::vector<double> subpath_coverages;

	// For each long read consider both the read and its reverse complement
	while(true) {

		entry = fastareader->next();

		// FastaReader returns empty id if there was nothing to read
		if(entry.id == "")
			break;

		// Reminder to self: The clean-up of old Tuples happens when either solver->solve is called or solver is destroyed
		// Do NOT start cleaning them here too, it'll cause pointer havoc!

		// When the last argument is "true", solver checks both forward and reverse complement
		ColinearChain bestchain = solver->solve(entry, gcsa_name, lcp_name, mem_length_threshold, true);

		// Check the score, it's -1 if no chain was found
		if(bestchain.coverageScore == -1)
			continue;

		// Convert the best chain into a subpath in the splicing graph
		std::vector<int> subpath = SGraph->convertChainToExons(bestchain, stringency);

		// convertChainToExons returns empty vector if it could not form a subpath that respects the structure of the graph with minimal adjustments
		// (See the function in SequenceGraph.ccp for details)
		if(subpath.size() == 0)
			continue;

		subpath_vector_list.push_back(subpath);

/*		if(subpath_list.size() == 0) {
			subpath_list.push_back(subpath);
			subpath_coverages.push_back(1.0);
		}
		else {
			for(unsigned i=0;i<subpath_list.size();i++) {
				if(subpath_list.at(i) == subpath) {
					subpath_coverages.at(i)++;
					break;
				}
				// We're at last element and string hasn't been found
				else if(i == subpath_list.size()-1) {
					subpath_list.push_back(subpath);
					subpath_coverages.push_back(1.0);
					// List size grew, break that the new insert isn't checked
					break;
				}


			}	

		}*/
	}

	// List holds a lot of duplicates, sort, count coverages, and convert the final products into strings
	sort(subpath_vector_list.begin(), subpath_vector_list.end());

	for(unsigned i=0;i<subpath_vector_list.size();i++) {
		std::vector<int> path_vector = subpath_vector_list.at(i);

		// Length 1 is not accepted
		if(path_vector.size() == 1)
			continue;

		for(unsigned j=0;j<path_vector.size()-1;j++)
			sstm << path_vector.at(j) << " ";
		sstm << path_vector.at(path_vector.size()-1);

		std::string path_string = sstm.str();
		sstm.str("");

		if(subpath_list.size() == 0) {
			subpath_list.push_back(path_string);
			subpath_coverages.push_back(1.0);
		}
		else if(subpath_list.at(subpath_list.size()-1) == path_string) {
			subpath_coverages.at(subpath_coverages.size()-1)++;
		}
		else {
			subpath_list.push_back(path_string);
			subpath_coverages.push_back(1.0);		
		}

	}


	// Divide the counts by the lengths

	// First read the nodes file back for the lengths of the exons
	nodes_in.open(nodesFile);

	std::vector<std::pair<long, long> > exon_list;

	while(getline(nodes_in, line)) {

		std::vector<std::string> parts = split(line, '\t');
		exon_list.push_back(std::make_pair(atol(parts.at(1).c_str()),atol(parts.at(2).c_str())));

	}


	for(unsigned i=0;i<subpath_list.size();i++) {

		long total_length = 0;
		std::vector<std::string> subpath = split(subpath_list.at(i), ' ');

		for(unsigned j=0;j<subpath.size();j++) {
			total_length += (exon_list.at(MPCUtil::getIntValue(subpath.at(j))).second-exon_list.at(MPCUtil::getIntValue(subpath.at(j))).first +1);
		}
		subpath_coverages.at(i) = subpath_coverages.at(i)/total_length;
	}


	// Write the subpaths into the original splicing graph file (append)
	std::ofstream output;
	output.open(graphFile.c_str(), std::ofstream::app);
	
	if(!output.is_open()) {
		std::cerr << "Failure to open the graph file for appending." << std::endl;
		return;
	}

	if(subpath_list.size() == 0)
		output << "0" << std::endl;
	else {
		output << subpath_list.size() << std::endl;
		for(unsigned i=0;i<subpath_list.size();i++)
			output << subpath_list.at(i) << std::endl;
		for(unsigned i=0;i<subpath_coverages.size();i++)
			output << subpath_coverages.at(i) << std::endl;
	}


	// Clean up
	fastareader->Close();
	delete fastareader;
	delete SGraph;
	delete solver;

}

void run_traphlor(std::string odir, int no_of_files, std::string samfile, bool debug) {

	std::stringstream sstm;

	// Run the flow engine
	for(int i=1;i<=no_of_files;i++) {

		sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << i;
		std::string graphFile = sstm.str();
		sstm.str("");
		sstm << odir << "/tmp/" << NODES_FILE_PREFIX << "_" << i;
		std::string nodeFile = sstm.str();
		sstm.str("");

		MPCGraph g = MPCGraph::createMPCGraph(graphFile, nodeFile);

		MPCPrintUtil::printGraph(g);

		std::string solution = g.solve();

		MPCPrintUtil::writeSolutionPath(solution, graphFile);
	}

	// Convert paths

	// Decide gene IDs in advance (preparation for parallel processing)

	int gene_id = 1;
	std::vector<int> gene_ids;

	for(int i=1;i<=no_of_files;i++) {
		bool isEmpty = false;
		std::ifstream temp_handle((odir + "/tmp/" + GRAPH_FILE_PREFIX + "_" + MPCUtil::getStringValue(i) + ".sol").c_str());
		std::string line;

		getline(temp_handle, line);
		if(line == "")
			isEmpty = true;

		temp_handle.close();
		if(!isEmpty) {
			gene_ids.push_back(gene_id);
			gene_id++;
		}
		else
			gene_ids.push_back(-1);
	}
	// Again preparation for parallel processing, forces order for print
	std::map<int, std::vector<std::string> > all_transcripts;

	for(int i=1;i<=no_of_files;i++) {
		sstm << odir << "/tmp/" << NODES_FILE_PREFIX << "_" << i;
		std::string nodes_file = sstm.str();
		sstm.str("");
		sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << i << ".sol";
		std::string solution_file = sstm.str();
		sstm.str("");

		if(gene_ids.at(i-1) != -1) {
			std::vector<std::string> transcripts = convertPaths::convert_path_no_annotation(nodes_file, solution_file, gene_ids.at(i-1), samfile, debug);

			all_transcripts[gene_ids.at(i-1)] = transcripts;

		}
	}

	std::ofstream out_handle((odir + "/" + TRANSCRIPT_FILE).c_str());

	for(std::map<int, std::vector<std::string> >::iterator it = all_transcripts.begin(); it != all_transcripts.end(); it++) {
		std::vector<std::string> transcripts = it->second;
		for(unsigned i=0;i<transcripts.size();i++)
			out_handle << transcripts.at(i);

	}

	out_handle.close();

}


int main(int argc, char **argv) {

	static struct option long_options[] =
	{
		{"short", required_argument, 0, 's'},
		{"long", required_argument, 0, 'l'},
		{"genome", required_argument, 0, 'g'},
		{"output", required_argument, 0, 'o'},
		{"minimum", required_argument, 0, 'm'},
		{"stringency", required_argument, 0, 't'},
		{"debug", no_argument, 0, 'd'},
		{0, 0, 0, 0}

	};



	std::string samfile = "";
	std::string fastafile = "";
	std::string genomefile = "";
	std::string odir = "";

	int mem_length_threshold = 5;
	bool debug = false;
	int stringency = 0;

	int option_index = 0;
	int c;
	while((c = getopt_long(argc, argv, "s:l:g:o:m:t:d", long_options, &option_index)) != -1){

		switch(c) {
			case 's':
				samfile = string(optarg); break;
			case 'l':
				fastafile = string(optarg); break;
			case 'g':
				genomefile = string(optarg); break;
			case 'o':
				odir = string(optarg); break;
			case 'm':
				mem_length_threshold = atoi(optarg); break;
			case 't':
				stringency = atoi(optarg); break;
			case 'd':
				debug = true; break;

		}
	}

	if(samfile == "" || fastafile == "" || genomefile == "" || odir == "") {
		print_help(argv[0]);
		return 1;
	}

	if(stringency < 0 || stringency > 5) {
		std::cerr << "Stringency value should be in range 0-5." << std::endl;
		print_help(argv[0]);
		return 1;
	}


	SamReader* sam_handle = new SamReader();

	bool opened = sam_handle->Open(samfile);	

	if(!opened) {
		std::cerr << "Failure to open short reads SAM file." << std::endl;
		return 1;
	}

	GenomeReader* genome_handle = new GenomeReader();

	opened = genome_handle->Open(genomefile);

	if(!opened) {
		std::cerr << "Failure to open the genome file." << std::endl;
		return 1;
	}

	// Test for the existence of the input file
	if (FILE *file = fopen(fastafile.c_str(), "r")) {
		fclose(file);
	} else {
		std::cerr << "Failure to open long reads FASTA file." << std::endl;
		return 1;
	}

	

	// Create the folders
	system(("mkdir " + odir).c_str());
	system(("mkdir " + odir + "/tmp/").c_str());

	std::vector<std::string> references = sam_handle->getReferences();

	std::vector<std::tuple<std::string, long, long> > ranges;

	// Collect the ranges (this is for future parallel support use)
	for(unsigned i=0;i<references.size();i++) {
		std::vector<SamAlignment> alignments = sam_handle->readAlignmentsRegion(references.at(i));

		// Saves the last currently read end position
		// Since the alignments are sorted by position, if the start of next alignment is farther than the farthest end, then there's a gap
		// -> process this area
		long last_end = 0;
		long start = 0;


		for(unsigned al=0;al<alignments.size();al++) {

			SamAlignment alignment = alignments.at(al);

			if(alignment.isSecondary())
				continue;

			if(alignment.getMappingQuality() == 0)
				continue;

			if(last_end == 0) {
				last_end = alignment.getReferenceEnd();
				start = alignment.getReferenceStart();
				continue;
			}
			if(alignment.getReferenceStart() > last_end) {
				ranges.push_back(make_tuple(references.at(i), start, last_end));
				start = alignment.getReferenceStart();
				last_end = alignment.getReferenceEnd();
			}
			else if(alignment.getReferenceEnd()> last_end) {
				last_end = alignment.getReferenceEnd();
			}
		}
		if(start != 0)
			ranges.push_back(make_tuple(references.at(i), start, last_end));
	}
	sam_handle->Close();

	sam_handle->Open(samfile);

	for(unsigned file_tally=0;file_tally<ranges.size();file_tally++) {

		std::string reference;
		long start;
		long end;
		tie(reference, start, end) = ranges.at(file_tally);

		process_region(sam_handle, genome_handle, reference, start, end, file_tally, mem_length_threshold, debug, fastafile, odir, stringency);
	}

	sam_handle->Close();
	delete sam_handle;

	genome_handle->Close();
	delete genome_handle;

	// Traphlor's flow engine and converting the paths to transcripts
	run_traphlor(odir, int(ranges.size()), samfile, debug);


}

