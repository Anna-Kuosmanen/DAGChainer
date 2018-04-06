/**
** Does all the stuff required to run Traphlor.
*
* Author:
* aekuosma
*
*/
#include <iostream>
#include <getopt.h>

#include "createGraphWithSubpathConstraints.h"
#include "convertPaths.h"
#include "MPCGraph.h"
#include "MPCUtil.h"
#include "MPCPrintUtil.h"

#ifdef PARALLEL_SUPPORT
#include <omp.h>
#endif



// Parameters for parsing command line
//enum parameter_t {long_opt_alt_threshold};

std::string GRAPH_FILE_PREFIX = "gene.graph";
std::string NODES_FILE_PREFIX = "gene.nodes";
std::string PATHS_FILE_SUFFIX = ".sol";
//string ANNOTATION_PATHS_PREFIX = "annotation.paths";
std::string TRANSCRIPT_FILE = "transcripts.gtf";

void print_help(char const *name) {

	std::cerr << "usage: " << name << " [options]" << std::endl
		<< std::endl
		<< "Required options:" << std::endl
		<< "    -i INPUT, --input INPUT        Input file" << std::endl
		<< "    -o OUTPUT, --output OUTPUT     Output directory" << std::endl << std::endl
		<< "Other options:" << std::endl
//         << "    --alternative-threshold THRESHOLD  " << endl
//         << "          Slope threshold for searching for alternative" << std::endl
//         << "          transcripts starts/ends. The smaller, the more" << std::endl
//         << "          sensitive, but too small value can cause false" << std::endl
//         << "          positives if coverage is very uneven. Default 0.3" << std::endl
		<< "    -d, --debug                    Debug mode" << std::endl;
#ifdef PARALLEL_SUPPORT
	std::cerr << "    -P <int>, --parallel <int>     Number of parallel threads to use (default is one," << std::endl
		 << "                                   give argument 0 to use all available cores)." << std::endl;
#else
	cerr << "    -P <int>, --parallel <int>     Parallel computation requires OpenMP, see README." << endl;
#endif
	std::cerr << std::endl;


}


// Create graphs for one area at a time
// Returns the number of graphs created
int create_graphs_no_annotation(std::string ifile, std::string odir, double alt_prime_threshold, bool debug, int parallel) { 

	SamReader* alignment_handle = new SamReader();
	alignment_handle->Open(ifile);
	if(!alignment_handle->isOpen())
		cerr << "Failure to open SAM file for counting ranges" << endl;

	std::vector<std::string> references = alignment_handle->getReferences();
    
	std::vector<std::tuple<std::string, long, long> > ranges;

	// Collect the ranges that can do a parallel for loop on them
	for(unsigned i=0;i<references.size();i++) {
 
		// Saves the last currently read end position
		// Since the alignments are sorted by position, if the start of next alignment is farther than the farthest end, then there's a gap
		// -> process this area
		long last_end = 0;
		long start = 0;

		std::vector<SamAlignment> alignments = alignment_handle->readAlignmentsRegion(references.at(i));

		for(unsigned al=0;al<alignments.size();al++) {

			SamAlignment alignment = alignments.at(al);

			if(alignment.isSecondary())
				continue;

			if(alignment.getMappingQuality() == 0)
				continue;

			if(last_end == 0) {
				last_end = alignment.getReferenceEnd();
				start = alignment.getReferenceStart();;
				continue;
			}
			if(alignment.getReferenceStart() > last_end) {
				ranges.push_back(make_tuple(references.at(i), start, last_end));
				start = alignment.getReferenceStart();
				last_end = alignment.getReferenceEnd();
			}
			else if(alignment.getReferenceEnd() > last_end) {
				last_end = alignment.getReferenceEnd();
			}

		}
		if(start != 0)
			ranges.push_back(make_tuple(references.at(i), start, last_end));

		alignments.clear();

	}
	alignment_handle->Close();

	delete alignment_handle;

#ifdef PARALLEL_SUPPORT
if(parallel != 0)
	omp_set_num_threads(parallel);
#pragma omp parallel
#endif
{
	// Different reader for each thread
	SamReader* reader = new SamReader();
	reader->Open(ifile);
	if(!reader->isOpen())
		std::cerr << "Failure to open SAM file" << std::endl;

#ifdef PARALLEL_SUPPORT
#pragma omp for
#endif
	for(unsigned file_tally=0;file_tally<ranges.size();file_tally++) {
		std::stringstream sstm;

		std::string reference;
		long start;
		long end;
		tie(reference, start, end) = ranges.at(file_tally);
		sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << (file_tally + 1);
		std::string graphfile = sstm.str();
		sstm.str("");
		sstm << odir << "/tmp/" << NODES_FILE_PREFIX << "_" << (file_tally + 1);
		std::string nodesfile = sstm.str();
		sstm.str("");

		createGraphWithSubpathConstraints::create_graph_with_subpath_constraints_without_annotation(reader, reference, start, end, graphfile, nodesfile, alt_prime_threshold, debug);
	}
	reader->Close();
	delete reader;

} // end pragma omp parallel
	return int(ranges.size());
}

// Calculates the average coverage over all the graphs
// Used for finding the threshold on which graphs to discard
double find_average_coverage(std::string odir, int no_of_files) {
	std::vector<double> coverages;

	std::stringstream sstm;

	for(int i=1;i<=no_of_files;i++) {
		sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << i;
		std::ifstream handle(sstm.str().c_str());
		sstm.str("");
		std::string line;

		getline(handle, line);
		int node_number = MPCUtil::getIntValue(line);
		for(int j=0;j<node_number;j++)
			getline(handle, line);
		getline(handle, line);
		std::vector<std::string> exon_covs;
		split(exon_covs, line, is_any_of(MPCUtil::SPACE_SEPARATOR));
		for(unsigned j=0;j<exon_covs.size();j++)
			coverages.push_back(MPCUtil::getDoubleValue(exon_covs.at(j)));

		handle.close();
	}

	double cov_sum = 0.0;

	for(unsigned i=0;i<coverages.size();i++)
		cov_sum += coverages.at(i);

	return cov_sum/coverages.size();
}

// Check that there's enough evidence that this graph isn't just created from mismappings
bool check_graph(std::string graph_file, std::string nodes_file, double avg_cov) {

	std::ifstream handle(graph_file.c_str());

	std::string line;

	getline(handle, line);

	int node_number = MPCUtil::getIntValue(line);

	double max_coverage = 0;

	// Read off the adjacencies
	for(int i=0;i<node_number;i++)
		getline(handle,line);

	std::vector<std::string> coverages;

	getline(handle, line);

	split(coverages, line, is_any_of(MPCUtil::SPACE_SEPARATOR));

	for(unsigned i=0;i<coverages.size();i++)
		if(MPCUtil::getDoubleValue(coverages.at(i)) > max_coverage)
			max_coverage = MPCUtil::getDoubleValue(coverages.at(i));

	handle.close();
	int max_length = 0;
	handle.open(nodes_file.c_str());
	while(true) {
		getline(handle, line);
		if(line == "")
			break;
		std::vector<std::string> parts;
		split(parts, line, is_any_of(MPCUtil::TAB_SEPARATOR));
		max_length += (MPCUtil::getIntValue(parts.at(2)) - MPCUtil::getIntValue(parts.at(1)));

	}

	handle.close();

	// If coverage of a graph is less than 1, it consists of a single read
	if(max_coverage > 1.0*avg_cov/10 and max_coverage > 1)
		return true;
/*    if(node_number >= 3)
		return true;
	if(max_length > 1000)
		return true;
*/
	return false;
}


// Called when there isn't enough evidence that the graph isn't just created from mismappings
// Writes empty solution file
void write_empty_solution_file(string odir, int i) {

	std::stringstream sstm;
	sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << i << ".sol";
	std::ofstream handle(sstm.str().c_str());
	handle.close();

}

// Run minimum path cover algorithm to find paths
void run_mpc(string odir, int no_of_files, int parallel) {

	double avg_cov = find_average_coverage(odir, no_of_files);

#ifdef PARALLEL_SUPPORT
if(parallel != 0)
	omp_set_num_threads(parallel);
#pragma omp parallel
#endif
{

	std::stringstream sstm;

#ifdef PARALLEL_SUPPORT
#pragma omp for
#endif
	for(int i=1;i<=no_of_files;i++) {
		sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << i;
		std::string graphFile=sstm.str();
		sstm.str("");
		sstm << odir << "/tmp/" << NODES_FILE_PREFIX << "_" << i;
		std::string nodeFile= sstm.str();
		sstm.str("");

		if(check_graph(graphFile, nodeFile, avg_cov)) {

			// Create graph for LEMON
			MPCGraph g=MPCGraph::createMPCGraph(graphFile, nodeFile);
			MPCPrintUtil::printGraph(g);
			// Solve the flow
			std::string solution=g.solve();
			// Print
			MPCPrintUtil::writeSolutionPath(solution, graphFile);

		}
		else
			write_empty_solution_file(odir, i);
		}
	}
}         


// Convert the paths into transcripts
void convert_paths_no_annotation(std::string odir, int no_of_files, std::string samfile, bool debug, int parallel) {

	// Map where gene_id is the key, forces the order of the transcripts for printing
	std::map<int, std::vector<std::string> > all_transcripts;

	int gene_id = 1;
	std::vector<int> gene_ids;

	// Decide gene IDs in advance, as otherwise in parallel processing they'll be in random order
	for(int i=1;i<=no_of_files;i++) {
		// Check the solution file isn't empty 
		bool isEmpty = false;
		std::ifstream temp_handle((odir + "/tmp/" + GRAPH_FILE_PREFIX + "_" + MPCUtil::getStringValue(i) + ".sol").c_str());
		std::string line;

		getline(temp_handle, line);
		if(line == "") {
			isEmpty = true;
		}

		temp_handle.close();
		if(!isEmpty) {
			gene_ids.push_back(gene_id);
			gene_id++;
		}
		else
			gene_ids.push_back(-1);
	}

#ifdef PARALLEL_SUPPORT
if(parallel != 0)
	omp_set_num_threads(parallel);
#pragma omp parallel
#endif
{

	std::stringstream sstm;

	std::string nodes_file;
	std::string solution_file;

#ifdef PARALLEL_SUPPORT
#pragma omp for
#endif

	for(int i=1;i<=no_of_files;i++) {

		sstm << odir << "/tmp/" << NODES_FILE_PREFIX << "_" << i;
		nodes_file = sstm.str();
		sstm.str("");
		sstm << odir << "/tmp/" << GRAPH_FILE_PREFIX << "_" << i << ".sol";
		solution_file = sstm.str();
		sstm.str("");

		if(gene_ids.at(i-1) != -1) {
			std::vector<std::string> transcripts = convertPaths::convert_path_no_annotation(nodes_file, solution_file, gene_ids.at(i-1), samfile, debug);

			all_transcripts[gene_ids.at(i-1)] = transcripts;
		}
        
	}
} //end pragma omp parallel

	// Print
	std::ofstream out_handle((odir + "/" + TRANSCRIPT_FILE).c_str());

	for(std::map<int, std::vector<std::string> >::iterator it = all_transcripts.begin(); it != all_transcripts.end();it++) {
		// This actually contains lines of transcripts and exons. Contains end lines.
		std::vector<std::string> transcripts = it->second;
		for(unsigned i=0;i<transcripts.size();i++)
			out_handle << transcripts.at(i);
	}

	out_handle.close();

}


int main(int argc, char **argv) {
 
	static struct option long_options[] =
	{
		{"input", required_argument, 0, 'i'},
		{"output", required_argument, 0, 'o'},
		//{"alternative-threshold", required_argument, 0, long_opt_alt_threshold},
		{"debug",     no_argument,       0, 'd'},
		{"parallel", required_argument, 0, 'P'},
		{0, 0, 0, 0}
	};

	std::string odir = "";
	std::string ifile = "";
	bool debug = false;
	double alt_threshold = 0.3;
	unsigned parallel = 1;


	int option_index = 0;
	int c;
	while ((c = getopt_long(argc, argv, "i:o:dP:", long_options, &option_index)) != -1){
		switch(c){
			case 'i':
				ifile = string(optarg); break;
			case 'o':
				odir = string(optarg); break;
			case 'd':
				debug = true; break;
//			case long_opt_alt_threshold:
//				alt_threshold = MPCUtil::getDoubleValue(string(optarg)); break;
			case 'P':
#ifdef PARALLEL_SUPPORT
				parallel = atoi(optarg); break;
#else
				std::cerr << "runTraphlor: Parallel processing not currently available!" << std::endl
				<< "Please recompile with parallel support; see README for more information." << std::endl;
				return 1;
#endif

		}
	}

	if(odir == "" || ifile == "") {
		print_help(argv[0]);
		return 1;
	}

	// Test for the existence of the input file
	if (FILE *file = fopen(ifile.c_str(), "r")) {
		fclose(file);
	} else {
		std::cerr << "Input file doesn't exist. Exiting." << std::endl;
		return 1;
	}

	// Test opening SAM file
	SamReader* reader = new SamReader();
	reader->Open(ifile);
	if(!reader->isOpen()) {
		std::cerr << "Failed to open SamReader. Exiting." << std::endl;
		return 1;
	}

	reader->Close();


	// Create folders
	bool dircreated = system(("mkdir " + odir).c_str());

	if(!dircreated) {
		std::cerr << "Failure to create the output folder." << std::endl;
	}

	dircreated = system(("mkdir " + odir + "/tmp/").c_str());

	if(!dircreated) {
		std::cerr << "Failure to create the output folder." << std::endl;
	}

	std::cout << "Creating graphs..." << std::endl;
	int no_of_files = create_graphs_no_annotation(ifile, odir, alt_threshold, debug, parallel);
	if(no_of_files == 0) {
		std::cerr << "Fatal error occurred, no graph files were created. Exiting." << std::endl;
		return 1;
	}
	std::cout << "Searching for paths..." << std::endl;
	run_mpc(odir, no_of_files, parallel);
	std::cout << "Converting paths..." << std::endl;
	convert_paths_no_annotation(odir, no_of_files, ifile, debug, parallel);

	if(!debug) {
		std::string cmd = "rm -r " + odir + "/tmp/";
		system(cmd.c_str());
	}

	std::cout << "Running Traphlor completed." << std::endl;

	delete reader;

}


