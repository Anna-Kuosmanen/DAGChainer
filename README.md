<h1>DAGChainer</h1>
<p>This is a preliminary implementation of the algorithms proposed in the article ”Using Minimum Path Cover to Boost Dynamic Programming in DAGs: Co-linear Chaining Extended” by A. Kuosmanen, T. Paavilainen, T. Gagie, R. Chikhi, A. Tomescu and V. Mäkinen (presented in <i>RECOMB2018</i>).</p>

<p>ColinearSolver class is the main class implementing the algorithms. SequenceGraph models a graph where each node represents a single base. Currently SequenceGraph class supports only creation from a splicing graph. There are two approaches available for finding maximal exact matches (to be used as anchors for the colinear chain algorithm) between a sequence graph and a sequence: slow naive approach based on graph traversal and fast approach based on GCSA2. Additionally the project includes various utility classes, such as a static Range Maximum Query Tree implementation, SAM format reader and FASTA format reader.</p>

<p>The folder MC-MPC contains a copy of Topi Paavilainen's Master's thesis project[1], for finding the minimum path cover in a graph. The folder traphlor contains a slightly modified copy of Traphlor software[2].</p>

<p>The project includes a pipeline (transcriptPipeline) that takes as input short read alignments and long read sequences, and predicts transcripts from them in the following steps:</p>
<ul>
<li>a splicing graph is created from the short reads</li>
<li>long reads are aligned to the splicing graph using colinear chaining</li>
<li>colinear chains are converted into subpath constraints (chains of exons)</li>
<li>Traphlor is used to predict transcripts based on both short and long reads</li>
</ul>

<h2>Branches</h2>
<p>Branch <code>master</code> implements the algorithms of the article exactly. Branch <code>heuristics</code> is testing ground for various heuristics to improve transcript prediction accuracy.</p>

<h1>Getting started</h1>
<h2>Prerequisites</h2>
<p>DAGChainer uses GCSA2[3], SDSL[4] and LEMON[5] libraries. Also the software VG[6] needs to be in the system path, as it is used to transform the splicing graphs into input format suitable for GCSA2 indexing.</p>

<h2>Building</h2>

<p>Change the variables for the library locations in Makefile to correspond to your installation locations. Then issue ”make” command in the main directory.</p>

<p>Note that ”divsufsort” and ”divsufsort64” libraries should have been installed when installing SDSL libraries, but they might be located in a different folder than SDSL libraries.</p>

<h2>Running the pipeline</h2>

<p>To run the transcript prediction pipeline, issue the command ”transcriptPipeline [options]” with the following options:</p>
<pre>
Required options:
    -s SHORT --short SHORT        Short reads alignment (SAM) file.
    -l LONG --long LONG           Long reads (FASTA) file.
    -g GENOME --genome GENOME     Genome (FASTA) file.
    -o OUTPUT --output OUTPUT     Output directory.
Other options:
    -m SEEDLEN --minimum SEEDLEN  Minimum seed length (default 5).
    -t INT --stringency INT       How strict the subpath reporting is.
                                  Higher values allow for more distant
                                  anchors to form a chain. (Range: 0-5. Default: 0)
    -d  --debug               	  Debug mode on.
</pre>    
   

<h1>Authors</h1>

<p>The main contributors are:</p>
<ul>
<li>Anna Kuosmanen (the modules in the top directory, parts of Traphlor)</li>
<li>Topi Paavilainen (the algorithm for finding the minimum path cover, MC-MPC)</li>
<li>Ahmed Sobih (parts of Traphlor)</li>
</ul>

<h1>Contact</h1>
aekuosma@cs.helsinki.fi

<h1>Citation</h1>

<p>If you use this project in an academic setting, please cite</p>

<p>Kuosmanen A., Paavilainen T., Gagie T., Chikhi R., Tomescu A., Mäkinen V. (2018) Using Minimum Path Cover to Boost Dynamic Programming on DAGs: Co-linear Chaining Extended. In: Raphael B. (eds) Research in Computational Molecular Biology. RECOMB 2018. Lecture Notes in Computer Science, vol 10812. Springer, Cham</p>

<p>The preprint of the paper is available <a href="https://arxiv.org/abs/1705.08754">here in arxiv.</a></p>

<h1>References</h1>
<p>[1] https://github.com/tobtobtob/MC-MPC</p>
<p>[2] https://sourceforge.net/projects/ilordmpc/</p>
<p>[3] https://github.com/jltsiren/gcsa2</p>
<p>[4] https://github.com/simongog/sdsl-lite</p>
<p>[5] http://lemon.cs.elte.hu/trac/lemon</p>
<p>[6] https://github.com/vgteam/vg</p>

