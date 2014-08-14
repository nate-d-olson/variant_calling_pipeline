//flash max overlap
MAXOVERLAP=150

//Legnth of kmer used for analysis
KMER=5

//Number of bootstrap replicates
REP=10

//Minimum quality
QMIN=30

//Minimum percent of bases below quality threshold
PQMIN=90

// Set location of executables
FLSH_BIN="/Users/nolson/Desktop/FLASH-1.2.11/"
FASTX_BIN="/Users/nolson/Desktop/bin"
//FFP_BIN ran sudo make install binaries are in path


// Log data from various tools will appear in here
LOG="pipeline.log"

// Number of threads
n=8 //bpipe running in parallel - only run individual stages as single threads

