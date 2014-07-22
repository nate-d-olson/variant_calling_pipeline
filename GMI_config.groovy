//Project directory
PRJ_HOME="/home/ubuntu/GMI_bioinf/variant_calling_pipeline"
// Set location of you reference files here (see below for the files required)
REFBASE="$PRJ_HOME/references"

// Set a good location for storing large temp files here (probably not /tmp)
TMPDIR="$PRJ_HOME/tmp"

// Set location of executables
BIN="~/bin"
BWA_BIN="~/bwa"

// Log data from various tools will appear in here
LOG="pipeline.log"

// Number of threads
n=8 //bpipe running in parallel - only run individual stages as single threads

// java heap size
JHEAP="-Xmx8g"
