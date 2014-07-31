//Project directory
<<<<<<< HEAD
PRJ_HOME="/home/nolson/Desktop/micro_rm_dev"
=======
PRJ_HOME="/media/nolson/second/mirror/micro_rm_dev"
>>>>>>> 2937bf88ff6840836e222b947c9a056ad4bfa317
// Set location of you reference files here (see below for the files required)
REFBASE="$PRJ_HOME/references"

// Set a good location for storing large temp files here (probably not /tmp)
TMPDIR="$PRJ_HOME/tmp"

// Set location of executables
BIN="$PRJ_HOME/bin"

// Log data from various tools will appear in here
LOG="pipeline.log"

// Number of threads
<<<<<<< HEAD
n=7
=======
n=1 //bpipe running in parallel - only run individual stages as single threads
>>>>>>> 2937bf88ff6840836e222b947c9a056ad4bfa317

// java heap size
JHEAP="-Xmx8g"
