#!/usr/bin/awk -f
BEGIN {
    if (max_length == 0) max_length = 100  # Default value if not set
    sequence = ""
}

# If line starts with '>', it's a new header
/^>/ {
    if (sequence != "") {
        process_sequence(header, sequence)
        sequence = ""
    }
    header = substr($0, 2)  # Remove '>'
    next
}

# Collect sequence lines (handling multi-line FASTA format)
{
    sequence = sequence $0
}

# Process last sequence when file ends
END {
    if (sequence != "") {
        process_sequence(header, sequence)
    }
}

# Function to split and print sequence
function process_sequence(header, sequence,    i, len, chunk) {
    len = length(sequence)
    
    # If sequence is within max_length, print it as is
    if (len <= max_length) {
        print ">" header
        print sequence
    } else {
        # Split the sequence into chunks
        for (i = 1; i <= len; i += max_length) {
            chunk = substr(sequence, i, max_length)
            print ">" header "_" int((i - 1) / max_length + 1)  # Add fragment index
            print chunk
        }
    }
}
