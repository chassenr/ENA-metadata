
    #!/bin/bash

    # Check if the correct number of arguments is provided
    if [ "$#" -ne 2 ]; then
        echo "Usage: $0 <input_file> <output_file>"
        exit 1
    fi

    input_file="$1"
    output_file="$2"

    # Check if the input file exists
    if [ ! -f "$input_file" ]; then
        echo "Input file not found: $input_file"
        exit 1
    fi

    # Read each line of the input file
    cat "$input_file" | while IFS= read -r line; do
        # Extract corpusid from the current line
        corpus_id=$(echo "$line" | jq -r '.corpusid')

        # Extract text content from the current line
        text_content=$(echo "$line" | jq -r '.content.text')

        # Use grep with the regular expression to find matches
        matches=$(echo "$text_content" | grep -oE '[A-Z]{1}RP[0-9]{6,8}')

        # Write matches to a text file if found
        if [ -n "$matches" ]; then
            for match in $matches; do
                echo "$corpus_id $match" >> "$output_file"
            done
        fi
        
        done

    echo "Output written to $output_file"
