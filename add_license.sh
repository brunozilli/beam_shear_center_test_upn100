#!/bin/bash
#=======================================================================
# add_license.sh - Add MIT license header to all Fortran source files
#=======================================================================
# MIT License - Copyright (c) 2024 Bruno Zilli & DeepSeek
#=======================================================================

LICENSE_HEADER='c=======================================================================
c FILENAME_PLACEHOLDER
c=======================================================================
c MIT License
c Copyright (c) 2024 Bruno Zilli & DeepSeek
c
c Permission is hereby granted, free of charge, to any person obtaining a copy
c of this software and associated documentation files (the "Software"), to deal
c in the Software without restriction, including without limitation the rights
c to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
c copies of the Software, and to permit persons to whom the Software is
c furnished to do so, subject to the following conditions:
c
c The above copyright notice and this permission notice shall be included in all
c copies or substantial portions of the Software.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
c IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
c FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
c AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
c LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
c SOFTWARE.
c======================================================================='

MAKEFILE_HEADER='#=======================================================================
# Makefile
#=======================================================================
# MIT License
# Copyright (c) 2024 Bruno Zilli & DeepSeek
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#======================================================================='

echo "Adding MIT license headers to Fortran source files..."

# Process all .f files in src/ and test/ directories
for dir in src test; do
    if [ -d "$dir" ]; then
        for file in "$dir"/*.f; do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                echo "  Processing: $file"
                
                # Create temporary file with header
                echo "$LICENSE_HEADER" | sed "s/FILENAME_PLACEHOLDER/$filename/" > "$file.tmp"
                echo "" >> "$file.tmp"
                
                # Check if file already has a license header (look for "MIT License")
                if grep -q "MIT License" "$file"; then
                    echo "    License already present, skipping..."
                    rm "$file.tmp"
                else
                    # Append original content (skip any existing header comments)
                    # Find the first line that is not a comment or blank
                    in_header=true
                    while IFS= read -r line; do
                        if [ "$in_header" = true ]; then
                            # Skip blank lines and comment lines at the beginning
                            if [[ ! "$line" =~ ^[[:space:]]*$ ]] && [[ ! "$line" =~ ^[cC*#] ]]; then
                                in_header=false
                                echo "$line" >> "$file.tmp"
                            fi
                        else
                            echo "$line" >> "$file.tmp"
                        fi
                    done < "$file"
                    
                    # Replace original file
                    mv "$file.tmp" "$file"
                    echo "    License added."
                fi
            fi
        done
    fi
done

# Process Makefile separately
if [ -f "Makefile" ]; then
    echo "  Processing: Makefile"
    if grep -q "MIT License" "Makefile"; then
        echo "    License already present, skipping..."
    else
        echo "$MAKEFILE_HEADER" > "Makefile.tmp"
        echo "" >> "Makefile.tmp"
        
        in_header=true
        while IFS= read -r line; do
            if [ "$in_header" = true ]; then
                if [[ ! "$line" =~ ^[[:space:]]*$ ]] && [[ ! "$line" =~ ^# ]]; then
                    in_header=false
                    echo "$line" >> "Makefile.tmp"
                fi
            else
                echo "$line" >> "Makefile.tmp"
            fi
        done < "Makefile"
        
        mv "Makefile.tmp" "Makefile"
        echo "    License added."
    fi
fi

echo ""
echo "Done! All files have been processed."
