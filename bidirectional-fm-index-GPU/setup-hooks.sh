#!/bin/bash

# Check if we are in a git repository
if [ ! -d ".git" ]; then
    echo "This script must be run from the root of a git repository."
    exit 1
fi

# Copy the pre-push hook to the .git/hooks directory
cp hooks/pre-push .git/hooks/pre-push
chmod +x .git/hooks/pre-push
