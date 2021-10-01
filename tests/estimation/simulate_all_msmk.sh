#!/bin/bash

# Remove "key3" and write results back to test.json (recreate it with result).
jq -c '.gen = "r6"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f6"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
