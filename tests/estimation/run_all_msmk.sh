#!/bin/bash

# Change key values and write results back to settings.json (recreate it with result).
jq -c '.gen = "r7"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f7"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate.R &
sleep 8s
jq -c '.gen = "r8"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f8"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate.R &
sleep 8s
jq -c '.gen = "r9"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f9"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate.R &
sleep 8s
jq -c '.gen = "r10"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f10"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate.R &
sleep 8s
jq -c '.gen = "r11"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f11"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate.R &
sleep 8s
jq -c '.gen = "r12"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f12"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate.R &
wait
