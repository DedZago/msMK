#!/bin/bash

# Change key values and write results back to settings.json (recreate it with result).

# DensArray=("6" "7" "8" "9" "10" "11" "12")
# nArray=(250 500 1000)
# pArray=(2 5 10)
# for dens in ${DensArray[*]}; do
#     for n in ${nArray[*]}; do
#         for p in ${pArray[*]}; do
#             jq -c ".gen = \"r${dens}\"" settings.json > tmp.$$.json && mv tmp.$$.json settings.json
#             jq -c ".df = \"f${dens}\"" settings.json > tmp.$$.json && mv tmp.$$.json settings.json
#             jq -c ".p = ${p}" settings.json > tmp.$$.json && mv tmp.$$.json settings.json
#             jq -c ".n = ${n}" settings.json > tmp.$$.json && mv tmp.$$.json settings.json
#             Rscript simulate_DP.R &
#             sleep 2s
#             pwait 6
#         done
#     done
# done
jq -c '.gen = "r7"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f7"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate_DP.R &
sleep 5s
jq -c '.gen = "r8"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f8"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate_DP.R &
sleep 5s
jq -c '.gen = "r9"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f9"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate_DP.R &
sleep 5s
jq -c '.gen = "r10"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f10"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate_DP.R &
sleep 5s
jq -c '.gen = "r11"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f11"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate_DP.R &
sleep 5s
jq -c '.gen = "r12"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
jq -c '.df = "f12"' settings.json > tmp.$$.json && mv tmp.$$.json settings.json
gnome-terminal -- Rscript simulate_DP.R &
wait
