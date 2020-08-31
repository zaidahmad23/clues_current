# clues_script
step.py creates .traj files for mssel 
case2.py runs helper functions for Ancient DNA stuff 
script_v8.py is the main simulation pipeline 


LAUNCH OPTIONS:
CASE 1 (--times, no --ancientSamps):
python3 script_v8.py -p0 0.2 -s 0.04 -n 4000 --nder 4 --ton 100 --toff 0 \
--converted-filename relate_input \
--path-to-relate-bin ./relate_v1.1.2_x86_64_static/bin/Relate \
--relate-output-filename relate_step_1 \
--relate-map-file-path clues/example/genetic_map.txt \
--sample-branch-length-coal-file-path clues/example/example.coal \
--path-to-sample-branch-length-script ./relate_v1.1.2_x86_64_static/scripts/SampleBranchLengths/SampleBranchLengths.sh \
--sample-branch-length-script-format b \
--sample-branch-length-script-output-filename clues_input_from_relate \
--sample-branch-length-script-n-samples 5 \
--inference-script-output-filename clues_output \
--output-directory output --inference-script-coalescence-times-file-path clues_input_from_relate

CASE 2 (no --times, --ancientSamps):
python3 script_v8.py -p0 0.2 -s 0.04 -n 4000 --nder 4 --ton 100 --toff 0 \
--converted-filename relate_input \
--path-to-relate-bin ./relate_v1.1.2_x86_64_static/bin/Relate \
--relate-output-filename relate_step_1 \
--relate-map-file-path clues/example/genetic_map.txt \
--sample-branch-length-coal-file-path clues/example/example.coal \
--path-to-sample-branch-length-script ./relate_v1.1.2_x86_64_static/scripts/SampleBranchLengths/SampleBranchLengths.sh \
--sample-branch-length-script-format b \
--sample-branch-length-script-output-filename clues_input_from_relate \
--sample-branch-length-script-n-samples 5 \
--inference-script-output-filename clues_output \
--output-directory output --create-ancient-samples --case2-script-ancient-samples-generation-gap 500 \
--case2-script-number-of-ancient-samples 100

CASE 3 (--times, --ancientSamps):
python3 script_v8.py -p0 0.2 -s 0.04 -n 4000 --nder 4 --ton 100 --toff 0 \
--converted-filename relate_input \
--path-to-relate-bin ./relate_v1.1.2_x86_64_static/bin/Relate \
--relate-output-filename relate_step_1 \
--relate-map-file-path clues/example/genetic_map.txt \
--sample-branch-length-coal-file-path clues/example/example.coal \
--path-to-sample-branch-length-script ./relate_v1.1.2_x86_64_static/scripts/SampleBranchLengths/SampleBranchLengths.sh \
--sample-branch-length-script-format b \
--sample-branch-length-script-output-filename clues_input_from_relate \
--sample-branch-length-script-n-samples 5 \
--inference-script-output-filename clues_output \
--output-directory output --inference-script-coalescence-times-file-path clues_input_from_relate \
--create-ancient-samples --case2-script-ancient-samples-generation-gap 500 \
--case2-script-number-of-ancient-samples 100


TO DO:
1. fix whack case1 bug (plot.pdf created, but no plot)
2. fix slight display bug for case 2, case 3 where a bit of the plot is graphically cut off
