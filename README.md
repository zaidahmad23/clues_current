# clues_script
step.py creates .traj files for mssel 

step2.py runs helper functions for Ancient DNA stuff 

case1.py runs inference with only --times specified 
case2.py runs inference with only --ancientSamps specified
case3.py runs inference with both --times and --ancientSamps specified

LAUNCH OPTIONS:

**CASE 2 (no --times, --ancientSamps):**
~~~
python3 case2.py -p0 0.2 -s 0 -n 10000 --ton 50 --toff 0 \
--converted-filename relate_input \
--inference-script-output-filename clues_output \
--output-directory output --mutation-rate 1.25e-8 \
--create-ancient-samples --step2-script-ancient-samples-generation-gap 600 --step2-script-number-of-ancient-samples 500
~~~

**TO DO:**
2. write better documentation lol
