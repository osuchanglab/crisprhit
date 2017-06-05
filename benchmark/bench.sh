#!/bin/bash
echo "pre"
../crisprhit.py --PAM CTT protospacers_rev.fas spacer.fna spacer_v_protospacer.tab --quiet --filters pre > benchmark_pre.crisprhit
bash summarize.sh benchmark_pre.crisprhit
echo "pre_first"
../crisprhit.py --PAM CTT protospacers_rev.fas spacer.fna spacer_v_protospacer.tab --quiet --filters pre --filters first > benchmark_pre_first.crisprhit
bash summarize.sh benchmark_pre_first.crisprhit
echo "pre_first_second"
../crisprhit.py --PAM CTT protospacers_rev.fas spacer.fna spacer_v_protospacer.tab --quiet --filters pre --filters first --filters second > benchmark_pre_first_second.crisprhit
bash summarize.sh benchmark_pre_first_second.crisprhit
echo "pre_first_second_third"
../crisprhit.py --PAM CTT protospacers_rev.fas spacer.fna spacer_v_protospacer.tab --quiet --filters pre --filters first --filters second --filters third > benchmark_pre_first_second_third.crisprhit
bash summarize.sh benchmark_pre_first_second_third.crisprhit
echo "pre_first_second_third_post"
../crisprhit.py --PAM CTT protospacers_rev.fas spacer.fna spacer_v_protospacer.tab --quiet --filters pre --filters first --filters second --filters third --filters post > benchmark_pre_first_second_third_post.crisprhit
bash summarize.sh benchmark_pre_first_second_third_post.crisprhit
