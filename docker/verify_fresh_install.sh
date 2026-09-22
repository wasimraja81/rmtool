#!/usr/bin/env bash
# Runs inside a container built from docker/dockerfile.docs-verify --
# walks through the exact commands BUILD.md/QUICKSTART.md/TUTORIAL.md/
# EXAMPLES.md tell a new user to run, in order, so a doc/build.sh/
# tests/run_tests.sh change that breaks the documented path fails here
# before it ships, not after a real user hits it.
set -e
cd /app/rmtool

echo "################ STEP 1: ./build.sh (fresh, no pre-existing venv) ################"
./build.sh

echo "################ STEP 2: full test suite (tests/run_tests.sh) ################"
bash tests/run_tests.sh

echo "################ STEP 3: casda_fetch.py via the auto-created venv ################"
~/venv/rmtool/bin/python3 scripts/casda_fetch.py --help >/dev/null
echo "casda_fetch.py --help: OK"

echo "################ STEP 4: TUTORIAL.md's exact documented rm_synthesis command ################"
bin/rm_synthesis_release_cpu_omp cfg/rmsynth-e2e-smalltest.cfg
echo "TUTORIAL.md rmsynth example: OK"

echo "################ STEP 5: EXAMPLES.md's exact documented multiband command ################"
bin/rm_synthesis_release_cpu_omp cfg/rmsynth-e2e-multiband-matched-smalltest.cfg
echo "EXAMPLES.md multiband rmsynth example: OK"

echo "################ STEP 6: standalone match_cubes (QUICKSTART section 7 syntax, real fixture data) ################"
# stages=reproject, not both: these fixture cubes carry no BEAMS table
# (see cfg/pipeline-e2e-multiband-smalltest.cfg's own header comment),
# so convolve has nothing to work from without an extra ASCII beamfile.
bin/match_cubes stages=reproject footprint_mode=reference \
  reffile=tests/data/TEST.Q.FITSCUBE infiles=tests/data/TEST.Q.FITSCUBE,tests/data/TEST_BAND2.Q.FITSCUBE \
  mem_frac_ram=0.25
echo "standalone match_cubes: OK"

echo "################ STEP 7: full match->rmsynth->rmclean pipeline, single-band (TUTORIAL.md's recommended cfg) ################"
scripts/run_pipeline.sh cfg/pipeline-e2e-smalltest.cfg
echo "pipeline (single-band): OK"

echo "################ STEP 8: full match->rmsynth->rmclean pipeline, multi-band ################"
scripts/run_pipeline.sh cfg/pipeline-e2e-multiband-smalltest.cfg
echo "pipeline (multi-band): OK"

echo "################ ALL FRESH-USER STEPS COMPLETED SUCCESSFULLY ################"
