#!/bin/bash

python3.4 run_archetype.py v51010_XAGN

# python3.4 run_archetype.py qso_BL
# python3.4 run_archetype.py XAGN
# python3.4 run_archetype.py ELG
# python3.4 run_archetype.py LRG

# python3.4 run_archetype.py XMMSL
# python3.4 run_archetype.py QSO
# python3.4 run_archetype.py qso_t2
# python3.4 run_archetype.py X_no_BL
# python3.4 run_archetype.py GAL_agn_all

run construct_archetypes.py 0.5 1.0 v51010_XAGN 3.
run construct_archetypes.py 1.0 1.5 v51010_XAGN 3.
run construct_archetypes.py 1.5 2.5 v51010_XAGN 3.
run construct_archetypes.py 2.5 4.5 v51010_XAGN 3.