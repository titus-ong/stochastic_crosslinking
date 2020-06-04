# stochastic_crosslinking
Stochastically modelling the crosslinking reaction between polymers

A script that allows for the input of different polymers and performs the crosslinking reaction stochastically. Users can define polymers with various hydroxyl groups and their amounts. With a certain number of crosslinks (obtained from modelling software e.g. gPROMS), users can input this number and run the script to obtain information about the crosslinking network (or a 'blob').

Functions include:
- accounting for both inter- and intra-molecular crosslinking
- tracking every blob formed, not just the major network
- tracking the number of polymer chains in each blob
- stochastic crosslinking based on randomly choosing the hydroxyl group rather than the polymer chain (more accurate)

Limitations include:
- assume perfect homogeneous mixing: each hydroxyl group has an equal chance to crosslink with any other hydroxyl group
- works within reasonable time for small numbers of polymer chains: unable to model actual scale crosslinking reactions due to huge Avogadro's number

Results from this script have shown that the crosslinking reaction follows the Boltzmann sigmoidal curve.

### A simplified non-technical description
This script aims to shed some light on the following scenario:

An alien party is hosted and lots of different species of aliens come. Each species has a different number of hands. The host decides to play a prank on its guests and sprays superglue on each of their hands when they enter, such that two aliens will have their hands joined together when they shake hands. With this script and given the number of handshakes that occur during the party, we can determine the different chains of aliens that are stuck together and how many aliens of each species there are in each chain.
