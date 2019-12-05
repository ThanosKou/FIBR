This repository contains files needed for a Matlab simulation regarding handover in a multi-connectivity framework.

To start the simulation, run the SimulationLOS.m found in the HO_Simulations folder. 

The simulation considers dynamic and self blockages and assumes user multi-connectivity, i.e. a user is able to connect to multiple Base Stations simultaneously, and can switch to a new BS (if available) in the case of blockage, such that the required degree of multi-connectivity is maintained.

The simulation's outputs include the probability of Radio Link Failure, the average blockage duration and the IDs of the serving BSs at all blockage times.
