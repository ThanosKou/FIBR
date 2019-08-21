This repository contains files needed for a Matlab simulation regarding handover in a multi-connectivity framework.

To start the simulation, run the SimulationLOS.m found in the HO_Simulations folder. 

The simulation considers dynamic and self blockages and assumes user multi-connectivity, i.e. a user is able to connect to multiple Base Stations simultaneously, and can switch to a new BS (if available) in the case of blockage, such that the required degree of multi-connectivity is maintained.

The simulation's outputs is the probability of Radio Link Failure, which occurs if the user is blocked by all its Base Stations for a prolonged period of time.
