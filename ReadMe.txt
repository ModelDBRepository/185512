%% MATLAB code to go with JCNS "Effect of Polarization Induced by Non-Weak Electric Fields on the Excitability of Elongated Neurons With Active Dendrite"
Robert I. Reznik  Ernest Barreto  Evelyn Sander  Paul So

RIR November 7, 2015

The paper explores the effects of polarization (VdsOut) on neural excitability. We use the two-compartment Pinsky-Rinzel model (1994)
In this work excitability is quantified by the TTFS of the neuron in response to stimulus.
The main stimulus used is ramp injected current (applied to soma). '
Later we use a model of synaptic AMPA applied to the dendrite instead of the ramp injection.
 Also we modify the PR neuron by adding an Ih current as modeled by Lippert and Booth (2009).

** I apologize for the extrememly long file names. I like to pack in a lot of metadata into the filenames themselves.

**************************************************************************

Main computational code

(1) Finding the equilibriums (if they exist)
	1.1 Use NumerEquilPR_db
		*input uses the PR structure
		*requires an initial guess but I found the guess can be quite off.
(2) Integrating the polarized PR neuron with ODE23(You are welcome to use other integrators I also used ODE45)
	2.1 In the folder Integration you will see a number of routines I would use
		2.1.1 SingIntegODE23PRWithInjCurr_db for use of the ramp protocol
		2.1.2 SingIntegODE23PRWithSyn_db for using the synaptic AMPA
		2.1.3 SingIntegODE23PRWithSynWithIntegParam_db is fairly versatile and can be used with a number of different protocols including syanptic and injected current,
		2.1.4 SingIntegODE23PRPlusHcurrWithSyn_db Synaptic AMPA and Ih current

	2.2 Also in the Integrate folder there is the initialize routine which uses a structure "aPR" and first populates it with default values. Later, you can set the values to whatever you want.
		2.2.1 IniPR_db

(3) The folder PR_WIh contains all of the additional routines needed to run code for the polarized PR with Ih current.	

***************************************************************	

General notes for the driver routines. All follow these general steps:

(1) Initialize parameters.
(2) Inside a loop over polarization (VdsOut)
	2.1 numerically search for an equilibrium (At this point there is no stimulus)
		(If it does not exist go to next increment in VdsOut loop)
	2.2 Intgrate ODE (with appropriate stimulus protocol)
		Default is to stop integration upon Vs passing up through 30 mV.
	2.3 optionally capture computational data such as currents and membrane potentials.
(3) There is the option, of course, to loop over additional parameters like Ek

************************************************
	
Some suggested Scripts and Figures to start with

(1) Fig4_ContIntmpt8 and Fig5_ContIntmpt3 reproduce figs 4 and 5 of paper.

(2) Plotting_ActiveDendCurr_AMPA_db
(3) TriipleLoopFixgAMPAVaryEkandgKAHP_db (This file shows a loop over 3 parameters for Polarized PR with AMPA
(4) PRwIhAMPALoopEkgKAHP.m (This one calls SingIntegODE23PRPlusHcurrWithSyn_db which uses both synaptic AMPA and the Ih currents) As currently set this might take awhile to run. The resulting dat is in computationaldata folder and figure 18 can be produced by running GridPlotLikeFig18PRwHJust3Grids.m
	