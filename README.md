# TestBeamProg

-Ntuples_wiki: basic Wiki for ntuples from T10 May 2018 Test Beam on crystal bars

-plotWF.C(Conf,nevento) plots the 6 wavwforms acquired in one event if given (string runname,int eventnumber), plots iteratively for all the events if given (string runname, int 0); expects a primitive to pass from one event to the next one

-ploWF_graph.C(Conf) plots max amplitude distribuitions  and creates a file histoConf.root contaiing such distributions

-plotHisto.C(histoConf1,histoconf2) plots max amplitude distributions overlapped

-plotWF_fit.C(Conf) plots max amplitude distributions and fits the distribuitions with a landau function

-plotWF_cut.C(Conf) same as plotWF_fit.C + cuts on data inbetween (0.8*MIP,3*MIP)

-plotWF_tamp.C(Conf) plots 2D histograms  t_r-t_MCP t_l-t_MCP t_ave-t_MCP vs max amplitude , fits on time averages of 1D histograms obtaining by summing 1D histograms at fixed max amplitude for bins of max amplitude. It uses time recorded with LED300 by default.

-plot_lsign.C(Conf) plots time distributions for uncorrected (no ampwalk, no tdiff) timestamps

-plotWF_su.C(Conf) plots max amplitude vs t_left-t_right

-plotWF_corr.C(Conf) plots same graphs as _tamp.C, redefines timestamps using fits done in _tamp.C and plots time distributions for the corrected timestamps

-plotWF_tdiff.C

-plotWF_time.C

-TResAmp.C

-RisRelVsMip

-RisVsMip

.RisRelVsTdiff

-RisVsTdiff

-DarkBkg














======T10 ntuples elementary wiki =======

This file contains information on the structure and main variable of the ntuples of the T10 May 2018 TesBeam (LYSO Crystal Bar time resolution).

The good ntuples(both for data and pedestals) are available on lxplus on /eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_T10_May2018/ntuples_v1/
The TestBeam logbook (info on configurations' different settings) is available at https://docs.google.com/spreadsheets/d/1ArGOxF1clg_I_9lgCssy9A57RwJ98hGXUgBIKj3aCq4/edit#gid=1831808220

Each ntuple contains 4 trees:
     --info
     --wf
     --h4
     --digi

==========TREES DESCRIPTION=============

The INFO tree contains:
_____index of reference for both run and event (Format run00000000event);
_____tableX, tableY: information on the position of the bar;
_____config: code of configuration;
_____Vbias_bar: bias potential on the bar left and right SIPMs;
_____SIPM_current bar: current in SIPMs (in ADC counts);
_____NINOthr_bar: threshold (in ADC counts) on the raw SIPM signal to trigger NINO acquisition;

The tree has also a section dedicated to voltage and current on crystalbar matrices, not used in these ntuples.




The WF tree  contains the information on the raw waveforms collected from each detector in the system:
____index of reference for both run and event (format run0000000event0)
____WF_samples: total number of samples collected for a single event;
____WF_ch(vector): each sample taken is associated to the channel it belongs by the value of this vector at given sample (e.g. wf_ch[45]=2 means that the 45th sample belongs to the second channel
waveform). There are 6 or 8 channels in the configurations studied since now:
	   	       	    	     	         --0 - MCP pulse signal
						 --1/2- NINO amplified signal x 2																					    --3/4- Amplified signal x 2
						 -- Unkown channel (almost always noise)*;																				    --6/7- trigger delayed signal x 2**



To be clarified:
* only in some configurations;
**in some configurations there is only one trigger signal


The h4 tree contains parameter and settings necessary for the first reco step. This reco is operated with the H4Analysis software on wf data, and the values obtained are stored in digi tree:

It contains all integer numbers:
___index: same as in info;
___start_time;
___time_stamp;
___run: run number;
___spill: number of spill;
___event: is the event number also present in the index, +1;


The DIGI tree contains partially reconstructed and organized data.The variables used to perform the time resolution analysis are amp_max and time.
-------amp_max is a vector containing the maximum recorded on different detectors;
-------time is the timestamps of the rising signals in different conditions;
Each different detector or detector+settings system is identified with an integer number: the system is flexible in terms of changes and upgrades of detector and/or electronics.

A stamp of the encoding is below:

There are also: --- index: to retrieve run and event number (format as in info);
      	     	---n_channel: the number of channel recorded (6/8);
			 ---n_timetypes: the number of different time acquisition recorded and encoded;

amp_max and time, are both n_timetypes long for each recorded event.

The first encoded values (MCP, NINOBAR1,NINOBAR2,AMPBAR1,AMPBAR2) identify the detectors: their index is to be used to retrieve max amplitudes: in fact the amp_max vectors are only full on their
first n_channel entries.

The time information is instead recorded in many different ways:
----------------timestamps on the MCP signal are always taken in CFD (Constant Fraction Discrimination, fixed at 50%) to take care of the detector noisy behaviour;
----------------timestamps on SIPMs are collected both on the amplified waveforms with CFD and on NINO with different LED (Leading Edge Discriminations): in this way lots of different timestamps
			      can be taken into account to identify the right combination of detector and detector settings to obtain the better resolution.
			      
# Sept2018TB
