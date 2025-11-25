# modfilterbank
 
Spectro-Temporal Modulation Filter Bank Model decompostiuon of natural sounds.

The model consistes of a cochler filterbank stage followed by a bank of 
spectro-temporal receptive field filters that model the auditory midbrain 
trasformations. These STRFs consist of non-separable Gabors function similar to 
those used to fit nueral data in IC by Qiu et al 2003. There are two versions of 
the code. 
  1) The first is used to generate the modualtion power spectrum (MPS) of a sound
     where the sound power can be plotted as a funciton of temporal and spectral 
     modulation frequency.
  2) The second approach genearates a high-dimensional model representation of the 
     auditory midbrain outputs for a sound. Unlike the MPS approach which collapses 
     the sound dimensions to two dimensional space (temporal and spectral modualtion)
     this model output represents the sound in a 4-dimensional space 
     (frequency, temporal modualtion, spectral modualation and time)
