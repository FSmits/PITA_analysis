# PITA_analysis
Pilot for tACS (PITA) study analysis scripts

Pre-registration: https://osf.io/tyfq5

This repository contains the data analysis scripts.


Study description:
This study has a counter-balanced crossover design with two conditions: active tACS and sham tACS (2-sec fade-in/fade-out). A washout period of minimally one day is applied to minimize carry-over effects between sessions. Experimental sessions last approximately 90 minutes.

Following screening, eligible participants are randomly assigned to receive sham or active tACS during a period of associative memory consolidation. To reduce individual state-dependent effects and between-subjects’ variability in attentional engagement, a video is shown during tACS intervention. 

To apply tACS, two 3 x 3 cm tACS electrodes in saline-soaked sponges covered by gel are positioned over 10-20 EEG electrode positions F3 and F4. This bifrontal montage was chosen as the frontal cortex is implicated in a neural network subserving associative memory processes. 

EEG is recorded from 30 scalp electrodes (F3 and F4 are missing due to tACS electrode placement, remaining electrodes: Fp1, Fp2, AF3, AF4, F7, F8, Fz, Cz, FC1, FC2, FC5, FC6, C3, C4, CP1, CP2, CP5, CP6, P7, P8, P3, P4, Pz, PO3, PO4, T3, T4, O1, O2 and Oz), according to standard procedure. EEG is recorded with a BioSemi ActiveTwo system at a sampling rate of 2048 Hz.

4 minutes of resting state EEG is recorded (alternating 1 min. eyes-open, 1 min. eyes-closed) at the start of the experimental session and right after cessation of the stimulation (tACS). 

Frontal tACS (5 Hz, 2mA peak-to-peak) is applied using an interleaved tACS/EEG protocol: twenty blocks of stimulation (1 minute) are alternated with blocks without stimulation (0.5 minute) to allow for EEG recording, adding up to a total stimulation time of 20 minutes spread over 30 minutes. 

Associative memory performance is assessed using a task adapted from a memory contextualization task (MCT) (van Ast et al., 2014; Sep et al., 2019). Individuals are instructed to memorize the combinations of context (scene stimuli) and face stimuli. Participants undergo 50 encoding trials. In each trial, context is presented on the screen, and 1 second after a face stimulus is presented on top of the context. The combined context-face presentation lasts 3 seconds. Trials are separated by 500 ms inter-trial interval. Total duration of the encoding phase is approximately 4 minutes and is followed by the tACS intervention. 

After the tACS intervention and second resting-state EEG (+/- 35 minutes after associative memory encoding phase), participants complete 150 self-paced memort retrieval trials: congruent trials (context-face combination shown during encoding), incongruent trials (context or face shown before, but not in this combination), and lure trials (new context and face stimuli) are presented in a randomized order. A total of 100 different faces depicting a neutral facial expression are divided into two gender-matched subsets: 50 stimuli are identical to encoding, while the other half is used as lure during the retrieval phase. The old faces are presented against either a different encoding context (incongruent non-target trials) or the original encoding context (congruent target trials), whereas lure faces are combined with new contexts (lure non-target trials).
Subjects are instructed to indicate whether the face-context combination was presented during the encoding phase or not (yes/no), and to what extent they are confident using a 3-point Likert Scale (1= not confident at all; 3= completely confident). Measures of interest are: hit rate (proportion of correctly recognizing targets as “old”), false alarms (FA, proportion of misjudging non-target trials as “old”), and d-prime sensitivity index (d’,
Z(hit) – Z(FA)).
