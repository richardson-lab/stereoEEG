# stereoEEG
Repository to house code for a collaborative project analyzing stereoEEG data during anesthesia induction, anesthesia emergence, and sleep for NINDS-R01NS113366: [The role of dynamical criticality in human perception](https://reporter.nih.gov/project-details/9969278)

## Downloading data from iEEG.org
Process and save .mat file for each session with:  
`process_ieeg_data;`  
If tone response data is available for the session, make sure the 'toneResponse_Induction' subfolder is in the specified 'data_dir/subject' path so that it is included in the .mat file.

## Sync behavioral-sEEG data
If tone response data is available for the session, run `toneResponseSync;`. 
This adds the tone on, tone off, button press, button release timestamped events to the session annotations. 

## Qualitative analysis: time domain
Load .mat file for session. Then produce scrolling time plot with:  
`timePlot(session.data,session.sample_rate,[],[],1,session.annotations,session.channel_labels(:,1));`.  
To also view the behavioral events on the scrolling plot, use `session.annotations2` rather than `session.annotations`. 

## Qualitative analysis: frequency domain
Load .mat file for session. Then produce power spectral density (PSD) plot with:  
`pwelchPlot(session.data,session.sample_rate,[0.5 500],5,1,500,session.channel_labels(:,1));`  
Identify bad channels by clicking on outlier PSD with Edit Plot tool and finding channel name:  
`get(gco,'Tag')`

## Stability analysis 
Replicating analysis of Solovey et al 2015 J Neurosci (Figures 2, 3, 4).
Load .mat file for session. Then perform analysis with:  
`seegStability(session);`.
All available sessions can be batch processed using:  
`seegStabilityBatch;` which will generate a plot for each session and save a .mat file for population analysis.  
After loading that population .mat file containing the structure `stability`, the population analysis can be run with:  
`seegStabilityPop(stability);`

## Time-frequency analysis
Computes both average power and average criticality in each time x frequency bin. 
Load .mat file for session. Then perform analysis with:  
`seegTF(session);`.

## Criticality-frequency analysis
Computes distribution of mode criticality vs mode frequency for a specified pair of time windows (awake vs anesthetized) in a session. 
Load .mat file for session. Then peform analyis with:
`seegCF(session,tawake,tanesth);`.
