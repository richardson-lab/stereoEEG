# stereoEEG
Repository to house code for a collaborative project analyzing stereoEEG data during anesthesia induction and sleep for NINDS-R01NS113366: [The role of dynamical criticality in human perception](https://reporter.nih.gov/project-details/9969278)

## Downloading data from iEEG.org
Process and save .mat file for each session with:  
`process_ieeg_data;`

## Sync behavioral-sEEG data
Run `syncToneResponse;`. 
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
