% Change working directory to GitHub folder
cd C:\Users\Victor\Documents\GitHub\PSC ;

% Clear all functions from RAM before going any further
clear functions ;

% Run PSC with each fitness function
%pscdemo('SEGMENTED',1,'HOLDOUT-20G',false) ;
%pscdemo('SEGMENTED',2,'HOLDOUT-20G',false) ;
%pscdemo('SEGMENTED',3,'HOLDOUT-20G',false) ;

% Alternative Run
pscdemo('SEGMENTED',2,'HOLDOUT',false) ;