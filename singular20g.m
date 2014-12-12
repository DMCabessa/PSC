% Change working directory to GitHub folder
cd C:\Users\Victor\Documents\GitHub\PSC ;

% Clear all functions from RAM before running any functions
clear functions ;

% Run PSC with each fitness function
pscdemo('SINGULAR',1,'HOLDOUT-20G') ;
pscdemo('SINGULAR',2,'HOLDOUT-20G') ;
pscdemo('SINGULAR',3,'HOLDOUT-20G') ;