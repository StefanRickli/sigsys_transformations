% This script saves the current state of the random number generator
% in 'seed.mat' to be used later in transformations_visualizations.m

seed = rng();
save('seed.mat', 'seed');
