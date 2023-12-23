# Radiator-Tools
Tools for modelling radiators, and aiding in the selection process

Radiator_rating.m uses equations and data from Compact Heat Exchangers Third Edition, by Kays and London to predict the heat dissipation rate of a radiator, given its dimensions and flow specification. Whilst it hasn't been experimentally verified, the outputs match what is expected. 

Rate.m is based on Radiator_rating.m, but incorporates further calculations to determine air flow rate, based on car speed and fans. Furthermore it has been rewritten to be a function of car speed and inlet temp, allowing it to be used in simulations. 
