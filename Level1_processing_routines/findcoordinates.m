function [lon,lat] = findcoordinates(time)

% Load the ship track
shippos = load('shiptrack.mat');

% Eliminate repeated valutes
[C,IA,~] = unique(shippos.dnum);

% Interpolate linearly for longitude and latitude
lon = interp1(C,shippos.lon(IA),time);
lat = interp1(C,shippos.lat(IA),time);
