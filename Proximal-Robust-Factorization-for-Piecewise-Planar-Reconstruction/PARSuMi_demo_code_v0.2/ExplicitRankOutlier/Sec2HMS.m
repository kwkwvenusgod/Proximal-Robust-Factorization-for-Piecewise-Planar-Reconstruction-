function [ hours, minutes,seconds ] = Sec2HMS( timeElapsed )
%SEC2HMS Summary of this function goes here
%   Detailed explanation goes here

hours = floor(timeElapsed/3600); 
timeMin = timeElapsed - (hours*3600);
minutes = floor(timeMin/60);
timeSec = timeElapsed - (hours*3600) - (minutes*60);
seconds = round(timeSec);

end

