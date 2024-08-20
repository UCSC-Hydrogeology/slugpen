%%%% ======================================================================
%%  Purpose
%     This function filters raw data based on the median
%
%%   Last edit
%      07/19/2023 by Kristin Dickerson, UCSC
%%% =======================================================================

function out = rbmmed(in,wl)

% median filtering of realtime deeptow mag
% fluxgate sensor
% Input array is set to: in
% Output array is set to: out
%
% Usage out=rbmmed(in,window_length)
%--------------------------------------------
% Maurice A. Tivey 12 Jan 1992
%--------------------------------------------
    
    wl2=wl/2;
    
    start=wl2+1;
    
    finish=length(in)-start;
    
% initialise array 'out'
    
    out=in;
    
    for i=start:finish
    
      out(i)=median(in(i-wl2:i+wl2));
    
    end

