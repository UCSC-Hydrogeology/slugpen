%%% =======================================================================
%%  Purpose: 
%     This function allows the user to turn OFF plots in SlugPen, including
%     all the axes children (various points, lines, etc.)
%%  Last edit:
%     07/19/2023 by Kristin Dickerson, UCSC
%%% =======================================================================


function axesChildrenOFF(C)

% Loop through all children for a given axes and turn visibility off
        for k = 1:numel(C)
            S(1).type = '()';
            S(1).subs = {k};
            S(2).type = '.';
            S(2).subs = 'Visible';
            C = subsasgn(C, S, 'off');
        end