%==========================================================================
%
%	This function create a new state-space object with parameters A,B,C,D 
%   from two series connected state-space model with paramters A1,B1,C1,D1 
%   and A2,B2,C2,D2 in the following form.
%
%                       sys1               sys2
%                 ---------------    ---------------
%       input --->| A1,B1,C1,D1 |--->| A2,B2,C2,D2 |---> output
%                 ---------------    ---------------
%   to
%                            -----------    
%                  input --->| A,B,C,D |---> output
%                            -----------
%   where
%
%           A 	=   [ A1     0  ]   B	=	[ B1    ]
%                   [ B2*C1  A2 ]           [ B2*D1 ]
%           C   =   [ D2*C1  C2 ]   D   =     D2*D1
%
%
%	INPUTS
%           sys1    =   state-space object of the left system
%           sys2    =   state-space object of the right system
%	OUTPUTS
%           sys     =   state-space object of the series connected system
%
%========================================================================== 
function [sys] = series_ss2ss(sys1,sys2)
A1 = sys1.A; A2 = sys2.A; C1 = sys1.C; C2 = sys2.C;
B1 = sys1.B; B2 = sys2.B; D1 = sys1.D; D2 = sys2.D;
[s1A1, ~] = size(A1); [~, s2A2] = size(A2);

A = [A1, zeros(s1A1, s2A2); B2*C1, A2];
B = [B1; B2*D1];
C = [D2*C1, C2];
D = D2*D1;

sys = ss(A,B,C,D,1);
end

