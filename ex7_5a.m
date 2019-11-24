%% Copyright (C) 2019   Eun Hwan Shin
%
% This file is the Octave implementation of the fist part of Example 7.5: 
%          Quadrilateral Adjustments by Directions.
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 151-155
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% 3. The name of the author may not be used to endorse or promote products
% derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
% NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

clear;

as2r = pi/180/3600; % arc-sec to rad 
r2as = 1/as2r;      % rad to arc-sec

%-- Adjustment using directions
fprintf('Adjustment using directions\n');

md_dms = [...
  0,  0,     0;  % 10 -  1
 22, 01, 42.51;  % 11 -  2
 38, 46, 13.71;  % 12 -  3
  0,  0,     0;  % 13 -  4
 57, 08, 57.10;  % 14 -  5
 76, 42, 11.23;  % 15 -  6
  0,  0,     0;  % 16 -  7
 86, 33, 13.45;  % 17 -  8
145, 19, 49.38;  % 18 -  9
  0,  0,     0;  % 19 - 10
 15, 06, 52.28;  % 20 - 11
 99, 11, 42.94]; % 21 - 12

% direction measurements in arc-sec 
md_as = md_dms(:,1)*3600 + md_dms(:,2)*60 + md_dms(:,3); 

% direction measurements in rad 
md_r = md_as * as2r;


A_d = [ ...
-1, 0, 1, -1,  1, 0,  0,  0, 0,  0, -1, 1;
 0, 0, 0,  0, -1, 1, -1,  0, 1, -1,  1, 0;
-1, 1, 0,  0,  0, 0,  0, -1, 1, -1,  0, 1;
zeros(1,12)];

% pi in arc-sec
pi_as = 180*3600;
f_d = -A_d * md_as + [pi_as; pi_as; pi_as; 0];

U = sin(md_r(2)-md_r(1)) * sin(md_r(11)-md_r(10)) * sin(md_r(8)-md_r(7)) * sin(md_r(5)-md_r(4));
V = sin(md_r(3)-md_r(2)) * sin(md_r(12)-md_r(11)) * sin(md_r(9)-md_r(8)) * sin(md_r(6)-md_r(5));

U_V = U/V;

f_d(4) = (1 - U_V) * r2as;

A_d(4, 1) = -U_V * cot(md_r(2)-md_r(1));
A_d(4, 2) =  U_V * (cot(md_r(2)-md_r(1)) + cot(md_r(3)-md_r(2)));
A_d(4, 3) = -U_V * cot(md_r(3)-md_r(2));
A_d(4, 4) = -U_V * cot(md_r(5)-md_r(4));
A_d(4, 5) =  U_V * (cot(md_r(5)-md_r(4)) + cot(md_r(6)-md_r(5)));
A_d(4, 6) = -U_V * cot(md_r(6)-md_r(5)); 
A_d(4, 7) = -U_V * cot(md_r(8)-md_r(7)); 
A_d(4, 8) =  U_V * (cot(md_r(8)-md_r(7)) + cot(md_r(9)-md_r(8)));
A_d(4, 9) = -U_V * cot(md_r(9)-md_r(8)); 
A_d(4,10) = -U_V * cot(md_r(11)-md_r(10)); 
A_d(4,11) =  U_V * (cot(md_r(11)-md_r(10)) + cot(md_r(12)-md_r(11)));
A_d(4,12) = -U_V * cot(md_r(12)-md_r(11)); 

% Let Q = eye(12)
Qe = A_d * A_d';
We = inv(Qe);

k_d = We * f_d;
v_d = A_d' * k_d

d_hat = (md_as + v_d)./3600;
d_hat_dms = deg2dms(d_hat);

% Adjusted eight angles of the quadrilateral
F = [ ...
-1,  1, 0,  0,  0, 0,  0,  0, 0,  0,  0, 0;
 0, -1, 1,  0,  0, 0,  0,  0, 0,  0,  0, 0;
 0,  0, 0, -1,  1, 0,  0,  0, 0,  0,  0, 0;
 0,  0, 0,  0, -1, 1,  0,  0, 0,  0,  0, 0; 
 0,  0, 0,  0,  0, 0, -1,  1, 0,  0,  0, 0;
 0,  0, 0,  0,  0, 0,  0, -1, 1,  0,  0, 0;
 0,  0, 0,  0,  0, 0,  0,  0, 0, -1,  1, 0; 
 0,  0, 0,  0,  0, 0,  0,  0, 0,  0, -1, 1]; 
l_hat = F * d_hat;
l_hat_dms =  deg2dms(l_hat);

% a-posteriori cofactor matrix of the directions
I_Qvv = eye(12) - A_d' * We * A_d;
Qdd = I_Qvv * I_Qvv;

% Weighted sum of residuals squared
fprintf('Weighted sum of residuals squared = %.4f (arc-sec^2)\n', k_d' * f_d);

Qll = F * Qdd * F';

