%% Copyright (C) 2019   Eun Hwan Shin
%
% This file is the Octave implementation of Example 7.1 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 140-141
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

% observed angles of a simple plane triangle.
obs = [ ...
40, 19, 02;
70, 30, 01;
69, 11, 00];

% angles in arc seconds
o_as = obs(:,1)*3600 + obs(:,2)*60 + obs(:,3);

A = [1, 1, 1];

f = 180*3600 - A * o_as;

Qe = A * A';
We = 1/Qe;
k = We * f;

% residuals
v = A' * k;

% adjusted angles
o_hat = (o_as + v)./3600;

fprintf('Adjusted angles:\n')
o_hat_dms = deg2dms(o_hat)

Qvv = A' * We * A;
Qll = eye(3) - Qvv

fprintf('rank(Qll) = %d\n',rank(Qll));
fprintf('Reference variance = %.3f (arc-sec^2)\n', f'*k);