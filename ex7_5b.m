%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of the second part of Example 7.5: 
%          Quadrilateral Adjustments by Angles.
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 155-159
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
%%

clear;

as2r = pi/180/3600; % arc-sec to rad 
r2as = 1/as2r;

% angle measurements in DMS
ma_dms = [ ...
 22, 01, 42.51;  % 1
 16, 44, 31.20;  % 2
 57, 08, 57.10;  % 3
 19, 33, 14.13;  % 4
 86, 33, 13.45;  % 5
 58, 46, 35.93;  % 6
 15, 06, 52.28;  % 7
 84, 04, 50.66]; % 8 
 
% angle measurements in arc-sec  
ma_as = ma_dms(:,1)*3600 + ma_dms(:,2)*60 + ma_dms(:,3); 

% angle measurements in rad 
ma_r  = ma_as * as2r;

F = [ ...
-1,  1, 0,  0,  0, 0,  0,  0, 0,  0,  0, 0;
 0, -1, 1,  0,  0, 0,  0,  0, 0,  0,  0, 0;
 0,  0, 0, -1,  1, 0,  0,  0, 0,  0,  0, 0;
 0,  0, 0,  0, -1, 1,  0,  0, 0,  0,  0, 0; 
 0,  0, 0,  0,  0, 0, -1,  1, 0,  0,  0, 0;
 0,  0, 0,  0,  0, 0,  0, -1, 1,  0,  0, 0;
 0,  0, 0,  0,  0, 0,  0,  0, 0, -1,  1, 0; 
 0,  0, 0,  0,  0, 0,  0,  0, 0,  0, -1, 1]; 
 
Qll = F * F';
 
A_a = [ ...
1, 1, 1, 0, 0, 0, 0, 1;
0, 0, 0, 1, 1, 1, 1, 0;
1, 0, 0, 0, 0, 1, 1, 1;
zeros(1,8)]; 

pi_as = 180*3600;

f = [pi_as; pi_as; pi_as; 0] - A_a * ma_as;

Y = sin(ma_r(1)) *  sin(ma_r(3)) * sin(ma_r(5)) * sin(ma_r(7));
Z = sin(ma_r(2)) *  sin(ma_r(4)) * sin(ma_r(6)) * sin(ma_r(8));
Y_Z = Y/Z;
f(4) = (1 - Y_Z)*r2as;

A_a(4,1) =  Y_Z * cot(ma_r(1));
A_a(4,2) = -Y_Z * cot(ma_r(2));
A_a(4,3) =  Y_Z * cot(ma_r(3));
A_a(4,4) = -Y_Z * cot(ma_r(4));
A_a(4,5) =  Y_Z * cot(ma_r(5));
A_a(4,6) = -Y_Z * cot(ma_r(6));
A_a(4,7) =  Y_Z * cot(ma_r(7));
A_a(4,8) = -Y_Z * cot(ma_r(8));

Qe = A_a * Qll * A_a';
We = inv(Qe);

% the vector of Lagrange multipliers
k = We * f;

% the vector of residuals
v = Qll * A_a' * k;

% Adjusted angles
a_hat = (ma_as + v)./3600;
a_hat_dms = deg2dms(a_hat);

% Weighted sum of residuals squared
vWv = k' * f;
fprintf('Weighted sum of residuals squared = %.4f (arc-sec^2)\n', vWv);
fprintf('a posteriori estimate of the reference variance = %.4f\n', vWv/4);

% Cofactor matrix of the estimated angles
I_Qvv = eye(8) - Qll * A_a' * We * A_a;
Qaa = I_Qvv * I_Qvv * Qll;