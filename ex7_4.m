%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 7.4: Adjustment of a
% Level Net
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 147-151
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

table7_1

Q = diag(0.1*dist); % (m^2)

A = [ ...
1,  1,  1, 0, 0,  0, 0,  0;
0,  0, -1, 0, 0, -1, 1,  0;
0, -1,  0, 1, 0,  0, 0, -1;
0,  0,  0, 0, 1,  1, 0,  1];

% Reduced normal equations coefficient matrix and its inverse
Qe = A * Q * A';
We = inv(Qe);

% Lagrange multipliers
f = -A * dE
k = We * f;

% vector of residuals
v = Q * A' * k;

% estimated observations
l_hat = dE + v;

% a posteriori cofactor matrix of the estimated observations
I_Qvv = eye(8) - Q * A' * We * A;
Qll = I_Qvv * I_Qvv * Q;

% Elevation of point A
EA = 800.0;

J = [...
1, 0,  0, 0, 0, 0,  0, 0;
0, 0, -1, 0, 0, 0,  0, 0;
1, 0,  0, 1, 0, 0,  0, 0;
0, 0,  0, 0, 0, 0, -1, 0];

% Elevation of points B, C, D, E
x = EA + J * l_hat;

% Cofactor matrix of x
Qxx = J * Qll * J'

ref_var = k' * f / 4;
fprintf('Reference variance = %.4f\n', ref_var)
 
% Covariance of x 
Cxx = ref_var * Qxx