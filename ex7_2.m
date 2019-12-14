%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 7.2 in
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 141-143
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

% Table 6.1: observation value and standard deviation 
table6_1 

Q = diag(obs(:,2).^2);

A = [-obs(5,1), -obs(4,1)-obs(5,1), obs(4,1), obs(3,1)-obs(2,1), -obs(1,1)-obs(2,1)];

f = obs(1,1)*obs(5,1) + obs(2,1)*obs(4,1) + obs(2,1)*obs(5,1) - obs(3,1)*obs(4,1); % (mm * m)

% Normal equation matrix
Qe = A * Q * A'; % (mm^2 * m^2)

% Lagrange multiplier
k = f/Qe; % 1/(mm * m)

% observation residual
v = Q * A' * k;

fprintf('Adjusted observations\n');
o_hat = obs(:,1) + v

% Adjusted coordinates of P
p = 100; % mm
x1_hat = o_hat(1)*o_hat(4)/(o_hat(1)+o_hat(2));
x2_hat = p * o_hat(4) / (o_hat(1)+o_hat(2));

fprintf('Adjusted coordinates of P = (%.3f, %.3f) (m)\n', x1_hat, x2_hat)

We = 1/Qe;
I_Qvv = eye(5) - Q * A' * We * A;

fprintf('a posteriori cofactor matrix\n')
Qll = I_Qvv * I_Qvv * Q

fprintf('Reference variance = %.4f\n', k*f);

% Compute cofactor matrix of P

J = zeros(2,3);
J(1,1) =  x1_hat/o_hat(1) - x1_hat/(o_hat(1)+o_hat(2));
J(1,2) = -x1_hat/(o_hat(1)*o_hat(2));
J(1,3) =  x1_hat/o_hat(4);
J(2,1) = -x2_hat/(o_hat(1)+o_hat(2));
J(2,2) = J(2,1);
J(2,3) = x2_hat / p; 

% Cofactor matrix of sub-vector [1,2,4]
S =[ ...
1, 0, 0, 0, 0;
0, 1, 0, 0, 0;
0, 0, 0, 1, 0];

Qs = S * Qll * S';

fprintf('Cofactor matrix of P\n')
Qx = J * Qs * J'
