%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 10.2: 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 263-266.
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

%% Part 1
% Anlge measurements of a triangle in (arc-min)
L = [ ...
  60*60+00;   % A
  60*60+02;   % A
  60*60+00;   % B
  60*60+01];  % C
Q = eye(4);

% condition equations
% L(1) - L(2) = 0
% L(1) + L(3) + L(4) - pi = 0

A = [ ...
  1, -1, 0, 0;
  1,  0, 1, 1]
f = [0;180*60] - A*L % in arc-minutes


Qe = A * Q * A'
We = inv(Qe);
k = We * f;
v = Q * A' * k;
Qvv = A' * We * A;
fprintf('k * 5 \n')
k * 5
fprintf('v * 5 \n')
v * 5
fprintf('Qvv * 5 \n')
Qvv * 5
Qll = Q - Qvv;
fprintf('Qll * 5 \n')
Qll * 5

%% Part 2: adjustment by steps
% 1. Adjust L(1) and L(2)

A1 = [1, -1, 0 , 0];
Qe1 = A1 * A1';
We1 = 1/Qe1;
k1 = We1 * f(1);
v1 = A1' * k1

L1 = L + v1

% cofactor of the measurements after first step.
Q1 = eye(4) - A1'*We1*A1

% Use derived measurements

D = [ ...
  0.5, 0.5, 0, 0;
  0,   0,   1, 0;
  0,   0,   0, 1];
  
Qmm = D * Q1 * D'

m = D * L1

% Apply 2nd condition on the derived measurements
Am = [1, 1, 1]
fm = 180*60 - Am * m
Qe2 = Am * Qmm * Am'
We2 = 1/Qe2
km = We2 * fm
vm = Qmm * Am' * km

% Redidual to the derived measurements
v2 = Q1 * D' * Am' * km

% Note that v and vL are identical.
vL = v1 + v2

% Check cofactor matrix
Qv2v2 = Q1 * D' * Am' * We2 * Am * D * Q1
Qll2 = Q1 - Qv2v2;

fprintf('Qll2 * 5 = \n')
Qll2 * 5

fprintf('Qll and Qll2 are identical.\n')

