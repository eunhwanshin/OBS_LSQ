%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 10.4: 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 270-273.
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

% Anlge measurements of two adjacent triangle in (degrees, minutes, seconds)

L_dms = [ ...
  60, 00, 03;
  60, 00, 02;
  60, 00, 01;
 120, 00, 05;
 120, 00, 05;
  60, 00, 01];

% Angles in arc-sec
L = L_dms(:,1)*3600 + L_dms(:,2)*60 + L_dms(:,3);
Q = eye(6);

% condition equations:
% L(1) + L(2) + L(3) - 180 = 0
% L(2) + L(4) + L(5) + L(6) - 360 = 0

%% Part 1: Direct Adjustment

A = [ ...
  1, 1, 1, 0, 0, 0;
  0, 1, 0, 1, 1, 1]
f = [180*3600;360*3600] - A*L % in arc-sec

Qe = A * Q * A'
We = inv(Qe);
k = We * f
v = Q * A' * k
Qvv = A' * We * A;

fprintf('Qvv * 11 \n')
Qvv * 11
Qll = Q - Qvv;
fprintf('Qll * 11 \n')
Qll * 11

%% Part 2: adjustment by steps
% Step 1: Adjust L(1) and L(2)

A1 = A(1,:);
Qe1 = A1 * Q * A1';
We1 = 1/Qe1;
k1 = We1 * f(1);
v1 = A1' * k1

L1 = L + v1

% cofactor of the measurements after first step.
Q1 = eye(6) - A1'*We1*A1;
fprintf('Q1 * 3 = \n')
Q1 * 3

% Step 2: apply 2nd condition
A2 = A(2,:);
f2 = 360*3600 - A2 * L1

Qe2 = A2 * Q1 * A2';

We2 = 1/Qe2;
k2 = We2 * f2
v2 = Q1 * A2' * k2

% Note that v and vL are identical.
vL = v1 + v2

% Check cofactor matrix
Qv2v2 = Q1 * A2' * We2 * A2 * Q1;
fprintf('Qv2v2 * 33 = \n')
Qv2v2 * 33

Qll2 = Q1 - Qv2v2;
fprintf('Qll2 * 11 = \n')
Qll2 * 11

fprintf('Qll and Qll2 are identical.\n')

