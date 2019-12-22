%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 6.2 in
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 124-127.
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

% Load x and y coordinates in Table 6.2:  
table6_2 

% Cofacror matrix
Q = 0.01 * eye(6);

% initial values of the transformation parameters
a = 1;
b = 2; 

T = [a -b; b a];

% construct linearized model: A v + B delta = f 
A = [ ...
  T, zeros(2,4); 
  zeros(2), T, zeros(2); 
  zeros(2,4), T];
 
B = [ ...
  x(1,1), -x(2,1); 
  x(2,1),  x(1,1);
  x(1,2), -x(2,2); 
  x(2,2),  x(1,2);
  x(1,3), -x(2,3); 
  x(2,3),  x(1,3)];

f = [ ...
  y(:,1) - T * x(:,1);
  y(:,2) - T * x(:,2);
  y(:,3) - T * x(:,3)];
  
% Run conditions-only LSQ

Qe = A*Q*A';
We = inv(Qe);

N = B'*We*B;
t = B'*We*f;
Qx = inv(N); 
delta = Qx * t;

a_hat = a + delta(1);
b_hat = b + delta(2);
fprintf('Adjusted parameters: a = %.6f, b = %.6f\n', a_hat, b_hat)

k = We*(-B*delta+f); % Lagrange multipliers
v = Q*A'*k;
Qv = Q*A'*(We - We*B*Qx*B'*We)*A*Q;

vWv = f'*We*f - delta'*t;
fprintf('vWv = %.6f\n', vWv);

dof = 4;

ref_var = vWv / dof;

fprintf('Reference variance = %.6f\n', ref_var);

Ql = Q - Qv;



