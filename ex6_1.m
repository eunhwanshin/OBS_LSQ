%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 6.1 in
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 120
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

L = obs(:,1);
Q = diag(obs(:,2).^2);

p = 100; % focal length (mm) 

% Initialize parameters
x = zeros(2,1);
x(2) = L(4)*p/(L(1)+L(2));
x(1) = L(1)*x(2)/p;

% Set up linearized condition model: A v + B delta = f
A = [ ...
x(2),    0,    0,  0,  0;
   0, x(2),    0, -p,  0;
   0,    0, x(2), -p, -p];

B = [-p, L(1); p, L(2); p, L(3)];

f = [ ...
 p*x(1)-L(1)*x(2); 
 p*(L(4)-x(1))-L(2)*x(2);
 p*(L(4)+L(5)-x(1))-L(3)*x(2)]

Qe = A*Q*A';
We = inv(Qe);

N = B'*We*B;
t = B'*We*f;
iN = inv(N); 

Del = iN * t;

x = x + Del;

fprintf('x1 = %.3f (m), x2 = %.3f (m)\n', x(1), x(2));

k = We*(-B*Del+f);
v = Q*A'*k;

Lhat = L + v;

% cofactor matrix of the residual
Qvv = Q*A'*(We - We*B*iN*B'*We)*A*Q;

% Note that the degree of freedom is 1.
ref_var = v'* inv(Q) * v;

fprintf('sigma_0^2 = %.5f\n', ref_var);


