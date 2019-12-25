%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 6.3 in
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 128-131.
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

% Unknowns:
% x_1 = x_1 coordinate of P
% x_2 = x_2 coordinate of P
% x_3 = x_1 coordinate of S_2
% x_4 = x_1 coordinate of S_3

% Condition equations
% L_1 x_2 - p x_1         = 0
% L_2 x_2 - p (L_4 - x_1) = 0
% L_3 x_2 - p (x_4 - x_1) = 0
% L_4 - x_3               = 0
% L_4 + L_5 - x_4         = 0 

% Initialize parameters
x = zeros(4,1);
x(2) = L(4)*p/(L(1)+L(2));
x(1) = L(1)*x(2)/p;
x(3) = L(4);
x(4) = L(4) + L(5);

% Set up linearized condition model: A v + B delta = f
A = [ ...
x(2),    0,    0,  0,  0;
   0, x(2),    0, -p,  0;
   0,    0, x(2),  0,  0;
   0,    0,    0,  1,  0;
   0,    0,    0,  1,  1];

B = [...
  -p, L(1),  0,  0; 
   p, L(2),  0,  0;
   p, L(3),  0, -p;
   0,    0, -1,  0;
   0,    0,  0, -1];

f = [ ...
 p*x(1)-L(1)*x(2); 
 p*(L(4)-x(1))-L(2)*x(2);
 p*(x(4)-x(1))-L(3)*x(2);
 x(3)-L(4);
 x(4)-L(4)-L(5)];

Qe = A*Q*A';
We = inv(Qe);

N = B'*We*B;
t = B'*We*f;
iN = inv(N); 

Del = iN * t;

x = x + Del;

fprintf('x1 = %.3f (m), x2 = %.3f (m), x3 = %.3f, x4 = %.3f\n', x(1), x(2), x(3), x(4));
