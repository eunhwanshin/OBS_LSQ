%% Eun Hwan Shin, 2019
%
% This file is an Octave implementation of Example 8.7.
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 186-199.
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

%------------------------------------------------------------------------------
% observations in Table 8.1 for the transformation:
% y = [a, -b; b, a] x

x = [ ...
  0, 1, 1, 1;   % x_1
  1, 0, 1, 2];  % x_2

y = [...
  -2.1, 1.0, -0.9;  % y_1
   1.1, 2.0,  2.8]; % y_2

% Cofactor matrix of each point with a-priori reference variance 0.01
Qi = diag([1, 1, 4, 4]);
ref_var = 0.01;

%------------------------------------------------------------------------------
% Tasks to do:
% 1. Compute the estimates of a and b and their covariance matrix.
% 2. Compute y coordinates of point 4 and its covariance.
% 3. Use point 1 to analyze the question of final coordinates and covariance 
%    matrix of the points in the y-coordinate system.

%% Task 1:

% Case 1: initialize parameters with zero
%a = 0; b = 0;

% Case 2:
% 
a = 1, b =2;

% Set up condition equations: A v + B \delta = f

iter = 0;

while iter < 5,
  iter = iter + 1;
  
  T = [ ...
      a, -b;
      b,  a];
	  
  Ai = [ ...
    a, -b, -1,  0;
    b,  a,  0, -1];

  Qei = Ai * Qi * Ai';
  Wei = inv(Qei);
  
  N = zeros(2);
  t = zeros(2,1);

  for i=1:3,
    Bi = [...
      x(1,i), -x(2,i);
      x(2,i),  x(1,i)];
	
    fi = y(:,i) - T * x(:,i);

    N = N + Bi' * Wei * Bi;
    t = t + Bi' * Wei * fi;  
  end

  Qd = inv(N);
  delta = Qd * t;

  a = a + delta(1);
  b = b + delta(2);

  fprintf('Iteration %d: a = %.3f, b = %.3f\n', iter, a, b);
end

fprintf('Covariance of parameters:\n');
Cd = ref_var * Qd

%% Task 2:
fprintf('Transformed coordinates of P4:\n'); 
y4 = T * x(:,4)

J = [...
  x(1,4), -x(2,4), a, -b;
  x(2,4),  x(1,4), b,  a];
  
fprintf('Covariane of y4\n');
Cy4 = J * [Cd, zeros(2); zeros(2), ref_var*eye(2)] * J'  

%% Task 3:
fprintf('Transformed coordindates of point 1\n');
y1 = T * x(:,1)

B1 = [...
  x(1,1), -x(2,1);
  x(2,1),  x(1,1)];

Qx1 = eye(2);	  
Qdx1 = -Qd * B1' * Wei * T * Qx1;

fprintf('Covariance of delta and x1\n');
Cdx1 = ref_var * Qdx1
Cx1 = ref_var * Qx1;

fprintf('Define an auxilisty vector t.\n')
t = [a; b; x(1,1); x(2,1)]

fprintf('Covariance matrix of vector t\n');

Ct = [ ...
     Cd, Cdx1;
  Cdx1', Cx1]
  
J = [...
  x(1,1), -x(2,1), a, -b;
  x(2,1),  x(1,1), b,  a];   

fprintf('Covariance of y1\n')
Cy1 = J * Ct * J' 

%% Best estimate of y1

ell = [y1; y(:,1)]

Qy1 =4*eye(2);

% cofactor between estimated y1 and observed y
Qypy1 = B1 * Qd * B1' * Wei * Qy1

% Cofactor of ell
Qell = [...
  Cy1./ref_var, Qypy1;
  Qypy1', Qy1]
  
Well = inv(Qell)

% To obtain the best estimate of y1 we do adjustment of
% indirect observations
%    v1 - y_hat = -y1;
%    v2 - y_hat = -y

By1 = -[eye(2); eye(2)];

Ny1 = By1' * Well * By1 

ty1 = -By1' * Well * ell

Qy1_hat = inv(Ny1)
y1_hat = Qy1_hat * ty1

Cy1_hat = ref_var * Qy1_hat

