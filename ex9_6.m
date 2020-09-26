%% Eun Hwan Shin, 2020
%
% This file is the Octave implementation of Example 9.6: 
% Constraints with added parameters. 
% 
% Mikhail, E. M. and Ackermann, F. (1976). Observations and Least Squares,
%   Thomas Y. Crowell Company, Inc., pp. 238-248.
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

% total number iterations
% set this to 1 to see the precision after 1st iteration
n_iter = 3; 
apply_constraints = 1; 
by_elimination = 0; % apply constraints by elimination

% Table 9-1:

p = 100;  % principal distance of the camera (mm) 
si = 0.1; % standard deviation of the image coordinates (mm)

% image coordinates on S1, S2, S3 (mm)
obs = [ ...
  14.1,  6.1, 22.1;  % point A
  16.6,  7.1, 26.3;  % point B
   6.5, 15.2, 32.5;  % point C
   8.3, 12.3, 28.7]; % point D
   
L4 = 10;   % distance measurement between S1 and S2 (m)
L5 =  8;   % distance measurement between S2 and S3 (m)
sd = 0.05; % standard deviation of distance measurement (m)

% approximate coordinates of A, B, C, D
x = [ ...
  8,  7,  4,  4;  % x-coordinates
 51, 41, 47, 48]; % y-coordinates

station = 'ABCD'; % station code

if apply_constraints
  fprintf('Exercise 9.6: part (b) by ');
  if by_elimination
    fprintf(' elimination\n');
  else
    fprintf(' direct solution\n');
  end
else
  fprintf('Exercise 9.6: part (a)\n');
end

% Part (a) Compute the coordinates of A, B, C, D
% Part (b) Compute the coordinates of A, B, C, D with the constraints
%     (x(1,i) - x_10)^2 + (x(2,i) - x_20)^2 - R^2 = 0
% where x_10, x_20 and R are three added parameters; thus the 4 points
% will be on a circle.

x_10 = 7; x_20 = 46; radius = 4;

% a priori cofactor matrix
Q = diag([ones(12,1)*si^2; ones(2,1)*sd^2]);

A = zeros(12,14);
B = zeros(12,8);
f = zeros(12,1);

iter = 0;

while iter < n_iter,
  iter = iter + 1;
  fprintf('\n== Iteration %d ==\n', iter);
  
  k = 1;
  for i=1:4,
    A(k,k) = x(2,i);
    A(k+1,k+1) = x(2,i);
    A(k+2,k+2) = x(2,i);
    A(k:k+2,13:14) = [ 0, 0; -p, 0; -p, -p];
  
    B(k:k+2,2*i-1:2*i) = [-p, obs(i,1); p, obs(i,2); p, obs(i,3)];
  
    f(k) = p * x(1,i) - obs(i,1) * x(2,i);
    f(k+1) = p * L4 - obs(i,2) * x(2,i) - p * x(1,i);
    f(k+2) = p * L4 + p * L5 - obs(i,3) * x(2,i) - p * x(1,i);  
    k = k + 3;
  end
  
  QAt = Q * A';
  Qe = A * QAt;
  We = inv(Qe);

  if ~apply_constraints || ~by_elimination
    N = B' * We * B;
    t = B' * We * f;
    invN = inv(N);
  
    delta = invN * t;
  end
  
  if apply_constraints
    % apply constraints with added parameters
    D1 = zeros(4,8);
    D2 = zeros(4,3);
    h = zeros(4,1);
    
    k = 1;
    for i=1:4,
      dx1 = (x(1,i)-x_10);
      dx2 = (x(2,i)-x_20);
      D1(i,k) = 2*dx1; D1(i,k+1) = 2*dx2;
      D2(i,1) = -2*dx1; D2(i,2) = -2*dx2; D2(i,3) = -2*radius;
      h(i) = radius^2 - dx1^2 - dx2^2;
      k = k + 2;
    end
    
    if by_elimination
      B1 = B(:,1); B2 = B(:,2:end);
      D11 = D1(1,:); D12 = D2(1,:);
      D21 = D1(2:end,:); D22 = D2(2:end,:);
      invD22 = inv(D22);
      C = D11 - D12 * invD22 * D21;
      g = h(1) - D12 * invD22 * h(2:end);
      C12 = C(2:end);
      
      invC11 = 1/C(1);
      B_bar = -B1 * invC11 * C12 + B2;
      f_bar = -B1 * invC11 * g + f;
      N_bar = B_bar' * B_bar;
      t_bar = B_bar' * f_bar;
      invN_bar = inv(N_bar);
      delta2 = invN_bar * t_bar;
      delta1 = invC11 * (g - C12 * delta2);
      delta = [delta1; delta2];
      dp = invD22 * (h(2:end) - D21 * delta);

    else
      P = D1 * invN * D1';
      invP = inv(P);
      R = D2' * invP * D2;
      invR = inv(R);
      h_D1d = h - D1 * delta;
      dp = invR * (D2' * invP * h_D1d);
      delta = delta + invN * D1' * invP * (h_D1d - D2 * dp);
    end

    % correction of the added parameters
    fprintf('dp = [%.5f, %.5f, %.5f]\n', dp(1), dp(2), dp(3));
    x_10 = x_10 + dp(1);
    x_20 = x_20 + dp(2);
    radius = radius + dp(3);
  end % if apply_constraints
  
  fprintf('delta = [');
  for i=1:8,
    fprintf('%.3f ', delta(i));
  end
  fprintf(']\n');
  
  % correction of the 4 coordinates
  k = 1;
  for i=1:4,
    x(:,i) = x(:,i) + delta(k:k+1,1);
    k = k + 2;
  end

end   % while iter < n_iter,

fprintf('\nSolution: ');
for i = 1:4,
  fprintf('%c(%.3f,%.3f) ', station(i), x(1,i), x(2,i));
end

fprintf('\nAdded Parameters: x_10 = %.3f, x_20 = %.3f, R = %.3f\n', x_10, x_20, radius);

r = 4; % degrees of freedom

if ~apply_constraints || ~by_elimination
  Qdd = invN;
end

if apply_constraints
  if by_elimination
    % Eq. (9.29)
    Qd1d2 =  -invC11 * C12 * invN_bar;
    Qdd = [invC11*C12*invN_bar*C12'*invC11', Qd1d2;
           Qd1d2', invN_bar];
    % Eq. (9.68)
    Qdpdp = invD22 * D21 * Qdd * D21' * invD22';
  else
    Qdpdp = invR;
    D1t_invP = D1' * invP;
    D1_invN = D1 * invN;
    Qdd = Qdd * (eye(8) - D1t_invP * D1_invN + D1t_invP * D2 * invR * D2' * invP * D1_invN);
  end
  fprintf('\nQdpdp * 1000\n');
  Qdpdp * 1000
  
   % added 4 constrains with 3 more parameters
   r = r + 1;
   
end

% precision estimation

QAtWe = Q * A' * We;
QAtWeB = QAtWe * B;
Qvv = QAtWe * QAt' - QAtWeB * Qdd * QAtWeB'; 
Qll = Q - Qvv;
fprintf('Qdd * 1000\n');
Qdd * 1000
fprintf('Qll * 1000\n');
Qll * 1000
  
% reference variance
k = We * (f - B * delta);
v = QAt * k;
ref_var = v' * inv(Q) * v / r;

fprintf('sigma_0^2 = %.4f\n', ref_var);